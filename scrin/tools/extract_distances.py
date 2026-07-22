import argparse
import csv
import gzip
import hashlib
import os
import pickle
import random
import re
import tempfile
from collections import defaultdict
from pathlib import Path

import joblib
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd
from tqdm import tqdm


DISTANCE_COLUMNS = ['gene_1', 'gene_2', 'pair', 'distance']
SUMMARY_COLUMNS = [
    'gene_1',
    'gene_2',
    'pair',
    'total_distance_count',
    'output_distance_count',
    'status',
]


def _positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('value must be a positive integer')
    return number


def _nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError('value must be a non-negative integer')
    return number


def _canonical_pair(gene_1, gene_2):
    gene_1 = str(gene_1)
    gene_2 = str(gene_2)
    if not gene_1 or not gene_2:
        raise ValueError('Gene names in a pair cannot be empty.')
    if gene_1 == gene_2:
        raise ValueError(f'Self-pair {gene_1}-{gene_2} is not supported.')
    return tuple(sorted((gene_1, gene_2)))


def _load_pairs(pair_entries, pair_path):
    if pair_entries is not None:
        input_pairs = pair_entries
    else:
        header = pd.read_csv(pair_path, nrows=0).columns
        if {'gene_A', 'gene_B'}.issubset(header):
            pair_columns = ['gene_A', 'gene_B']
        elif {'gene_1', 'gene_2'}.issubset(header):
            pair_columns = ['gene_1', 'gene_2']
        else:
            raise ValueError(
                'The pair file must contain gene_A/gene_B or gene_1/gene_2 columns.'
            )
        pair_df = pd.read_csv(
            pair_path,
            usecols=pair_columns,
            dtype={column: str for column in pair_columns},
            keep_default_na=False,
        )
        if pair_df.empty:
            raise ValueError('The pair file does not contain any pairs.')
        input_pairs = pair_df.itertuples(index=False, name=None)

    pairs = []
    seen = set()
    for gene_1, gene_2 in input_pairs:
        pair = _canonical_pair(gene_1, gene_2)
        if pair not in seen:
            pairs.append(pair)
            seen.add(pair)
    return pairs


def _natural_sort_key(path):
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r'(\d+)', str(path))
    ]


def _detect_background(intermediate_dir, requested_background):
    all_available = (
        (intermediate_dir / 'Distance_split').is_dir()
        and (intermediate_dir / 'molecule_distribution_ID.pkl').is_file()
    )
    cooccurrence_merge_dir = Path(f'{intermediate_dir}_gene_id_merge')
    cooccurrence_available = (
        (intermediate_dir / 'gene_rank_dict.pkl').is_file()
        and cooccurrence_merge_dir.is_dir()
    )

    if requested_background == 'all':
        if not all_available:
            raise ValueError('Global-background distance intermediates were not found.')
        return 'all'
    if requested_background == 'cooccurrence':
        if not cooccurrence_available:
            raise ValueError('Co-occurrence distance intermediates were not found.')
        return 'cooccurrence'

    available = []
    if all_available:
        available.append('all')
    if cooccurrence_available:
        available.append('cooccurrence')
    if len(available) == 1:
        return available[0]
    if not available:
        raise ValueError(
            'No supported distance intermediates were found. The original SCRIN run must '
            'enable --distribution_analysis and retain its intermediate directories.'
        )
    raise ValueError('Both distance layouts were found; specify --background explicitly.')


def _load_global_file(path):
    with open(path, 'rb') as handle:
        return pickle.load(handle)


def _load_cooccurrence_file(path):
    return joblib.load(path)


def _build_file_requests(intermediate_dir, background, pairs):
    file_requests = defaultdict(list)
    source_keys = {}
    unavailable = set()

    if background == 'all':
        with open(intermediate_dir / 'molecule_distribution_ID.pkl', 'rb') as handle:
            task_gene_dict = pickle.load(handle)
        gene_to_task = {
            str(gene): task_id
            for task_id, genes in task_gene_dict.items()
            for gene in genes
        }
        for pair in pairs:
            gene_1, gene_2 = pair
            source_gene, target_gene = gene_1, gene_2
            task_id = gene_to_task.get(source_gene)
            if task_id is None and gene_2 in gene_to_task:
                source_gene, target_gene = gene_2, gene_1
                task_id = gene_to_task[gene_2]
            if task_id is None:
                unavailable.add(pair)
                continue
            file_path = intermediate_dir / 'Distance_split' / f'task_results{task_id}.pkl'
            if not file_path.is_file():
                unavailable.add(pair)
                continue
            source_keys[pair] = (source_gene, target_gene)
            file_requests[file_path].append(pair)
        loader = _load_global_file
    else:
        gene_rank_dict = joblib.load(intermediate_dir / 'gene_rank_dict.pkl')
        gene_to_rank = {
            str(gene): rank_id
            for rank_id, genes in gene_rank_dict.items()
            for gene in genes
        }
        merge_root = Path(f'{intermediate_dir}_gene_id_merge')
        for pair in pairs:
            gene_1, gene_2 = pair
            source_gene, target_gene = gene_1, gene_2
            rank_id = gene_to_rank.get(source_gene)
            if rank_id is None and gene_2 in gene_to_rank:
                source_gene, target_gene = gene_2, gene_1
                rank_id = gene_to_rank[gene_2]
            if rank_id is None:
                unavailable.add(pair)
                continue
            rank_dir = merge_root / f'rank_{rank_id}'
            rank_files = sorted(rank_dir.glob('*.pkl'), key=_natural_sort_key)
            if not rank_files:
                unavailable.add(pair)
                continue
            source_keys[pair] = (source_gene, target_gene)
            for file_path in rank_files:
                file_requests[file_path].append(pair)
        loader = _load_cooccurrence_file

    ordered_requests = dict(sorted(file_requests.items(), key=lambda item: _natural_sort_key(item[0])))
    return ordered_requests, source_keys, unavailable, loader


def _get_distances(file_data, source_key):
    source_gene, target_gene = source_key
    return file_data.get(source_gene, {}).get(target_gene, [])


def _count_distances(file_requests, source_keys, loader, pairs):
    counts = {pair: 0 for pair in pairs}
    for file_path, requested_pairs in tqdm(file_requests.items(), desc='Counting distances'):
        file_data = loader(file_path)
        for pair in requested_pairs:
            counts[pair] += len(_get_distances(file_data, source_keys[pair]))
        del file_data
    return counts


def _pair_seed(random_seed, pair):
    payload = f'{random_seed}\0{pair[0]}\0{pair[1]}'.encode('utf-8')
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], 'big')


def _prepare_selection(pairs, counts, unavailable, min_distance_count, max_distance_count, random_seed):
    selected_indices = {}
    summary = {}

    for pair in pairs:
        total = counts[pair]
        if pair in unavailable or total == 0:
            status = 'not_found'
            output_count = 0
            selected_indices[pair] = set()
        elif total < min_distance_count:
            status = 'below_minimum'
            output_count = 0
            selected_indices[pair] = set()
        elif max_distance_count is not None and total > max_distance_count:
            status = 'sampled'
            output_count = max_distance_count
            rng = random.Random(_pair_seed(random_seed, pair))
            selected_indices[pair] = set(rng.sample(range(total), max_distance_count))
        else:
            status = 'complete'
            output_count = total
            selected_indices[pair] = None

        summary[pair] = {
            'gene_1': pair[0],
            'gene_2': pair[1],
            'pair': f'{pair[0]}_{pair[1]}',
            'total_distance_count': total,
            'output_distance_count': output_count,
            'status': status,
        }
    return selected_indices, summary


class _CSVWriter:
    def __init__(self, path, compressed=False):
        if compressed:
            self.handle = gzip.open(path, 'wt', newline='', encoding='utf-8')
        else:
            self.handle = open(path, 'w', newline='', encoding='utf-8')
        self.writer = csv.writer(self.handle)
        self.writer.writerow(DISTANCE_COLUMNS)

    def write(self, pair, distances):
        pair_name = f'{pair[0]}_{pair[1]}'
        self.writer.writerows((pair[0], pair[1], pair_name, distance) for distance in distances)

    def close(self):
        self.handle.close()


class _ParquetWriter:
    def __init__(self, path, batch_size=100_000):
        self.schema = pa.schema([
            ('gene_1', pa.string()),
            ('gene_2', pa.string()),
            ('pair', pa.string()),
            ('distance', pa.float64()),
        ])
        self.writer = pq.ParquetWriter(path, self.schema, compression='snappy')
        self.batch_size = batch_size
        self.buffer = {column: [] for column in DISTANCE_COLUMNS}

    def write(self, pair, distances):
        pair_name = f'{pair[0]}_{pair[1]}'
        for distance in distances:
            self.buffer['gene_1'].append(pair[0])
            self.buffer['gene_2'].append(pair[1])
            self.buffer['pair'].append(pair_name)
            self.buffer['distance'].append(float(distance))
            if len(self.buffer['distance']) >= self.batch_size:
                self._flush()

    def _flush(self):
        if not self.buffer['distance']:
            return
        table = pa.Table.from_pydict(self.buffer, schema=self.schema)
        self.writer.write_table(table)
        self.buffer = {column: [] for column in DISTANCE_COLUMNS}

    def close(self):
        self._flush()
        self.writer.close()


def _resolve_output_format(save_path, output_format):
    if output_format != 'auto':
        return output_format
    name = save_path.name.lower()
    if name.endswith('.parquet'):
        return 'parquet'
    if name.endswith('.csv') or name.endswith('.csv.gz'):
        return 'csv'
    raise ValueError('Cannot infer output format; use a .csv, .csv.gz, or .parquet path.')


def _write_distances(file_requests, source_keys, loader, selected_indices, save_path,
                     output_format):
    save_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        prefix=f'.{save_path.name}.', suffix='.tmp', dir=save_path.parent, delete=False
    )
    temporary.close()
    temp_path = Path(temporary.name)
    compressed = save_path.name.lower().endswith('.csv.gz')
    writer = None
    offsets = {pair: 0 for pair in selected_indices}

    try:
        if output_format == 'parquet':
            writer = _ParquetWriter(temp_path)
        else:
            writer = _CSVWriter(temp_path, compressed=compressed)

        for file_path, requested_pairs in tqdm(file_requests.items(), desc='Writing distances'):
            file_data = loader(file_path)
            for pair in requested_pairs:
                selection = selected_indices[pair]
                distances = _get_distances(file_data, source_keys[pair])
                start = offsets[pair]
                if selection is None:
                    writer.write(pair, distances)
                elif selection:
                    chosen = (
                        distance for index, distance in enumerate(distances, start=start)
                        if index in selection
                    )
                    writer.write(pair, chosen)
                offsets[pair] += len(distances)
            del file_data

        writer.close()
        writer = None
        os.replace(temp_path, save_path)
    except BaseException:
        if writer is not None:
            writer.close()
        if temp_path.exists():
            temp_path.unlink()
        raise


def _default_summary_path(save_path):
    name = save_path.name
    if name.lower().endswith('.csv.gz'):
        stem = name[:-7]
    else:
        stem = save_path.stem
    return save_path.with_name(f'{stem}_summary.csv')


def _write_summary(summary, pairs, summary_path):
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        mode='w', newline='', encoding='utf-8', prefix=f'.{summary_path.name}.',
        suffix='.tmp', dir=summary_path.parent, delete=False
    )
    temp_path = Path(temporary.name)
    try:
        writer = csv.DictWriter(temporary, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for pair in pairs:
            writer.writerow(summary[pair])
        temporary.flush()
        os.fsync(temporary.fileno())
        temporary.close()
        os.replace(temp_path, summary_path)
    except BaseException:
        temporary.close()
        if temp_path.exists():
            temp_path.unlink()
        raise


def _print_summary(summary, pairs, background, save_path, summary_path):
    print(f'Distance intermediate layout: {background}')
    print('Pair extraction summary:')
    for pair in pairs:
        row = summary[pair]
        print(
            f"  {row['pair']}: total={row['total_distance_count']}, "
            f"output={row['output_distance_count']}, status={row['status']}"
        )
    print(f'Distances saved to {save_path}')
    print(f'Summary saved to {summary_path}')


def main():
    parser = argparse.ArgumentParser(
        description='Extract distance observations for undirected gene pairs from SCRIN intermediates.'
    )
    parser.add_argument('--intermediate_dir', required=True)
    pair_group = parser.add_mutually_exclusive_group(required=True)
    pair_group.add_argument(
        '--pair', action='append', nargs=2, metavar=('GENE_1', 'GENE_2'),
        help='An undirected pair; repeat --pair to extract multiple pairs.',
    )
    pair_group.add_argument('--pair_path', help='CSV containing gene_A/gene_B or gene_1/gene_2.')
    parser.add_argument('--save_path', required=True)
    parser.add_argument(
        '--output_format', choices=['auto', 'csv', 'parquet'], default='auto',
        help='Default: infer from .csv, .csv.gz, or .parquet.',
    )
    parser.add_argument(
        '--background', choices=['auto', 'all', 'cooccurrence'], default='auto',
        help='Default: detect the intermediate layout automatically.',
    )
    parser.add_argument(
        '--min_distance_count', type=_positive_int, default=1,
        help='Discard pairs with fewer distance observations. Default: 1.',
    )
    parser.add_argument(
        '--max_distance_count', type=_positive_int, default=None,
        help='Uniformly sample pairs with more distance observations. Default: no maximum.',
    )
    parser.add_argument('--random_seed', type=_nonnegative_int, default=42)
    parser.add_argument('--summary_path', default=None)
    args = parser.parse_args()

    if args.max_distance_count is not None and args.max_distance_count < args.min_distance_count:
        parser.error('--max_distance_count must be greater than or equal to --min_distance_count.')

    try:
        pairs = _load_pairs(args.pair, args.pair_path)
    except ValueError as error:
        parser.error(str(error))

    intermediate_dir = Path(args.intermediate_dir).resolve()
    save_path = Path(args.save_path).resolve()
    summary_path = (
        Path(args.summary_path).resolve()
        if args.summary_path is not None else _default_summary_path(save_path)
    )
    if save_path == summary_path:
        parser.error('--summary_path must be different from --save_path.')

    try:
        output_format = _resolve_output_format(save_path, args.output_format)
        background = _detect_background(intermediate_dir, args.background)
        file_requests, source_keys, unavailable, loader = _build_file_requests(
            intermediate_dir, background, pairs
        )
        counts = _count_distances(file_requests, source_keys, loader, pairs)
        selected_indices, summary = _prepare_selection(
            pairs,
            counts,
            unavailable,
            args.min_distance_count,
            args.max_distance_count,
            args.random_seed,
        )
        _write_distances(
            file_requests,
            source_keys,
            loader,
            selected_indices,
            save_path,
            output_format,
        )
        _write_summary(summary, pairs, summary_path)
    except ValueError as error:
        parser.error(str(error))

    _print_summary(summary, pairs, background, save_path, summary_path)


if __name__ == '__main__':
    main()
