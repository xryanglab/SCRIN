import argparse
import csv
import math
import multiprocessing
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm


OUTPUT_COLUMNS = [
    'cell',
    'cell_x',
    'cell_y',
    'gene_A',
    'gene_B',
    'r_check',
    'gene_A_count',
    'gene_B_count',
    'colocalization_count',
]

_WORKER_PAIRS = None
_WORKER_RADIUS = None
_WORKER_Z_MODE = None


def _parse_column_names(column_name):
    column_names = [name.strip() for name in column_name.split(',')]
    if len(column_names) not in {4, 5} or any(not name for name in column_names):
        raise ValueError(
            "--column_name must contain x,y,geneID,cell or x,y,z,geneID,cell columns."
        )
    if len(set(column_names)) != len(column_names):
        raise ValueError("Duplicate names are not allowed in --column_name.")
    return column_names


def _canonical_pair(gene_a, gene_b):
    gene_a = str(gene_a)
    gene_b = str(gene_b)
    if not gene_a or not gene_b:
        raise ValueError("Gene names in a pair cannot be empty.")
    if gene_a == gene_b:
        raise ValueError(
            f"Self-pair {gene_a}-{gene_b} is not supported. Please provide two different genes."
        )
    return tuple(sorted((gene_a, gene_b)))


def _load_pairs(pair, pair_path):
    if pair is not None:
        input_pairs = [pair]
    else:
        pair_df = pd.read_csv(
            pair_path,
            usecols=['gene_A', 'gene_B'],
            dtype={'gene_A': str, 'gene_B': str},
            keep_default_na=False,
        )
        if pair_df.empty:
            raise ValueError("The pair file does not contain any gene pairs.")
        input_pairs = pair_df[['gene_A', 'gene_B']].itertuples(index=False, name=None)

    output_pairs = []
    seen = set()
    for gene_a, gene_b in input_pairs:
        canonical = _canonical_pair(gene_a, gene_b)
        if canonical not in seen:
            output_pairs.append(canonical)
            seen.add(canonical)
    return output_pairs


def _load_transcripts(data_path, column_name, target_genes):
    column_names = _parse_column_names(column_name)
    dtype = {column_names[0]: float, column_names[1]: float}
    if len(column_names) == 5:
        dtype.update({column_names[2]: float, column_names[3]: str, column_names[4]: str})
    else:
        dtype.update({column_names[2]: str, column_names[3]: str})

    df = pd.read_csv(
        data_path,
        usecols=column_names,
        dtype=dtype,
        na_values='',
        keep_default_na=False,
    )

    if len(column_names) == 5:
        rename_map = dict(zip(column_names, ['x', 'y', 'z', 'geneID', 'cell']))
    else:
        rename_map = dict(zip(column_names, ['x', 'y', 'geneID', 'cell']))
    df.rename(columns=rename_map, inplace=True)
    if len(column_names) == 4:
        df['z'] = 0.0

    if df.empty:
        raise ValueError("The input data file does not contain any transcripts.")
    if df[['x', 'y', 'z', 'geneID', 'cell']].isnull().any().any():
        raise ValueError("Input data contains null values in required columns.")
    if (df['geneID'].str.len() == 0).any() or (df['cell'].str.len() == 0).any():
        raise ValueError("Input data contains empty gene or cell identifiers.")
    if not np.isfinite(df[['x', 'y', 'z']].to_numpy()).all():
        raise ValueError("Input data contains non-finite coordinates.")

    present_genes = set(df['geneID'].unique())
    missing_genes = sorted(target_genes - present_genes)
    if missing_genes:
        raise ValueError(
            "Requested genes were not found in the input data: " + ', '.join(missing_genes)
        )

    cell_centers = (
        df.groupby('cell', sort=False, observed=True)[['x', 'y']]
        .median()
        .rename(columns={'x': 'cell_x', 'y': 'cell_y'})
    )
    target_df = df.loc[df['geneID'].isin(target_genes), ['cell', 'geneID', 'x', 'y', 'z']].copy()
    return cell_centers, target_df


def _initialize_worker(pairs, radius, z_mode):
    global _WORKER_PAIRS, _WORKER_RADIUS, _WORKER_Z_MODE
    _WORKER_PAIRS = pairs
    _WORKER_RADIUS = radius
    _WORKER_Z_MODE = z_mode


def _build_gene_trees(gene_ids, coordinates):
    gene_counts = {}
    gene_trees = {}

    for gene in set(gene_ids):
        gene_coordinates = coordinates[gene_ids == gene]
        gene_counts[gene] = len(gene_coordinates)
        if _WORKER_Z_MODE == 'continuous':
            gene_trees[gene] = cKDTree(gene_coordinates)
        else:
            z_trees = {}
            for z_value in np.unique(gene_coordinates[:, 2]):
                layer_coordinates = gene_coordinates[gene_coordinates[:, 2] == z_value, :2]
                z_trees[z_value] = cKDTree(layer_coordinates)
            gene_trees[gene] = z_trees

    return gene_counts, gene_trees


def _count_from_trees(gene_a, gene_b, gene_trees):
    if gene_a not in gene_trees or gene_b not in gene_trees:
        return 0

    if _WORKER_Z_MODE == 'continuous':
        return int(gene_trees[gene_a].count_neighbors(gene_trees[gene_b], _WORKER_RADIUS))

    count = 0
    trees_a = gene_trees[gene_a]
    trees_b = gene_trees[gene_b]
    for z_value in trees_a.keys() & trees_b.keys():
        count += int(trees_a[z_value].count_neighbors(trees_b[z_value], _WORKER_RADIUS))
    return count


def _process_cell(task):
    cell, cell_x, cell_y, gene_ids, coordinates = task
    if len(gene_ids):
        gene_counts, gene_trees = _build_gene_trees(gene_ids, coordinates)
    else:
        gene_counts, gene_trees = {}, {}

    rows = []
    for gene_a, gene_b in _WORKER_PAIRS:
        rows.append([
            cell,
            cell_x,
            cell_y,
            gene_a,
            gene_b,
            _WORKER_RADIUS,
            gene_counts.get(gene_a, 0),
            gene_counts.get(gene_b, 0),
            _count_from_trees(gene_a, gene_b, gene_trees),
        ])
    return rows


def _cell_task_iterator(cell_centers, target_df):
    grouped_indices = target_df.groupby('cell', sort=False, observed=True).indices

    for cell, cell_x, cell_y in cell_centers.itertuples(name=None):
        indices = grouped_indices.get(cell)
        if indices is None:
            gene_ids = np.empty(0, dtype=object)
            coordinates = np.empty((0, 3), dtype=float)
        else:
            cell_df = target_df.iloc[indices]
            gene_ids = cell_df['geneID'].to_numpy()
            coordinates = cell_df[['x', 'y', 'z']].to_numpy()
        yield cell, cell_x, cell_y, gene_ids, coordinates


def _write_results(cell_centers, target_df, pairs, radius, z_mode, save_path,
                   n_jobs, task_chunk_size):
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    temp_handle = tempfile.NamedTemporaryFile(
        mode='w',
        newline='',
        encoding='utf-8',
        prefix=f'.{save_path.name}.',
        suffix='.tmp',
        dir=save_path.parent,
        delete=False,
    )
    temp_path = Path(temp_handle.name)

    try:
        writer = csv.writer(temp_handle)
        writer.writerow(OUTPUT_COLUMNS)
        tasks = _cell_task_iterator(cell_centers, target_df)

        if n_jobs == 1:
            _initialize_worker(pairs, radius, z_mode)
            results = map(_process_cell, tasks)
            for rows in tqdm(results, total=len(cell_centers), desc='Processing cells'):
                writer.writerows(rows)
        else:
            context = multiprocessing.get_context()
            with context.Pool(
                processes=n_jobs,
                initializer=_initialize_worker,
                initargs=(pairs, radius, z_mode),
            ) as pool:
                results = pool.imap(_process_cell, tasks, chunksize=task_chunk_size)
                for rows in tqdm(results, total=len(cell_centers), desc='Processing cells'):
                    writer.writerows(rows)

        temp_handle.flush()
        os.fsync(temp_handle.fileno())
        temp_handle.close()
        os.replace(temp_path, save_path)
    except BaseException:
        temp_handle.close()
        if temp_path.exists():
            temp_path.unlink()
        raise


def _positive_finite_float(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("value must be a finite positive number")
    return number


def _positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return number


def main():
    parser = argparse.ArgumentParser(
        description="Count transcript-level colocalization events for gene pairs in each cell."
    )
    parser.add_argument('--data_path', required=True, help='Path to the transcript CSV file.')
    parser.add_argument('--save_path', required=True, help='Path for the output CSV file.')
    parser.add_argument(
        '--column_name',
        default='x,y,z,geneID,cell',
        help='Input columns ordered as x,y,z,geneID,cell or x,y,geneID,cell.',
    )
    pair_group = parser.add_mutually_exclusive_group(required=True)
    pair_group.add_argument('--pair', nargs=2, metavar=('GENE_A', 'GENE_B'))
    pair_group.add_argument(
        '--pair_path',
        help='CSV containing gene_A and gene_B columns; SCRIN result files are accepted.',
    )
    parser.add_argument('--r_check', required=True, type=_positive_finite_float)
    parser.add_argument(
        '--z_mode', choices=['discrete', 'continuous'], default='discrete',
        help="Use same-layer 2D distances or continuous 3D distances. Default: discrete.",
    )
    parser.add_argument(
        '--n_jobs', type=_positive_int, default=1,
        help='Number of local worker processes. Default: 1.',
    )
    parser.add_argument(
        '--task_chunk_size', type=_positive_int, default=10,
        help='Number of cell tasks sent to each worker at a time. Default: 10.',
    )
    args = parser.parse_args()

    save_path_resolved = Path(args.save_path).resolve()
    source_paths = [Path(args.data_path)]
    if args.pair_path is not None:
        source_paths.append(Path(args.pair_path))
    if any(path.resolve() == save_path_resolved for path in source_paths):
        parser.error("--save_path must be different from all input file paths.")

    pairs = _load_pairs(args.pair, args.pair_path)
    target_genes = {gene for pair in pairs for gene in pair}
    cell_centers, target_df = _load_transcripts(
        args.data_path, args.column_name, target_genes
    )

    print(f'Loaded {len(cell_centers)} cells and {len(pairs)} unique gene pair(s).')
    print(f'Using {args.n_jobs} local worker process(es).')
    _write_results(
        cell_centers,
        target_df,
        pairs,
        args.r_check,
        args.z_mode,
        args.save_path,
        args.n_jobs,
        args.task_chunk_size,
    )
    print(f'Results saved to {args.save_path}')


if __name__ == '__main__':
    main()
