import csv
import heapq
import math
import os
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RESULT_CHUNK_ROWS = 250_000
DEFAULT_READ_CHUNK_ROWS = 100_000
TARGET_BUCKET_BYTES = 32 * 1024 * 1024
MAX_BUCKETS = 256
GENE_ID_DTYPES = {'gene_A': str, 'gene_B': str}


def prepare_hyper_result(df):
    """Add the existing SCRIN post-processing columns without row-wise apply."""
    if df is None or df.empty:
        return df

    output = df.copy()
    gene_a = output['gene_A'].astype(str)
    gene_b = output['gene_B'].astype(str)
    output['pair'] = np.where(gene_a <= gene_b, gene_a + '_' + gene_b, gene_b + '_' + gene_a)
    output['enrichment_ratio'] = (
        (output['gene_B_around'] / output['gene_around']) /
        (output['gene_B_slice'] / output['gene_slice'])
    )
    output['support_ratio'] = (
        output['gene_B_around'] /
        np.minimum(output['gene_A_N'], output['gene_B_N'])
    )
    return output


def natural_sort_key(path):
    """Sort numbered intermediate files by their numeric components."""
    return [int(token) if token.isdigit() else token.lower()
            for token in re.split(r'(\d+)', str(path))]


def collect_qvalue_pairs(file_list, qvalue_threshold, chunksize=DEFAULT_READ_CHUNK_ROWS):
    """Collect directed gene pairs eligible for distance-shape calculation."""
    eligible_pairs = set()
    for file_path in sorted((str(path) for path in file_list), key=natural_sort_key):
        for chunk in pd.read_csv(
            file_path,
            usecols=['gene_A', 'gene_B', 'qvalue_BH'],
            dtype={'gene_A': str, 'gene_B': str},
            chunksize=chunksize,
            float_precision='round_trip',
        ):
            selected = chunk[chunk['qvalue_BH'] < qvalue_threshold]
            eligible_pairs.update(zip(selected['gene_A'], selected['gene_B']))
    return eligible_pairs

class ChunkedCSVWriter:
    """Write incoming DataFrames immediately, rotating at a bounded row count."""

    def __init__(self, output_dir, prefix, max_rows=DEFAULT_RESULT_CHUNK_ROWS, transform=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.prefix = prefix
        self.max_rows = max_rows
        self.transform = transform
        self.file_index = 0
        self.rows_in_file = 0
        self.paths = []

    def write(self, df):
        if df is None or not isinstance(df, pd.DataFrame) or df.empty:
            return

        if self.transform is not None:
            df = self.transform(df)

        start = 0
        while start < len(df):
            if self.rows_in_file == self.max_rows:
                self.file_index += 1
                self.rows_in_file = 0

            row_count = min(self.max_rows - self.rows_in_file, len(df) - start)
            part = df.iloc[start:start + row_count]
            path = self.output_dir / f'{self.prefix}_{self.file_index}.csv'
            new_file = self.rows_in_file == 0
            part.to_csv(
                path, mode='w' if new_file else 'a', header=new_file, index=False,
                lineterminator='\n',
            )
            if new_file:
                self.paths.append(str(path))
            self.rows_in_file += row_count
            start += row_count


def combine_csv_files(file_list, save_path):
    """Concatenate CSV files without constructing a combined DataFrame."""
    files = sorted((str(path) for path in file_list), key=natural_sort_key)
    if not files:
        return

    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    with open(save_path, 'w', newline='', encoding='utf-8') as output_file:
        wrote_header = False
        for file_path in files:
            with open(file_path, 'r', newline='', encoding='utf-8') as input_file:
                header = input_file.readline()
                if not header:
                    continue
                if not wrote_header:
                    output_file.write(header)
                    wrote_header = True
                shutil.copyfileobj(input_file, output_file, length=1024 * 1024)


def _bucket_count(file_list):
    total_size = sum(os.path.getsize(path) for path in file_list if os.path.exists(path))
    return max(1, min(MAX_BUCKETS, math.ceil(total_size / TARGET_BUCKET_BYTES)))


def _append_frame(path, df):
    exists = os.path.exists(path)
    df.to_csv(
        path, mode='a' if exists else 'w', header=not exists, index=False,
        lineterminator='\n',
    )


def _partition_files(file_list, bucket_dir, bucket_count, add_pair=False):
    Path(bucket_dir).mkdir(parents=True, exist_ok=True)
    columns = None

    for file_path in sorted(file_list, key=natural_sort_key):
        for chunk in pd.read_csv(
            file_path,
            dtype=GENE_ID_DTYPES,
            chunksize=DEFAULT_READ_CHUNK_ROWS,
            float_precision='round_trip',
        ):
            if add_pair and 'pair' not in chunk.columns:
                gene_a = chunk['gene_A'].astype(str)
                gene_b = chunk['gene_B'].astype(str)
                chunk['pair'] = np.where(gene_a <= gene_b, gene_a + '_' + gene_b, gene_b + '_' + gene_a)
            if columns is None:
                columns = list(chunk.columns)
            buckets = pd.util.hash_pandas_object(chunk['pair'], index=False).to_numpy() % bucket_count
            for bucket_id in np.unique(buckets):
                bucket_frame = chunk.loc[buckets == bucket_id]
                _append_frame(os.path.join(bucket_dir, f'bucket_{int(bucket_id)}.csv'), bucket_frame)

    return columns


def _numeric_sort_key(row):
    return float(row['qvalue_BH']), -float(row['enrichment_ratio'])


def _merge_sorted_files(file_list, output_path, columns):
    """K-way merge bucket outputs already sorted by q-value and enrichment."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    readers = []
    handles = []
    heap = []
    serial = 0

    try:
        for file_path in sorted(file_list, key=natural_sort_key):
            handle = open(file_path, 'r', newline='', encoding='utf-8')
            handles.append(handle)
            reader = csv.DictReader(handle)
            readers.append(reader)
            row = next(reader, None)
            if row is not None:
                heapq.heappush(heap, (*_numeric_sort_key(row), serial, len(readers) - 1, row))
                serial += 1

        with open(output_path, 'w', newline='', encoding='utf-8') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=columns, lineterminator='\n')
            writer.writeheader()
            while heap:
                _, _, _, reader_index, row = heapq.heappop(heap)
                writer.writerow(row)
                next_row = next(readers[reader_index], None)
                if next_row is not None:
                    heapq.heappush(
                        heap,
                        (*_numeric_sort_key(next_row), serial, reader_index, next_row),
                    )
                    serial += 1
    finally:
        for handle in handles:
            handle.close()


def stream_postprocess(hyper_files, output_path, final_output_path, filter_threshold,
                       keep='last', work_dir=None, shape_files=None, shape_output_path=None):
    """Create raw and post-processed outputs with memory bounded by one hash bucket."""
    hyper_files = sorted((str(path) for path in hyper_files), key=natural_sort_key)
    shape_files = sorted((str(path) for path in (shape_files or [])), key=natural_sort_key)
    if not hyper_files:
        raise ValueError('No hyper test results found. Please check the input data and parameters.')

    combine_csv_files(hyper_files, output_path)
    if shape_files and shape_output_path is not None:
        combine_csv_files(shape_files, shape_output_path)

    work_dir = Path(work_dir or f'{output_path}.postprocess')
    if work_dir.exists():
        shutil.rmtree(work_dir)
    hyper_bucket_dir = work_dir / 'hyper_buckets'
    shape_bucket_dir = work_dir / 'shape_buckets'
    selected_dir = work_dir / 'selected_buckets'
    selected_dir.mkdir(parents=True, exist_ok=True)

    bucket_count = _bucket_count(hyper_files)
    hyper_columns = _partition_files(hyper_files, hyper_bucket_dir, bucket_count)
    shape_columns = None
    if shape_files:
        shape_columns = _partition_files(shape_files, shape_bucket_dir, bucket_count, add_pair=True)

    output_columns = list(hyper_columns)
    if shape_columns:
        output_columns.extend(column for column in shape_columns
                              if column not in {'gene_A', 'gene_B', 'pair'} and column not in output_columns)

    selected_files = []
    for bucket_id in range(bucket_count):
        hyper_path = hyper_bucket_dir / f'bucket_{bucket_id}.csv'
        if not hyper_path.exists():
            continue

        bucket_df = pd.read_csv(
            hyper_path, dtype=GENE_ID_DTYPES, float_precision='round_trip'
        )
        if shape_columns:
            shape_path = shape_bucket_dir / f'bucket_{bucket_id}.csv'
            if shape_path.exists():
                shape_df = pd.read_csv(
                    shape_path, dtype=GENE_ID_DTYPES, float_precision='round_trip'
                ).drop(columns=['pair'])
                bucket_df = pd.merge(bucket_df, shape_df, on=['gene_A', 'gene_B'], how='left')
            else:
                for column in output_columns:
                    if column not in bucket_df.columns:
                        bucket_df[column] = np.nan

        bucket_df = bucket_df.sort_values(
            by=['qvalue_BH', 'enrichment_ratio'], ascending=[True, False]
        )
        bucket_df = bucket_df.drop_duplicates(subset='pair', keep=keep).reset_index(drop=True)
        bucket_df = bucket_df[bucket_df['qvalue_BH'] < filter_threshold]
        if bucket_df.empty:
            continue

        bucket_df = bucket_df.sort_values(
            by=['qvalue_BH', 'enrichment_ratio'], ascending=[True, False]
        )
        selected_path = selected_dir / f'selected_{bucket_id}.csv'
        bucket_df.to_csv(
            selected_path, index=False, columns=output_columns, lineterminator='\n'
        )
        selected_files.append(str(selected_path))

    _merge_sorted_files(selected_files, final_output_path, output_columns)

