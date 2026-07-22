import argparse
import csv
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


DEFAULT_QVALUE_THRESHOLDS = [0.05, 0.01, 0.001, 0.00001]
DEFAULT_READ_CHUNK_ROWS = 100_000
TARGET_BUCKET_BYTES = 32 * 1024 * 1024
MAX_BUCKETS = 256
REQUIRED_COLUMNS = ['gene_A', 'gene_B', 'pair', 'qvalue_BH', 'enrichment_ratio']
MATCH_MODES = {'exact', 'prefix', 'suffix', 'contains'}
SUMMARY_COLUMNS = [
    'qvalue_threshold',
    'total_significant_pairs',
    'control_group',
    'control_associated_pair_count',
    'control_associated_pair_fraction',
    'control_associated_percentage',
    'matched_control_feature_count',
]


def _positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('value must be a positive integer')
    return number


def _qvalue_threshold(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0 or number > 1:
        raise argparse.ArgumentTypeError('q-value thresholds must be in the interval (0, 1]')
    return number


def _parse_control_groups(entries):
    groups = {}
    for group_name, match_mode, pattern in entries:
        if not group_name:
            raise ValueError('Control group names cannot be empty.')
        if group_name == 'any_control':
            raise ValueError("The control group name 'any_control' is reserved.")
        if match_mode not in MATCH_MODES:
            raise ValueError(
                f"Invalid match mode '{match_mode}' for group '{group_name}'. "
                f"Choose from: {', '.join(sorted(MATCH_MODES))}."
            )
        if not pattern:
            raise ValueError(f"The match pattern for group '{group_name}' cannot be empty.")
        groups.setdefault(group_name, []).append((match_mode, pattern))
    return groups


def _match_name(name, rules, ignore_case=False):
    value = str(name)
    if ignore_case:
        value = value.casefold()

    for match_mode, pattern in rules:
        candidate = pattern.casefold() if ignore_case else pattern
        if match_mode == 'exact' and value == candidate:
            return True
        if match_mode == 'prefix' and value.startswith(candidate):
            return True
        if match_mode == 'suffix' and value.endswith(candidate):
            return True
        if match_mode == 'contains' and candidate in value:
            return True
    return False


def _match_series(series, rules, ignore_case=False):
    values = series.astype(str)
    if ignore_case:
        values = values.str.casefold()

    output = pd.Series(False, index=series.index)
    for match_mode, pattern in rules:
        candidate = pattern.casefold() if ignore_case else pattern
        if match_mode == 'exact':
            output |= values.eq(candidate)
        elif match_mode == 'prefix':
            output |= values.str.startswith(candidate)
        elif match_mode == 'suffix':
            output |= values.str.endswith(candidate)
        else:
            output |= values.str.contains(candidate, regex=False)
    return output


def _bucket_count(result_path):
    file_size = os.path.getsize(result_path)
    return max(1, min(MAX_BUCKETS, math.ceil(file_size / TARGET_BUCKET_BYTES)))


def _validate_chunk(chunk):
    if chunk[REQUIRED_COLUMNS].isnull().any().any():
        raise ValueError('The raw result contains null values in required columns.')
    if (
        (chunk['gene_A'].str.len() == 0).any()
        or (chunk['gene_B'].str.len() == 0).any()
        or (chunk['pair'].str.len() == 0).any()
    ):
        raise ValueError('The raw result contains empty gene or pair identifiers.')
    numeric = chunk[['qvalue_BH', 'enrichment_ratio']].to_numpy()
    if not np.isfinite(numeric).all():
        raise ValueError('The raw result contains non-finite q-values or enrichment ratios.')


def _partition_results(result_path, bucket_dir, bucket_count, control_groups,
                       ignore_case, read_chunk_size):
    bucket_dir.mkdir(parents=True, exist_ok=True)
    matched_features = {group_name: set() for group_name in control_groups}
    row_count = 0
    dtype = {'gene_A': str, 'gene_B': str, 'pair': str}

    try:
        chunks = pd.read_csv(
            result_path,
            usecols=REQUIRED_COLUMNS,
            dtype=dtype,
            chunksize=read_chunk_size,
        )
        for chunk in tqdm(chunks, desc='Partitioning raw results'):
            _validate_chunk(chunk)
            row_count += len(chunk)

            observed_features = set(chunk['gene_A']) | set(chunk['gene_B'])
            for group_name, rules in control_groups.items():
                matched_features[group_name].update(
                    feature for feature in observed_features
                    if _match_name(feature, rules, ignore_case)
                )

            bucket_ids = (
                pd.util.hash_pandas_object(chunk['pair'], index=False).to_numpy()
                % bucket_count
            )
            for bucket_id in np.unique(bucket_ids):
                bucket_frame = chunk.loc[bucket_ids == bucket_id]
                bucket_path = bucket_dir / f'bucket_{int(bucket_id)}.csv'
                exists = bucket_path.exists()
                bucket_frame.to_csv(
                    bucket_path,
                    mode='a' if exists else 'w',
                    header=not exists,
                    index=False,
                )
    except ValueError as error:
        if 'Usecols do not match columns' in str(error):
            raise ValueError(
                'The input must be a raw SCRIN result containing: '
                + ', '.join(REQUIRED_COLUMNS)
            ) from error
        raise

    if row_count == 0:
        raise ValueError('The raw SCRIN result does not contain any rows.')
    return matched_features


def _analyze_buckets(bucket_dir, bucket_count, control_groups, ignore_case,
                     thresholds, pair_keep, matched_features):
    group_names = list(control_groups)
    include_union = len(group_names) > 1
    output_groups = group_names + (['any_control'] if include_union else [])
    total_counts = {threshold: 0 for threshold in thresholds}
    control_counts = {
        threshold: {group_name: 0 for group_name in output_groups}
        for threshold in thresholds
    }

    for bucket_id in tqdm(range(bucket_count), desc='Analyzing pair buckets'):
        bucket_path = bucket_dir / f'bucket_{bucket_id}.csv'
        if not bucket_path.exists():
            continue

        bucket_df = pd.read_csv(
            bucket_path,
            dtype={'gene_A': str, 'gene_B': str, 'pair': str},
        )
        bucket_df = bucket_df.sort_values(
            by=['qvalue_BH', 'enrichment_ratio'], ascending=[True, False]
        )
        bucket_df = bucket_df.drop_duplicates(subset='pair', keep=pair_keep)

        group_masks = {}
        union_mask = pd.Series(False, index=bucket_df.index)
        for group_name, rules in control_groups.items():
            group_mask = (
                _match_series(bucket_df['gene_A'], rules, ignore_case)
                | _match_series(bucket_df['gene_B'], rules, ignore_case)
            )
            group_masks[group_name] = group_mask
            union_mask |= group_mask
        if include_union:
            group_masks['any_control'] = union_mask

        for threshold in thresholds:
            significant_mask = bucket_df['qvalue_BH'] < threshold
            total_counts[threshold] += int(significant_mask.sum())
            for group_name in output_groups:
                control_counts[threshold][group_name] += int(
                    (significant_mask & group_masks[group_name]).sum()
                )

    union_features = set().union(*matched_features.values()) if matched_features else set()
    summary_rows = []
    for threshold in thresholds:
        total = total_counts[threshold]
        for group_name in output_groups:
            control_count = control_counts[threshold][group_name]
            fraction = control_count / total if total else None
            if group_name == 'any_control':
                matched_count = len(union_features)
            else:
                matched_count = len(matched_features[group_name])
            summary_rows.append({
                'qvalue_threshold': threshold,
                'total_significant_pairs': total,
                'control_group': group_name,
                'control_associated_pair_count': control_count,
                'control_associated_pair_fraction': fraction,
                'control_associated_percentage': fraction * 100 if fraction is not None else None,
                'matched_control_feature_count': matched_count,
            })
    return summary_rows


def _print_summary(summary_rows, matched_features, pair_keep):
    print('--- SCRIN empirical false-positive baseline analysis ---')
    print(f'Pair deduplication: keep={pair_keep}')
    print('Control features observed in the raw result:')
    for group_name, features in matched_features.items():
        feature_text = ', '.join(sorted(features)) if features else '(none)'
        print(f'  {group_name} ({len(features)}): {feature_text}')
        if not features:
            print(f"  Warning: no features matched control group '{group_name}'.")
    print()

    headers = ['Threshold', 'Total pairs', 'Control group', 'Control pairs', 'Fraction', 'Percentage']
    display_rows = []
    for row in summary_rows:
        fraction = row['control_associated_pair_fraction']
        percentage = row['control_associated_percentage']
        display_rows.append([
            f"{row['qvalue_threshold']:.10g}",
            str(row['total_significant_pairs']),
            row['control_group'],
            str(row['control_associated_pair_count']),
            'NA' if fraction is None else f'{fraction:.6g}',
            'NA' if percentage is None else f'{percentage:.4f}%',
        ])

    widths = [
        max(len(headers[index]), *(len(row[index]) for row in display_rows))
        for index in range(len(headers))
    ]
    print('  '.join(header.ljust(widths[index]) for index, header in enumerate(headers)))
    print('  '.join('-' * width for width in widths))
    for row in display_rows:
        print('  '.join(value.ljust(widths[index]) for index, value in enumerate(row)))


def _write_summary(summary_rows, save_path):
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
        writer = csv.DictWriter(temp_handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(summary_rows)
        temp_handle.flush()
        os.fsync(temp_handle.fileno())
        temp_handle.close()
        os.replace(temp_path, save_path)
    except BaseException:
        temp_handle.close()
        if temp_path.exists():
            temp_path.unlink()
        raise


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Estimate empirical false-positive baselines from control-associated pairs '
            'in a raw SCRIN result.'
        )
    )
    parser.add_argument('--result_path', required=True, help='Path to a raw SCRIN result CSV.')
    parser.add_argument(
        '--control_group',
        action='append',
        nargs=3,
        required=True,
        metavar=('GROUP_NAME', 'MATCH_MODE', 'PATTERN'),
        help='Repeatable control rule. MATCH_MODE: exact, prefix, suffix, or contains.',
    )
    parser.add_argument(
        '--qvalue_thresholds',
        nargs='+',
        type=_qvalue_threshold,
        default=DEFAULT_QVALUE_THRESHOLDS,
        help='One or more q-value cutoffs. Default: 0.05 0.01 0.001 1e-5.',
    )
    parser.add_argument(
        '--pair_keep', choices=['first', 'last'], default='last',
        help='Bidirectional pair deduplication behavior. Default: last.',
    )
    parser.add_argument('--ignore_case', action='store_true', help='Match control names case-insensitively.')
    parser.add_argument('--save_path', default=None, help='Optional path for the summary CSV.')
    parser.add_argument('--work_dir', default=None, help='Optional parent directory for temporary pair buckets.')
    parser.add_argument(
        '--read_chunk_size', type=_positive_int, default=DEFAULT_READ_CHUNK_ROWS,
        help='Rows read from the raw result at a time. Default: 100000.',
    )
    args = parser.parse_args()

    thresholds = []
    for threshold in args.qvalue_thresholds:
        if threshold not in thresholds:
            thresholds.append(threshold)
    try:
        control_groups = _parse_control_groups(args.control_group)
    except ValueError as error:
        parser.error(str(error))

    result_path = Path(args.result_path).resolve()
    if args.save_path is not None and Path(args.save_path).resolve() == result_path:
        parser.error('--save_path must be different from --result_path.')

    temporary_parent = None
    if args.work_dir is not None:
        temporary_parent = Path(args.work_dir)
        temporary_parent.mkdir(parents=True, exist_ok=True)

    bucket_count = _bucket_count(result_path)
    with tempfile.TemporaryDirectory(
        prefix='scrin_false_positive_', dir=temporary_parent
    ) as temporary_dir:
        bucket_dir = Path(temporary_dir) / 'pair_buckets'
        matched_features = _partition_results(
            result_path,
            bucket_dir,
            bucket_count,
            control_groups,
            args.ignore_case,
            args.read_chunk_size,
        )
        summary_rows = _analyze_buckets(
            bucket_dir,
            bucket_count,
            control_groups,
            args.ignore_case,
            thresholds,
            args.pair_keep,
            matched_features,
        )

    _print_summary(summary_rows, matched_features, args.pair_keep)
    if args.save_path is not None:
        _write_summary(summary_rows, args.save_path)
        print(f'\nSummary saved to {args.save_path}')


if __name__ == '__main__':
    main()
