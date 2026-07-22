import heapq

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, skew, kurtosis


def compute_shape_para(distances, bw=0.05):
    shape_count = len(distances)

    kde = gaussian_kde(distances, bw_method=bw)
    x = np.linspace(min(distances), max(distances), 1000)
    density = kde.evaluate(x)
    mode = x[np.argmax(density)]

    skewness = skew(distances)
    kurt = kurtosis(distances)

    median = np.median(distances)

    q75, q25 = np.percentile(distances, [75, 25])
    iqr = q75 - q25

    return shape_count, mode, skewness, kurt, median, q25, q75, iqr


def _sample_shape_distances(distances, shape_max_distance_count):
    """Return the distances used for shape calculation with reproducible sampling."""
    shape_count = len(distances)
    if shape_max_distance_count is None or shape_count <= shape_max_distance_count:
        return distances

    rng = np.random.default_rng(0)
    sample_indices = rng.choice(shape_count, size=shape_max_distance_count, replace=False)
    return np.asarray(distances)[sample_indices]


def build_distribution_shape_tasks(distribution_dict, around_count_threshold=100,
                                   shape_max_distance_count=None, target_task_count=1, eligible_pairs=None):
    """Flatten gene-level distance dictionaries into approximately balanced pair batches."""
    pair_entries = []
    for gene_id, gene_dist_dict in distribution_dict.items():
        for gene_id_around, distances in gene_dist_dict.items():
            if eligible_pairs is not None and (str(gene_id), str(gene_id_around)) not in eligible_pairs:
                continue

            shape_count = len(distances)
            if shape_count < around_count_threshold:
                continue

            shape_distances = _sample_shape_distances(distances, shape_max_distance_count)
            pair_entries.append((gene_id, gene_id_around, shape_count, shape_distances))

    if not pair_entries:
        return []

    batch_count = min(len(pair_entries), max(1, int(target_task_count)))
    # KDE cost is dominated by the number of retained distance observations. Assign
    # the heaviest pairs first to the currently lightest batch to limit stragglers.
    pair_entries.sort(key=lambda entry: len(entry[3]), reverse=True)
    batch_heap = [(0, batch_index, []) for batch_index in range(batch_count)]
    heapq.heapify(batch_heap)

    for entry in pair_entries:
        batch_load, batch_index, batch = heapq.heappop(batch_heap)
        batch.append(entry)
        heapq.heappush(batch_heap, (batch_load + len(entry[3]), batch_index, batch))

    batches = [batch for _, _, batch in batch_heap if batch]
    batches.sort(key=lambda batch: sum(len(entry[3]) for entry in batch), reverse=True)
    return batches


def distribution_shape_calculate_batch(pair_batch):
    """Calculate distribution-shape descriptors for a balanced batch of gene pairs."""
    records = []
    for gene_id, gene_id_around, shape_count, shape_distances in pair_batch:
        _, mode, skewness, kurt, median, q25, q75, iqr = compute_shape_para(shape_distances)
        records.append({
            'gene_A': gene_id,
            'gene_B': gene_id_around,
            'shape_count': shape_count,
            'mode': mode,
            'skewness': skewness,
            'kurtosis': kurt,
            'median': median,
            'q25': q25,
            'q75': q75,
            'iqr': iqr,
        })

    if not records:
        return None
    return pd.DataFrame.from_records(records)


def distribution_shape_calculate(tuple_proc, around_count_threshold=100, shape_max_distance_count=None):
    gene_id, gene_dist_dict = tuple_proc
    tasks = build_distribution_shape_tasks(
        {gene_id: gene_dist_dict},
        around_count_threshold=around_count_threshold,
        shape_max_distance_count=shape_max_distance_count,
        target_task_count=1,
    )
    if not tasks:
        return None
    return distribution_shape_calculate_batch(tasks[0])


def shape_correct(df_result, opt):
    df_result.loc[:, 'skewness_adjusted'] = 1 / (1 + np.exp(df_result['skewness']))
    df_result.loc[:, 'skewness_relative'] = df_result['skewness'] + 0.566
    df_result.loc[:, 'kurtosis_relative'] = df_result['kurtosis'] + 0.6
    df_result.loc[:, 'median_relative'] = df_result['median'] - (0.707 * opt.r_dist)
    df_result.loc[:, 'q25_relative'] = df_result['q25'] - (0.5 * opt.r_dist)
    df_result.loc[:, 'q75_relative'] = df_result['q75'] - (0.866 * opt.r_dist)
    return df_result


