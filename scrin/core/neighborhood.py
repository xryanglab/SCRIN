import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import itertools
import os
import joblib
from collections import defaultdict


def _core_neighbor_search(df_source, query_coords, r_check, dims=['x', 'y']):
    """
    Core function: Use KDTree to search for neighbor indices.
    dims: Specify the dimensions involved in distance calculation. ['x', 'y'] for 2D, ['x', 'y', 'z'] for 3D.
    """
    if df_source.empty or len(query_coords) == 0:
        return [[] for _ in range(len(query_coords))]

    # Build spatial index
    tree = cKDTree(df_source[dims].values)
    # Query points within the radius in batches
    return tree.query_ball_point(query_coords, r=r_check)


def region_function_discrete_z(df_cell, r_check=0.5):
    """
    Function for discrete Z-axis: Calculates neighbors in 2D (x, y) for each Z layer.

    Parameters:
    - df_cell: DataFrame containing spatial and gene information.
    - r_check: Radius for neighbor search in 2D space.

    Returns:
    - dict_out: A dictionary where keys are gene IDs and values are lists containing:
        1. The number of neighboring points within the specified radius.
        2. A dictionary of neighboring gene counts (excluding the current gene).
    """
    df_cell = df_cell.reset_index(drop=True)

    dict_out = {}
    gene_list = df_cell['geneID'].unique()

    # Group data by Z-axis
    z_groups = {z: group for z, group in df_cell.groupby('z')}

    for gene in gene_list:
        df_gene = df_cell[df_cell['geneID'] == gene]
        all_neighbor_global_indices = set()

        for z, sub_group in df_gene.groupby('z'):
            if z not in z_groups:
                continue

            target_layer = z_groups[z]

            # Search for neighbors in the current Z layer
            indices = _core_neighbor_search(
                target_layer,
                sub_group[['x', 'y']].values,
                r_check,
                dims=['x', 'y']
            )

            # Collect global indices of neighbors
            for idx_list in indices:
                all_neighbor_global_indices.update(target_layer.iloc[idx_list].index)

        # Extract neighboring points from the DataFrame
        neighbor_df = df_cell.loc[list(all_neighbor_global_indices)]

        # Count the number of neighboring points
        around_num = len(neighbor_df)

        # Count occurrences of other genes in the neighboring points
        gene_counts = (
            neighbor_df[neighbor_df['geneID'] != gene]
            .groupby('geneID')
            .size()
            .to_dict()
        )

        dict_out[gene] = [around_num, gene_counts]

    return dict_out


def region_function_continuous_z(df_cell, r_check=0.5):
    """
    Function for continuous Z-axis: Calculates 3D spatial distances in XYZ space.
    """
    df_cell = df_cell.reset_index(drop=True)

    dict_out = {}
    gene_list = df_cell['geneID'].unique()

    # Build a global 3D spatial index
    tree = cKDTree(df_cell[['x', 'y', 'z']].values)

    for gene in gene_list:
        df_gene = df_cell[df_cell['geneID'] == gene]

        # Retrieve all neighbors for the gene in 3D space
        indices_list = tree.query_ball_point(df_gene[['x', 'y', 'z']].values, r=r_check)

        # Flatten and deduplicate neighbor indices
        flat_indices = list(set([i for sublist in indices_list for i in sublist]))

        # Extract neighboring gene information
        neighbor_genes = df_cell.iloc[flat_indices]['geneID']

        around_num = len(neighbor_genes)
        gene_counts = neighbor_genes[neighbor_genes != gene].value_counts().to_dict()

        dict_out[gene] = [around_num, gene_counts]

    return dict_out


def region_function_cell_nine_grid(df_cell, r_check=1):
    """
    returns a dictionary with geneID as keys and a list as values.
    The list contains two elements:
    1. The number of points around the gene within the specified radius.
    2. A dictionary where keys are geneIDs and values are the number of points around the gene within the specified radius.
    {geneID: [around_num, {geneID1: num1, geneID2: num2, ...}]}
    """
    df_cell = df_cell.reset_index(drop=True)

    dict_out = {}
    gene_list = df_cell['geneID'].unique()

    # pre z dict
    df_cell_z_group = df_cell.groupby('z')
    df_cell_z_dict = {}
    for z, group in df_cell_z_group:
        df_cell_z_dict[z] = group

    for gene in gene_list:
        dict_out[gene] = []
        df_cell_gene = df_cell[df_cell['geneID'] == gene]

        gene_around_df_list = []
        for x, y, z in zip(df_cell_gene['x'], df_cell_gene['y'], df_cell_gene['z']):
            # df_cell_z = df_cell[df_cell['z'] == z]
            df_cell_z = df_cell_z_dict[z]
            df_cell_gene_around = df_cell_z[
                (abs(x - df_cell_z['x']) <= r_check) & (abs(y - df_cell_z['y']) <= r_check)]
            gene_around_df_list.append(df_cell_gene_around)

        gene_around_df = pd.concat(gene_around_df_list)
        gene_around_num = len(gene_around_df.drop_duplicates(subset=['x', 'y', 'z']))
        dict_out[gene].append(gene_around_num)

        gene_around_df_flt = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
        gene_cell_num_dict = gene_around_df_flt.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        dict_out[gene].append(gene_cell_num_dict)

    return dict_out


# For glb_distribution; TODO: optimize with KDTree, support continuous z
def region_function_cell_glb_distribution(df_cell, r_check=0.5, r_dist=0.5):
    """
    For a cell, return statistics for all gene types:
    - total number of surrounding points
    - number of surrounding points for each gene type
    Output format: {geneID: [around_num, {geneID1: num1, geneID2: num2, ...}]}
    """
    df_cell = df_cell.reset_index(drop=True)

    dict_out, dict_dist = {}, {}
    gene_list = df_cell['geneID'].unique()
    for gene in gene_list:
        dict_out[gene] = []
        df_cell_gene = df_cell[df_cell['geneID'] == gene]

        gene_around_df_list, gene_around_dist_df_list = [], []
        for x, y, z in zip(df_cell_gene['x'], df_cell_gene['y'], df_cell_gene['z']):
            df_cell_z = df_cell[df_cell['z'] == z]
            distances_xy = np.sqrt((x - df_cell_z['x']) ** 2 + (y - df_cell_z['y']) ** 2)
            df_cell_z_copy = df_cell_z.copy()
            df_cell_z_copy['distances_xy'] = distances_xy

            df_cell_gene_around = df_cell_z_copy[distances_xy <= r_check]
            gene_around_df_list.append(df_cell_gene_around)
            df_cell_gene_around_dist = df_cell_z_copy[distances_xy <= r_dist]
            gene_around_dist_df_list.append(df_cell_gene_around_dist)

        gene_around_df = pd.concat(gene_around_df_list)
        gene_around_df = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
        gene_around_num = len(gene_around_df)
        dict_out[gene].append(gene_around_num)

        gene_around_dist_df = pd.concat(gene_around_dist_df_list)

        gene_cell_num_dict = gene_around_df.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        dict_out[gene].append(gene_cell_num_dict)

        gene_cell_dist_dict = {}
        for gene_id, group in gene_around_dist_df.groupby('geneID'):
            gene_cell_dist_dict[gene_id] = group['distances_xy'].tolist()
        if gene in gene_cell_dist_dict:
            del gene_cell_dist_dict[gene]
        dict_dist[gene] = gene_cell_dist_dict

    return [dict_out, dict_dist]


# For clb_distribution
def region_function_cell_multi_proc(df_cell, r_check=0.5, gene_rank_id2rank=None, z_continuous=False):
    """
    Return local updates from a cell dataframe.
    Local update 1: gene B around
    Local update 2: gene around, gene B, gene all, gene A
    """
    df_cell = df_cell.reset_index(drop=True)

    local_updates1 = defaultdict(list)  # Store updates of gene B around
    local_updates2 = defaultdict(list)  # Store updates of gene around, gene B, gene all, gene A

    gene_list = df_cell['geneID'].unique()

    df_cell_z_dict = {}
    kdtree_dict = {}  # Store KDTree for each z layer
    kdtree_3d = None

    if not z_continuous:
        df_cell_z_group = df_cell.groupby('z')
        for z, group in df_cell_z_group:
            df_cell_z_dict[z] = group
            # Build KDTree from x,y numpy array for faster neighbor queries
            kdtree_dict[z] = cKDTree(group[['x', 'y']].values)
    else:
        # For continuous z, we build a single KDTree for the entire dataset in 3D space
        kdtree_3d = cKDTree(df_cell[['x', 'y', 'z']].values)

    gene_cell_num = len(df_cell)
    gene_cell_num_all_dict = df_cell.groupby('geneID').size().to_dict()
    all_gene_keys = list(gene_cell_num_all_dict.keys())

    for gene in gene_list:
        df_cell_gene = df_cell[df_cell['geneID'] == gene]

        gene_around_df_list = []

        if not z_continuous:
            # Process by z layer instead of iterating each transcript
            for z in df_cell_gene['z'].unique():
                gene_points_in_z = df_cell_gene[df_cell_gene['z'] == z][['x', 'y']].values
                if len(gene_points_in_z) == 0:
                    continue

                # Batch query r_check neighborhood for all target points on this z layer
                # Returns nested list of neighbor row indices in the original group
                idx_list_of_lists = kdtree_dict[z].query_ball_point(gene_points_in_z, r=r_check)

                # Flatten nested list with itertools and deduplicate using set
                # This operates on a small set of integers in memory, avoiding expensive DataFrame copies
                unique_indices = set(itertools.chain.from_iterable(idx_list_of_lists))

                # Use deduplicated absolute indices to slice the DataFrame for this z layer at once
                if unique_indices:
                    gene_around_df_list.append(df_cell_z_dict[z].iloc[list(unique_indices)])
        else:
            # For continuous z, we can directly query the 3D KDTree with all points of the gene at once
            gene_points_3d = df_cell_gene[['x', 'y', 'z']].values
            if len(gene_points_3d) > 0:
                # Batch query r_check neighborhood for all target points in 3D space
                idx_list_of_lists = kdtree_3d.query_ball_point(gene_points_3d, r=r_check)
                unique_indices = set(itertools.chain.from_iterable(idx_list_of_lists))
                if unique_indices:
                    # Since the KDTree was built on the entire df_cell, we can directly slice df_cell with the unique indices
                    gene_around_df_list.append(df_cell.iloc[list(unique_indices)])

        if len(gene_around_df_list) > 0:
            gene_around_df = pd.concat(gene_around_df_list)
            # Since deduplication in each z layer has been applied, merged result is nearly clean
            # Keep a final deduplication as a safe fallback (minimal overhead)
            gene_around_df = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
            gene_around_num = len(gene_around_df)
        else:
            gene_around_num = 0
            gene_around_df = pd.DataFrame(columns=df_cell.columns)

        gene_cell_num_dict = gene_around_df.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        around_gene_keys = list(gene_cell_num_dict.keys())

        gene_rank_belong = gene_rank_id2rank[gene]  # Get the rank to which the current gene belongs

        for gene_around in around_gene_keys:
            gene_B_around_num = gene_cell_num_dict[gene_around]
            seq_B_around = [gene, gene_around, gene_B_around_num]  # gene_B_around
            local_updates1[gene_rank_belong].append(seq_B_around)

        for g in all_gene_keys:
            if g == gene:
                continue

            #  Last 4 values are: gene_around, gene_B, gene_all, gene_A
            seq_gene_around = [gene, g, gene_around_num,
                               gene_cell_num_all_dict[g], gene_cell_num, gene_cell_num_all_dict[gene]]
            local_updates2[gene_rank_belong].append(seq_gene_around)

    return [local_updates1, local_updates2]


def region_function_cell_multi_proc_nine_grid(df_cell, r_check=0.5, gene_rank_id2rank=None):
    """
    Return local updates from a cell dataframe.
    Local update 1: gene B around
    Local update 2: gene around, gene B, gene all, gene A
    """
    df_cell = df_cell.reset_index(drop=True)

    local_updates1 = defaultdict(list)  # Store updates of gene B around
    local_updates2 = defaultdict(list)  # Store updates of gene around, gene B, gene all, gene A

    gene_list = df_cell['geneID'].unique()

    # Precompute dictionary by z
    df_cell_z_group = df_cell.groupby('z')
    df_cell_z_dict = {}
    for z, group in df_cell_z_group:
        df_cell_z_dict[z] = group

    # gene_cell_num = len(df_cell)
    gene_cell_num = len(df_cell[['x', 'y', 'z']].drop_duplicates())
    df_cell_pixel = df_cell.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
    gene_cell_num_all_dict = df_cell_pixel.groupby('geneID').size().to_dict()
    all_gene_keys = list(gene_cell_num_all_dict.keys())
    for gene in gene_list:
        df_cell_gene = df_cell[df_cell['geneID'] == gene]

        gene_around_df_list = []
        for x, y, z in zip(df_cell_gene['x'], df_cell_gene['y'], df_cell_gene['z']):
            # df_cell_z = df_cell[df_cell['z'] == z]
            df_cell_z = df_cell_z_dict[z]
            df_cell_gene_around = df_cell_z[
                (abs(x - df_cell_z['x']) <= r_check) & (abs(y - df_cell_z['y']) <= r_check)]
            gene_around_df_list.append(df_cell_gene_around)

        gene_around_df = pd.concat(gene_around_df_list)
        gene_around_num = len(gene_around_df.drop_duplicates(subset=['x', 'y', 'z']))

        gene_around_df_flt = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
        gene_cell_num_dict = gene_around_df_flt.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        around_gene_keys = list(gene_cell_num_dict.keys())

        gene_rank_belong = gene_rank_id2rank[gene]  # Get the rank to which the current gene belongs

        for gene_around in around_gene_keys:
            gene_B_around_num = gene_cell_num_dict[gene_around]
            seq_B_around = [gene, gene_around, gene_B_around_num]  # gene_B_around
            local_updates1[gene_rank_belong].append(seq_B_around)

        for g in all_gene_keys:
            if g == gene:
                continue

            #  Last 4 values are: gene_around, gene_B, gene_all, gene_A
            seq_gene_around = [gene, g, gene_around_num,
                               gene_cell_num_all_dict[g], gene_cell_num, gene_cell_num_all_dict[gene]]
            local_updates2[gene_rank_belong].append(seq_gene_around)

    return [local_updates1, local_updates2]


def region_function_cell_multi_proc_distribution(df_cell, r_check=0.5, r_dist=0.5, distribution_save_interval=10, local_distribution_temp=None, local_save_path=None, gene_rank_id2rank=None):
    """
    Return local updates from a cell dataframe, store gene distribution information.
    Local update 1: gene B around
    Local update 2: gene around, gene B, gene all, gene A
    """
    df_cell = df_cell.reset_index(drop=True)

    local_updates1 = defaultdict(list)  # Store updates of gene B around
    local_updates2 = defaultdict(list)  # Store updates of gene around, gene B, gene all, gene A
    distribution_dict = {}  # Store gene distribution information

    gene_list = df_cell['geneID'].unique()

    # Precompute dictionary by z
    df_cell_z_group = df_cell.groupby('z')
    df_cell_z_dict = {}
    for z, group in df_cell_z_group:
        df_cell_z_dict[z] = group

    gene_cell_num = len(df_cell)
    gene_cell_num_all_dict = df_cell.groupby('geneID').size().to_dict()
    all_gene_keys = list(gene_cell_num_all_dict.keys())
    for gene in gene_list:
        df_cell_gene = df_cell[df_cell['geneID'] == gene]

        gene_around_df_list, gene_distribution_df_list = [], []
        for x, y, z in zip(df_cell_gene['x'], df_cell_gene['y'], df_cell_gene['z']):
            # df_cell_z = df_cell[df_cell['z'] == z]
            df_cell_z = df_cell_z_dict[z]
            distances_xy = np.sqrt((x - df_cell_z['x']) ** 2 + (y - df_cell_z['y']) ** 2)

            df_cell_z_copy = df_cell_z.copy()
            df_cell_z_copy['distances_xy'] = distances_xy

            df_cell_gene_around = df_cell_z_copy[df_cell_z_copy['distances_xy'] <= r_check]
            df_cell_gene_distribution = df_cell_z_copy[df_cell_z_copy['distances_xy'] <= r_dist]

            gene_around_df_list.append(df_cell_gene_around)
            gene_distribution_df_list.append(df_cell_gene_distribution)

        gene_around_df = pd.concat(gene_around_df_list)
        gene_distribution_df = pd.concat(gene_distribution_df_list)

        gene_around_df = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
        gene_around_num = len(gene_around_df)

        gene_cell_num_dict = gene_around_df.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        around_gene_keys = list(gene_cell_num_dict.keys())

        gene_rank_belong = gene_rank_id2rank[gene]  # Get the rank to which the current gene belongs

        for gene_around in around_gene_keys:
            gene_B_around_num = gene_cell_num_dict[gene_around]
            seq_B_around = [gene, gene_around, gene_B_around_num]  # gene_B_around
            local_updates1[gene_rank_belong].append(seq_B_around)

        for g in all_gene_keys:
            if g == gene:
                continue

            #  Last 4 values are: gene_around, gene_B, gene_all, gene_A
            seq_gene_around = [gene, g, gene_around_num,
                               gene_cell_num_all_dict[g], gene_cell_num, gene_cell_num_all_dict[gene]]
            local_updates2[gene_rank_belong].append(seq_gene_around)

        # process gene distribution
        gene_cell_distribution_dict = {}
        for gene_id, group in gene_distribution_df.groupby('geneID'):
            if gene_id == gene:
                continue

            gene_cell_distribution_dict[gene_id] = group['distances_xy'].tolist()

        if len(gene_cell_distribution_dict) > 0:
            distribution_dict[gene] = gene_cell_distribution_dict

    # process temp distribution list
    if len(local_distribution_temp) >= distribution_save_interval:
        local_save_path_file_count = len(os.listdir(local_save_path))
        distribution_save_path = os.path.join(local_save_path, f"distribution_{local_save_path_file_count}.pkl")
        joblib.dump(local_distribution_temp, distribution_save_path, compress=3)
        # clear the temp list after saving
        local_distribution_temp.clear()

    if len(distribution_dict) > 0:
        local_distribution_temp.append(distribution_dict)  # structrue: [{gene: {gene_id: [distances_xy, ...]}, ...}, ...], 1 dict to 1 cell

    return [local_updates1, local_updates2]


# For glb_nocell; TODO: optimize with KDTree, support continuous z
def region_function_without_cell(point_comb_list, r_check=1, detection_method='radius'):
    """
    For a cell-free dataset, return statistics for all gene types:
    - Total number of surrounding points
    - Number of surrounding points for each gene type
    Output format: {geneID: [around_num, {geneID1: num1, geneID2: num2, ...}]}
    """
    point_list = point_comb_list[0]
    point_list_ex = point_comb_list[1]
    df_cell = pd.DataFrame(point_list, columns=['x', 'y', 'geneID'])
    df_cell = df_cell.reset_index(drop=True)
    df_cell_ex = pd.DataFrame(point_list_ex, columns=['x', 'y', 'geneID'])
    df_cell_ex = df_cell_ex.reset_index(drop=True)

    dict_out = {}
    gene_list = df_cell['geneID'].unique()
    for gene in gene_list:
        dict_out[gene] = []
        df_cell_gene = df_cell[df_cell['geneID'] == gene]

        gene_around_df_list = []
        for x, y in zip(df_cell_gene['x'], df_cell_gene['y']):
            if detection_method == 'radius':
                distances_xy = np.sqrt((x - df_cell_ex['x']) ** 2 +
                                       (y - df_cell_ex['y']) ** 2)
                df_cell_gene_around = df_cell_ex[distances_xy <= r_check]
            else:
                df_cell_gene_around = df_cell_ex[
                    (abs(x - df_cell_ex['x']) <= r_check) & (abs(y - df_cell_ex['y']) <= r_check)]
            gene_around_df_list.append(df_cell_gene_around)

        gene_around_df = pd.concat(gene_around_df_list)
        if detection_method == 'radius':
            gene_around_df = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y'])
            gene_around_num = len(gene_around_df)
        else:
            gene_around_df = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y'])
            gene_around_df_total = gene_around_df.drop_duplicates(subset=['x', 'y'])
            gene_around_num = len(gene_around_df_total)

        dict_out[gene].append(gene_around_num)

        gene_cell_num_dict = gene_around_df.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        dict_out[gene].append(gene_cell_num_dict)

    return dict_out