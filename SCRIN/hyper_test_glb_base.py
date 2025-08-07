import argparse
import os
import csv
import pickle
import time
import pandas as pd
import numpy as np
from mpi4py import MPI
from mpi4py.util import pkl5
from tqdm import tqdm
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from functools import partial
from collections import defaultdict
from SCRIN.tools.result_proc import add_pair_column, test_result_df_filter, test_result_df_ratio_proc


def large_bcast(data, comm, rank, size, root=0):
    if rank == root:
        for i in range(size):
            if i != root:
                comm.send(data, dest=i, tag=77)
    else:
        data = comm.recv(source=root, tag=77)
    return data


def robust_send(data, dest, tag, comm, rank, max_retries=5, retry_interval=1):
    retries = 0
    while retries < max_retries:
        try:
            comm.send(data, dest=dest, tag=tag)
            return
        except Exception as e:
            print(f"Rank {rank}: Send failed with exception: {e}. Retrying...")
            time.sleep(retry_interval)
            retries += 1
    print("Failed to send data after retries. Raising exception.")
    raise


def merge_dicts(list_of_dicts):
    merged_list = []  # creating a list to store the merged results

    merged_dict = defaultdict(lambda: [0, defaultdict(int)])  # Create a default dictionary to temporarily store the merge results

    for gene_data in list_of_dicts:
        for gene_id, (around_num, sub_dict) in gene_data.items():
            merged_dict[gene_id][0] += around_num  # Merge around_num
            for sub_gene_id, sub_num in sub_dict.items():
                merged_dict[gene_id][1][sub_gene_id] += sub_num  # Merge values from sub-dictionaries

    # Convert the temporarily stored result into a list
    for gene_id, (around_num, sub_dict) in merged_dict.items():
        merged_list.append((gene_id, around_num, dict(sub_dict)))  # Add the results to a list as a tuple

    return merged_list


def distribute_tasks_dynamic(comm, rank, size, tasks, task_processor, task_tag):
    """
    Dynamically assign tasks to processes and collect results.
    :param comm: MPI communicator.
    :param rank: rank of the process.
    :param size: total number of processes.
    :param tasks: list of tasks, each task is a processable unit.
    :param task_processor: a function that receives a task as input and returns the processing result.
    :param task_tag: tag of the task.
    :return: If rank 0, return a list of processing results of all tasks; otherwise return None.
    """

    task_results = []

    if rank == 0:
        task_index = 0
        active_workers = 0
        task_completed = 0

        # send initial tasks to workers
        for i in range(1, size):
            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=i, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=i, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1

        # collect results until all tasks are processed, sending new tasks as workers become available
        while active_workers > 0:
            status = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            task_results.append(result)
            active_workers -= 1
            task_completed += 1

            print(f"{task_tag} processing: Task {task_completed} completed, {len(tasks) - task_completed} remaining.")

            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=status.source, tag=0, comm=comm, rank=rank, max_retries=10,
                            retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=status.source, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        return task_results

    else:
        while True:
            task = comm.recv(source=0, tag=MPI.ANY_TAG)
            if task is None:
                break
            result = task_processor(task)
            robust_send(result, dest=0, tag=1, comm=comm, rank=rank, max_retries=10, retry_interval=60)


def region_function_cell_radius(df_cell, r_check=0.5):
    """
    returns a dictionary with geneID as keys and a list as values.
    The list contains two elements:
    1. The number of points around the gene within the specified radius.
    2. A dictionary where keys are geneIDs and values are the number of points around the gene within the specified radius.
    {geneID: [around_num, {geneID1: num1, geneID2: num2, ...}]}
    """

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
            distances_xy = np.sqrt((x - df_cell_z['x']) ** 2 + (y - df_cell_z['y']) ** 2)
            df_cell_gene_around = df_cell_z[distances_xy <= r_check]
            gene_around_df_list.append(df_cell_gene_around)

        gene_around_df = pd.concat(gene_around_df_list)
        gene_around_df = gene_around_df.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])
        gene_around_num = len(gene_around_df)
        dict_out[gene].append(gene_around_num)

        gene_cell_num_dict = gene_around_df.groupby('geneID').size().to_dict()
        if gene in gene_cell_num_dict:
            del gene_cell_num_dict[gene]
        dict_out[gene].append(gene_cell_num_dict)

    return dict_out


def region_function_cell_nine_grid(df_cell, r_check=1):
    """
    returns a dictionary with geneID as keys and a list as values.
    The list contains two elements:
    1. The number of points around the gene within the specified radius.
    2. A dictionary where keys are geneIDs and values are the number of points around the gene within the specified radius.
    {geneID: [around_num, {geneID1: num1, geneID2: num2, ...}]}
    """

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
                (abs(x - df_cell['x']) <= r_check) & (abs(y - df_cell['y']) <= r_check)]
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


def test_function(gene_B_around_num, gene_A_N, gene_B_N,
                  gene_slice_num, gene_around_num, minGenenumber, expressionLevel):

    if (gene_around_num < minGenenumber or
            gene_A_N / gene_B_N > expressionLevel or gene_B_N / gene_A_N > expressionLevel):
        return None

    a = hypergeom.sf(k=gene_B_around_num - 1, M=gene_slice_num, n=gene_B_N,
                     N=gene_around_num)
    para_list = [gene_B_around_num, gene_B_N, gene_around_num, gene_slice_num, gene_A_N]

    return [a, para_list]


def hyper_test(tuple_proc, gene_num_dict, pixel_num_slice_all, min_gene_number, expression_level):
    gene_id = tuple_proc[0]
    pixel_num_around = tuple_proc[1]
    gene_around_dict = tuple_proc[2]

    gene_id_around_list, p_list, para_list_all = [], [], []

    for gene_id_around, around_num in gene_around_dict.items():
        tf_list = test_function(gene_B_around_num=around_num,
                                gene_A_N=gene_num_dict[gene_id],
                                gene_B_N=gene_num_dict[gene_id_around],
                                gene_slice_num=pixel_num_slice_all,
                                gene_around_num=pixel_num_around,
                                minGenenumber=min_gene_number,
                                expressionLevel=expression_level)
        if tf_list is not None:
            gene_id_around_list.append(gene_id_around)
            p_list.append(tf_list[0])
            para_list_all.append(tf_list[1])

    if len(p_list) == 0:
        return None

    gene_B_around_num_list, gene_B_N_list, gene_around_num_list, gene_slice_num_list, gene_A_N_list = zip(*para_list_all)

    p_bh = multipletests(p_list, method='fdr_bh')[1]
    p_bo = multipletests(p_list, method='bonferroni')[1]

    results_df = pd.DataFrame({
        'gene_A': [gene_id] * len(p_list),
        'gene_B': gene_id_around_list,
        'pvalue': p_list,
        'qvalue_BH': p_bh,
        'qvalue_BO': p_bo,
        'gene_B_around': gene_B_around_num_list,
        'gene_B_slice': gene_B_N_list,
        'gene_around': gene_around_num_list,
        'gene_slice': gene_slice_num_list,
        'gene_A_N': gene_A_N_list,
        'gene_B_N': gene_B_N_list
    })

    return results_df


def hyper_test_glb_base(opt):
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--detection_method", type=str, choices=['Radius', 'Nine_grid'],
    #                     default='Radius', help="Method for neighbor detection, can be 'Radius' or 'Nine_grid'")
    # parser.add_argument("--r_check", type=float, default=None, help="radius of checking")
    # parser.add_argument("--grid_check", type=int, default=None,
    #                     help="grid size for Nine_grid detection method, default is 1")
    # parser.add_argument("--column_name", type=str, default="x,y,z,geneID,cell",
    #                     help="column name used in data")
    # parser.add_argument("--min_gene_number", type=int, default=5,
    #                     help="minimum number of transcripts for a gene to be considered")
    # parser.add_argument("--min_neighbor_number", type=int, default=1,
    #                     help="minimum number of neighbors for a pair to be considered")
    # parser.add_argument("--expression_level", type=float, default=100,
    #                     help="For gene A and gene B in the pair, the maximum ratio of their expression count.")
    # parser.add_argument("--filter_threshold", type=float, default=0.00001,
    #                     help="filter threshold for qvalue_BH in post processing")
    # parser.add_argument("--pair_keep", type=str, default='last',
    #                     help="keep method for pair post processing, can be 'first' or 'last'")
    # parser.add_argument("--data_path", type=str,
    #                     default="/data2/yangxr009/ST_STA_stereo/data_ori/E16.5_E1S3_WholeBrain_GEM_CellBin.tsv.gz",
    #                     help="path of data")
    # parser.add_argument("--save_path", type=str,
    #                     default="/data2/yangxr009/ST_STA_stereo/hyperTest_re/E16.5_E1S3_WholeBrain_GEM_CellBin_hyper_test.csv",
    #                     help="path of result saving")
    # parser.add_argument("--num_nodes", type=int, default=6, help="Number of nodes")
    # parser.add_argument("--cores_per_node", type=int, default=16, help="Number of cores per node")
    # opt = parser.parse_args()

    # comm = pkl5.Intracomm(MPI.COMM_WORLD)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        # Print the options
        print("\n--- Options ---")
        for arg_name, arg_value in vars(opt).items():
            print(f"{arg_name}: {arg_value}")
        print("--------------------")

        # Check the options of detection method
        if opt.detection_method == 'Radius':
            if opt.r_check is None:
                raise ValueError("Detection method 'Radius' requires --r_check to be set.")

        if opt.detection_method == 'Nine_grid':
            if opt.grid_check is None:
                raise ValueError("Detection method 'Nine_grid' requires --grid_check to be set.")

        column_names = opt.column_name.split(',')

        if len(column_names) == 5:
            df_flt = pd.read_csv(opt.data_path, sep=',', header=0, usecols=column_names,
                                 dtype={column_names[0]: float, column_names[1]: float,
                                        column_names[2]: float, column_names[3]: str, column_names[4]: str},
                                 na_values='', keep_default_na=False)
        elif len(column_names) == 4:
            df_flt = pd.read_csv(opt.data_path, sep=',', header=0, usecols=column_names,
                                 dtype={column_names[0]: float, column_names[1]: float,
                                        column_names[2]: str, column_names[3]: str},
                                 na_values='', keep_default_na=False)
        else:
            raise ValueError("Invalid number of columns specified in --column_name")

        if df_flt.isnull().any().any():
            raise ValueError("Data contains null values. Process aborted.")

        # rename columns based on the provided column names
        if len(column_names) == 5:
            df_flt.rename(columns={column_names[0]: 'x',
                                   column_names[1]: 'y',
                                   column_names[2]: 'z',
                                   column_names[3]: 'geneID',
                                   column_names[4]: 'cell'}, inplace=True)
            print(f"The column {column_names[0]} is used as x, {column_names[1]} as y, "
                  f"{column_names[2]} as z, {column_names[3]} as geneID, {column_names[4]} as cell.")
        else:
            df_flt.rename(columns={column_names[0]: 'x',
                                   column_names[1]: 'y',
                                   column_names[2]: 'geneID',
                                   column_names[3]: 'cell'}, inplace=True)
            print(f"The column {column_names[0]} is used as x, {column_names[1]} as y, "
                  f"{column_names[2]} as geneID, {column_names[3]} as cell.")
            # add a dummy z column
            df_flt['z'] = 0.0

        # Calculate the number of transcripts for each geneID
        gene_counts = df_flt.groupby('geneID').size().reset_index(name='transcript_number')

        # Filter genes with more than 5 transcripts
        gene_flt = gene_counts[gene_counts['transcript_number'] >= opt.min_gene_number]
        print("Gene with less than 5 transcripts are filtered out.")
        print(f"Total number of genes after filtering: {len(gene_flt)}")

        # Filter the main DataFrame based on the filtered genes
        df_flt_region = df_flt[df_flt['geneID'].isin(gene_flt['geneID'])]
        df_flt = None

        # If nine_grid detection method is used, drop duplicate rows
        if opt.detection_method == 'Nine_grid':
            df_flt_region = df_flt_region.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])

        # Calculate the number of unique genes in DataFrame
        gene_num_dict = df_flt_region.groupby('geneID').size().to_dict()
        print("len of gene_num_dict: ", len(gene_num_dict))

        # Calculate the total number of pixels in the slice
        if opt.detection_method == 'Nine_grid':
            # For nine_grid, pixel_num_slice_all is the number of unique (x, y, z) coordinates
            pixel_num_slice_all = len(df_flt_region[['x', 'y', 'z']].drop_duplicates())
        else:
            # For radius detection, pixel_num_slice_all is the total number of rows in the DataFrame
            pixel_num_slice_all = len(df_flt_region)

        # Group the DataFrame by 'cell' to prepare for parallel processing
        cell_group = df_flt_region.groupby('cell')
        cell_group_list = [cell_group.get_group(cell_id) for cell_id in cell_group.groups]
        print("cell num: ", len(cell_group_list))

        df_flt_region = None
        cell_group = None

    else:
        gene_num_dict = None
        pixel_num_slice_all = None
        cell_group_list = None

    gene_num_dict = comm.bcast(gene_num_dict, root=0)
    pixel_num_slice_all = comm.bcast(pixel_num_slice_all, root=0)

    min_gene_number_local = opt.min_neighbor_number
    expression_level_local = opt.expression_level

    if opt.detection_method == 'Nine_grid':
        region_function_cell = region_function_cell_nine_grid
        print("Using Nine_grid detection method.")
    else:
        region_function_cell = region_function_cell_radius
        print("Using Radius detection method.")

    # Task 1: Detect neighbors in each cell
    partial_func = partial(region_function_cell, r_check=opt.r_check)

    if rank == 0:
        dict_around_list = distribute_tasks_dynamic(comm, rank, size, cell_group_list, partial_func, 'NeighborDetection')
    else:
        distribute_tasks_dynamic(comm, rank, size, cell_group_list, partial_func, 'NeighborDetection')
        dict_around_list = None

    comm.Barrier()

    if rank == 0:
        tuple_around = merge_dicts(dict_around_list)
        print("len of tuple_around: ", len(tuple_around))
        dict_around_list = None
    else:
        tuple_around = None

    # Task 2: Perform hypergeometric test on the detected neighbors
    partial_func2 = partial(hyper_test, gene_num_dict=gene_num_dict,
                            pixel_num_slice_all=pixel_num_slice_all, min_gene_number=min_gene_number_local,
                            expression_level=expression_level_local)

    if rank == 0:
        result_output = distribute_tasks_dynamic(comm, rank, size, tuple_around, partial_func2, 'HyperTest')
    else:
        distribute_tasks_dynamic(comm, rank, size, tuple_around, partial_func2, 'HyperTest')
        result_output = None

    comm.Barrier()

    if rank == 0:
        df_result = pd.concat(result_output)
        df_result = add_pair_column(df_result)
        df_result = test_result_df_ratio_proc(df_result)
        df_result.to_csv(opt.save_path, index=False)

        df_flt, df_end = test_result_df_filter(df_result, filter_threshold=opt.filter_threshold, keep=opt.pair_keep)
        save_path = opt.save_path.replace('.csv', f'_dedup_{opt.filter_threshold}_post_proc.csv')
        df_end.to_csv(save_path, index=False)

