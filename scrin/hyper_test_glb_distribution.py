import argparse
import os
import csv
import pickle
import time
import shutil
import copy
import pandas as pd
import numpy as np
from mpi4py import MPI
from mpi4py.util import pkl5
from scipy.stats import hypergeom, gaussian_kde, skew, kurtosis
from statsmodels.stats.multitest import multipletests
from functools import partial
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from scrin.core.neighborhood import region_function_cell_glb_distribution
from scrin.tools.result_proc import add_pair_column, test_result_df_filter, test_result_df_ratio_proc
from scrin.tools.stream_results import ChunkedCSVWriter, collect_qvalue_pairs, prepare_hyper_result, stream_postprocess
from scrin.core.parallel import split_list_into_sublists, robust_send
from scrin.core.statistics import test_function, hyper_test
from scrin.core.distribution import build_distribution_shape_tasks, distribution_shape_calculate_batch, shape_correct
from scrin.core.intermediate import merge_dicts


def large_bcast(data, comm, rank, size, root=0):
    # If root process, send data to all other processes
    if rank == root:
        for i in range(size):
            if i != root:
                comm.send(data, dest=i, tag=77)
    else:
        # Non-root processes receive data from root
        data = comm.recv(source=root, tag=77)
    return data


def merge_dists(list_of_dicts, output_mode='list'):
    merged_list = []  # Store merged results

    merged_dict = defaultdict(lambda: defaultdict(list))  # Temporary storage for merging

    for gene_data in list_of_dicts:
        for gene_id, sub_dict in gene_data.items():
            for sub_gene_id, sub_dist_list in sub_dict.items():
                merged_dict[gene_id][sub_gene_id].extend(sub_dist_list)  # Merge lists in sub-dict

    if output_mode != 'list':
        return merged_dict

    # Convert temporary results to list, [(gene_id, {sub_gene_id: [dist1, dist2, ...]})]
    for gene_id, sub_dict in merged_dict.items():
        merged_list.append((gene_id, dict(sub_dict)))

    return merged_list


def update_pkl_file(task_result_sub, gene_id_list_split, opt, task_tag, monitor_list=None):
    task_re_num_sub = [_[0] for _ in task_result_sub]
    task_re_dist_sub = [_[1] for _ in task_result_sub]

    merge_results_num = merge_dicts(task_re_num_sub, output_mode='dict')
    merge_results_dist = merge_dists(task_re_dist_sub, output_mode='dict')

    merge_gene_id_list_num = list(merge_results_num.keys())
    merge_gene_id_list_dist = list(merge_results_dist.keys())

    # Update pkl files for NeighborDetection results
    for i, gene_id_list in enumerate(gene_id_list_split):
        common_gene_id = list(set(gene_id_list) & set(merge_gene_id_list_num))
        if len(common_gene_id) != 0:
            with open(f'{opt.intermediate_dir}/NeighborDetection/task_results{i}.pkl', 'rb') as f:
                gene_dict_pkl = pickle.load(f)
                # Update intersection gene_ids
                for gene_id in common_gene_id:
                    gene_dict_pkl[gene_id][0] += merge_results_num[gene_id][0]
                    for sub_gene_id, sub_num in merge_results_num[gene_id][1].items():
                        gene_dict_pkl[gene_id][1][sub_gene_id] += sub_num

                with open(f'{opt.intermediate_dir}/NeighborDetection/task_results{i}.pkl', 'wb') as f:
                    pickle.dump(gene_dict_pkl, f)
        else:
            pass

    # Update pkl files for distance results
    for i, gene_id_list in enumerate(gene_id_list_split):
        common_gene_id = list(set(gene_id_list) & set(merge_gene_id_list_dist))
        if not common_gene_id:
            continue

        main_path = f'{opt.intermediate_dir}/Distance/task_results{i}.pkl'
        chunk_size = 100 * 1024 * 1024  # 100MB

        if os.path.getsize(main_path) > chunk_size:
            # Handle chunked files
            extra_path = f'{opt.intermediate_dir}/extra_im_file_index_dict.pkl'
            with open(extra_path, 'rb') as f:
                extra_dict = pickle.load(f)

            current_chunk, size = get_current_chunk(i, extra_dict, opt.intermediate_dir)
            if not current_chunk or size > chunk_size:
                new_index = (extra_dict[i][-1] + 1) if extra_dict[i] else 1
                extra_dict[i].append(new_index)
                current_chunk = f'{opt.intermediate_dir}/Distance/task_results{i}_extra{new_index}.pkl'

                with open(extra_path, 'wb') as f:
                    pickle.dump(extra_dict, f)

            update_chunk_file(current_chunk, gene_id_list, common_gene_id, merge_results_dist)
        else:
            update_chunk_file(main_path, gene_id_list, common_gene_id, merge_results_dist)

    if monitor_list is not None:
        monitor_list.pop(0)


def get_current_chunk(i, extra_dict, intermediate_dir):
    """Get current chunk info and maintain index."""
    indices = extra_dict.get(i, [])
    if not indices:
        return None, 0  # No extra chunk

    latest_extra = indices[-1]
    chunk_path = f"{intermediate_dir}/Distance/task_results{i}_extra{latest_extra}.pkl"
    return chunk_path, os.path.getsize(chunk_path)


def update_chunk_file(file_path, gene_id_list, common_gene_id, merge_results_dist):
    """General function to update pkl files."""
    try:
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
    except FileNotFoundError:
        data = {gi: defaultdict(list) for gi in gene_id_list}

    for gene_id in common_gene_id:
        for sub_id, dist_list in merge_results_dist[gene_id].items():
            data[gene_id][sub_id].extend(dist_list)

    with open(file_path, 'wb') as f:
        pickle.dump(data, f)
    return os.path.getsize(file_path)


def distribute_tasks_dynamic(comm, rank, size, tasks, task_processor, opt, task_tag, intermediate_split=100,
                             intermediate_save=False, gene_id_list=None, result_handler=None):
    """
    Dynamically distribute tasks to processes and collect results.
    :param comm: MPI communicator.
    :param rank: Process rank.
    :param size: Total number of processes.
    :param tasks: List of tasks.
    :param task_processor: Function to process a task.
    :param opt: Argument container.
    :param task_tag: Task label.
    :param intermediate_split: Split size for saving intermediate results.
    :param intermediate_save: Whether to save intermediate results.
    :param gene_id_list: List of gene IDs.
    :return: If rank 0, returns all results; else None.
    """

    executor = ThreadPoolExecutor(max_workers=1)  # Thread for updating pkl files

    task_results = []
    monitor_list = []

    if rank == 0:
        task_index = 0
        active_workers = 0
        task_completed = 0

        task_result_proc = 0

        intermediate_size = len(tasks) // intermediate_split
        gene_id_list_split = None

        if gene_id_list is not None and task_tag == 'NeighborDetection':
            if len(gene_id_list) < intermediate_split:
                intermediate_split = len(gene_id_list)
                intermediate_size = len(tasks) // intermediate_split
                gene_id_list_split = split_list_into_sublists(gene_id_list, intermediate_split)
            else:
                gene_id_list_split = split_list_into_sublists(gene_id_list, intermediate_split)

            # Create empty dicts for each split and save as pkl
            extra_im_file_index_dict = defaultdict(list)
            with open(f'{opt.intermediate_dir}/extra_im_file_index_dict.pkl', 'wb') as f:
                pickle.dump(extra_im_file_index_dict, f)

            for i, gene_id_list in enumerate(gene_id_list_split):
                dict_pkl_num, dict_pkl_distance = {}, {}

                for gene_id in gene_id_list:
                    dict_pkl_num[gene_id] = [0, defaultdict(int)]
                    dict_pkl_distance[gene_id] = defaultdict(list)

                with open(f'{opt.intermediate_dir}/NeighborDetection/task_results{i}.pkl', 'wb') as f:
                    pickle.dump(dict_pkl_num, f)
                with open(f'{opt.intermediate_dir}/Distance/task_results{i}.pkl', 'wb') as f:
                    pickle.dump(dict_pkl_distance, f)

        # Initial task distribution
        for i in range(1, size):
            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=i, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=i, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1

        # Collect results and continue distributing remaining tasks
        while active_workers > 0:
            status = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            if result_handler is not None and not intermediate_save:
                result_handler(result)
            else:
                task_results.append(result)
            active_workers -= 1
            task_completed += 1

            if intermediate_save:
                print(f"asyn proc num: {len(monitor_list)}")
                # Wait for pkl update thread if monitor_list is too large
                while len(monitor_list) >= 3:
                    print(f"asyn proc num: {len(monitor_list)}, waiting...")
                    time.sleep(10)

                if len(task_results) == intermediate_size or task_completed == len(tasks):
                    task_result_sub = copy.copy(task_results)
                    task_results = []
                    monitor_list.append(0)
                    future = executor.submit(update_pkl_file, task_result_sub, gene_id_list_split, opt, task_tag, monitor_list)

            print(f"{task_tag} processing: Task {task_completed} completed, {len(tasks) - task_completed} remaining.")

            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=status.source, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=status.source, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        executor.shutdown()

        return task_results

    else:
        while True:
            task = comm.recv(source=0, tag=MPI.ANY_TAG)
            if task is None:
                break
            result = task_processor(task)
            robust_send(result, dest=0, tag=1, comm=comm, rank=rank, max_retries=10, retry_interval=60)


def result_combine(file_list, save_path):
    with open(save_path, 'w', newline='') as f:
        output_writer = csv.writer(f)
        # Initialize a flag to skip the header row of subsequent files
        skip_header = False

        for file in file_list:
            with open(file, 'r', newline='') as input_csvfile:
                input_reader = csv.reader(input_csvfile)

                # If it's the first file, keep the header row
                if not skip_header:
                    for row in input_reader:
                        output_writer.writerow(row)
                    skip_header = True
                else:
                    # Skip the header row for subsequent files and append data to the output file
                    next(input_reader)  # Skip header row
                    for row in input_reader:
                        output_writer.writerow(row)


def split_dictionary(input_dict, num_splits):
    """
    Split a dictionary into a specified number of sub-dictionaries (not necessarily even).
    Args:
        input_dict (dict): Original dictionary to split.
        num_splits (int): Number of sub-dictionaries.
    Returns:
        list: List containing num_splits sub-dictionaries.
    """
    dict_len = len(input_dict)
    if num_splits > dict_len:
        raise ValueError("Number of splits cannot exceed dictionary length")

    keys = list(input_dict.keys())
    new_dicts = [{} for _ in range(num_splits)]
    current_dict_index = 0

    for key in keys:
        new_dicts[current_dict_index][key] = input_dict[key]
        current_dict_index = (current_dict_index + 1) % num_splits

    return new_dicts


def hyper_test_glb_distribution(opt):
    # parser = argparse.ArgumentParser()
    # # parser.add_argument("--detection_method", type=str, choices=['radius', 'nine_grid'],
    # #                     default='radius', help="Method for neighbor detection, can be 'radius' or 'nine_grid'")
    # parser.add_argument("--r_check", type=float, default=0.5, help="radius of checking")
    # parser.add_argument("--r_dist", type=float, default=1.0, help="radius of distribution")
    # parser.add_argument("--around_count_threshold", type=int, default=100, help="around count threshold")
    # parser.add_argument("--column_name", type=str, default="x,y,z,geneID,cell", help="column name used in data")
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
    # parser.add_argument("--molecule_distribution_id_path", type=str,
    #                     default="/data2/yangxr009/ST_STA_stereo/hyperTest_re/E16.5_E1S3_WholeBrain_GEM_CellBin_molecule_distribution.pkl",
    #                     help="path of molecule distribution saving")
    # parser.add_argument("--intermediate_dir", type=str,
    #                     default="/data2/yangxr009/ST_STA_stereo/hyperTest_re/E16.5_E1S3_WholeBrain_GEM_CellBin_hyper_test_intermediate",
    #                     help="path of intermediate result saving")
    # parser.add_argument("--intermediate_split", type=int, default=100, help="interval of intermediate result saving")
    # parser.add_argument("--num_nodes", type=int, default=6, help="Number of nodes")
    # parser.add_argument("--cores_per_node", type=int, default=16, help="Number of cores per node")
    # opt = parser.parse_args()

    comm = pkl5.Intracomm(MPI.COMM_WORLD)
    # comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        # Print the options
        print("--- Options ---")
        for arg_name, arg_value in vars(opt).items():
            print(f"{arg_name}: {arg_value}")
        print("--------------------")

        # Remove existing folders and their contents
        shutil.rmtree(f'{opt.intermediate_dir}/NeighborDetection', ignore_errors=True)
        shutil.rmtree(f'{opt.intermediate_dir}/Distance', ignore_errors=True)
        shutil.rmtree(f'{opt.intermediate_dir}/HyperTest', ignore_errors=True)
        shutil.rmtree(f'{opt.intermediate_dir}/DistanceShape', ignore_errors=True)
        shutil.rmtree(f'{opt.intermediate_dir}/Distance_split', ignore_errors=True)

        # Create folders
        os.makedirs(f'{opt.intermediate_dir}/NeighborDetection')
        os.makedirs(f'{opt.intermediate_dir}/Distance')
        os.makedirs(f'{opt.intermediate_dir}/HyperTest')
        os.makedirs(f'{opt.intermediate_dir}/DistanceShape')
        os.makedirs(f'{opt.intermediate_dir}/Distance_split')

    if rank == 0:
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

        # Calculate the number of transcripts for each geneID in the filtered DataFrame
        gene_expression_num = pd.DataFrame(df_flt_region['geneID'].value_counts()).reset_index()
        gene_expression_num.columns = ['geneName', 'Number']
        # Distribute genes evenly across multiple sublists
        sorted_genes = gene_expression_num.sort_values('Number', ascending=False)['geneName'].tolist()
        gene_expression_num.set_index('geneName', inplace=True)

        # Split the sorted genes into sublists
        num_chunks = opt.intermediate_split
        chunks = [[] for _ in range(num_chunks)]
        # Distribute genes evenly across chunks
        for idx, gene in enumerate(sorted_genes):
            chunk_index = idx % num_chunks
            chunks[chunk_index].append(gene)

        # flatten the list of lists into a single list
        gene_id_list = [gene for chunk in chunks for gene in chunk]

        # Calculate the number of unique genes in DataFrame
        gene_num_dict = df_flt_region.groupby('geneID').size().to_dict()
        print("len of gene_num_dict: ", len(gene_num_dict))

        # Calculate the total number of pixels in the slice
        pixel_num_slice_all = len(df_flt_region)

        # Group the DataFrame by 'cell' to prepare for parallel processing
        cell_group = df_flt_region.groupby('cell')
        cell_group_list = [cell_group.get_group(cell_id) for cell_id in cell_group.groups]
        print("cell num: ", len(cell_group_list))

        df_flt_region = None
        cell_group = None

    else:
        gene_id_list = None
        gene_num_dict = None
        pixel_num_slice_all = None
        cell_group_list = None

    gene_num_dict = comm.bcast(gene_num_dict, root=0)
    pixel_num_slice_all = comm.bcast(pixel_num_slice_all, root=0)

    min_gene_number_local = opt.min_neighbor_number
    expression_level_local = opt.expression_level

    # Task 1: Region Function
    partial_func = partial(region_function_cell_glb_distribution, r_check=opt.r_check, r_dist=opt.r_dist)

    if rank == 0:
        dict_around_list = distribute_tasks_dynamic(comm, rank, size, cell_group_list, partial_func, opt, 'NeighborDetection', intermediate_split=opt.intermediate_split, intermediate_save=True, gene_id_list=gene_id_list)
    else:
        distribute_tasks_dynamic(comm, rank, size, cell_group_list, partial_func, opt, 'NeighborDetection', intermediate_split=opt.intermediate_split, intermediate_save=True, gene_id_list=gene_id_list)
        dict_around_list = None

    comm.Barrier()

    cell_group_list = None

    around_file_list = [f for f in os.listdir(f'{opt.intermediate_dir}/NeighborDetection') if f.startswith('task_results')]
    around_file_list.sort()
    distance_file_list = [f for f in os.listdir(f'{opt.intermediate_dir}/Distance') if f.startswith('task_results')]
    distance_file_list.sort()

    if rank == 0:
        # Merge results from all processes
        print("Processing extra distance files...")
        distance_inter_new_path = f'{opt.intermediate_dir}/Distance_split'
        os.makedirs(distance_inter_new_path, exist_ok=True)
        end_index = opt.intermediate_split
        for i_ex in range(opt.intermediate_split):
            im_temps = []
            pre_name = f'task_results{i_ex}_extra'
            for file in distance_file_list:
                if pre_name in file:
                    im_temps.append(file)

            if len(im_temps) > 0:
                re_ori = f'task_results{i_ex}.pkl'
                print(f"Cancatenating {re_ori}...")
                with open(f'{opt.intermediate_dir}/Distance/{re_ori}', 'rb') as f:
                    dict_pkl_distance = pickle.load(f)

                for file in im_temps:
                    with open(f'{opt.intermediate_dir}/Distance/{file}', 'rb') as f:
                        dict_pkl_distance_extra = pickle.load(f)
                    for gene_id, dist_dict in dict_pkl_distance_extra.items():
                        for gene_id_around, dist_list in dist_dict.items():
                            dict_pkl_distance[gene_id][gene_id_around].extend(dist_list)

                    # Remove the extra file after concatenation
                    # os.remove(f'{opt.intermediate_dir}/Distance/{file}')

                with open(f'{distance_inter_new_path}/{re_ori}', 'wb') as f:
                    pickle.dump(dict_pkl_distance, f)

                # Check the size of the concatenated file
                re_ori_size = os.path.getsize(f'{distance_inter_new_path}/{re_ori}')
                # If the size exceeds 500MB, split it into smaller files
                if re_ori_size > 500*1024*1024:
                    print(f"Cancatenated {re_ori} is too large, splitting...")
                    sp_num = (re_ori_size // (500*1024*1024)) + 1
                    if len(dict_pkl_distance) < sp_num:
                        sp_num = len(dict_pkl_distance)
                    # Split the dictionary into smaller dictionaries
                    dict_pkl_distance_split = split_dictionary(dict_pkl_distance, sp_num)
                    for i_sp, dict_pkl_distance_sub in enumerate(dict_pkl_distance_split):
                        if i_sp == 0:
                            with open(f'{distance_inter_new_path}/{re_ori}', 'wb') as f:
                                pickle.dump(dict_pkl_distance_sub, f)
                        else:
                            with open(f'{distance_inter_new_path}/task_results{end_index}.pkl', 'wb') as f:
                                pickle.dump(dict_pkl_distance_sub, f)
                            end_index += 1
            else:
                shutil.copy(f'{opt.intermediate_dir}/Distance/task_results{i_ex}.pkl', f'{distance_inter_new_path}/task_results{i_ex}.pkl')

        # Memory cleanup
        dict_pkl_distance = None
        dict_pkl_distance_split = None

    # Synchronize all processes before proceeding
    comm.Barrier()

    # Update distance_file_list
    distance_file_list = [f for f in os.listdir(f'{opt.intermediate_dir}/Distance_split') if f.startswith('task_results')]
    distance_file_list.sort()

    if rank == 0:
        hyper_writer = ChunkedCSVWriter(
            f'{opt.intermediate_dir}/HyperTest', 'HyperTest', transform=prepare_hyper_result
        )
        shape_writer = ChunkedCSVWriter(
            f'{opt.intermediate_dir}/DistanceShape', 'DistanceShape',
            transform=lambda df: shape_correct(df, opt),
        )
    else:
        hyper_writer = None
        shape_writer = None

    for i, file in enumerate(around_file_list):
        if rank == 0:
            print(f"Hyper test: Processing intermediate file {file}...")
            with open(f'{opt.intermediate_dir}/NeighborDetection/{file}', 'rb') as f:
                dict_around_list = pickle.load(f)
            tuple_around = [(gene_id, pixel_num_around, gene_around_dict) for gene_id, (pixel_num_around, gene_around_dict) in dict_around_list.items()]
            print("len of tuple_around: ", len(tuple_around))
            dict_around_list = None
        else:
            tuple_around = None

        # Task 2: Hyper Test
        # Use partial function to pass parameters to dynamic task distributor
        partial_func2 = partial(hyper_test, gene_num_dict=gene_num_dict,
                                pixel_num_slice_all=pixel_num_slice_all, min_gene_number=min_gene_number_local,
                                expression_level=expression_level_local)

        if rank == 0:
            distribute_tasks_dynamic(
                comm, rank, size, tuple_around, partial_func2, opt, 'HyperTest',
                result_handler=hyper_writer.write,
            )
        else:
            distribute_tasks_dynamic(comm, rank, size, tuple_around, partial_func2, opt, 'HyperTest')

        comm.Barrier()

    if rank == 0:
        shape_eligible_pairs = collect_qvalue_pairs(
            hyper_writer.paths, opt.shape_qvalue_threshold
        )
        print(f"Directed pairs eligible for shape calculation (qvalue_BH < "
              f"{opt.shape_qvalue_threshold}): {len(shape_eligible_pairs)}")
    else:
        shape_eligible_pairs = None

    if rank == 0:
        task_gene_dict = {}
    else:
        task_gene_dict = None

    for i, file in enumerate(distance_file_list):
        if rank == 0:
            print(f"Distance Shape: Processing intermediate file {file}...")
            task_id = int(file.replace('task_results', '').replace('.pkl', ''))
            with open(f'{opt.intermediate_dir}/Distance_split/{file}', 'rb') as f:
                dict_distance_list = pickle.load(f)
            tuple_dist = build_distribution_shape_tasks(
                dict_distance_list,
                around_count_threshold=opt.around_count_threshold,
                shape_max_distance_count=opt.shape_max_distance_count,
                target_task_count=max(1, (size - 1) * 4),
                eligible_pairs=shape_eligible_pairs,
            )
            shape_pair_count = sum(len(batch) for batch in tuple_dist)
            print(f"Prepared {shape_pair_count} q-value-eligible pairs in "
                  f"{len(tuple_dist)} balanced distribution-shape tasks.")
            gene_keys = list(dict_distance_list.keys())
            task_gene_dict[task_id] = gene_keys

            dict_distance_list = None
        else:
            tuple_dist = None

        # Task 3: Distribution Shape Calculation
        partial_func3 = distribution_shape_calculate_batch

        if rank == 0:
            distribute_tasks_dynamic(
                comm, rank, size, tuple_dist, partial_func3, opt, 'ShapeCalculate',
                result_handler=shape_writer.write,
            )
        else:
            distribute_tasks_dynamic(comm, rank, size, tuple_dist, partial_func3, opt, 'ShapeCalculate')

        comm.Barrier()

    if rank == 0:
        # Save the gene IDs assigned to each distance task.
        with open(f'{opt.intermediate_dir}/molecule_distribution_ID.pkl', 'wb') as f:
            pickle.dump(task_gene_dict, f)

        shape_output_path = opt.save_path.replace('.csv', '_shape.csv')
        stream_postprocess(
            hyper_writer.paths,
            opt.save_path,
            opt.save_path.replace('.csv', f'_dedup_{opt.filter_threshold}_post_proc_shape.csv'),
            filter_threshold=opt.filter_threshold,
            keep=opt.pair_keep,
            work_dir=f'{opt.intermediate_dir}/PostProcess',
            shape_files=shape_writer.paths,
            shape_output_path=shape_output_path if shape_writer.paths else None,
        )


