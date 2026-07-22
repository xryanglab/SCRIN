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
from tqdm import tqdm
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from functools import partial
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from rtree import index
import time
from scrin.core.neighborhood import region_function_without_cell
from scrin.tools.result_proc import add_pair_column, test_result_df_filter, test_result_df_ratio_proc
from scrin.tools.stream_results import ChunkedCSVWriter, prepare_hyper_result, stream_postprocess
from scrin.core.parallel import split_list_into_sublists, robust_send
from scrin.core.statistics import test_function, hyper_test
from scrin.core.intermediate import merge_dicts, update_pkl_file


def large_bcast(data, comm, rank, size, root=0):
    # If root process, broadcast data to all other processes
    if rank == root:
        for i in range(size):
            if i != root:
                comm.send(data, dest=i, tag=77)
    else:
        # Non-root processes receive data from the root process
        data = comm.recv(source=root, tag=77)
    return data


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

        if gene_id_list is not None:
            if len(gene_id_list) < intermediate_split:
                intermediate_split = len(gene_id_list)
                intermediate_size = len(tasks) // intermediate_split
                gene_id_list_split = split_list_into_sublists(gene_id_list, intermediate_split)
            else:
                gene_id_list_split = split_list_into_sublists(gene_id_list, intermediate_split)

            # Create empty dictionaries for each split and save as pkl
            for i, gene_id_list in enumerate(gene_id_list_split):
                dict_pkl = {}
                for gene_id in gene_id_list:
                    dict_pkl[gene_id] = [0, defaultdict(int)]
                with open(f'{opt.intermediate_dir}/{task_tag}/task_results{i}.pkl', 'wb') as f:
                    pickle.dump(dict_pkl, f)

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
    """
    Combine multiple CSV files into a single file.
    """
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


def hyper_test_glb_nocell(opt):
    # # Parse command-line arguments
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--detection_method", type=str, choices=['radius', 'nine_grid'],
    #                     default='radius', help="Method for neighbor detection, can be 'radius' or 'nine_grid'")
    # parser.add_argument("--r_check", type=float, default=None, help="radius of checking")
    # parser.add_argument("--grid_check", type=int, default=None,
    #                     help="grid size for nine_grid detection method, default is 1")
    # parser.add_argument("--rect_length", type=float, default=20,
    #                     help="length of the rectangle for splitting the data, recommended value is the cell diameter")
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
    # parser.add_argument("--rtree_path", type=str, default=None, help="path of rtree index", required=True)
    # parser.add_argument("--intermediate_dir", type=str,
    #                     default="/data2/yangxr009/ST_STA_stereo/hyperTest_re/E16.5_E1S3_WholeBrain_GEM_CellBin_hyper_test_intermediate",
    #                     help="path of intermediate result saving")
    # parser.add_argument("--intermediate_split", type=int, default=100, help="interval of intermediate result saving")
    # parser.add_argument("--num_nodes", type=int, default=6, help="Number of nodes")
    # parser.add_argument("--cores_per_node", type=int, default=16, help="Number of cores per node")
    # opt = parser.parse_args()

    # comm = pkl5.Intracomm(MPI.COMM_WORLD)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        # Print the options
        print("--- Options ---")
        for arg_name, arg_value in vars(opt).items():
            print(f"{arg_name}: {arg_value}")
        print("--------------------")

        # Check the options of detection method
        # if opt.detection_method == 'radius':
        #     if opt.r_check is None:
        #         raise ValueError("Detection method 'radius' requires --r_check to be set.")
        #
        # if opt.detection_method == 'nine_grid':
        #     if opt.grid_check is None:
        #         raise ValueError("Detection method 'nine_grid' requires --grid_check to be set.")

        # Remove existing folders and their contents
        shutil.rmtree(f'{opt.intermediate_dir}/NeighborDetection', ignore_errors=True)
        shutil.rmtree(f'{opt.intermediate_dir}/HyperTest', ignore_errors=True)

        # Create folders
        os.makedirs(f'{opt.intermediate_dir}/NeighborDetection')
        os.makedirs(f'{opt.intermediate_dir}/HyperTest')

    if rank == 0:
        column_names = opt.column_name.split(',')
        if len(column_names) == 3:
            df_flt = pd.read_csv(opt.data_path, sep=',', header=0, usecols=column_names,
                                 dtype={column_names[0]: float, column_names[1]: float,
                                        column_names[2]: str},
                                 na_values='', keep_default_na=False)

            df_flt.rename(columns={column_names[0]: 'x',
                                   column_names[1]: 'y',
                                   column_names[2]: 'geneID'
                                   }, inplace=True)
            print(f"The column {column_names[0]} is used as x, {column_names[1]} as y, "
                  f"{column_names[2]} as geneID.")
            df_flt['z'] = 0.0

        else:
            raise ValueError("Invalid number of columns specified in --column_name. Please provide three columns: x, y, and geneID."
                             " Note: This detection method used for cell-free ID data is not applicable to multiple z-axes.")

        if df_flt.isnull().any().any():
            raise ValueError("Data contains null values. Process aborted.")

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
        if opt.detection_method == 'nine_grid':
            df_flt_region = df_flt_region.drop_duplicates(subset=['geneID', 'x', 'y', 'z'])

        # Get unique gene IDs from the filtered DataFrame
        gene_id_list = df_flt_region['geneID'].unique()

        # Calculate the number of unique genes in DataFrame
        gene_num_dict = df_flt_region.groupby('geneID').size().to_dict()
        print("len of gene_num_dict: ", len(gene_num_dict))

        # Calculate the total number of pixels in the slice
        if opt.detection_method == 'nine_grid':
            # For nine_grid, pixel_num_slice_all is the number of unique (x, y, z) coordinates
            pixel_num_slice_all = len(df_flt_region[['x', 'y', 'z']].drop_duplicates())
        else:
            # For radius detection, pixel_num_slice_all is the total number of rows in the DataFrame
            pixel_num_slice_all = len(df_flt_region)

        # 建立r树索引
        min_x = min(df_flt_region['x'])
        min_y = min(df_flt_region['y'])
        max_x = max(df_flt_region['x'])
        max_y = max(df_flt_region['y'])

        point_list = []

        # Load or build rtree index
        if opt.rtree_path is not None and os.path.exists(opt.rtree_path + '.dat'):
            idx = index.Index(opt.rtree_path)
            print("Load rtree index successfully.")
            x_list = df_flt_region['x'].tolist()
            y_list = df_flt_region['y'].tolist()
            gene_id_region_list = df_flt_region['geneID'].tolist()
            point_list = list(zip(x_list, y_list, gene_id_region_list))
        else:
            if opt.rtree_path is None:
                idx = index.Index()
                print("Building index in memory...")
            else:
                idx = index.Index(opt.rtree_path)
                print(f"Start to build index and save to {opt.rtree_path}.")

            x_vals = df_flt_region['x'].values
            y_vals = df_flt_region['y'].values
            gene_ids = df_flt_region['geneID'].values

            for i in tqdm(range(len(df_flt_region))):
                x, y, label = x_vals[i], y_vals[i], gene_ids[i]
                idx.insert(i, (x, y, x, y), obj=label)
                point_list.append((x, y, label))

            print("Index built successfully.")

        # Generate rectangles for detection
        rect_list, point_comb_list = [], []
        detect_length = opt.rect_length
        for min_y_rect in np.arange(min_y, max_y, detect_length):
            for min_x_rect in np.arange(min_x, max_x, detect_length):
                if min_y_rect + detect_length > max_y and min_x_rect + detect_length > max_x:
                    rect = (min_x_rect, min_y_rect, max_x, max_y)
                elif min_y_rect + detect_length > max_y and min_x_rect + detect_length <= max_x:
                    rect = (min_x_rect, min_y_rect, min_x_rect + detect_length, max_y)
                elif min_y_rect + detect_length <= max_y and min_x_rect + detect_length > max_x:
                    rect = (min_x_rect, min_y_rect, max_x, min_y_rect + detect_length)
                else:
                    rect = (min_x_rect, min_y_rect, min_x_rect + detect_length, min_y_rect + detect_length)
                rect_list.append(rect)
        print("Rect list generated.")

        print("Start to build point comb list.")
        if opt.detection_method == 'radius':
            additional_check = opt.r_check
        else:
            additional_check = opt.grid_check

        for rect in tqdm(rect_list):
            hits = list(idx.intersection(rect))
            if len(hits) != 0:
                point_re = [point_list[i] for i in hits]

                hits_ex = list(idx.intersection(
                    (rect[0] - additional_check, rect[1] - additional_check, rect[2] + additional_check, rect[3] + additional_check)))
                point_re_ex = [point_list[i] for i in hits_ex]

                point_comb_list.append([point_re, point_re_ex])
        print("Point comb list generated.")

        idx.close()

        df_flt_region = None
        point_list = None
        idx = None

    else:
        gene_id_list = None
        gene_num_dict = None
        pixel_num_slice_all = None
        point_comb_list = None

    gene_num_dict = comm.bcast(gene_num_dict, root=0)
    pixel_num_slice_all = comm.bcast(pixel_num_slice_all, root=0)

    min_gene_number_local = opt.min_neighbor_number
    expression_level_local = opt.expression_level

    if opt.detection_method == 'radius':
        partial_func = partial(region_function_without_cell, r_check=opt.r_check, detection_method=opt.detection_method)
    else:
        partial_func = partial(region_function_without_cell, r_check=opt.grid_check, detection_method=opt.detection_method)

    if rank == 0:
        dict_around_list = distribute_tasks_dynamic(comm, rank, size, point_comb_list, partial_func, opt, 'NeighborDetection', intermediate_split=opt.intermediate_split, intermediate_save=True, gene_id_list=gene_id_list)
    else:
        distribute_tasks_dynamic(comm, rank, size, point_comb_list, partial_func, opt, 'NeighborDetection', intermediate_split=opt.intermediate_split, intermediate_save=True, gene_id_list=gene_id_list)
        dict_around_list = None

    comm.Barrier()

    around_file_list = [f for f in os.listdir(f'{opt.intermediate_dir}/NeighborDetection') if f.startswith('task_results')]
    around_file_list.sort()

    if rank == 0:
        hyper_writer = ChunkedCSVWriter(
            f'{opt.intermediate_dir}/HyperTest', 'HyperTest', transform=prepare_hyper_result
        )
    else:
        hyper_writer = None

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

        # Task 2: Perform hypergeometric test on the detected neighbors
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
        save_path = opt.save_path.replace('.csv', f'_dedup_{opt.filter_threshold}_post_proc.csv')
        stream_postprocess(
            hyper_writer.paths,
            opt.save_path,
            save_path,
            filter_threshold=opt.filter_threshold,
            keep=opt.pair_keep,
            work_dir=f'{opt.intermediate_dir}/PostProcess',
        )
