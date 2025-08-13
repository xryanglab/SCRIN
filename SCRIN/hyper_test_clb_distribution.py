import argparse
import os
import csv
import pickle
import joblib
import time
import pandas as pd
import numpy as np
import shutil
import msgpack
import zlib
from mpi4py import MPI
from mpi4py.util import pkl5
from tqdm import tqdm
from scipy.stats import hypergeom, gaussian_kde, skew, kurtosis
from statsmodels.stats.multitest import multipletests
from functools import partial
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from threading import Event
from SCRIN.tools.result_proc import add_pair_column, test_result_df_filter, test_result_df_ratio_proc


def split_list_into_sublists(input_list, num_sublists):
    avg = len(input_list) // num_sublists
    remainder = len(input_list) % num_sublists
    sublists = []

    start = 0
    for i in range(num_sublists):
        sublist_size = avg + (1 if i < remainder else 0)
        sublist = input_list[start:start + sublist_size]
        sublists.append(sublist)
        start += sublist_size

    return sublists


def large_bcast(data, comm, rank, size, root=0):
    # If it is the root process, send data to all other processes
    if rank == root:
        for i in range(size):
            if i != root:
                comm.send(data, dest=i, tag=77)
    else:
        # Non-root process receives data from root process
        data = comm.recv(source=root, tag=77)
    return data


# Function for receiving large data from all non-root processes
def large_recv(data, comm, rank, size, root=0):
    # If it is the root process, send data to all other processes
    recv_list = []
    if rank == root:
        for i in range(size):
            if i != root:
                recv_data = comm.recv(source=i, tag=88)
                recv_list.append(recv_data)
        return recv_list
    else:
        # Non-root process sends data to root process
        comm.send(data, dest=root, tag=88)
        return None  # Non-root processes do not need to return data


def compress_dict(data_dict):
    packed = msgpack.packb(data_dict, use_bin_type=True)
    compressed = zlib.compress(packed)
    return compressed

def decompress_dict(compressed_data):
    packed = zlib.decompress(compressed_data)
    return msgpack.unpackb(packed, raw=False)


def robust_send(data, dest, tag, comm, rank, max_retries=5, retry_interval=1):
    retries = 0
    while retries < max_retries:
        try:
            comm.send(data, dest=dest, tag=tag)
            # print("Send successful")
            return  # Successfully sent, exit function
        except Exception as e:  # Catch all exceptions
            print(f"Rank {rank}: Send failed with exception: {e}. Retrying...")
            time.sleep(retry_interval)
            retries += 1
    print("Failed to send data after retries. Raising exception.")
    raise  # Rethrows the last caught exception


def send_gene_list(gene_list, comm, rank, size):
    gene_rank_dict = {}

    sub_lists = split_list_into_sublists(gene_list, size-1)
    # Distribute the sublist to each non-root process in turn, and record the gene ID that each root process is responsible for
    for i in range(1, size):
        send_data = sub_lists[i-1]
        robust_send(send_data, dest=i, tag=99, comm=comm, rank=rank, max_retries=10, retry_interval=60)
        gene_rank_dict[i] = send_data

    return gene_rank_dict


def get_gene_rank(gene_list, size):
    gene_rank_dict = {}

    sub_lists = split_list_into_sublists(gene_list, size-1)
    # Distribute the sublist to each non-root process in turn, and record the gene ID that each root process is responsible for
    for i in range(1, size):
        send_data = sub_lists[i-1]
        gene_rank_dict[i] = send_data

    return gene_rank_dict  # structure: {rank_id: [gene_id1, gene_id2, ...]}


def merge_dicts(count_matrix, gene2idx):
    merged_list = []  # Create an empty list to store the merge results

    merged_dict = defaultdict(lambda: defaultdict(list))  # Create a default dictionary to temporarily store the merge results

    idx2gene = {v: k for k, v in gene2idx.items()}  # Reverse a dictionary for indexing

    b_around_matrix = count_matrix[:, :, 0]  # gene_B_around
    row_indices, col_indices = np.nonzero(b_around_matrix)  # Get the row and column indices of nonzero elements
    non_zero_indices = list(zip(row_indices, col_indices))

    for r, c in tqdm(non_zero_indices, desc="Getting dictionaries"):
        merged_dict[idx2gene[r]][idx2gene[c]].extend([
            count_matrix[r, c, 0],
            count_matrix[r, c, 1],
            count_matrix[r, c, 2],
            count_matrix[r, c, 3],
            count_matrix[r, c, 4]
        ])

    # Convert the temporarily stored result into a list
    for gene_id, sub_dict in merged_dict.items():
        merged_list.append((gene_id, sub_dict))  # Append the results to a list as a tuple

    return merged_list


def distribute_tasks_dynamic(comm, rank, size, tasks, task_processor, task_tag, result_recv=True):
    """
    Dynamically assign tasks to processes and collect results
    :param comm: MPI communicator
    :param rank: Rank of the process
    :param size: Total number of processes
    :param tasks: A task list, where each task is a processable unit
    :param task_processor: A function that receives a task as input and returns the result of the processing
    :param task_tag: The tag of the task
    :param result_recv: If True, receive the processing result; otherwise, do not
    :return: If it is rank 0, return a list of processing results of all tasks; otherwise return None
    """

    task_results = []

    if rank == 0:
        task_index = 0
        active_workers = 0
        task_completed = 0

        # Initial distribution task
        for i in range(1, size):
            if task_index < len(tasks):
                # comm.send(tasks[task_index], dest=i, tag=0)
                robust_send(tasks[task_index], dest=i, tag=133, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=i, tag=133, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1

        # Collect the results and continue distributing the remaining tasks
        while active_workers > 0:
            status = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=155, status=status)

            if result_recv:
                task_results.append(result)

            active_workers -= 1
            task_completed += 1

            print(f"{task_tag} processing: Task {task_completed} completed, {len(tasks) - task_completed} remaining.")

            if task_index < len(tasks):
                # comm.send(tasks[task_index], dest=status.source, tag=0)
                robust_send(tasks[task_index], dest=status.source, tag=133, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                # comm.send(None, dest=status.source, tag=2)
                robust_send(None, dest=status.source, tag=133, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        return task_results

    else:
        while True:
            task = comm.recv(source=0, tag=133)
            if task is None:
                break
            result = task_processor(task)
            # comm.send(result, dest=0)
            robust_send(result, dest=0, tag=155, comm=comm, rank=rank, max_retries=10, retry_interval=60)


def update_recv_asyn(comm, matrix, gene2idx, local_gene2idx, is_idle, local_pool=None, local_pool_lock=None, stop_event=None):
    while not stop_event.is_set():
        while comm.iprobe(source=MPI.ANY_SOURCE, tag=177):
            is_idle.clear()
            # recv_update = comm.recv(source=MPI.ANY_SOURCE, tag=177)
            recv_update_compressed = comm.recv(source=MPI.ANY_SOURCE, tag=177)
            recv_update = decompress_dict(recv_update_compressed)  # Decompress the received update data
            update_local_matrix(matrix, recv_update, gene2idx, local_gene2idx)
            is_idle.set()

        with local_pool_lock:
            if local_pool is not None and len(local_pool) > 0:
                # Check if there are any updates in the local pool that need to be processed
                for local_update in local_pool:
                    is_idle.clear()
                    update_local_matrix(matrix, local_update, gene2idx, local_gene2idx)
                    is_idle.set()
                local_pool.clear()  # Clear the local pool
        time.sleep(0.01)  # Reduce CPU usage
    print("Async update thread exiting.")


def distribute_tasks_dynamic_region(comm, rank, size, tasks, task_processor, task_tag, matrix, gene2idx, local_gene2idx):
    """
    Dynamically assign tasks to processes and collect results
    :param comm: MPI communicator.
    :param rank: Rank of the process.
    :param size: Total number of processes.
    :param tasks: A task list, where each task is a processable unit.
    :param task_processor: A function that receives a task as input and returns the result of the processing.
    :param task_tag: The tag of the task.
    :param matrix: The matrix to be updated.
    :param gene2idx: A dictionary mapping gene IDs to their indices in the matrix.
    :param local_gene2idx: A dictionary mapping local gene IDs to their indices in the matrix.
    :return: If it is rank 0, return a list of processing results of all tasks; otherwise return None.
    """

    if rank == 0:
        task_index = 0
        active_workers = 0
        task_completed = 0

        # Initial distribution task
        for i in range(1, size):
            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=i, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=i, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1

        # Collect the results and continue distributing the remaining tasks
        while active_workers > 0:
            status = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=1, status=status)

            active_workers -= 1
            task_completed += 1

            print(f"{task_tag} processing: Task {task_completed} completed, {len(tasks) - task_completed} remaining.")

            if task_index < len(tasks):
                # time_send_start = time.time()
                robust_send(tasks[task_index], dest=status.source, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                # time_send_end = time.time()
                # print(f"[Rank {rank}] Task {task_index} sent to worker {status.source} in {time_send_end - time_send_start:.4f} seconds.")
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=status.source, tag=0, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        # Wait for all workers to finish and send their status
        done_flags = [False] * size
        done_flags[0] = True  # Rank 0 is always done
        time_waiter = 0  # Wait time counter

        while not (all(done_flags) and time_waiter >= 30):
            while comm.iprobe(source=MPI.ANY_SOURCE, tag=3):
                status = MPI.Status()
                msg = comm.recv(source=MPI.ANY_SOURCE, tag=3, status=status)
                src = status.Get_source()
                if msg == "done":
                    # print(f"[Rank {rank}] Received 'done' signal from process {src}. Marking as completed.")
                    done_flags[src] = True
                elif msg == "not_done":
                    # print(f"[Rank {rank}] Received 'not_done' signal from process {src}. Marking as incomplete.")
                    done_flags[src] = False

            time.sleep(0.05)  # Reduce CPU usage

            if all(done_flags):
                time_waiter += 0.05
                print(f"[Rank {rank}] All processes marked as done. Waiting to confirm stability: {time_waiter:.2f} seconds")
            else:
                # if time_waiter > 0:
                #     print(f"[Rank {rank}] Status changed. Resetting stability wait timer.")
                time_waiter = 0
                # print(f"[Rank {rank}] Current done_flags status: {done_flags}")

        # Notify all processes that all tasks are completed
        print(f"All tasks completed. Notifying all processes to exit. Rank {rank}")
        for i in range(1, size):
            robust_send("all_done", dest=i, tag=4, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        return None

    else:
        executor = ThreadPoolExecutor(max_workers=1)  # Create a separate thread for updating the matrix
        is_idle = Event()
        is_idle.set()  # Initially set to idle
        stop_event = Event()  # Asynchronous exit flag

        task_finished = False  # Task finished flag for non-root processes

        local_matrix_pool = []  # Local update pool for storing local updates
        local_pool_lock = Lock()  # Lock to protect the local update pool

        future = executor.submit(update_recv_asyn, comm, matrix, gene2idx, local_gene2idx, is_idle, local_matrix_pool, local_pool_lock, stop_event)

        while True:
            if comm.iprobe(source=0, tag=0):
                task = comm.recv(source=0, tag=0)
                if task is None:
                    # print(f"Rank {rank} received None from root rank, no more tasks.")
                    task_finished = True
            else:
                task = None

            now_message = comm.iprobe(source=MPI.ANY_SOURCE, tag=177)

            while task is None and is_idle.is_set() and len(local_matrix_pool) == 0 and not now_message and task_finished:
                # If there are no tasks and the current process is idle, wait for new tasks or state changes
                # print(f"Rank {rank} is idle and has no local updates. Waiting for new tasks or region processing state change.")
                time.sleep(1.0)  # Send an idle signal every second
                # Send idle signal to root process
                info_done = comm.isend("done", dest=0, tag=3)
                info_done.wait()  # Wait for the send to complete

                # Check if there are any update requests from other processes, if so, stop waiting and send not_done signal
                now_message = comm.iprobe(source=MPI.ANY_SOURCE, tag=177)
                task_from_root = comm.iprobe(source=0, tag=0)
                if now_message:
                    # print(f"Rank {rank} received an update request while waiting. Processing updates.")
                    info_not_done = comm.isend("not_done", dest=0, tag=3)
                    info_not_done.wait()
                    break

                if task_from_root:
                    # print(f"Rank {rank} received a new task from root while waiting.")
                    info_not_done = comm.isend("not_done", dest=0, tag=3)
                    info_not_done.wait()
                    break

                if comm.iprobe(source=0, tag=4):  # Check if the root process has sent an exit signal
                    break

            if task is not None:
                result = task_processor(task)
            else:
                result = None

            if result is not None:
                local_update_dict1 = result[0]  # local_updates1
                local_update_dict2 = result[1]  # local_updates2

                des_ranks1 = list(local_update_dict1.keys())
                des_ranks2 = list(local_update_dict2.keys())
                des_ranks = list(set(des_ranks1 + des_ranks2))

                if len(des_ranks) == 0:
                    pass
                else:
                    merge_update_dict = {_: [[], []] for _ in des_ranks}  # Create a dictionary to merge updates

                    for des_rank in des_ranks:
                        if des_rank in local_update_dict1:
                            merge_update_dict[des_rank][0] = local_update_dict1[des_rank]
                        if des_rank in local_update_dict2:
                            merge_update_dict[des_rank][1] = local_update_dict2[des_rank]

                    for des_rank, updates in merge_update_dict.items():
                        if des_rank == rank:
                            with local_pool_lock:
                                local_matrix_pool.append(updates)  # Add updates to the local pool
                        else:
                            updates_compressed = compress_dict(updates)  # Compress the updates before sending
                            req = comm.isend(updates_compressed, dest=des_rank, tag=177)
                            req.wait()  # Wait for the send to complete

            if result is not None:
                robust_send("yes", dest=0, tag=1, comm=comm, rank=rank, max_retries=10, retry_interval=60)

            # Check if all processes have completed and wait for the root process to notify
            status2 = MPI.Status()
            if comm.iprobe(source=0, tag=4, status=status2):
                msg = comm.recv(source=0, tag=4)
                if msg == "all_done":
                    print(f"Rank {rank} received all_done signal. Exiting.")
                    stop_event.set()
                    break

            time.sleep(0.01)  # Reduce CPU usage

        # executor.shutdown()


def asyn_update_output_dict(output_dict, receive_update):
    gene_id, sub_gene_distribution_dict = receive_update[0], receive_update[1]  # 接收基因id和分布字典
    sub_dict = output_dict.setdefault(gene_id, {})
    for sub_gene_id, distribution_list in sub_gene_distribution_dict.items():
        sub_dict.setdefault(sub_gene_id, []).extend(distribution_list)


def distribution_proc_asyn(comm, output_path, is_idle, tag, local_pool=None, local_pool_lock=None, stop_event=None):
    output_dict = {}
    while not stop_event.is_set():
        while comm.iprobe(source=MPI.ANY_SOURCE, tag=tag):
            is_idle.clear()

            recv_update = comm.recv(source=MPI.ANY_SOURCE, tag=tag)
            asyn_update_output_dict(output_dict, recv_update)  # Update the output dictionary with the received data

            is_idle.set()

        with local_pool_lock:
            if local_pool is not None and len(local_pool) > 0:
                # Check if there are any updates in the local pool that need to be processed
                is_idle.clear()
                for local_update in local_pool:
                    asyn_update_output_dict(output_dict, local_update)
                local_pool.clear()  # Clear the local pool after processing
                is_idle.set()
        time.sleep(0.01)  # Reduce CPU usage
    # Save the output dictionary to a file
    is_idle.clear()
    output_path_file_count = len(os.listdir(output_path))
    output_path_i = output_path + f"/gene_distribution_dict_{output_path_file_count}.pkl"
    joblib.dump(output_dict, output_path_i, compress=3)
    is_idle.set()
    print("Async update thread exiting.")


def distribute_tasks_dynamic_multi_asyn(comm, rank, size, tasks, task_processor, task_tag, asyn_func, distribution_output_path=None):
    """
    Dynamically assign tasks to processes and collect results
    :param comm: MPI communicator。
    :param rank: Rank of the process。
    :param size: Total number of processes。
    :param tasks: A task list, where each task is a processable unit。
    :param task_processor: A function that receives a task as input and returns the result of the processing。
    :param task_tag: The tag of the task。
    :param asyn_func: An asynchronous function to handle the results.
    :param distribution_output_path: The output path for storing results.
    :return: If it is rank 0, return a list of processing results of all tasks; otherwise return None.
    """

    if rank == 0:
        task_index = 0
        active_workers = 0
        task_completed = 0

        # Initial distribution task
        for i in range(1, size):
            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=i, tag=700, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=i, tag=700, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1

        # Collect the results and continue distributing the remaining tasks
        while active_workers > 0:
            status = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=701, status=status)

            active_workers -= 1
            task_completed += 1

            print(f"{task_tag} processing: Task {task_completed} completed, {len(tasks) - task_completed} remaining.")

            if task_index < len(tasks):
                robust_send(tasks[task_index], dest=status.source, tag=700, comm=comm, rank=rank, max_retries=10, retry_interval=60)
                task_index += 1
                active_workers += 1
            else:
                robust_send(None, dest=status.source, tag=700, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        # Wait for all workers to finish and send their status
        done_flags = [False] * size
        done_flags[0] = True  # Rank 0 is always done
        time_waiter = 0  # Wait time counter

        while not (all(done_flags) and time_waiter >= 5):
            while comm.iprobe(source=MPI.ANY_SOURCE, tag=703):
                status = MPI.Status()
                msg = comm.recv(source=MPI.ANY_SOURCE, tag=703, status=status)
                src = status.Get_source()
                if msg == "done":
                    # print(f"[Rank {rank}] Received 'done' signal from process {src}. Marking as completed.")
                    done_flags[src] = True
                elif msg == "not_done":
                    # print(f"[Rank {rank}] Received 'not_done' signal from process {src}. Marking as incomplete.")
                    done_flags[src] = False

            time.sleep(0.05)  # Reduce CPU usage

            if all(done_flags):
                time_waiter += 0.05
                print(f"[Rank {rank}] All processes marked as done. Waiting to confirm stability: {time_waiter:.2f} seconds")
            else:
                # if time_waiter > 0:
                #     print(f"[Rank {rank}] Status changed. Resetting stability wait timer.")
                time_waiter = 0
                # print(f"[Rank {rank}] Current done_flags status: {done_flags}")

        # Notify all processes that all tasks are completed
        print(f"All tasks completed. Notifying all processes to exit. Rank {rank}")
        for i in range(1, size):
            robust_send("all_done", dest=i, tag=704, comm=comm, rank=rank, max_retries=10, retry_interval=60)

        return None

    else:
        executor = ThreadPoolExecutor(max_workers=1)  # Create a separate thread for processing results
        is_idle = Event()
        is_idle.set()  # Initially set to idle
        # Asynchronous exit flag
        stop_event = Event()

        task_finished = False  # Task finished flag for non-root processes

        local_matrix_pool = []  # Local update pool for storing local updates
        local_pool_lock = Lock()  # Lock to protect the local update pool

        output_path = distribution_output_path + f"/rank_{rank}"
        future = executor.submit(asyn_func, comm, output_path, is_idle, 777, local_matrix_pool, local_pool_lock, stop_event)

        while True:
            if comm.iprobe(source=0, tag=700):
                task = comm.recv(source=0, tag=700)
                if task is None:
                    # print(f"Rank {rank} received None from root rank, no more tasks.")
                    task_finished = True
            else:
                task = None

            now_message = comm.iprobe(source=MPI.ANY_SOURCE, tag=777)

            while task is None and is_idle.is_set() and len(local_matrix_pool) == 0 and not now_message and task_finished:
                # If there are no tasks and the current process is idle, wait for new tasks or state changes
                # print(f"Rank {rank} is idle and has no local updates. Waiting for new tasks or region processing state change.")
                time.sleep(1.0)  # Send an idle signal every second
                # Send idle signal to root process
                info_done = comm.isend("done", dest=0, tag=703)
                info_done.wait()

                # Check if there are any update requests from other processes, if so, stop waiting and send not_done signal
                now_message = comm.iprobe(source=MPI.ANY_SOURCE, tag=777)
                task_from_root = comm.iprobe(source=0, tag=700)
                if now_message:
                    # print(f"Rank {rank} received an update request while waiting. Processing updates.")
                    info_not_done = comm.isend("not_done", dest=0, tag=703)
                    info_not_done.wait()
                    break

                if task_from_root:
                    # print(f"Rank {rank} received a new task from root while waiting.")
                    info_not_done = comm.isend("not_done", dest=0, tag=703)
                    info_not_done.wait()
                    break

                if comm.iprobe(source=0, tag=704):  # Check if the root process has sent an exit signal
                    break

            if task is not None:
                result = task_processor(task)
            else:
                result = None

            if result is not None:
                for re_i in result:
                    des_rank = re_i[0]  # Get the destination rank
                    updates = re_i[1]  # Get the updates to be sent

                    if des_rank == rank:
                        with local_pool_lock:
                            local_matrix_pool.append(updates)  # Add updates to the local pool
                    else:
                        req = comm.isend(updates, dest=des_rank, tag=777)
                        req.wait()

            if result is not None:
                robust_send("yes", dest=0, tag=701, comm=comm, rank=rank, max_retries=10, retry_interval=60)

            # Check if all processes have completed and wait for the root process to notify
            status2 = MPI.Status()
            if comm.iprobe(source=0, tag=704, status=status2):
                msg = comm.recv(source=0, tag=704)
                if msg == "all_done":
                    print(f"Rank {rank} received all_done signal. Exiting.")
                    stop_event.set()
                    break

            time.sleep(0.01)  # Reduce CPU usage

        while not is_idle.is_set():
            time.sleep(0.05)

        # executor.shutdown()
        return None


def region_function_cell_multi_proc(df_cell, r_check=0.5, gene_rank_id2rank=None):
    """
    Return local updates from a cell dataframe.
    Local update 1: gene B around
    Local update 2: gene around, gene B, gene all, gene A
    """

    local_updates1 = defaultdict(list)  # Store updates of gene B around
    local_updates2 = defaultdict(list)  # Store updates of gene around, gene B, gene all, gene A

    gene_list = df_cell['geneID'].unique()

    # pre z dict
    df_cell_z_group = df_cell.groupby('z')
    df_cell_z_dict = {}
    for z, group in df_cell_z_group:
        df_cell_z_dict[z] = group

    gene_cell_num = len(df_cell)
    gene_cell_num_all_dict = df_cell.groupby('geneID').size().to_dict()
    all_gene_keys = list(gene_cell_num_all_dict.keys())
    for gene in gene_list:
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

    local_updates1 = defaultdict(list)  # Store updates of gene B around
    local_updates2 = defaultdict(list)  # Store updates of gene around, gene B, gene all, gene A

    gene_list = df_cell['geneID'].unique()

    # pre z dict
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
                (abs(x - df_cell['x']) <= r_check) & (abs(y - df_cell['y']) <= r_check)]
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

    local_updates1 = defaultdict(list)  # Store updates of gene B around
    local_updates2 = defaultdict(list)  # Store updates of gene around, gene B, gene all, gene A
    distribution_dict = {}  # Store gene distribution information

    gene_list = df_cell['geneID'].unique()

    # pre z dict
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


def distribution_file_proc(joblib_file_path, gene_rank_id2rank):
    rank_data = joblib.load(joblib_file_path)

    output_list = []
    combine_gene_distribution_dict = {}

    for cell_gene_dict in rank_data:
        for gene_id, gene_around_distribution_dict in cell_gene_dict.items():
            gene_level = combine_gene_distribution_dict.setdefault(gene_id, {})
            for gene_around_id, distribution_list in gene_around_distribution_dict.items():
                gene_level.setdefault(gene_around_id, []).extend(distribution_list)

    # Combine the gene distribution dictionary with gene rank information
    for gene_id, gene_around_distribution_dict in combine_gene_distribution_dict.items():
        gene_rank_belong = gene_rank_id2rank[gene_id]
        output_list.append([gene_rank_belong, [gene_id, gene_around_distribution_dict]])

    return output_list


def update_local_matrix(matrix, local_update_seqs_list, gene2idx, local_gene2idx):
    local_update_seqs1 = local_update_seqs_list[0]  # local_updates1
    local_update_seqs2 = local_update_seqs_list[1]  # local_updates2

    if len(local_update_seqs1) != 0:
        for seq in local_update_seqs1:
            idx1 = local_gene2idx[seq[0]]
            idx2 = gene2idx[seq[1]]
            matrix[idx1, idx2, 0] += seq[2]  # gene_B_around

    if len(local_update_seqs2) != 0:
        for seq in local_update_seqs2:
            idx1 = local_gene2idx[seq[0]]
            idx2 = gene2idx[seq[1]]
            matrix[idx1, idx2, 1] += seq[2]  # gene_around
            matrix[idx1, idx2, 2] += seq[3]  # gene_B
            matrix[idx1, idx2, 3] += seq[4]  # gene_all
            matrix[idx1, idx2, 4] += seq[5]  # gene_A


def merge_distribution_dicts_from_file(intermediate_dir):
    """
    Merge multiple distribution dictionaries from joblib files into a single dictionary.
    """
    merged_distribution_dict = {}
    for file in intermediate_dir:
        top_data = joblib.load(file)
        for gene_id, gene_around_distribution_dict in top_data.items():
            gene_level = merged_distribution_dict.setdefault(gene_id, {})
            for gene_around_id, distribution_list in gene_around_distribution_dict.items():
                gene_level.setdefault(gene_around_id, []).extend(distribution_list)
    return merged_distribution_dict


def test_function(gene_B_around_num, gene_A_N, gene_B_N,
                  gene_slice_num, gene_around_num, minGenenumber, expressionLevel):

    if (gene_around_num < minGenenumber or
            gene_A_N / gene_B_N > expressionLevel or gene_B_N / gene_A_N > expressionLevel):
        return None

    a = hypergeom.sf(k=gene_B_around_num - 1, M=gene_slice_num, n=gene_B_N,
                     N=gene_around_num)
    para_list = [gene_B_around_num, gene_B_N, gene_around_num, gene_slice_num, gene_A_N]

    return [a, para_list]


def hyper_test(tuple_proc, min_gene_number, expression_level):

    gene_id = tuple_proc[0]
    gene_around_dict = tuple_proc[1]

    gene_id_around_list, p_list, para_list_all = [], [], []

    for gene_id_around, around_num_list in gene_around_dict.items():
        tf_list = test_function(gene_B_around_num=around_num_list[0], gene_A_N=around_num_list[4],
                                gene_B_N=around_num_list[2], gene_slice_num=around_num_list[3],
                                gene_around_num=around_num_list[1], minGenenumber=min_gene_number,
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


def distribution_shape_calculate(tuple_proc, around_count_threshold=100):
    gene_id = tuple_proc[0]
    gene_dist_dict = tuple_proc[1]

    gene_id_around_list, dist_list = [], []

    for gene_id_around, dist_list_around in gene_dist_dict.items():
        if len(dist_list_around) >= around_count_threshold:
            gene_id_around_list.append(gene_id_around)
            shape_count, mode, skewness, kurt, median, q25, q75, iqr = compute_shape_para(dist_list_around)
            dist_list.append([shape_count, mode, skewness, kurt, median, q25, q75, iqr])
        else:
            continue

    if len(dist_list) == 0:
        return None

    results_df = pd.DataFrame({
        'gene_A': [gene_id] * len(dist_list),
        'gene_B': gene_id_around_list,
        'shape_count': [_[0] for _ in dist_list],
        'mode': [_[1] for _ in dist_list],
        'skewness': [_[2] for _ in dist_list],
        'kurtosis': [_[3] for _ in dist_list],
        'median': [_[4] for _ in dist_list],
        'q25': [_[5] for _ in dist_list],
        'q75': [_[6] for _ in dist_list],
        'iqr': [_[7] for _ in dist_list]
    })

    return results_df


def shape_correct(df_result, opt):
    df_result.loc[:, 'skewness_adjusted'] = 1 / (1 + np.exp(df_result['skewness']))
    df_result.loc[:, 'skewness_relative'] = df_result['skewness'] + 0.566
    df_result.loc[:, 'kurtosis_relative'] = df_result['kurtosis'] + 0.6
    df_result.loc[:, 'median_relative'] = df_result['median'] - (0.707 * opt.r_dist)
    df_result.loc[:, 'q25_relative'] = df_result['q25'] - (0.5 * opt.r_dist)
    df_result.loc[:, 'q75_relative'] = df_result['q75'] - (0.866 * opt.r_dist)
    return df_result


def hyper_test_clb_distribution(opt):
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--detection_method", type=str, choices=['radius', 'nine_grid'],
    #                     default='radius', help="Method for neighbor detection, can be 'radius' or 'nine_grid'")
    # parser.add_argument("--r_check", type=float, default=None, help="radius of checking")
    # parser.add_argument("--grid_check", type=int, default=None, help="grid size for nine_grid detection method, default is 1")
    # parser.add_argument("--r_dist", type=float, default=None, help="if not None, enable co-localization distribution saving")
    # parser.add_argument("--around_count_threshold", type=int, default=100,
    #                     help="threshold for the number of points around a gene to consider it for distribution analysis")
    # parser.add_argument("--column_name", type=str, default="x,y,z,geneID,cell", help="column name used in data")
    # parser.add_argument("--min_gene_number", type=int, default=5,
    #                     help="minimum number of transcripts for a gene to be considered")
    # parser.add_argument("--min_neighbor_number", type=int, default=1,
    #                     help="minimum number of neighbors for a pair to be considered")
    # parser.add_argument("--expression_level", type=float, default=100,
    #                     help="For gene A and gene B in the pair, the maximum ratio of their expression count.")
    # parser.add_argument("--filter_threshold", type=float, default=0.00001,
    #                     help="filter threshold for qvalue_BH in post processing")
    # parser.add_argument("--distribution_save_interval", type=int, default=10,
    #                     help="interval for saving distribution data to file")
    # parser.add_argument("--intermediate_dir", type=str, default=None,
    #                     help="File path to save co-localization distribution data")
    # parser.add_argument("--pair_keep", type=str, default='last', help="keep method for pair post processing, can be 'first' or 'last'")
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

        # Check the options of distribution
        # if opt.r_dist is not None:
        #     if opt.intermediate_dir is None:
        #         raise ValueError("If r_dist is set, intermediate_dir must be provided.")
        #     if opt.detection_method == 'nine_grid':
        #         raise ValueError("Nine_grid detection method is not compatible with co-localization distribution analysis."
        #                          "Please use 'radius' detection method or close co-localization distribution analysis.")

        # if opt.r_dist is not None:
        if opt.distribution_analysis:
            if os.path.exists(opt.intermediate_dir):
                print(f"Removing existing distribution file path: {opt.intermediate_dir}")
                shutil.rmtree(opt.intermediate_dir)
            os.makedirs(opt.intermediate_dir)
            # make n=size-1 directories, no root rank directory
            for i in range(1, size):
                os.makedirs(os.path.join(opt.intermediate_dir, f"rank_{i}"))

            if os.path.exists(opt.intermediate_dir + '_gene_id_merge'):
                print(f"Removing existing distribution file path for gene ID merge: {opt.intermediate_dir + '_gene_id_merge'}")
                shutil.rmtree(opt.intermediate_dir + '_gene_id_merge')
            os.makedirs(opt.intermediate_dir + '_gene_id_merge')
            for i in range(1, size):
                os.makedirs(os.path.join(opt.intermediate_dir + '_gene_id_merge', f"rank_{i}"))

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

        # Get the unique gene list and create a mapping from gene to index
        gene_list_all = df_flt_region['geneID'].unique()
        gene2idx = {gene: idx for idx, gene in enumerate(gene_list_all)}
        matrix_size = len(gene_list_all)

        # Group the DataFrame by 'cell' to prepare for parallel processing
        cell_group = df_flt_region.groupby('cell')
        cell_group_list = [cell_group.get_group(cell_id) for cell_id in cell_group.groups]
        print("cell num: ", len(cell_group_list))

        df_flt_region = None
        cell_group = None

        # Split the gene list into chunks for each rank
        gene_rank_dict = get_gene_rank(gene_list_all, size)

    else:
        matrix_size = None
        gene_rank_dict = None
        start_time = None
        cell_group_list = None
        gene2idx = None

    # Broadcast the necessary variables to all ranks
    matrix_size = comm.bcast(matrix_size, root=0)
    gene_rank_dict = comm.bcast(gene_rank_dict, root=0)
    gene2idx = comm.bcast(gene2idx, root=0)

    if rank == 0:
        region_proc_state = False
        np_zero_matrix = None
        local_gene2idx = None
        gene_rank_len = None
        gene_rank_id2rank = None
        local_distribution_temp = None
        local_save_path = None
    else:
        region_proc_state = False
        gene_rank_len = len(gene_rank_dict[rank])
        # Initialize a zero matrix for the count matrix, only for the worker ranks
        np_zero_matrix = np.zeros((gene_rank_len, matrix_size, 5), dtype=int)  # 5 for gene_B_around, gene_around, gene_B, gene_all, gene_A
        local_gene2idx = {gene: idx for idx, gene in enumerate(gene_rank_dict[rank])}
        gene_rank_id2rank = {v: k for k, l in gene_rank_dict.items() for v in l}

        # For distribution saving
        local_distribution_temp = []
        # if opt.r_dist is not None and opt.intermediate_dir is not None:
        if opt.distribution_analysis:
            local_save_path = os.path.join(opt.intermediate_dir, f"rank_{rank}")
        else:
            local_save_path = None

    min_gene_number_local = opt.min_neighbor_number
    expression_level_local = opt.expression_level

    # Task 1: Detecting neighbors in each cell
    if opt.distribution_analysis:
        partial_func = partial(region_function_cell_multi_proc_distribution, r_check=opt.r_check, r_dist=opt.r_dist,
                                 distribution_save_interval=opt.distribution_save_interval,
                               local_distribution_temp=local_distribution_temp, local_save_path=local_save_path,
                                    gene_rank_id2rank=gene_rank_id2rank)
    else:
        if opt.detection_method == 'nine_grid':
            partial_func = partial(region_function_cell_multi_proc_nine_grid, r_check=opt.grid_check,
                                   gene_rank_id2rank=gene_rank_id2rank)
        else:
            # Default to radius method
            partial_func = partial(region_function_cell_multi_proc, r_check=opt.r_check, gene_rank_id2rank=gene_rank_id2rank)

    if rank == 0:
        distribute_tasks_dynamic_region(comm, rank, size, cell_group_list, partial_func, 'NeighborDetection', matrix=np_zero_matrix, gene2idx=gene2idx, local_gene2idx=local_gene2idx)
    else:
        distribute_tasks_dynamic_region(comm, rank, size, cell_group_list, partial_func, 'NeighborDetection', matrix=np_zero_matrix, gene2idx=gene2idx, local_gene2idx=local_gene2idx)

    comm.Barrier()

    # Clear local_distribution_temp which not meet the save condition
    # if opt.r_dist is not None and opt.intermediate_dir is not None:
    if opt.distribution_analysis:
        if rank != 0 and len(local_distribution_temp) > 0:
            local_save_path_file_count = len(os.listdir(local_save_path))
            distribution_save_path = os.path.join(local_save_path, f"distribution_{local_save_path_file_count}.pkl")
            joblib.dump(local_distribution_temp, distribution_save_path, compress=3)
            # clear the temp list after saving
            local_distribution_temp.clear()

    # Wait for all processes to finish the region detection
    comm.Barrier()

    # Gather the count matrices from all processes
    count_matrix_list = large_recv(np_zero_matrix, comm, rank, size, root=0)
    # After gathering, the np_zero_matrix is no longer needed
    np_zero_matrix = None
    if rank == 0:
        print("Root process received count matrix from all work processes. Total matrices received: ", len(count_matrix_list))
        print("Combining count matrices...")
        global_count_matrix = np.concatenate(count_matrix_list, axis=0)
        print("Global count matrix shape: ", global_count_matrix.shape)
        # Memory cleanup
        count_matrix_list = None
    else:
        global_count_matrix = None

    if rank == 0:
        tuple_around = merge_dicts(global_count_matrix, gene2idx)
        print("len of tuple_around: ", len(tuple_around))
        dict_around_list = None
    else:
        tuple_around = None

    # Task 2: Hypergeometric test
    partial_func2 = partial(hyper_test, min_gene_number=min_gene_number_local,
                            expression_level=expression_level_local)

    if rank == 0:
        result_output = distribute_tasks_dynamic(comm, rank, size, tuple_around, partial_func2, 'HyperTest')
    else:
        distribute_tasks_dynamic(comm, rank, size, tuple_around, partial_func2, 'HyperTest')
        result_output = None

    # Wait for all processes to finish the hypergeometric test
    comm.Barrier()

    # if opt.r_dist is not None and opt.intermediate_dir is not None:
    if opt.distribution_analysis:
        # Task 3: Process distribution files
        partial_func3 = partial(distribution_file_proc, gene_rank_id2rank=gene_rank_id2rank)
        for i in range(1, size):
            input_path = os.path.join(opt.intermediate_dir, f"rank_{i}")
            distribution_file_list = os.listdir(input_path)
            distribution_file_list = [os.path.join(input_path, file) for file in distribution_file_list if file.endswith('.pkl')]
            merge_output_path = opt.intermediate_dir + '_gene_id_merge'

            if len(distribution_file_list) > 0:
                if rank == 0:
                    print(f"Distribution file from cell is processing. Path: {input_path}")
                    distribute_tasks_dynamic_multi_asyn(comm, rank, size, distribution_file_list, partial_func3,
                                                        'DistributionFileProc', distribution_proc_asyn,
                                                        merge_output_path)
                    print(f"Distribution file from cell processing completed. Path: {input_path}")
                else:
                    distribute_tasks_dynamic_multi_asyn(comm, rank, size, distribution_file_list, partial_func3,
                                                        'DistributionFileProc', distribution_proc_asyn,
                                                        merge_output_path)

                comm.Barrier()

    if rank == 0:
        result_output_shape_list = []
    else:
        result_output_shape_list = None

    # if opt.r_dist is not None and opt.intermediate_dir is not None:
    if opt.distribution_analysis:
        for i in range(1, size):
            if rank == 0:
                print(f"Processing distribution data for rank {i}...")
                merged_output_path = os.path.join((opt.intermediate_dir + "_gene_id_merge"), f"rank_{i}")
                pkl_list = os.listdir(merged_output_path)

                while len(pkl_list) < (size - 1):
                    print(f"Waiting for more distribution files for rank {i}, current: {len(pkl_list)} / {size - 1}")
                    time.sleep(3)
                    pkl_list = os.listdir(merged_output_path)

                pkl_list = [os.path.join(merged_output_path, file) for file in pkl_list if file.endswith('.pkl')]
                print(f"Found {len(pkl_list)} distribution files to merge for rank {i}.")
                merged_distribution_dict = merge_distribution_dicts_from_file(pkl_list)

                tuple_dist = [(gene_id, gene_distance_dict) for gene_id, gene_distance_dict in
                              merged_distribution_dict.items()]

                merged_distribution_dict = None
            else:
                tuple_dist = None
                merged_output_path = None

            # Task 4: Calculate distribution shape
            partial_func4 = partial(distribution_shape_calculate, around_count_threshold=opt.around_count_threshold)

            if rank == 0:
                print(f"Distribution shape calculation is processing. Path: {merged_output_path}")
                result_output_shape = distribute_tasks_dynamic(comm, rank, size, tuple_dist, partial_func4, 'ShapeCalculate')
                result_output_shape_list.append(result_output_shape)

            else:
                distribute_tasks_dynamic(comm, rank, size, tuple_dist, partial_func4, 'ShapeCalculate')
                result_output_shape = None

            comm.Barrier()


    if rank == 0:
        if opt.distribution_analysis:
            result_output_shape_list = [item for sublist in result_output_shape_list for item in sublist if item is not None]

        # joblib.dump(gene_rank_dict, os.path.join(opt.save_path.replace('.csv', '_gene_rank_dict.pkl')))
        if opt.distribution_analysis:
            joblib.dump(gene_rank_dict, f"{opt.intermediate_dir}/gene_rank_dict.pkl")

        df_result = pd.concat(result_output)
        df_result = add_pair_column(df_result)
        df_result = test_result_df_ratio_proc(df_result)
        df_result.to_csv(opt.save_path, index=False)

        if len(result_output_shape_list) == 0:
            # print("No distribution shape data. Around count threshold may be too high.")
            df_merge = df_result
        else:
            df_shape = pd.concat(result_output_shape_list)
            df_shape = shape_correct(df_shape, opt)
            df_shape.to_csv(opt.save_path.replace('.csv', '_shape.csv'), index=False)

            df_merge = pd.merge(df_result, df_shape, on=['gene_A', 'gene_B'], how='left')

        df_merge_pair_dedup, df_merge_pair_filter = test_result_df_filter(df_merge, filter_threshold=opt.filter_threshold,
                                                                keep=opt.pair_keep)

        if opt.distribution_analysis:
            save_path = opt.save_path.replace('.csv', f'_dedup_{opt.filter_threshold}_post_proc_shape.csv')
        else:
            save_path = opt.save_path.replace('.csv', f'_dedup_{opt.filter_threshold}_post_proc.csv')
        df_merge_pair_filter.to_csv(save_path, index=False)

