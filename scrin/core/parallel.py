import time


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
