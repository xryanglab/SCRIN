import pickle
from collections import defaultdict


def merge_dicts(list_of_dicts, output_mode='list'):
    merged_list = []
    merged_dict = defaultdict(lambda: [0, defaultdict(int)])

    for gene_data in list_of_dicts:
        for gene_id, (around_num, sub_dict) in gene_data.items():
            merged_dict[gene_id][0] += around_num
            for sub_gene_id, sub_num in sub_dict.items():
                merged_dict[gene_id][1][sub_gene_id] += sub_num

    if output_mode != 'list':
        return merged_dict

    for gene_id, (around_num, sub_dict) in merged_dict.items():
        merged_list.append((gene_id, around_num, dict(sub_dict)))

    return merged_list


def update_pkl_file(task_result_sub, gene_id_list_split, opt, task_tag, monitor_list=None):
    merge_results = merge_dicts(task_result_sub, output_mode='dict')
    merge_gene_id_list = list(merge_results.keys())

    for i, gene_id_list in enumerate(gene_id_list_split):
        common_gene_id = list(set(gene_id_list) & set(merge_gene_id_list))
        if len(common_gene_id) != 0:
            with open(f'{opt.intermediate_dir}/{task_tag}/task_results{i}.pkl', 'rb') as f:
                gene_dict_pkl = pickle.load(f)
                for gene_id in common_gene_id:
                    gene_dict_pkl[gene_id][0] += merge_results[gene_id][0]
                    for sub_gene_id, sub_num in merge_results[gene_id][1].items():
                        gene_dict_pkl[gene_id][1][sub_gene_id] += sub_num
                with open(f'{opt.intermediate_dir}/{task_tag}/task_results{i}.pkl', 'wb') as f:
                    pickle.dump(gene_dict_pkl, f)

    if monitor_list is not None:
        monitor_list.pop(0)
