import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def test_function(gene_B_around_num, gene_A_N, gene_B_N,
                  gene_slice_num, gene_around_num, minGenenumber, expressionLevel):

    if (gene_B_around_num < minGenenumber or
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
