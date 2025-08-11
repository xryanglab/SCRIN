import pandas as pd


def merge_columns(row):
    values = [row['gene_A'], row['gene_B']]
    values.sort()
    return '_'.join(values)


def add_pair_column(df):
    if 'pair' in df.columns:
        return df

    pair = df.apply(merge_columns, axis=1)
    df = df.assign(pair=pair)

    return df


def test_result_df_filter(file, qvalue_column='qvalue_BH', filter_threshold=0.01, keep='last'):
    if type(file) == str:
        df = pd.read_csv(file)
    elif type(file) == pd.DataFrame:
        df = file
    else:
        raise ValueError("file must be a path or a DataFrame")

    # df = df.sort_values(by=qvalue_column)
    df = df.sort_values(by=[qvalue_column, 'enrichment_ratio'], ascending=[True, False])

    filtered_df = df.drop_duplicates(subset='pair', keep=keep)
    filtered_df = filtered_df.reset_index(drop=True)

    filtered_df_bh = filtered_df[filtered_df[qvalue_column] < filter_threshold]

    return filtered_df, filtered_df_bh


def test_result_df_filter_bi(file, qvalue_column='qvalue_BH', filter_threshold=0.01):
    if type(file) == str:
        df = pd.read_csv(file)
    elif type(file) == pd.DataFrame:
        df = file
    else:
        raise ValueError("file must be a path or a DataFrame")

    df = df.sort_values(by=qvalue_column)
    df = df.sort_values(by=[qvalue_column, 'enrichment_ratio'], ascending=[True, False])

    filtered_df = df
    pair_retain_list = filtered_df[filtered_df[qvalue_column] < filter_threshold]['pair'].tolist()

    filtered_df_bh = filtered_df[filtered_df['pair'].isin(pair_retain_list)]
    filtered_df_bh = filtered_df_bh.reset_index(drop=True)

    return filtered_df, filtered_df_bh


def test_result_df_ratio_proc(file):
    if type(file) == str:
        df = pd.read_csv(file)
    elif type(file) == pd.DataFrame:
        df = file
    else:
        raise ValueError("file must be a path or a DataFrame")

    df_output = df.copy()

    # df_output['around_sample_ratio'] = df['gene_around'] / df['gene_slice']

    #  k: gene_B_around, n: gene_around, N: gene_slice, K: gene_B_slice
    #  odds_ratio = (k * (N - n - K + k)) / (n - k) * (K - k)
    # df_output['odds_ratio'] = (df['gene_B_around'] * (
    #         df['gene_slice'] - df['gene_around'] - df['gene_B_slice'] + df['gene_B_around'])) / (
    #                            (df['gene_around'] - df['gene_B_around']) * (df['gene_B_slice'] - df['gene_B_around']))

    #  risk_ratio = k * (N - n) / n * (K - k)
    # df_output['risk_ratio'] = (df['gene_B_around'] * (df['gene_slice'] - df['gene_around'])) / (
    #         df['gene_around'] * (df['gene_B_slice'] - df['gene_B_around']))

    df_output['enrichment_ratio'] = (df['gene_B_around'] / df['gene_around']) / (df['gene_B_slice'] / df['gene_slice'])

    return df_output

