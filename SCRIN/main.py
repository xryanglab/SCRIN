import argparse
from SCRIN.hyper_test_glb_base import hyper_test_glb_base
from SCRIN.hyper_test_glb import hyper_test_glb
from SCRIN.hyper_test_glb_distribution import hyper_test_glb_distribution
from SCRIN.hyper_test_glb_nocell import hyper_test_glb_nocell
from SCRIN.hyper_test_clb_distribution import hyper_test_clb_distribution


def main():
    parser = argparse.ArgumentParser(description="SCRIN")
    parser.add_argument('--version', required=True, choices=['v1', 'v2', 'v3', 'v4', 'v5'], help="Select the version")
    parser.add_argument("--detection_method", type=str, choices=['Radius', 'Nine_grid'],
                        default='Radius', help="Method for neighbor detection, can be 'Radius' or 'Nine_grid'")
    parser.add_argument("--r_check", type=float, default=None, help="radius of checking")
    parser.add_argument("--grid_check", type=int, default=None, help="grid size for Nine_grid detection method, default is 1")
    parser.add_argument("--rect_length", type=float, default=20,
                        help="length of the rectangle for splitting the data, recommended value is the cell diameter")
    parser.add_argument("--r_dist", type=float, default=None, help="if not None, enable co-localization distribution saving")
    parser.add_argument("--around_count_threshold", type=int, default=100,
                        help="threshold for the number of points around a gene to consider it for distribution analysis")
    parser.add_argument("--column_name", type=str, default="x,y,z,geneID,cell", help="column name used in data")
    parser.add_argument("--min_gene_number", type=int, default=5,
                        help="minimum number of transcripts for a gene to be considered")
    parser.add_argument("--min_neighbor_number", type=int, default=1,
                        help="minimum number of neighbors for a pair to be considered")
    parser.add_argument("--expression_level", type=float, default=100,
                        help="For gene A and gene B in the pair, the maximum ratio of their expression count.")
    parser.add_argument("--filter_threshold", type=float, default=0.00001,
                        help="filter threshold for qvalue_BH in post processing")
    parser.add_argument("--distribution_save_interval", type=int, default=10,
                        help="interval for saving distribution data to file")
    parser.add_argument("--pair_keep", type=str, default='last', help="keep method for pair post processing, can be 'first' or 'last'")
    parser.add_argument("--data_path", type=str,
                        default="st_data.csv",
                        help="path of data")
    parser.add_argument("--save_path", type=str,
                        default="st_hypertest_result.csv",
                        help="path of result saving")
    parser.add_argument("--molecule_distribution_id_path", type=str,
                        default=None,
                        help="path of molecule distribution saving")
    parser.add_argument("--intermediate_path", type=str,
                        default=None,
                        help="path of intermediate result saving")
    parser.add_argument("--distribution_file_path", type=str, default=None,
                        help="File path to save co-localization distribution data")
    parser.add_argument("--rtree_path", type=str, default=None, help="path of rtree index")
    parser.add_argument("--save_split_size", type=int, default=100, help="interval of intermediate result saving")

    # 其他通用参数
    args = parser.parse_args()

    if args.version == 'v1':
        hyper_test_glb_base(args)
    elif args.version == 'v2':
        hyper_test_glb(args)
    elif args.version == 'v3':
        hyper_test_glb_distribution(args)
    elif args.version == 'v4':
        hyper_test_glb_nocell(args)
    elif args.version == 'v5':
        hyper_test_clb_distribution(args)


if __name__ == '__main__':
    main()