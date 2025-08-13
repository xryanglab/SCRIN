import argparse
from SCRIN.hyper_test_glb_base import hyper_test_glb_base
from SCRIN.hyper_test_glb import hyper_test_glb
from SCRIN.hyper_test_glb_distribution import hyper_test_glb_distribution
from SCRIN.hyper_test_glb_nocell import hyper_test_glb_nocell
from SCRIN.hyper_test_clb_distribution import hyper_test_clb_distribution


def main():
    parser = argparse.ArgumentParser(description="SCRIN")

    mode_group = parser.add_argument_group("Mode Options")
    mode_group.add_argument("--detection_method", type=str, choices=['radius', 'nine_grid'],
                        required=True, help="Method for neighbor detection, can be 'radius' or 'nine_grid'")
    mode_group.add_argument("--background", type=str, choices=['all', 'cooccurrence'],
                        required=True, help="Background calculation method, all cell background or gene a/b co-occurrence cell background")
    mode_group.add_argument("--mode", type=str, choices=['robust', 'fast'],
                        required=True, help="Method for neighbor detection, can be 'radius' or 'nine_grid'")

    base_params_group = parser.add_argument_group("Base Parameters")
    base_params_group.add_argument("--data_path", type=str,
                        default="st_data.csv",
                        help="path of data")
    base_params_group.add_argument("--save_path", type=str,
                        default="st_hypertest_result.csv",
                        help="path of result saving")
    base_params_group.add_argument("--column_name", type=str, default="x,y,z,geneID,cell", help="column name used in data")
    base_params_group.add_argument("--r_check", type=float, default=None, help="radius of checking")
    base_params_group.add_argument("--grid_check", type=int, default=None, help="grid size for nine_grid detection method, default is 1")
    base_params_group.add_argument("--min_gene_number", type=int, default=5,
                        help="minimum number of transcripts for a gene to be considered")
    base_params_group.add_argument("--min_neighbor_number", type=int, default=1,
                        help="minimum number of neighbors for a pair to be considered")
    base_params_group.add_argument("--expression_level", type=float, default=100,
                        help="For gene A and gene B in the pair, the maximum ratio of their expression count.")
    base_params_group.add_argument("--filter_threshold", type=float, default=0.00001,
                        help="filter threshold for qvalue_BH in post processing")
    base_params_group.add_argument("--pair_keep", type=str, default='last', help="keep method for pair post processing, can be 'first' or 'last'")

    intermediate_group = parser.add_argument_group("Intermediate Result Options: for large datasets, save intermediate results to avoid memory overflow")
    intermediate_group.add_argument("--intermediate_dir", type=str,
                        default=None,
                        help="path of intermediate result saving")
    intermediate_group.add_argument("--intermediate_split", type=int, default=100, help="interval of intermediate result saving")

    distribution_group = parser.add_argument_group("Distribution Options: for co-localization distribution analysis")
    distribution_group.add_argument("--distribution_analysis", action='store_true',
                        help="enable co-localization distribution analysis")
    distribution_group.add_argument("--r_dist", type=float, default=None, help="if not None, enable co-localization distribution saving")
    distribution_group.add_argument("--around_count_threshold", type=int, default=100,
                        help="threshold for the number of points around a gene to consider it for distribution analysis")
    distribution_group.add_argument("--distribution_save_interval", type=int, default=10,
                        help="interval for saving distribution data to file")

    unsegmented_group = parser.add_argument_group("Unsegmented Options: for data without cell segmentation")
    unsegmented_group.add_argument("--unsegmented", action='store_true',
                        help="enable unsegmented data processing")
    unsegmented_group.add_argument("--rect_length", type=float, default=20,
                        help="length of the rectangle for splitting the data, recommended value is the cell diameter")
    unsegmented_group.add_argument("--rtree_path", type=str, default=None, help="path of rtree index")

    args = parser.parse_args()

    function_run = None

    # Check
    if args.unsegmented:
        if args.mode == "robust":
            raise ValueError("If unsegmented data is used, only 'fast' mode is allowed.")
        if args.background == "cooccurrence":
            raise ValueError("Co-occurrence background is not supported for unsegmented data.")
        if args.distribution_analysis:
            raise ValueError("Distribution analysis is not supported for unsegmented data.")
    else:
        if args.mode == "robust" and args.distribution_analysis:
            raise ValueError("Distribution analysis is only allowed in 'fast' mode.")

    # Additional checks
    if args.mode == "fast":
        if args.intermediate_dir is None:
            raise ValueError("Intermediate directory must be specified in 'fast' mode. "
                             "Large datasets may cause memory overflow without it. Please set --intermediate_dir.")

    if args.detection_method == 'radius':
        if args.r_check is None:
            raise ValueError("Detection method 'radius' requires --r_check to be set.")

    if args.detection_method == 'nine_grid':
        if args.grid_check is None:
            raise ValueError("Detection method 'nine_grid' requires --grid_check to be set.")

    if args.distribution_analysis and args.background == "all":
        if args.r_dist is None:
            raise ValueError("Co-localization distribution analysis requires --r_dist to be set.")
        if args.detection_method == "nine_grid":
            raise ValueError("Nine_grid detection method is not compatible with co-localization distribution analysis. "
                             "Please use 'radius' detection method or close co-localization distribution analysis.")

    if args.distribution_analysis and args.background == "cooccurrence":
        if args.r_dist is None:
            raise ValueError("Co-localization distribution analysis requires --r_dist to be set.")
        if args.detection_method == "nine_grid":
            raise ValueError("Nine_grid detection method is not compatible with co-localization distribution analysis. "
                             "Please use 'radius' detection method or close co-localization distribution analysis.")

    # Select function to run
    func_map = {
        (True, "fast", "all", False): hyper_test_glb_nocell,
        (False, "robust", "all", False): hyper_test_glb_base,
        (False, "robust", "cooccurrence", False): hyper_test_clb_distribution,
        (False, "fast", "all", False): hyper_test_glb,
        (False, "fast", "all", True): hyper_test_glb_distribution,
        (False, "fast", "cooccurrence", False): hyper_test_clb_distribution,
        (False, "fast", "cooccurrence", True): hyper_test_clb_distribution,
    }

    func_key = (args.unsegmented, args.mode, args.background, args.distribution_analysis)
    function_run = func_map[func_key]


if __name__ == '__main__':
    main()