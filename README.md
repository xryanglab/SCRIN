# SCRIN: Subcellular Colocalized RNA Interaction Network

SCRIN is a tool for identifying RNA co-localization networks within subcellular spatial transcriptomics data.

## Requirements & Compatibility

The following dependencies are required to run this project.
The versions listed have been tested thoroughly and confirmed to be compatible.

In most cases, other versions also work, as the project relies mainly on stable and widely supported APIs.
If you encounter issues, we recommend reverting to the specified versions.

> **Note**:
> This project has been tested and is currently supported **only on Linux**.
> Support for Windows and macOS may be added in the future, but compatibility is not guaranteed at this time.

**Tested Environment:**

  - Python 3.9
  - Linux

### System Dependencies

  - **MPICH** (== 4.2.1): Required for parallel computing. SCRIN utilizes high-speed parallel processing to efficiently handle large-scale spatial transcriptomics data.

### Python Dependencies

  - mpi4py==3.1.5
  - msgpack==1.1.1
  - numpy==2.0.2
  - pandas==2.3.1
  - pyarrow==21.0.0
  - rtree==1.4.0
  - scikit-learn==1.6.1
  - scipy==1.13.1
  - statsmodels==0.14.5
  - tools==1.0.2
  - tqdm==4.67.1

## Installation

### Step 1: Set up Python Environment

We recommend using **Anaconda** to manage your environment. Create and activate a new environment:

```bash
conda create -n scrin_env python=3.9
conda activate scrin_env
```

### Step 2: Install MPICH

SCRIN leverages `mpi4py` for high-speed parallel computing to tackle the challenges of large-scale spatial transcriptomics data. This requires a functional MPI (Message Passing Interface) implementation on your system, such as MPICH.

Please install MPICH using one of the following methods before proceeding.

#### Method 1: Install via Conda

The easiest way to ensure compatibility is to let Conda install `mpi4py` and its required MPI implementation (MPICH) together.

```bash
conda install -c conda-forge mpi4py=3.1.5 mpich=4.2.1
```

> **Note on Version Availability**:
> If the command above fails because the specified versions cannot be found for your system, you can try installing without specifying the versions:
>
> ```bash
> conda install -c conda-forge mpi4py mpich
> ```
>
> Please be aware that this will install the latest available packages, which have not been officially tested by us and may lead to unexpected behavior.

#### Method 2: Install via System Package Manager

For Debian-based systems like Ubuntu, you can use `apt`:

```bash
sudo apt update
sudo apt install mpich=4.2.1
```

#### Method 3: Install from Source

For advanced users or specific system configurations, you can compile and install MPICH from the official source. Please refer to the [official MPICH installation guide](https://www.mpich.org/documentation/guides/) for detailed instructions.

### Step 3: Install Python Dependencies

Before installing SCRIN, install the dependencies listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

> **Note**: If you did not use the Conda method to install MPICH in Step 1, `pip` will attempt to compile `mpi4py` using the system's MPI compiler (`mpicc`). Ensure your MPICH installation is correctly configured in your system's PATH.

### Step 4: Install SCRIN

Once dependencies are installed, SCRIN can be installed in **two ways**:

#### 1. Install from PyPI

```bash
pip install scrin
```

#### 2. Install from local clone

```bash
git clone https://github.com/xryanglab/SCRIN
cd SCRIN
pip install .
```

## Usage

The basic command structure to run SCRIN is as follows:

```bash
mpirun -n <number_of_processes> scrin [OPTIONS]
```

### Examples

#### Input Data Format

SCRIN expects a CSV (Comma-Separated Values) file as input. The file should contain columns for spatial coordinates (x, y, and optionally z), a gene identifier, and a cell identifier. The header names can be arbitrary, as they will be mapped using the ``--column_name`` parameter.

Here are the first few lines of the example file, ``Mouse_brain_CosMX_1000cells.csv``:

```csv
x_global_px,y_global_px,z,target,cell
-494295.483333333,7129.73,-1,Prickle2,1_23
-494312.553333333,7164.95,-1,Meg3,1_23
...
```

**Column requirements:**
* **Coordinates (``x``, ``y``, ``z``)**: At least ``x`` and ``y`` columns are required. The ``z`` column is optional.
* **Gene ID (``target`` in the example)**: A column containing the names or identifiers of the RNA species.
* **Cell ID (``cell`` in the example)**: A column indicating which cell each transcript belongs to. This is highly recommended for standard analysis. For data without pre-existing cell segmentation, please refer to the ``Unsegmented Data Options``.

This section provides an example to demonstrate a typical workflow for using SCRIN. We will use a sample dataset derived from the [CosMx SMI Mouse Brain FFPE dataset](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-smi-mouse-brain-ffpe-dataset/) by NanoString. For demonstration purposes, we have randomly sampled 1000 cells from the original data.

**Download the example dataset here:** [https://zenodo.org/records/17019789] 

```bash
# Launch SCRIN on 16 parallel processes. Adjust the value of -n as needed.
mpirun -n 16 scrin \
	--detection_method "radius" \
	--background "cooccurrence" \
	--mode "fast" \
	--data_path "Mouse_brain_CosMX_1000cells.csv" \
	--save_path "Mouse_brain_CosMX_1000cells_hyper_test_cb.csv" \
	--column_name "x_global_px,y_global_px,z,target,cell" \
	--r_check 4.16 \
	--filter_threshold 0.00001 \
	--min_gene_number 5 \
	--min_neighbor_number 1 \
	--expression_level 100 \
	--intermediate_dir "Mouse_brain_CosMX_1000cells_hyper_test_cb"
```

**Explanation of parameters:**

* `--detection_method "radius"`: Since CosMx data provides continuous spatial coordinates, we use the `radius` method to define neighbors based on their straight-line distance.
* `--background "cooccurrence"`: The mouse brain is a highly heterogeneous tissue containing many different cell types. Using `cooccurrence` is recommended here, as it calculates the statistical background for a gene pair (A-B) using only the cells where both A and B are expressed. This provides a more specific and relevant context compared to the `all` option (which would be suitable for more homogeneous samples like single cell type).
* `--mode "fast"`: We use the `fast` mode to enable high-speed, low-memory parallel processing, which is essential for large datasets. This requires an intermediate directory (`--intermediate_dir`) to store temporary files.
* `--column_name`: This parameter maps the column names in our input CSV (`x_global_px`, `y_global_px`, etc.) to the fields SCRIN expects (`x`, `y`, `z`, `geneID`, `cell`).
* `--r_check 4.16`: Sets the search radius. For this dataset, this value corresponds to approximately 0.5 µm.
* `--filter_threshold 0.00001`: Sets the q-value cutoff for the final results, ensuring that only statistically significant interactions are reported.
* `--min_gene_number 5`: A pre-filtering step to improve efficiency by excluding sparsely expressed genes (those with fewer than 5 total transcripts in the dataset) from the analysis.
* `--min_neighbor_number 1`: Skips significance testing for gene pairs with zero observed co-localization events, as they cannot be statistically significant.
* `--expression_level 100`: Filters out gene pairs with highly imbalanced expression levels (where one gene's transcript count is over 100 times that of the other) to avoid potential artifacts.

## Command-line Options

### Mode Options

-   **`--detection_method`** `[radius|nine_grid]` (Required): Method for defining neighboring transcripts.
    -   `radius`: Defines neighbors based on the straight-line distance between transcripts. Any transcript within the distance specified by `--r_check` is considered a neighbor. This is suitable for continuous-coordinate data, such as from MERFISH.
    -   `nine_grid`: Defines neighbors as any transcripts located within the same grid square or its eight adjacent squares. This is suitable for array-based data with orderly coordinates, such as from Stereo-seq.

-   **`--background`** `[all|cooccurrence]` (Required): Define the statistical scope used to calculate the parameters for the hypergeometric test.
    -   `all`: All cells in the dataset are used to calculate the background parameters (`n`, `M`, `N`). For each gene, background parameters are computed only once, which enables more consistent comparison of co-localization strength across gene pairs and provides higher computational efficiency. Recommended for homogeneous data (e.g., single cell lines or types) or when using a global background is needed to find weak co-localization signals.
    -   `cooccurrence`: For a given gene pair A-B, only cells where both A and B are present are used to calculate the background parameters. Recommended for heterogeneous data with mixed or highly specific cell types.
    -   *Note: The value `k` (observed co-localizations) is calculated the same way in both modes, but the background parameters `n`, `M`, and `N` will differ.*

-   **`--mode`** `[robust|fast]` (Required): The running mode for the program.
    -   `fast`: Designed for large-scale spatial transcriptomics datasets, this mode employs complex asynchronous threading to enable low-memory, high-speed parallel processing, but requires higher network bandwidth for inter-process communication.
    -   `robust`: A more stable running mode but requires higher memory. Can be used for simple tests or if issues are encountered with the `fast` mode.

### Base Parameters

-   **`--data_path`** `[str]` (Required): Path to the input data file. The file must contain transcript spatial coordinates and gene IDs at a minimum. Including cell IDs is recommended. Please refer to `Mouse_brain_CosMX_1000cells.csv` for the standard input format.
-   **`--save_path`** `[str]` (Required): Path for saving the results. 
-   **`--column_name`** `[str]` (Required): A comma-separated string specifying which columns from the input file to use. The provided names are mapped sequentially to the expected fields: `x` (x-coordinate), `y` (y-coordinate), `z` (z-coordinate, optional), `geneID` (gene ID), and `cell` (cell ID, optional). If an optional field like `z` is not present in your data, simply omit it from the string while maintaining the order of the remaining fields. For example, if your file provides columns for x, y, geneID, and cell (but no z), and their names are `pos_x, pos_y, gene_name, cell_label`, your input should be `"pos_x,pos_y,gene_name,cell_label"`. The minimum required fields correspond to `x`, `y`, and `geneID`. Default: `"x,y,z,geneID,cell"`.
-   **`--r_check`** `[float]`: The search radius for the `'radius'` detection method. Transcripts with a distance between them less than this value are considered neighbors.
-   **`--grid_check`** `[int]`: Sets the search window size for the `'nine_grid'` method. It defines a square area of `(2 * grid_check + 1) x (2 * grid_check + 1)` grid cells around a central transcript. For example, `grid_check=1` defines a 3x3 grid (9 cells), while `grid_check=2` defines a 5x5 grid (25 cells). Transcripts within this area are considered neighbors.
-   **`--min_gene_number`** `[int]`: A pre-filtering step to remove sparsely expressed genes. Any gene whose total transcript count across the entire dataset is below this value will be excluded from the analysis. Default: `5`.
-   **`--min_neighbor_number`** `[int]`: Filters out gene pairs with insufficient co-localization events. For a given pair A-B, if the number of times transcripts of gene B are detected as neighbors of transcripts of gene A is below this threshold, that pair will be skipped during the significance calculation. Default: `1`.
-   **`--expression_level`** `[float]`: A filter to exclude gene pairs with highly imbalanced expression. This value sets the maximum allowable fold-difference in total transcript counts between two genes. For example, with the default of `100`, any pair where one gene is over 100 times more abundant than the other will be ignored. Default: `100`.
-   **`--filter_threshold`** `[float]`: The q-value (Benjamini-Hochberg adjusted p-value) threshold for filtering results in post-processing. Default: `0.00001`.
-   **`--pair_keep`** `[first|last]`: Method for deduplicating bidirectional pairs (e.g., A-B and B-A) during post-processing. Pairs are first sorted by their q-value in ascending order. `first` keeps the pair with the smaller q-value, while `last` keeps the one with the larger q-value. Default: `'last'`.

### Intermediate Result Options
For large datasets, use these options to save intermediate results and prevent memory overflow.

-   **`--intermediate_dir`** `[str]`: Directory path to save intermediate results. This parameter is required when using `fast` mode.
-   **`--intermediate_split`** `[int]`: Controls the chunk size for processing. A larger value reduces memory usage but may decrease computational efficiency. It is not recommended to set this value higher than the total number of genes or `1000`, as excessive partitioning can lead to issues. Default: `100`.

### Distribution Options
Options for analyzing the distance distribution of co-localized gene pairs.

-   **`--distribution_analysis`**: A flag to enable the analysis. This will save the distance distribution for each neighboring pair and calculate its statistical features. **Warning:** This can generate very large files and significantly increase runtime. Ensure you have sufficient disk space before enabling.
-   **`--r_dist`** `[float]`: Defines the maximum radius for the distance distribution analysis. For a pair A-B, all observed distances between their transcripts that are less than this value will be recorded.
-   **`--around_count_threshold`** `[int]`: A filter to ensure the statistical reliability of the distance distribution. For a gene pair, the analysis is performed only if the total number of observed co-localization events (i.e., distances less than `--r_dist`) exceeds this threshold. This prevents analyzing pairs with too few data points to be meaningful. Default: `100`.
-   **`--distribution_save_interval`** `[int]`: Controls how often collected distance data is written to intermediate files to manage memory. A smaller value decreases memory usage. For whole transcriptome datasets, a value no higher than `100` is recommended. Default: `10`.

### Unsegmented Data Options
For data without prior cell segmentation.

-   **`--unsegmented`**: A flag to enable processing of unsegmented data.
-   **`--rect_length`** `[float]`: The side length of the rectangle used to partition the data. The recommended value is the approximate cell diameter. Default: `20`.
-   **`--rtree_path`** `[str]`: Path to an R-tree index file for accelerating spatial queries. If the file does not exist, a new index will be built and saved to this path. If the file already exists, it will be loaded to save time.

## Output

SCRIN generates two main types of output files: a raw results file and a post-processed file for downstream analysis.

### Final Processed Output

This is the primary result file you will typically use. It is generated by performing several post-processing steps on the raw output. For the example command, this file would be named `Mouse_brain_CosMX_1000cells_hyper_test_cb_dedup_1e-05_post_proc.csv`.

The post-processing includes:
1.  **Adding a `pair` column**: A standardized, sorted identifier for each gene pair (e.g., `GeneA_GeneB`) is added to facilitate deduplication.
2.  **Calculating the Enrichment Ratio**: An enrichment ratio is calculated to prioritize pairs with a higher degree of enrichment.
3.  **Sorting**: The results are sorted by the Benjamini-Hochberg adjusted q-value (`qvalue_BH`) in ascending order.
4.  **Deduplication**: Bidirectional pairs (e.g., A-B and B-A) are deduplicated based on the `--pair_keep` parameter.
5.  **Filtering**: The results are filtered to keep only the interactions with a `qvalue_BH` below the `--filter_threshold`.

An example snippet from the final output file is shown below:

```
gene_A,gene_B,pvalue,qvalue_BH,qvalue_BO,gene_B_around,gene_B_slice,gene_around,gene_slice,gene_A_N,gene_B_N,pair,enrichment_ratio
Scd2,Plp1,9.545615285730711e-303,8.447869527871679e-300,8.44786952787168e-300,1175,13953,19063,803815,6850,13953,Plp1_Scd2,3.550872927582488
Meg3,Malat1,5.393831772255401e-178,4.967719062247224e-175,4.967719062247224e-175,2561,68433,26756,1287763,10284,68433,Malat1_Meg3,1.8011867965562922
...
```

### Raw Output

The program also saves the raw, unfiltered results from the hypergeometric test. Using the example command, this file would be `Mouse_brain_CosMX_1000cells_hyper_test_cb.csv`. This file contains all calculated gene pairs before any sorting, deduplication, or q-value filtering, and can be useful for custom analyses or diagnostics.

### Column Descriptions

The columns in the output files correspond to the standard parameters of a hypergeometric test and additional calculated metrics.

| Column           | Description                                                                                             |
|------------------|---------------------------------------------------------------------------------------------------------|
| `gene_A`         | The **central gene** in the pair.                                                                       |
| `gene_B`         | The **target gene** found in the neighborhood of the central gene (`gene_A`).                           |
| `pvalue`         | The raw p-value calculated from the hypergeometric test.                                                |
| `qvalue_BH`      | The p-value adjusted for multiple testing using the Benjamini-Hochberg (BH) method. This is the primary value for assessing significance. |
| `qvalue_BO`      | The p-value adjusted using the more conservative Bonferroni method.                                     |
| `gene_B_around`  | **(k)** The number of observed co-localizations: how many times transcripts of the target gene (`gene_B`) were found in the neighborhood of the central gene (`gene_A`). |
| `gene_B_slice`   | **(M)** The total number of target gene (`gene_B`) transcripts within the defined **statistical background**.     |
| `gene_around`    | **(n)** The total number of all transcripts (from any gene) found in the neighborhood of the central gene (`gene_A`). |
| `gene_slice`     | **(N)** The total number of all transcripts in the defined **statistical background**.              |
| `gene_A_N`       | The total transcript count for the central gene (`gene_A`) in the defined **statistical background**.                   |
| `gene_B_N`       | The total transcript count for the target gene (`gene_B`) in the defined **statistical background**; same as `gene_B_slice`.                    |
| `pair`           | A standardized identifier for the gene pair (e.g., alphabetically sorted), used for post-processing.    |
| `enrichment_ratio`| **($\delta$)** A metric to prioritize pairs that are statistically significant but differ in their degree of enrichment. It is calculated as $\delta = \frac{kN}{nM}$ |

## Citation

If you use SCRIN in your research, please cite our publication.

## For any questions, please contact:

Xuerui Yang (yangxuerui@tsinghua.edu.cn); Xu Chen (chenxu22@mails.tsinghua.edu.cn)