# SCRIN: Subcellular Colocalized RNA Interaction Network

Systematic dissections of the subcellular RNA colocalization landscapes in high-resolution spatial transcriptomics.

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

```bash
mpirun -n 16 scrin \
	--detection_method "radius" \
	--background "cooccurrence" \
	--mode "fast" \
	--data_path "Run1000_S1_Half_tx_file_cell_1000cells.csv" \
	--save_path "halfbrain_1000cells_hyper_test_cb.csv" \
	--column_name "x_global_px,y_global_px,z,target,cell" \
	--r_check 4.16 \
	--filter_threshold 0.00001 \
	--min_gene_number 6 \
	--min_neighbor_number 1 \
	--expression_level 100 \
	--intermediate_dir "halfbrain_1000cells_hyper_test_cb"
```

## Command-line Options

### Mode Options

-   **`--detection_method`** `[radius|nine_grid]` (Required): Method for defining neighboring transcripts.
    -   `radius`: Defines neighbors based on the straight-line distance between transcripts. Any transcript within the distance specified by `--r_check` is considered a neighbor. This is suitable for continuous-coordinate data, such as from MERFISH.
    -   `nine_grid`: Defines neighbors as any transcripts located within the same grid square or its eight adjacent squares. This is suitable for array-based data with orderly coordinates, such as from Stereo-seq.

-   **`--background`** `[all|cooccurrence]` (Required): Define the statistical scope used to calculate the parameters for the hypergeometric test.
    -   `all`: All cells in the dataset are used to calculate the background parameters (`n`, `M`, `N`). Recommended for homogeneous data (e.g., single cell lines or types) or when using a global background is necessary to find weak co-localization signals.
    -   `cooccurrence`: For a given gene pair A-B, only cells where both A and B are present are used to calculate the background parameters. Recommended for heterogeneous data with mixed or highly specific cell types.
    -   *Note: The value `k` (observed co-localizations) is calculated the same way in both modes, but the background parameters `n`, `M`, and `N` will differ.*

-   **`--mode`** `[robust|fast]` (Required): The running mode for the program.
    -   `fast`: Designed for large-scale spatial transcriptomics datasets, this mode employs complex asynchronous threading to enable low-memory, high-speed parallel processing, but requires higher network bandwidth for inter-process communication.
    -   `robust`: A more stable running mode but requires higher memory. Can be used for simple tests or if issues are encountered with the `fast` mode.

### Base Parameters

-   **`--data_path`** `[str]`: Path to the input data file. Default: `"st_data.csv"`.
-   **`--save_path`** `[str]`: Path for saving the results. Default: `"st_hypertest_result.csv"`.
-   **`--column_name`** `[str]`: A comma-separated string of column names used in the data. Default: `"x,y,z,geneID,cell"`.
-   **`--r_check`** `[float]`: The radius for neighbor searching when `detection_method` is set to `'radius'`. Default: `None`.
-   **`--grid_check`** `[int]`: The grid size for the `'nine_grid'` detection method. Default: `None`.
-   **`--min_gene_number`** `[int]`: Minimum number of transcripts for a gene to be included in the analysis. Default: `5`.
-   **`--min_neighbor_number`** `[int]`: Minimum number of neighbors for a gene pair to be considered. Default: `1`.
-   **`--expression_level`** `[float]`: The maximum allowed ratio of expression counts between the two genes in a pair. Default: `100`.
-   **`--filter_threshold`** `[float]`: The q-value (Benjamini-Hochberg adjusted p-value) threshold for filtering results in post-processing. Default: `0.00001`.
-   **`--pair_keep`** `[first|last]`: Method for handling duplicate pairs during post-processing. Default: `'last'`.

### Intermediate Result Options
For large datasets, use these options to save intermediate results and prevent memory overflow.

-   **`--intermediate_dir`** `[str]`: Directory path to save intermediate results. Enables this feature if set. Default: `None`.
-   **`--intermediate_split`** `[int]`: The interval at which to save intermediate results. Default: `100`.

### Distribution Options
For co-localization distribution analysis.

-   **`--distribution_analysis`**: A flag to enable co-localization distribution analysis.
-   **`--r_dist`** `[float]`: The radius for co-localization distribution analysis. If set, distribution data will be saved. Default: `None`.
-   **`--around_count_threshold`** `[int]`: The threshold for the number of points around a gene to consider it for distribution analysis. Default: `100`.
-   **`--distribution_save_interval`** `[int]`: The interval for saving distribution data to a file. Default: `10`.

### Unsegmented Data Options
For data without prior cell segmentation.

-   **`--unsegmented`**: A flag to enable processing of unsegmented data.
-   **`--rect_length`** `[float]`: The side length of the rectangle used to partition the data. The recommended value is the approximate cell diameter. Default: `20`.
-   **`--rtree_path`** `[str]`: Path to an R-tree index file to accelerate spatial queries. Can be used to load or save an index. Default: `None`.

## Output

An example output file (`test.csv`) is shown below:

```
gene_A,gene_B,pvalue,qvalue_BH,qvalue_BO,gene_B_around,gene_B_slice,gene_around,gene_slice,gene_A_N,gene_B_N,pair,enrichment_ratio
Kcnj8,Vtn,0.0,0.0,0.0,566,51096,3603,4885769,2241,51096,Kcnj8_Vtn,17.820363130077027,15.17803020428641,15.020977497702212
Bgn,Vtn,0.0,0.0,0.0,1728,51096,15717,4885769,9904,51096,Bgn_Vtn,12.062024287226444,10.84581394375585,10.512841372618968
...
```

## Note

Confirm that all required dependencies are installed before running SCRIN.

## References

(Add references here)

## For any questions, please contact:
Xuerui Yang (yangxuerui@tsinghua.edu.cn); Xu Chen (chenxu22@mails.tsinghua.edu.cn)