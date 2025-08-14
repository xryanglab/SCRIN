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
mpirun -n 16 SCRIN \
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

- **`--save_path`** : Path for saving the results.

## Output

An example output file (`test.csv`) is shown below:

```
gene_A,gene_B,pvalue,qvalue_BH,qvalue_BO,gene_B_around,gene_B_slice,gene_around,gene_slice,gene_A_N,gene_B_N,pair,enrichment_ratio
Kcnj8,Vtn,0.0,0.0,0.0,566,51096,3603,4885769,2241,51096,Kcnj8_Vtn,17.820363130077027,15.17803020428641,15.020977497702212
Bgn,Vtn,0.0,0.0,0.0,1728,51096,15717,4885769,9904,51096,Bgn_Vtn,12.062024287226444,10.84581394375585,10.512841372618968
Cspg4,Vtn,0.0,0.0,0.0,2074,51096,20569,4885769,12918,51096,Cspg4_Vtn,11.017065217111972,10.007031026811507,9.641433282377363
Cspg4,Pdgfra,0.0,0.0,0.0,658,16515,20569,4885769,12918,16515,Cspg4_Pdgfra,10.106358468927027,9.815047084194957,9.463832976934873
Bgn,Igf2,0.0,0.0,0.0,695,30296,15717,4885769,9904,30296,Bgn_Igf2,7.565480080472584,7.27515694909074,7.131202827107043
Vtn,Igf2,0.0,0.0,0.0,2309,30296,70250,4885769,51096,30296,Igf2_Vtn,5.813629355563021,5.655413409911846,5.300602558199196
Igf2,Timp3,0.0,0.0,0.0,1554,38561,41507,4885769,30296,38561,Igf2_Timp3,5.052599780665756,4.900872600692388,4.743668274521491
Vtn,Rgs5,0.0,0.0,0.0,2969,43823,70250,4885769,51096,43823,Rgs5_Vtn,5.157343774794808,4.981640519743338,4.711885124103652
Flt1,Cldn5,0.0,0.0,0.0,2066,46402,50865,4885769,41692,46402,Cldn5_Flt1,4.574565975489097,4.42937668412253,4.276687312341203
Unc5b,Sox10,5.78605e-319,1.1745686e-316,1.1745686e-316,542,15859,18330,4885769,13112,15859,Sox10_Unc5b,9.65228563955043,9.396446096907969,9.109487664186856
```

## Note

Confirm that all required dependencies are installed before running SCRIN.

## References

(Add references here)

## For any questions, please contact:
Xuerui Yang (yangxuerui@tsinghua.edu.cn); Xu Chen (chenxu22@mails.tsinghua.edu.cn)