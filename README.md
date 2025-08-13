# SCRIN: Subcellular Colocalized RNA Interaction Network

Systematic dissections of the subcellular RNA colocalization landscapes in high-resolution spatial transcriptomics.

## Requirements & Compatibility

The following dependencies are required to run this project.  
The versions listed have been tested thoroughly and confirmed to be compatible.  

In most cases, newer versions may also work, as the project relies mainly on stable and widely supported APIs.  
If you encounter issues, we recommend reverting to the specified versions.

> **Note**  
> This project has been tested and is currently supported **only on Linux**.  
> Support for Windows and macOS may be added in the future, but compatibility is not guaranteed at this time.

**Tested Environment:**
- Python 3.9
- Linux

**Required dependencies:**
- mpi4py==4.1.0  
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

Before installing SCRIN, install the dependencies listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

Once dependencies are installed, SCRIN can be installed in **two ways**:

### 1. Install from PyPI
```bash
pip install scrin
```

### 2. Install from local clone
```bash
git clone https://github.com/xryanglab/SCRIN
cd SCRIN
pip install .
```

## Usage

```bash
python scrin.py --save_path /path/to/save
```

## Flags

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