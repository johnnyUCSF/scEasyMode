# ScEasyMode: Wrappers for automating single cell workflows in Python

## Features
- Multiseq correction in python using a within-barcode zscore correction
- Plotting for stacked barplots in your dataset
- Mouse cell filtering/separation from mixed dataset
- Scanpy wrapper that simplifies the workflow

## Installation
Install the package using the Github extension of Pip:
```sh
pip3 install git+https://github.com/johnnyUCSF/scEasyMode#egg='scEasyMode'
```

## Usage
### Load the modules

```sh
from scEasyMode import mousefilter
from scEasyMode import clusterplot
from scEasyMode import pymulti
from scEasyMode import sceasy
```

### Demultiplex your samples

```python
import pandas as pd
from scEasyMode import pymulti

# Define parameters

len_10x=16 # Number of bases in cell barcodes
len_umi=12 # Length of UMI
len_multi=15, # Number of bases in the HTO barcodes / HashTag O

fastq_r1 = 'path/to/file'
fastq_r2 = 'path/to/file'
sample_name = 'test_demultiplexing'

cell_BC_file = 'path/to/cell_barcodes' # Counts Matrix after alignment and pre-processing
cell_bcs = pd.read_csv(cell_BC_file, sep='\t', header=None)[0].tolist()

multi_BC_file = 'path/to/barcodes' # Barcodes TSV file from 10x or Illumina
bcsmulti = pd.read_csv(multi_BC_file,sep=',',index_col=1,header=None) 
bcsmulti.columns = ['multi']
bcsmulti = bcsmulti['multi'].tolist()

pymulti.pymulti(fastq_r1, fastq_r2, bcsmulti=bcsmulti, bcs10x=cell_bcs, len_10x=len_10x, len_multi=len_multi, 
                len_umi=len_umi, split=True, hamming=True, median_only=True, sampname= sample_name)
```
