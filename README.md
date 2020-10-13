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

```
pymulti.pymulti(R1,R2,bcsmulti,bcs10x,split=True)
```
