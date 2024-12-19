# SingleCellArchetype
A wrapper of the `ulfaslak/py_pcha` python package that implements the PCHA algorithm for Archetypal Analysis by Mørup et. al.

## Table of Contents
 1. [Getting started](#Getting-started)
 2. [Setting up locally](#Setting-up-locally)

## Getting started
Go directly [here](https://github.com/FangmingXie/SingleCellArchetype/blob/main/sca/tutorial_minimum.ipynb) or
[![Open in colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/FangmingXie/SingleCellArchetype/blob/main/sca/tutorial_minimum.ipynb) for a short tutorial. You can run through this example application of Archetypal Analysis in your web browser in ~1 minute. No need to set up anything else.

Alternatively, you can also check out [this](https://github.com/FangmingXie/SingleCellArchetype/blob/main/sca/tutorial_complete.ipynb) complete tutorial.

## Setting up locally
### Step 1. get the data
Download the sample data [here](https://raw.githubusercontent.com/FangmingXie/SingleCellArchetype/main/data/data_snrna_v1.h5ad)
or with the following command.
```
wget 'https://raw.githubusercontent.com/FangmingXie/SingleCellArchetype/main/data/data_snrna_v1.h5ad'
```

### Step 2. get the code
```
pip install anndata
pip install py_pcha
git clone git@github.com:FangmingXie/SingleCellArchetype.git
```

### Step 3. follow the tutorials
Run through the tutorials in your own jupyter notebook or jupyter lab.
- `sca/tutorial_minimum.ipynb`
- `sca/tutorial_complete.ipynb`
