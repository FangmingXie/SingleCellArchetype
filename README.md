# SingleCellArchetype
## Table of Contents
1. [Project overview](#Project-overview)
2. [Repository structure](#Repository-structure)
3. [Getting started](#Getting-started)
4. [Setting up locally](#Setting-up-locally)
3. [Cite](#Cite)

## Project overview
This repository contains code for Archetypal Analysis of single-cell RNA-seq data as described in [Xie et al. 2024](https://www.biorxiv.org/content/10.1101/2023.12.18.572244v2), an application and extension of Pareto multi-tasking theory as described in [Adler et al. 2019](https://doi.org/10.1016/j.celrep.2023.112412). 

The code base is a wrapper and extension of the `ulfaslak/py_pcha` package ([link](https://github.com/ulfaslak/py_pcha)) that implements the PCHA algorithm.

![img](https://github.com/FangmingXie/SingleCellArchetype/blob/main/data/img.png)
*Figure: Archetypal analysis of the transcriptomic continuum of L2/3 excitatory neurons (Xie et al. 2024).*

## Repository structure
```
SingleCellArchetype/
├── README.md
├── data
│   └── data_snrna_v1.h5ad        # sample data
└── sca
    ├── utils.py                  # utility functions
    ├── sca.py                    # the SingleCellArchetype class
    ├── tutorial_complete.ipynb   # a short tutorial (~1 min)
    └── tutorial_minimum.ipynb    # a long tutorial  (~10 min)
```

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

## Cite
- Fangming Xie, Saumya Jain, Runzhe Xu et al. 2024 *bioRxiv* "**Spatial profiling of the interplay between cell type- and vision-dependent transcriptomic programs in the visual cortex**" ([link](https://www.biorxiv.org/content/10.1101/2023.12.18.572244v2))
- Miri Adler et al 2019 *Cell Reports* "**Emergence of division of labor in tissues through cell interactions and spatial cues**" ([link](https://doi.org/10.1016/j.celrep.2023.112412))
- The `py_pcha` package: https://github.com/ulfaslak/py_pcha
