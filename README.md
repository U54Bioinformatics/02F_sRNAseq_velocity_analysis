# RNA velocity of single cells
Analyze RNA velocity of 10x genomics scRNA data.

## Pre-requisites and installation 
Install velocyto and seurat-wrappers with Conda

```
# install velocyto
export PYTHONPATH=/PATH/envs/scRNA_velocyto/lib/python3.6/site-packages/
conda env create -f velocyto.yml -n scRNA_velocyto
# install additional tools for velocyto
conda install -n scRNA_velocyto -c conda-forge python-igraph
conda install -n scRNA_velocyto -c conda-forge numba pytables louvain

# install seurat wrapper
export R_LIBS=/PATH/envs/scRNA_seurat_wrapper/lib/R/library
conda create -n scRNA_seurat_wrapper -c r -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.6
conda install -n scRNA_seurat_wrapper -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr bioconductor-batchelor
# in R
/PATH/envs/scRNA_seurat_wrapper/bin/R
install.packages("devtools")
library(devtools)
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
install.packages("dplyr")
devtools::install_github('satijalab/seurat-wrappers')
```

## Workflow
1. Count spliced and unspliced reads
2. Analyze RNA velocity
3. Project RNA velocity on UMAP/TSNE

### Options/Arguments
Add parameters different from default

## Updates   
Add dataset-specific updates 

### Cluster parameters and runtime  
Cluster: N/A
Desktop: 1-2 hours 

## References and documentation 
Reference: [Manno et al. Nature, 2018](https://www.nature.com/articles/s41586-018-0414-6)

Vignette: [Velocyto](http://velocyto.org/velocyto.py/tutorial/index.html), [Estimating RNA Velocity using Seurat](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html)
