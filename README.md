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

```
# in bash
echo "Analyze spliced and unspliced reads and genrate a loom file for velocityR implemented in seurat-wrappers"
export PATH=$PATH:/PATH/envs/scRNA_velocyto/bin/
export PYTHONPATH=/PATh/envs/scRNA_velocyto/lib/python3.6/site-packages/
module load samtools

velocyto run10x --samtools-threads $CPU --samtools-memory 500 $TenXfolder $GTF
```

2. Analyze RNA velocity

```
# in R (scRNA_RNA_velocity_run_velocityR.R)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

args <- commandArgs(trailingOnly = TRUE)
fq_dir=args[1]
sample=args[2]
loomfile=paste0(fq_dir, '/veloctyto', sample, '.loom')
#loomfile="/net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/CellLines/Vince_resistence/preprocess/VG_CAMA1_count03/VG_CAMA1/velocyto/VG_CAMA1.loom"
ldat <- ReadVelocity(file =loomfile)
bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
```

3. Project RNA velocity on UMAP/TSNE

```
# in R (scRNA_RNA_velocity_run_velocityR.R)
pdf(paste0(sample, ".scRNA_velocityR.pdf"))
#plot umap labeled by clusters
DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1, group_by = 'ident')
DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1, group_by = 'Phase')

#plot umap with RNA velocity
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm,
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
    do.par = FALSE, cell.border.alpha = 0.1)
dev.off()
```

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
