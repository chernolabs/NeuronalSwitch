# NeuronalSwitch

Computational scripts for the paper: *Transcriptional dynamics
orchestrating the development and integration of neurons born in the
adult hippocampus*


#### Install

```
    devtools::install_github("chernolabs/NeuronalSwitch")
```

Note: the following `R` packages should also be installed in your system:

`SingleCellExperiment`, `BiocParallel`, `scater`, `scran`, `robustbase`, `batchelor`, `slingshot`,`Seurat`,`igraph`, `ggplot`, `patchwork`,`dendsort`,`pheatmap`,`ComplexHeatmap`




#### Get raw data from GEO
In order to run this
computational pipeline you should download `01_DS1x_raw_ssce.Rdata` and
`02_DS2_raw_ssce.Rdata` from `GEOXXXX` and put them in the `data_raw`
folder.

#### Pipeline

Directory structure
  + `data_raw` : initial raw data (you should put the count matrices downloaded from GEO here) 
  + `data_aux` : provided pre-calculated intermediate steps to alleviate processing time
  + `data`     : working directory 
  + `data_figs`: output figures of the paper

You should sequentially apply the following R scripts.

| script                   |                                          |
|:-------------------------|:-----------------------------------------|
| 01_sce_DS1x.R            | Preprocessing and QC of DS1 raw data     |
| 02_sce_DS2.R             | Preprocessing and QC of DS2 raw data     |
| 03_DS1x_DS2_comparison.R | Data integration, replicability analysis |
| 04_sce_DS1.R             | De-novo reprocessing and QC of dataset1  |
| 05_DS2_labelTransfer.R   | Seurat cluster label's transfer          |
| 06_DS1_MKNN_fig.R        | MKNN graph cluster and density figures   |
| 07_DS1_pseudotemp.R      | Pseudotime estimation for DS1            |
| 07_DS2_pseudotemp.R      | Pseudotime estimation for DS2            |
| 08_density_pseudot.R     | Pseudotime density plots                 |
| 09_cascades_TFs.R        | Heatmap of TF expression profiles        |
