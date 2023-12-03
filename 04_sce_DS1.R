# Process de-novo Dataset1 discarding non-reproducible nuclei detected in DS1x
# Generate the k=40 MKNN graph
# Louvain clustering + refinement 
#
# input:
#   data_raw/01_DS1x_raw_ssce.Rdata   :
#   data_aux/01_DS1x_DblFinder.Rdata  : pre-calculate doublet detection results
#   data_aux/04_DS1_clusterIDs.Rdata  : cluster annotation data
#output:
#   data/04_DS1.ssce.MKNN40.Rdata : igraph K=40 MKNN graph
#   data/04_DS1.ssce.Rdata        : QC filtered, logNormalized dataset1, clusters 
#   data/04_DS1.MKNN40.adjlist.csv: adjacency list of K=40 MKNN graph (to import into Gephi)
#   data/04_DS1.MKNN40.nodes.csv  : node metadata of K=40 MKNN graph (to import into Gephi)
#   data_figs/04_DS1_refinedLouvain.pdf: Visual comparison Louvain vs Curated partition

require(SingleCellExperiment)
library(scater)
library(scran)
library(BiocParallel)
require(ggplot2)
require(patchwork)

source("aux_fun.R")
options(bitmapType = 'cairo')


bsave  <- FALSE
bcheck <- FALSE

# (0) Load DS1x raw data ----
(load('data_raw/01_DS1x_raw_ssce.Rdata'))

# (1) QC Filtering ----
# (1.1) QC stat ----
is.mito   <-  which(rowData(ssce)$chromosome_name == "MT")

bpp  <- MulticoreParam(param$ncores)
ssce <- addPerCellQC(ssce, subsets=list(Mt=is.mito), BPPARAM=bpp)

# (1.2) Discard nuclei with less than minDetected features ----
iin  <- ssce$detected > param$minDetected
ssce <- ssce[,iin]

# (1.3) Discard 'outlier' ----
out <- read.table('data/03_DS1x_discardedNuclei.txt')[,1]

ssce <- ssce[,!colnames(ssce)%in%out]
dim(ssce)
# [1] 44847 15895

# (1.3) Discard 636 doublets ----
(load('data_aux/01_DS1x_DblFinder.Rdata'))
nucDbl <- rownames(dblf.clus)[dblf.clus$scDblFinder.class=='doublet' & dblf.rnd$scDblFinder.class=='doublet']
ssce   <- ssce[,!colnames(ssce)%in%nucDbl]
dim(ssce)
# [1] 44847 15259

# (2) Analisis de-novo ----
# (2.1) reset params ----
param <- c()
param$bfiltered        <- TRUE   # TRUE: load Solo.GenFull/filtered, FALSE: load Solo.GeneFull/raw
param$pathSTAR         <-"/data1/Schinder/DG002-005-Run1n2/03_STAR/"
param$setName          <- 'DS1'

param$minDetected    <- 1000
param$Mt_percent_max <- 1                            #max % of reads mapped to Mt 
param$cellFiltMode   <- c('Mt','robustbase','scater')[2] #cell filtering mode
param$nmads          <- 3                                #cell outlier detection (deviations larger than nmads MAD from median)


param$features.minMolecules   <- 20 
param$features.minExpressedPct<- 1
param$features.maxExpressedPct<- 80
param$seed                    <- 123457
param$numHVG                  <- 3000     

param$ncores   <- 30



if(param$cellFiltMode=='robustbase'){
  library(robustbase)
  df <- colData(ssce)[,c('SampleID','sum','detected','subsets_Mt_percent')]
  df[,2] <- log10(df[,2])
  df[,3] <- log10(df[,3])
  if(all(df[,4]==0)) df <- df[,-4]
  
  #Per SampleID outlier detection
  usid <- unique(df$SampleID)
  multio <- c()
  for(i in seq_along(usid)){
    set.seed(123457)
    outlying <- adjOutlyingness(df[df$SampleID%in%usid[i],-1], only.outlyingness = TRUE)
    multio   <- rbind(multio,
                      cbind(multi_outlier= isOutlier(outlying, type = "higher"),
                            Mt_percent_max=colData(ssce)[ssce$SampleID%in%usid[i],'subsets_Mt_percent']>param$Mt_percent_max))
  }
  multio <- cbind(multio, discard=multio[,1] | multio[,2])
}


# (2.2) Nuclei filtering ----
if(param$cellFiltMode=='robustbase'){
  filt <- multio[,'discard']
}
if(param$cellFiltMode=='scater'){
  filt <- qc[,'discard']
}

ssce <- ssce[,!filt]
dim(ssce)
# [1] 44847 14441


## (2.3) Feature filtering ----
featureQC   <- perFeatureQCMetrics(ssce)
ssce        <- addPerFeatureQC(ssce)

out1<-rowData(ssce)$mean*ncol(ssce) < param$features.minMolecules
out2<-rowData(ssce)$detected       > param$features.maxExpressedPct
out3<-rowData(ssce)$detected       < param$features.minExpressedPct
discardGF <- out1 | out2 | out3

## Discard features
ssce <- ssce[!discardGF,]
dim(ssce)
# [1] 13455 14441

# (3) sce refinement ----
# (3.1) Filter Sex, Stress, Mitochondria genes ----
sexStress <- c("Ehd2","Espl1","Jarid1d","Pnpla4","Rps4y1","Xist","Tsix",
               "Eif2s3y", "Ddx3y", "Uty","Kdm5d","Rpl26","Gstp1","Rpl35a",
               "Erh","Slc25a5","Pgk1","Eno1","Tubb2a","Emc4", "Scg5" )

#descarto tambien los genes de mitocondria
mito <- unique(c(which(rowData(ssce)$chromosome_name == "MT"), 
                 grep("mt-",rowData(ssce)$mgi_symbol),
                 grep("Mrps",rowData(ssce)$mgi_symbol),
                 grep("Mrpl",rowData(ssce)$mgi_symbol)))
mito     <- rowData(ssce)$mgi_symbol[mito]
sexStress <- c(sexStress,mito)

ssce <- ssce[!rowData(ssce)$mgi_symbol%in%sexStress,]
dim(ssce)
# [1] 13354 14441

# (3.2) Normalizacion ----
bpp  <- MulticoreParam(param$ncores)

cat("Normalization:\n")
set.seed(param$seed)
clust<-quickCluster(ssce,graph.fun=igraph::infomap.community,block=factor(ssce$SampleID),BPPARAM=bpp)

# Cells in each cluster are normalized separately and the size factors are 
# rescaled to be comparable across clusters. Then cell-specific normalization factors are estimated
ssce <- computeSumFactors(ssce,cluster=clust,min.mean=0.1,BPPARAM=bpp)
ssce <- logNormCounts(ssce)

# (3.3) HVG ----
#Modeling  mean variance trend
mgv          <- modelGeneVar(ssce,span=.8,block=factor(ssce$SampleID))
rowData(ssce)<-cbind(rowData(ssce),hvg.mvBio=mgv$bio)

# Set isHVG column 
chosen.hvgs <- rowData(ssce)$hvg.mvBio > 0
chosen.hvgs <- rank(-rowData(ssce)$hvg.mvBio) <= param$numHVG & rowData(ssce)$hvg.mvBio>0
table(chosen.hvgs)
rowData(ssce)$isHVG <- chosen.hvgs

# Exclude duplicated mgi_symbol features
iout <- which(duplicated(rowData(ssce)[,'mgi_symbol']))
if(length(iout)>0) ssce <- ssce[-iout,]
rownames(ssce) <-   rowData(ssce)[,'mgi_symbol']
dim(ssce)
# [1] 13353 14441

# (4) runPCA, runTSNE, runUMAP  ----
require(BiocParallel)
pparam <- MulticoreParam(workers = 30)

set.seed(12534)
ssce <- runPCA(ssce, subset_row=rowData(ssce)$isHVG,
               BSPARAM=BiocSingular::RandomParam(),
               BPPARAM=pparam)

set.seed(1111001)
ssce <- runTSNE(ssce, dimred="PCA",BPPARAM=pparam)
set.seed(123457)
ssce <- runTSNE(ssce, dimred="PCA",name="TSNE3D",ncomponents=3,BPPARAM=pparam)

set.seed(1111001)
ssce <- runUMAP(ssce, dimred="PCA",BPPARAM=pparam)
set.seed(123457)
ssce <- runUMAP(ssce, dimred="PCA",name="UMAP3D",ncomponents=3,BPPARAM=pparam)

# (*) save point ----
bsave <- TRUE
if(bsave){
  readme<-'DS1x #detected>1000 -> outlierMKNN20->out doublets-> robustbased->qcFeatures->outSexStress->Symbol'
  save(ssce,param,readme,file=paste0('data/04_',param$setName,'.ssce.Rdata'))
}


# (5) MKNN manifold computation ----
library(igraph)

# (5.1) Distance matrix in PCA space
dimPCA <- 20  # PCA to keep
bcalcul<- !(file.exists('data/04_DS1.ssce.MKNN40.Rdata') & file.exists('data/04_DS1.ssce.Rdata'))
if(bcalcul){
  dm     <- buildDistMatrix(ssce,use.reducedDim="PCA",num.PCA=dimPCA,mode="pearson",bverbose = TRUE)
  
  # (5.2) MKNN graph
  mk  <- 40  # number of mutual neighbors
  gg  <- buildMKNN(ssce,dist.mat = dm,mutualK=mk,bOutGaph = TRUE,bWeighted = TRUE,bverbose = TRUE)

  # (5.3) Discard nuclei clustered in small components
  memb <- clusters(gg)$membership
  (tt<-table(memb))
  # memb
  #     1     2     3     4     5     6     7     8     9 
  # 14133    58   131    45    21    17     2     2     2 
  nn <- names(memb)[memb%in%c(1,2,3,4)]
  
  gg   <- induced_subgraph(gg,nn)
  ssce <- ssce[,nn]
  dim(ssce)
  # [1] 13353 14367
  
  bsave <- TRUE
  if(bsave){
    save(gg,mk,file=paste0('data/04_DS1.ssce.MKNN',mk,'.Rdata'))
    readme<-'DS1x #detected>1000 -> outlierMKNN2->out doublets-> robustbased->qcFeatures->outSexStress->Symbol'
    save(ssce,param,readme,file=paste0('data/04_',param$setName,'.ssce.Rdata'))
  }
}else{
  (load('data/04_DS1.ssce.MKNN40.Rdata'))
  (load('data/04_DS1.ssce.Rdata'))
}

# add Louvain clusterization slot
set.seed(104519)
ssce    <- doClustering(ssce,gg,list(louvain=igraph::cluster_louvain),prefix=paste0('MKNN',mk,'_'))

# (6) Pre-calculated partition refinement ----
(load('data_aux/04_DS1_clusterIDs.Rdata'))
ssce$clusP <- clusp


# (*) savepoint ----
if(bsave) {
  save(ssce,file='data/04_DS1.ssce.Rdata')
  exportGraphToGephi(ssce,gg,file=paste0('data/04_DS1.MKNN',mk))
}

foi <- c('louvain','cluster')
df0  <- data.frame(SampleID=ssce$SampleID,
                   louvain=factor(ssce$MKNN40_louvain),
                   cluster=ssce$clusP)
saux <- apply(df0[,foi],1,function(x){paste(foi,x,collapse="\n")})

poi <- foi[1]

# (6.1) Louvain and Refined partition pdf ----
lfig<-list()
smin <- 30
for(poi in c('louvain','cluster')){
  
  df <- data.frame(reducedDim(ssce,"TSNE")[,1:2],df0)
  tt <- table(df[,poi])
  df <- df[df[,poi]%in%names(tt)[tt>smin],]
  
  colnames(df)[1:2]<-c("x","y")
  
  colorize <- poi
  
  library(RColorBrewer)
  colourCount = length(levels(df[,poi]))
  getPalette  = colorRampPalette(brewer.pal(9, "Set1"))
  df$col <- getPalette(colourCount)[df[,poi]]
  fig <- ggplot(df,aes(x=x,y=y,col=col)) +
    geom_point(size=0.5)  +
    theme_bw() + 
    theme (legend.position="none") + ggtitle(poi) 
  
  lfig[[poi]]<-fig
}
lfig[['louvain']] + lfig[['cluster']]
#(*) savepoint ----
if(bsave) ggsave(filename = 'data_figs/04_DS1_refinedLouvain.pdf',width=8,height = 5)



