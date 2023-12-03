# Preprocessing and QC of DS1 raw data
#  input : 
#   data_raw/01_DS1x_raw_ssce.Rdata 
#  output:
#   data/01_DS1x_DblFinder.Rdata   : computed scDblFinder doublet detection results
#   data/01_DS1x.ssce.Rdata        : QC filtered, logNormalized dataset
#   data/01_DS1x.ssce.MKNN20.Rdata : igraph k=20 MKNN graph used to identified nuclei communities not observed in DS2


require(SingleCellExperiment)
library(scater)
library(scran)
library(BiocParallel)
library(igraph)
source('aux_fun.R')
options(bitmapType = 'cairo')


# (1) Load raw dataset 1 consolidated into a SingleCellExperiment object:ssce ----
(load('data_raw/01_DS1x_raw_ssce.Rdata'))


# (2) QC Filtering ----
# param$cellFiltModer:
#  a) scater::quickPerCellQC (outlier detection based on mads)
#  b) robustbase::adjOutlyingness (multivariate outlier detection)
#     additional filtering: max Mt%


# (2.1) QC stat ----
is.mito   <-  which(rowData(ssce)$chromosome_name == "MT")

bpp  <- MulticoreParam(param$ncores)
ssce <- addPerCellQC(ssce, subsets=list(Mt=is.mito), BPPARAM=bpp)


# Discard nuclei with less than minDetected features
iin  <- ssce$detected > param$minDetected
table(iin)
# FALSE  TRUE 
#    80 16966
ssce <- ssce[,iin]

# (2.2) DoubletFinder ----
# Attention: For the sake of reproducibility  we included scDblFinder results in dblFinder_data folder
#
# We performed doublet identification with and without cluster assignments
# A doublet call was produced for nuclei identified as doublet using both methodologies
bcalcul <- !file.exists('data_aux/01_DS1x_DblFinder.Rdata')
if(bcalcul){
  library(scDblFinder)
  
  pparam <- MulticoreParam(workers = 30, RNGseed=123457)

  set.seed(123457)
  fcluster <- fastcluster(ssce,nstart=10,ndims=20,nfeatures=1000)
  names(fcluster)<-colnames(ssce)
  
  ssce$fcluster <- factor(fcluster)
  
  a  <- scDblFinder(ssce,clusters='fcluster',samples = ssce$SampleID,multiSampleMode = 'split',returnType = 'table',verbose=TRUE)
  dblf.clus <- colData(a)
  b  <- scDblFinder(ssce,samples = ssce$SampleID, multiSampleMode = 'split',returnType = 'table',verbose=TRUE)
  dblf.rnd <- colData(b)
  bsave<-TRUE
  if(bsave)save(dblf.clus,dblf.rnd,file='data/01_DS1x_DblFinder.Rdata')
}else{
  (load('data_aux/01_DS1x_DblFinder.Rdata'))
}
a <- dblf.clus
b <- dblf.rnd
table(clus=a$scDblFinder.class,rnd=b$scDblFinder.class)
#                rnd
# clus      singlet doublet
#   singlet   15684     406
#   doublet     240     636

# Discard doublets  (intersection)
nucDbl <- rownames(a)[a$scDblFinder.class=='doublet' & b$scDblFinder.class=='doublet']
ssce   <- ssce[,!colnames(ssce)%in%nucDbl]


# (2.3) QC ----
qc   <- quickPerCellQC(colData(ssce), percent_subsets ='subsets_Mt_percent',
                      batch=ssce$SampleID,nmads=param$nmads)

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
    # ojo!
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
  filt <- multio
}else{
  if(param$cellFiltMode=='scater'){
    filt <- qc
  }else{
    warning('ohoh\n')
  }
}  
ssce <- ssce[,!filt[,'discard']]
dim(ssce)
# [1] 44847 15491

## (2.3) Feature filtering ----
featureQC   <- perFeatureQCMetrics(ssce)
ssce        <- addPerFeatureQC(ssce)

out1<-rowData(ssce)$mean*ncol(ssce) < param$features.minMolecules
out2<-rowData(ssce)$detected       > param$features.maxExpressedPct
out3<-rowData(ssce)$detected       < param$features.minExpressedPct
discardGF <- out1 | out2 | out3
table(discardGF)
# discardGF
# FALSE  TRUE 
# 15113 29734

## Discard features
ssce <- ssce[!discardGF,]

# (3) sce refinement ----
# (3.1) Filter Sex, Stress, Mitochondria genes ----
sexStress <- c("Ehd2","Espl1","Jarid1d","Pnpla4","Rps4y1","Xist","Tsix",
               "Eif2s3y", "Ddx3y", "Uty","Kdm5d","Rpl26","Gstp1","Rpl35a",
               "Erh","Slc25a5","Pgk1","Eno1","Tubb2a","Emc4", "Scg5" )

mito <- unique(c(which(rowData(ssce)$chromosome_name == "MT"), 
                 grep("mt-",rowData(ssce)$mgi_symbol),
                 grep("Mrps",rowData(ssce)$mgi_symbol),
                 grep("Mrpl",rowData(ssce)$mgi_symbol)))
mito     <- rowData(ssce)$mgi_symbol[mito]
sexStress <- c(sexStress,mito)

ssce <- ssce[!rowData(ssce)$mgi_symbol%in%sexStress,]
dim(ssce)
# > dim(ssce)
# [1] 15009 15491

# (3.2) Normalizacion ----
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


# (4) runPCA, runTSNE  ----
set.seed(12534)
require(BiocParallel)
pparam <- MulticoreParam(workers = 30)

ssce <- runPCA(ssce, subset_row=rowData(ssce)$isHVG,
               BSPARAM=BiocSingular::RandomParam(),
               BPPARAM=pparam)
set.seed(1111001)
ssce <- runTSNE(ssce, dimred="PCA",BPPARAM=pparam)

# (5) Include Louvain community information ----
# (5.1) Distance matrix in PCA space
dimPCA <- 20  # PCA to keep
bcalcul<- !file.exists('data_aux/01_DS1x.ssce.MKNN20.Rdata')
if(bcalcul){
  dm     <- buildDistMatrix(ssce,use.reducedDim="PCA",num.PCA=dimPCA,mode="pearson",bverbose = TRUE)
  
  # (5.2) MKNN graph
  mk  <- 20  # number of mutual neighbors: This resolution served to
             # clearly identify nuclei structures not reproduced in DS2.
  gg  <- buildMKNN(ssce,dist.mat = dm,mutualK=mk,bOutGaph = TRUE,bWeighted = TRUE,bverbose = TRUE)
  bsave <- TRUE
  if(bsave) save(gg,mk,file=paste0('data/01_DS1x.ssce.MKNN',mk,'.Rdata'))
}else{
  (load('data_aux/01_DS1x.ssce.MKNN20.Rdata'))
}
set.seed(123457)
ssce    <- doClustering(ssce,gg,list(louvain=igraph::cluster_louvain),prefix=paste0('MKNN',mk,'_'))

bsave <- TRUE
if(bsave) save(ssce,param,file=paste0('data/01_',param$setName,'.ssce.Rdata'))


#---------------------End ----

