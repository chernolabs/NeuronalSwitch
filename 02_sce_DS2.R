# Preprocessing and QC of DS2 raw data
#  input : 
#   data_raw/02_DS2_raw_ssce.Rdata 
#  output:
#   data/02_DS2_DblFinder.Rdata   : computed scDblFinder doublet detection results
#   data/02_DS2.ssce.Rdata        : QC filtered, logNormalized dataset
#   data/02_DS2.ssce.MKNN40.Rdata : igraph k=40 MKNN graph
#   data/02_DS2.MKNN40.adjlist.csv: adjacency list of K=40 MKNN graph (to import into Gephi)
#   data/02_DS2.MKNN40.nodes.csv  : node metadata of K=40 MKNN graph (to import into Gephi)

require(DropletUtils)
require(SingleCellExperiment)
library(scater)
library(scran)
library(BiocParallel)
library(igraph)
source("aux_fun.R")
options(bitmapType = 'cairo')


# (1) Load raw dataset 2 consolidated into a SingleCellExperiment object:ssce ----
(load('data_raw/02_DS2_raw_ssce.Rdata'))



# (2) QC Filtering ----
# param$cellFiltModer:
#  a) scater::quickPerCellQC (outlier detection based on mads)
#  b) robustbase::adjOutlyingness (multivariate outlier detection)
#     additional filtering: max Mt%

# (2.1) QC stat ----
is.mito   <-  which(rowData(ssce)$chromosome_name == "MT")

bpp  <- MulticoreParam(param$ncores)
ssce <- addPerCellQC(ssce, subsets=list(Mt=is.mito), BPPARAM=bpp)
# dim: 44847 32621 

# Discard nuclei with less than minDetected features
iin  <- ssce$detected > param$minDetected
ssce <- ssce[,iin]
# dim: 44847 30295 

# (2.2) DoubletFinder ----
# Attention: For the sake of reproducibility  we included scDblFinder results in dblFinder_data folder
#
# We performed doublet identification with and without cluster assignments
# A doublet call was produced for nuclei identified as doublet using both methodologies
bcalcul <- !file.exists('data_aux/02_DS2_DblFinder.Rdata')
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
  
  bsave <- TRUE
  if(bsave)save(dblf.clus,dblf.rnd,file='data/02_DS2_DblFinder.Rdata')
}else{
  (load('data_aux/02_DS2_DblFinder.Rdata'))
}
a <- dblf.clus
b <- dblf.rnd
table(clus=a$scDblFinder.class,rnd=b$scDblFinder.class)
#             rnd
# clus      singlet doublet
# singlet     28264     714
# doublet       302    1015

# Discard doublets (interection)
nucDbl <- rownames(a)[a$scDblFinder.class=='doublet' & b$scDblFinder.class=='doublet']
ssce   <- ssce[,!colnames(ssce)%in%nucDbl]
dim(ssce)
# [1] 44847 29280

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
# [1] 44847 26798

## (2.3) Feature filtering ----
featureQC   <- perFeatureQCMetrics(ssce)
ssce        <- addPerFeatureQC(ssce)

out1<-rowData(ssce)$mean*ncol(ssce) < param$features.minMolecules
out2<-rowData(ssce)$detected       > param$features.maxExpressedPct
out3<-rowData(ssce)$detected       < param$features.minExpressedPct
discardGF <- out1 | out2 | out3

## Descarto features
ssce <- ssce[!discardGF,]
dim(ssce)
# [1] 14324 26798

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
# [1] 14221 26798

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


# (4) runPCA, runTSNE, runUMAP  ----
set.seed(12534)
require(BiocParallel)
pparam <- MulticoreParam(workers = 30)

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

# (*) savepoint (temporary) ----
bsave <- TRUE
if(bsave) save(ssce,param,file=paste0('data/02_',param$setName,'.ssce.Rdata'))


# (5) MKNN manifold computation ----
# (5.1) MKNN graph ----
dimPCA <- 20  # PCA to keep
bcalcul<- !file.exists('data/02_DS2.ssce.MKNN40.Rdata')
if(bcalcul){
  dm     <- buildDistMatrix(ssce,use.reducedDim="PCA",num.PCA=dimPCA,mode="pearson",bverbose = TRUE)
  
  mk  <- 40  # number of mutual neighbors
  gg  <- buildMKNN(ssce,dist.mat = dm,mutualK=mk,bOutGaph = TRUE,bWeighted = TRUE,bverbose = TRUE)
  if(bsave) save(gg,mk,file=paste0('data/02_DS2.ssce.MKNN',mk,'.Rdata'))
}else{
  (load('data/02_DS2.ssce.MKNN40.Rdata'))
}


# (5.1) Discard nuclei clustered in small components of the graph
memb <- clusters(gg)$membership
(tt<-table(memb))
# memb
#     1     2     3     4     5     6     7     8     9 
# 26144   101   369   102     2     3     2     2     2 
nn <- names(memb)[memb%in%names(tt)[tt>30]]

ssce <- ssce[,nn]
dim(ssce)

# (*) save point ----
bsave <- TRUE
if(bsave) save(ssce,param,file=paste0('data/02_DS2.ssce.Rdata'))

# (*) export to gephi ----
bsave <- TRUE
if(bsave) exportGraphToGephi(ssce,gg,file=paste0('data/02_DS2.MKNN',mk))



