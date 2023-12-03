# input:
#  data/04_DS1.ssce.Rdata
#
# output:
#  data/07_DS1_pseudot_umap.Rdata: pseudotime assignment

options(bitmapType = 'cairo')

library(SingleCellExperiment)
library(igraph)
library(scater)
library(slingshot)

bsave <- TRUE


# (1) Loading data ----
(load("data/04_DS1.ssce.Rdata")) #OJO con NB1/NB2/GCimm1
sce <- ssce;rm(ssce);gc()


# (2) ventrality signal ----
# We detected a subset of GCmat1 nuclei (mainly from cluster-Louvain 6), 
# that occupied a characteristic location on the embedding space,
# and exhibited high ventrality signature levels. We decided to disregard those nuclei in order to get
# a robust pseudotime estimation.
lfields <- list()
lfields[['ventrality']]   <- c('Atp2b4','Scn9a','Stk32b','Thsd7b','Trhr','Pde1a')
for(i in seq_along(lfields)){
  genes <- lfields[[i]]
  genes <- genes[genes%in%rownames(sce)]
  if(length(genes)>0){
    aux <- apply(logcounts(sce)[genes,],2,mean)
    colData(sce)[,paste0(names(lfields)[i],'_score')]  <- aux
  }
}
ventralityScore <- sce$ventrality_score
hist(ventralityScore);abline(v=seq(min(ventralityScore),max(ventralityScore),length=5),col='gray',lty=2)
a <- cut(ventralityScore,seq(min(ventralityScore),max(ventralityScore),length=10))
a[which(is.na(a))] <- levels(a)[1] 

# These nuclei mainly belong to cluster 6 from Louvain partition
table(a,sce$MKNN40_louvain,useNA='ifany')
# a                  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
#   (0,0.217]      577 1600 2119  624 1085  265  485  184  199  110  145    9   41  673   45  414 1024
#   (0.217,0.433]   78  198  327  104  291  220  217   59  142    8   40   21   47  289    0  339  573
#   (0.433,0.65]     9   27   47   13   92  148  140   25   74    4   27   24   28   64    0  164  239
#   (0.65,0.866]     1    5    7    1   15  126   31    1    6    1    9    4   13   15    0   46   96
#   (0.866,1.08]     0    0    1    0    2  101    0    0    0    0    7    0    2    1    0    6   40
#   (1.08,1.3]       0    0    0    0    1   57    0    0    0    1    3    0    0    0    0    4   20
#   (1.3,1.52]       0    0    0    0    0   35    0    0    0    0    0    0    0    0    0    1    2
#   (1.52,1.73]      0    0    0    0    0   21    0    0    0    0    0    0    0    0    0    0    0
#   (1.73,1.95]      0    0    0    0    0    7    0    0    0    0    0    0    0    0    0    0    1


# Criterio: ventrality > 0.866 & louvain 6
nVentral <- names(which(sce$ventrality_score>0.866 & sce$MKNN40_louvain==6))


# (2.1) Discard ventral nuclei ----
iin <- !colnames(sce)%in%nVentral
table(iin)
sce <- sce[,iin]


# (3) Trayectoria UMAP3D ----
# clusters of interest
coi <- c("RGL","NPC","NB1","NB2","GCimm1","GCimm2","GCyoung","GCmat1")
goi <- rowData(sce)$isHVG

# estimating the principal curves
breverse <- TRUE
if(breverse){
  slingUMAP3D <- slingshot(sce[goi,sce$clusP%in%coi],reducedDim='UMAP3D',
                           clusterLabels=sce$clusP[sce$clusP%in%coi],
                           end.clus='RGL',start.clus='GCmat1',omega=Inf)
}else{
  slingUMAP3D <- slingshot(sce[goi,sce$clusP%in%coi],reducedDim='UMAP3D',
                           clusterLabels=sce$clusP[sce$clusP%in%coi],
                           start.clus='RGL',end.clus='GCmat1',omega=Inf)
}


# (4) Pseudotime assignemnt ----

sce.sling         <- slingUMAP3D ; newDimRed='UMAP3D'
pseudoPaths       <- slingPseudotime(sce.sling)
curve.assignments <- slingBranchID(sce.sling)

shared.pseudo <- rowMeans(pseudoPaths,na.rm=TRUE)
ptime <- rep(NA,ncol(sce))
names(ptime) <- colnames(sce)

if(breverse) shared.pseudo <- max(shared.pseudo)-shared.pseudo

ptime[names(shared.pseudo)] <- shared.pseudo
sce$pseudot_umap <- ptime


if(bsave) save(ptime,file='data/07_DS1_pseudot_umap.Rdata')


