# input:
#  data/05_DS2_seuratTL.Rdata
#
# output:
#  data/07_DS2_pseudot_umap.Rdata: pseudotime assignment


options(bitmapType = 'cairo')

library(SingleCellExperiment)
library(igraph)
library(scater)
library(slingshot)

bsave <- TRUE  # save intermediate and final results

# (1) Load data ----
(load("data/05_DS2_seuratLT.Rdata")) #OJO con NB1/NB2/GCimm1
sce <- ssce;rm(ssce);gc()

# (2) ventrality signal ----
# We repeat the ventrality signal analysis performed on DS1
lfields <- list()
lfields[['ventrality']]    <- c('Atp2b4','Scn9a','Stk32b','Thsd7b','Trhr','Pde1a')
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

# clusP2: seuratLT + refinement
table(a,sce$clusP2,useNA='ifany')
# RGL  NPC  NB1  NB2 GCimm1 GCimm2 GCyoung GCmat1 GCmat2 Astro Oligo Peri  OPC <NA>
#   (0,0.278]     2049  459 4738    5   4140   1547    1151   1312     38   917   106   42   86   50
#   (0.278,0.557] 1062  179  890    2   1133    987     890    841     14   120     6   25  112   66
#   (0.557,0.835]  642  111  169    0    296    460     373    304      6    37     0   44   82   46
# (0.835,1.11]    74   15   26    0     74    161     151    144      8     6     0   17   65   35
# (1.11,1.39]      8    1   11    0     10     52      46     87      2     0     0    3   14   30
# (1.39,1.67]      0    0    4    0      1     13      21     38      0     0     0    0    0   13
# (1.67,1.95]      0    0    0    0      0      3       5     22      0     0     0    0    0    6
# (1.95,2.23]      0    0    0    0      0      1       0      7      0     0     0    1    0    0
# (2.23,2.5]       0    0    0    0      0      0       0      4      0     0     0    0    0    0

table(sce$ventrality_score>1.11)/length(sce$ventrality_score)
#      FALSE       TRUE 
# 0.98476568 0.01523432

(tt<-table(sce$clusP2[sce$ventrality_score>1.11]))
signif(tt/sum(tt),2)
#  clusP2
#    RGL     NPC     NB1     NB2  GCimm1  GCimm2 GCyoung  GCmat1  GCmat2   Astro   Oligo    Peri     OPC 
# 0.0220  0.0028  0.0450  0.0000  0.0310  0.1900  0.2000  0.4400  0.0056  0.0000  0.0000  0.0110  0.0450



# (2.1) Discard ventral nuclei ----
nVentral <- names(which(sce$ventrality_score>1.11 & sce$clusP2%in%'GCmat1'))

iin <- !colnames(sce)%in%nVentral
table(iin)
sce <- sce[,iin]


# (3.1) Trayectoria UMAP3D ----
# clusters de interes
coi <- c("RGL","NPC","NB1","GCimm1","GCimm2","GCyoung","GCmat1")
goi <- rowData(sce)$isHVG

# use time-reversed fit for stability reasons
breverse <- TRUE
if(breverse){
  slingUMAP3D <- slingshot(sce[goi,sce$clusP2%in%coi],
                           clusterLabels=sce$clusP2[sce$clusP2%in%coi],reducedDim='UMAP3D',
                           end.clus='RGL',start.clus='GCmat1',omega=Inf)
}else{
  slingUMAP3D <- slingshot(sce[goi,sce$clusP2%in%coi],
                           clusterLabels=sce$clusP2[sce$clusP2%in%coi],reducedDim='UMAP3D',
                           start.clus='RGL',end.clus='GCmat1',omega=Inf)
}




# (4) Pseudotime assignemnt ----
sce.sling         <- slingUMAP3D ; newDimRed='UMAP3D'
pseudoPaths       <- slingPseudotime(sce.sling)
curve.assignments <- slingBranchID(sce.sling)
table(curve.assignments)

shared.pseudo <- rowMeans(pseudoPaths,na.rm=TRUE)
ptime <- rep(NA,ncol(sce))
names(ptime) <- colnames(sce)

if(breverse) shared.pseudo <- max(shared.pseudo)-shared.pseudo

ptime[names(shared.pseudo)]               <- shared.pseudo
colData(sce)[names(ptime),'pseudot_umap'] <- ptime

colData(sce.sling)[names(shared.pseudo),'pseudot_umap'] <- shared.pseudo


if(bsave) save(ptime,file='data/07_DS2_pseudot_umap.Rdata')


