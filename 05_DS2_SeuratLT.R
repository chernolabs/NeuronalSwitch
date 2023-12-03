# Transfer cluster labels from DS1 to DS2 using Seurat
# Input:
#   data/04_DS1.ssce.Rdata: Dataset1
#   data/02_DS2.ssce.Rdata: Dataset2
#
# Output:
#   data/05_DS1DS2_seuratLT_anchors.Rdata : Anchors estimated bty FindTransferAnchors function
#   data/05_DS2_seuratLT.Rdata            : Dataset2 SingleCellExpression object with transferred cluster labels (clusP)
#   data_figs/05_DS2_refinedSeuratLT.pdf

library(SingleCellExperiment)
library(Seurat)
library(Matrix)
options(bitmapType = 'cairo')

# (0) Load data ----
(load('data/04_DS1.ssce.Rdata'))
sceRef <- ssce

(load('data/02_DS2.ssce.Rdata'))
sceQuery <- ssce


# (1) common genes ----
genes <- intersect(rownames(sceQuery),rownames(sceRef))
genes <- genes[genes!='']
length(genes)
#[1] 13194
sceQuery <- sceQuery[genes,]
sceRef <- sceRef[genes,]
  
# (2) SCE to Seurat: ----
dsQuery                   <- as.Seurat(sceQuery)
smatdata                  <- Matrix(GetAssayData(object = dsQuery, slot = "data"), sparse = TRUE)
dsQuery                   <- SetAssayData(dsQuery, slot = "data", new.data = smatdata)
VariableFeatures(dsQuery) <- rownames(sceQuery)[rowData(sceQuery)$isHVG]
  
dsRef                   <- as.Seurat(sceRef)
smatdata                <- Matrix(GetAssayData(object = dsRef, slot = "data"), sparse = TRUE)
dsRef                   <- SetAssayData(dsRef, slot = "data", new.data = smatdata)
VariableFeatures(dsRef) <- rownames(sceRef)[rowData(sceRef)$isHVG]
  
# (3) Anchors ----
bcalcul <- !file.exists('data_aux/05_DS1DS2_seuratLT_anchors.Rdata')
if(bcalcul){
    anchors <- FindTransferAnchors(reference = dsRef, query = dsQuery, dims = 1:30)
    if(TRUE) save(anchors,file='data/05_DS1DS2_seuratLT_anchors.Rdata')
}else{
    (load('data_aux/05_DS1DS2_seuratLT_anchors.Rdata'))
}

# (4) Transfer labels ----
pred        <- TransferData(anchorset = anchors, refdata = dsRef$clusP)$predicted.id
names(pred) <- colnames(dsQuery) 

dds1  <- dsQuery
dds1  <- AddMetaData(object = dds1, metadata = pred,col.name = 'clusP')

ssceAux       <- as.SingleCellExperiment(dds1)
ssceAux$clusP <- factor(ssceAux$clusP,levels=levels(sceRef$clusP))

table(week=ssceAux$week,clus=ssceAux$clusP)
#    clus
# week  RGL  NPC  NB1  NB2 GCimm1 GCimm2 GCyoung GCmat1 GCmat2 Astro Oligo Peri  OPC
#   2w 1580  291 2726    6   2441    409       1     21      2    88    26    8   86
#   3w  657   79 1006    0   1248   1431     342     24      1    83     7   28   69
#   4w 1267  170 1416    1   1355    996    1706    321      9   122    28   10   97
#   5w  547   66  440    0    343    280     419    836     26   100    15   19   48
#   8w  382   47  363    0    267    111     175   1645     43   268    36   23   59

# (4.1) Seurat LT refinement ----
# Refinement to smooth away heterogeneities in UMAP space. They mainly involved
# RGL, NB, NPC reassignations to align labels  with UMAP large nuclei aggregation structures
# Astro and Peri groups were also affected (eventhough these groups were not considered for pseudo time estimation)

df <- reducedDim(ssceAux,"UMAP")[,1:2]
df <- data.frame(df,clusP=ssceAux$clusP)
df$clusP2 <- df$clusP
colnames(df)[1:2]<-c('x','y')


# (a) NB1/Pericytes  -> Pericytes
#i <- which(df$x > 3  & df$x < 4.5  & df$y< -0.5 & df$y > -1.5)
i <- which(df$x > 1.85  & df$x < 4.55  & df$y< -0.4 & df$y > -1.6)
tt<-table(df$clusP2[i])
tt/sum(tt)
df$clusP2[i] <- 'Peri'


# (b) RGL/Astrocytes -> Astrocytes
#i <- which(df$x > 4  & df$x < 7  & df$y< -6.5 & df$y > -9.5)
i <- which(df$x > 4.5  & df$x < 7  & df$y< -6 & df$y > -9)
tt<-table(df$clusP2[i])
tt/sum(tt)
df$clusP2[i] <- 'Astro'

# (c) isolated RGL -> NA (non-conclusive)
#i <- which(df$x > 2.2  & df$x < 2.8  & df$y< -3.3 & df$y > -3.8)
i <- which(df$x > 1  & df$x < 1.4  & df$y< -7.7 & df$y > -7.9)
df$clusP2[i] <- NA

# (d) NPC+RGL+GCMat  -> NA (heterogeneous aggregation)
#i <- which(df$x > -1  & df$x < 0  & df$y< -3 & df$y > -3.8)
i <- which(df$x > -1.4  & df$x < -0.8  & df$y< -3.4 & df$y > -4.2)
table(df$clusP2[i])
# RGL     NPC     NB1     NB2  GCimm1  GCimm2 GCyoung  GCmat1  GCmat2   Astro   Oligo    Peri     OPC 
#   0      24      11       0       0       3       6      88      13       0       0       0       0 
df$clusP2[i] <- NA

# (e) NPC+RGL -> RGL
#i <- which(df$x > 7  & df$x < 9  & df$y< -3 & df$y > -5)
i <- which(df$x > 8  & df$x < 9  & df$y< -3.5 & df$y > -5)
tt<-table(df$clusP2[i])
tt/sum(tt)  
df$clusP2[i] <- 'RGL'

# (f) NB1+RGL+NPC -> NPC
#i <- which(df$x > 6.5  & df$x < 8.5  & df$y< 1 & df$y > -2 )
i <- which(df$x > 7  & df$x < 9  & df$y< 1 & df$y > -3 )
tt<-table(df$clusP2[i])
tt/sum(tt)
df$clusP2[i] <- 'NPC'

# (g) Oligo -> OPC
i <- which(df$x > 0.6  & df$x < 0.8  & df$y< -5.1 & df$y > -5.3 )
tt<-table(df$clusP2[i])
tt/sum(tt)
df$clusP2[i] <- 'OPC'


# (4.2) Visualization ----
library(ggplot2)
library(patchwork)

fig1 <- ggplot(df,aes(x=x,y=y,col=clusP)) +
  geom_point() + 
  ggtitle('DS2',subtitle='Seurat LabelTransfer partition') + 
  xlab('UMAP 1') + ylab('UMAP 2') #+
  #theme(legend.position="none")

fig2 <- ggplot(df,aes(x=x,y=y,col=clusP2)) +
  geom_point() + 
  ggtitle('DS2',subtitle='refined partition') + 
  xlab('UMAP 1') + ylab('UMAP 2') #+
  #theme(legend.position="none")


rect1 <- data.frame(xmin=1.85, xmax=4.55, ymin=-1.6, ymax=-0.4)
rect2 <- data.frame(xmin=4.5, xmax=7, ymin=-9, ymax=-6)
rect3 <- data.frame(xmin=1, xmax=1.4, ymin=-7.7, ymax=-7.9)
rect4 <- data.frame(xmin=-1.4, xmax=-0.8, ymin=-3.4, ymax=-4.2)
rect5 <- data.frame(xmin=8, xmax=9, ymin=-5, ymax=-3.5)
rect6 <- data.frame(xmin=7, xmax=9, ymin=-3, ymax=1)
rect7 <- data.frame(xmin=0.6, xmax=0.8, ymin=-5.3, ymax=-5.1)


fig1<-fig1 + 
  geom_rect(data=rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) + 
  geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=rect4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=rect5, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=rect6, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=rect7, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) 


fig1 + fig2
ggsave(filename = 'data_figs/05_DS2_refinedSeuratLT.pdf',width=11,height=4)

# Summary of changes
(tt<-table(seuratLT=df$clusP,refined=df$clusP2,useNA='ifany'))
# refined
# seuratLT   RGL  NPC  NB1  NB2 GCimm1 GCimm2 GCyoung GCmat1 GCmat2 Astro Oligo Peri  OPC <NA>
#   RGL     3830   84    0    0      0      0       0      0      0   418     0    0    0  101
#   NPC        9  619    0    0      0      0       0      0      0     1     0    0    0   24
#   NB1        1   59 5819    0      0      0       0      0      0     1     0   60    0   11
#   NB2        0    0    0    7      0      0       0      0      0     0     0    0    0    0
#   GCimm1     0    0    0    0   5654      0       0      0      0     0     0    0    0    0
#   GCimm2     0    0    0    0      0   3224       0      0      0     0     0    0    0    3
#   GCyoung    0    0    0    0      0      0    2637      0      0     0     0    0    0    6
#   GCmat1     0    0    0    0      0      0       0   2759      0     0     0    0    0   88
#   GCmat2     0    0    0    0      0      0       0      0     68     0     0    0    0   13
#   Astro      0    0    0    0      0      0       0      0      0   661     0    0    0    0
#   Oligo      0    0    0    0      0      0       0      0      0     0   102    0   10    0
#   Peri       0    0    0    0      0      0       0      0      0     0     0   88    0    0
#   OPC        0    0    0    0      0      0       0      0      0     0     0    0  359    0

nna      <- sum(tt[,14])
nrelabel <- sum(tt[,-14])-sum(diag(tt))

cat(' Number of disregarded nuclei:',nna,signif(nna/sum(tt),2),'\n',
    'Number of relabelled nuclei  :',nrelabel,signif(nrelabel/sum(tt),2),'\n')


ssceAux$clusP2 <- df$clusP2


colData(ssce)  <- colData(ssceAux)

if(ncol(rowData(ssce))==0)rowData(ssce) <- rowData(sceQuery)[rownames(ssce),]


# (*) savepoint ----
bsave<-TRUE
if(bsave) save(ssce,param,file='data/05_DS2_seuratLT.Rdata')


