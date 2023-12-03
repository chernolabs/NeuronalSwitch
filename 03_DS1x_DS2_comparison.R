# Data integration for replicability analysis
#
# input:
#   data/01_DS1x.ssce.Rdata
#   data/02_DS2x.ssce.Rdata
# output:
#   data/03_DS1x_discardedNuclei.txt  : 1071 discarded DS1x nuclei
#   data_figs/03_DS1x_DS2_umap.pdf


require(batchelor)
library(BiocParallel)
options(bitmapType = 'cairo')

# (0) Load data ----
(load('data/01_DS1x.ssce.Rdata'))
sce1 <- ssce
(load('data/02_DS2.ssce.Rdata'))
sce2 <- ssce
rm(ssce);gc()

# (1) Dataset quick integration ----
bpp  <- MulticoreParam(30)
quick.corrected <- quickCorrect(sce1, sce2, 
                                PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam(),
                                                   BPPARAM=bpp))
sce  <- quick.corrected$corrected

# (1.1) Add metadata ----
sce$batch   <- factor(sce$batch)

cohort        <- rep(NA,ncol(sce))
names(cohort) <- colnames(sce)
cohort[colnames(sce1)] <- sce1$week
cohort[colnames(sce2)] <- sce2$week
sce$cohort <- factor(cohort)
table(sce$batch,sce$cohort,useNA='ifany')
#     1w   2w   3w   4w   5w   8w
# 1 2893 4255    0 4906    0 3437
# 2    0 7685 4975 7498 3139 3419

# (1.2) DS1x clusters
cluster                 <- rep(NA,ncol(sce))
names(cluster)          <- colnames(sce)
cluster[colnames(sce1)] <- sce1$MKNN20_louvain
sce$clusDS1             <- factor(cluster)

# (2) TSNE
library(scater)
set.seed(123457)
sce <- runTSNE(sce, dimred="corrected")
sce <- runUMAP(sce, dimred="corrected")

# (*) savepoint ----
bsave <- TRUE
if(bsave) save(sce,file='data/03_DS1xDS2.sce.Rdata')

# (2) Visualization ----
# Community 1 
library(patchwork)

sce$color <- rep(NA,ncol(sce))
df <- data.frame(x=reducedDim(sce,'UMAP')[,1],y=reducedDim(sce,'UMAP')[,2],col=sce$color,batch=sce$batch,cohort=sce$cohort)
df$col[df$batch==1 & sce$clusDS1%in%c(25,27)] <- 'red'
fig1 <- ggplot(df[df$batch==1,],aes(x=x,y=y,col=col)) +
          geom_point() + 
          ggtitle('DS1x',subtitle='red: non-reproducible clusters') + 
          xlab('UMAP 1') + ylab('UMAP 2') +
          theme(legend.position="none")
fig2 <- ggplot(df[df$batch==2,],aes(x=x,y=y)) +
          geom_point(col='gray') + 
          ggtitle('DS2',subtitle='.') + 
          xlab('UMAP 1') + ylab('UMAP 2') +
          theme(legend.position="none")
fig1 + fig2
ggsave(filename = 'data_figs/03_DS1x_DS2_umap.pdf',width=8,height=4)



# The offending communities were 4-week old clusDS1-25 and clusDS1-27
# We also detected some under-representation of clusDS1x-1 nuclei, but as it was mainly composed
# of 1-week old nuclei (age not present in Dataset2) we decided to reteain it
table(sce1$MKNN20_louvain[sce1$MKNN20_louvain%in%c(1,25,27)],sce1$week[sce1$MKNN20_louvain%in%c(1,25,27)])
#     1w  2w  4w
# 1  546 270  22
# 25   0   0 167
# 27   0   0 904

# (3) Export discarded nuclei list ----
bsave <- TRUE
outNuclei <- colnames(sce1)[sce1$MKNN20_louvain%in%c(25,27)]
if(bsave) write.table(outNuclei,'data/03_DS1x_discardedNuclei.txt',row.names = FALSE,quote=FALSE,col.names=FALSE)



 