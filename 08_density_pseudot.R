# in
#  data/04_DS1.ssce.Rdata
#  data/05_DS2_seuratLT.Rdata
#  data/07_DS1_pseudot_umap.Rdata
#  data/07_DS2_pseudot_umap.Rdata
#
# out
#  data_figs/08_density_pseudot.pdf

library(ggplot2)
library(SingleCellExperiment)
options(bitmapType = 'cairo')

# (0) Load Datasets ----
# (0.1) Load DS1 ----
(load('data/04_DS1.ssce.Rdata'))
sce_ds1 <- ssce;rm(ssce);gc()

(load('data/07_DS1_pseudot_umap.Rdata'))
ptime0        <- rep(NA,ncol(sce_ds1))
names(ptime0) <- colnames(sce_ds1)
ptime0[names(ptime)] <- ptime
colData(sce_ds1)[names(ptime0),'pseudot_umap'] <- ptime0

# (0.2) Load DS2 ----
(load('data/05_DS2_seuratLT.Rdata'))
sce_ds2 <- ssce;rm(ssce);gc()

(load('data/07_DS2_pseudot_umap.Rdata'))
ptime0        <- rep(NA,ncol(sce_ds2))
names(ptime0) <- colnames(sce_ds2)
ptime0[names(ptime)] <- ptime
colData(sce_ds2)[names(ptime0),'pseudot_umap'] <- ptime0

# (1.0 ) Density plot ----

#DS1
df_ds1 <- data.frame(week=factor(sce_ds1$week, levels = c('1w','2w', '4w', '8w')),
                    pseudotime=sce_ds1$pseudot_umap, type = sce_ds1$clusP)
df_ds1 <- df_ds1[!is.na(df_ds1$pseudotime),]
df_ds1$dataset <- "Dataset 1"

#DS2
df_ds2 <- data.frame(week=factor(c(sce_ds2$week), levels = c('2w', '3w', '4w', '5w', '8w')),
                    pseudotime=sce_ds2$pseudot_umap, type = sce_ds2$clusP2)
df_ds2 <- df_ds2[!is.na(df_ds2$pseudotime),]
df_ds2$dataset <- "Dataset 2"

#Plot
df <- rbind(df_ds1, df_ds2)

df$type <- factor(df$type, levels = c("RGL", "NPC", "NB1", "NB2", "GCimm1", "GCimm2", "GCyoung", "GCmat1"))
df$week <- factor(df$week, levels = c("1w", "2w", "3w", "4w", "5w", "8w"))
df$dataset <- factor(df$dataset, levels =c("Dataset 1", "Dataset 2")) 
#df$type <- factor(df$type, levels = c("RGL", "NPC", "NB1", "GCimm1", "GCimm2", "GCyoung", "GCmat1"))

colClus <- c("#F5BA6C","#EB7604","#C6ADCF","#663E8E",
             "#A75528","#DC2A85","#2E9737","#DEA709")
names(colClus) <- c("RGL","NPC","NB1","NB2",
                    "GCimm1","GCimm2","GCyoung","GCmat1")

colWeeks <- c("#31c7ba", "#b5361c", "#369cc9", "#3a507f", "#1c9d7c", "#e35e28")
names(colWeeks) <- c("1w", "2w", "3w", "4w", "5w", "8w")

p22 <- ggplot(data=df,aes(x=pseudotime,fill=week)) +
  geom_density(alpha=0.6, adjust=1.2) +
  geom_rug(aes(color = type), length = unit(0.05, "npc"), sides = 'b', outside = F) +
  facet_grid(week ~ dataset, scales=c('free','fixed')[1]) + theme_light() +#rows=vars(week), cols=vars(dataset)
  theme(panel.spacing = unit(0.5, "lines"), strip.text.y = element_text(angle = 0))
 p22 + 
   scale_color_manual(values = colClus) +
   scale_fill_manual(values=rep(colWeeks, 2))

bpdf <- TRUE
if(bpdf){
  ggsave(filename = 'data_figs/08_density_pseudot.pdf',width=6,height=8)
}

