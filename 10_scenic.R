#
# Regulon analysis for DS2
#
# This script launches an entire SCENIC analysis of DS2
# In order to reproduce the paper results you should consider the following cisTarget DBs:
#   mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather
#   mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather
# If not alreadyavailable at 'cisTargetPath' directory they will be downloaded from https://resources.aertslab.org/cistarget/databases
# 
# input: 
#   data/05_DS2_seuratLT.Rdata         : DS2
#   data/07_DS2_pseudot_umap.Rdata     : pseudo time
# output:
#  scenic_DS2_Rgeni3                   : SCENIC temporary directory
#  data_figs/10_DS2_RSS_regulons.pdf   : RSS plot for peaks 3 and 4 (Fig 7A)

require(SingleCellExperiment)
require(SCENIC)
require(ggplot2)
require(AUCell)
require(reticulate)
options(bitmapType='cairo')

# (0) Initialize settings ----
bpdf            <- TRUE                       # export pdf figure?
sname           <- paste0('DS2_RSS_regulons') # figure id
nameCOI         <- 'DS2'                      # analysis name 
label           <- 'clusP2'                   # clusterization label in the SingleCellObject

numcores        <- 30                         # number of cores
bOnlyPosCorr    <- FALSE                      # if FALSE, also consider negative correlations
method          <- 'Rgenie3'                  # scenic computation method
resumeFrom      <- NULL                       # resume from previous scenic calculation?   
cisTargetPath   <- 'scenic_cisTarget'         # folder of downloaded cisTarget dbs 


scenicPath           <- paste('scenic',nameCOI,method,sep='_')

# (0.1) cisTarget DBs for scenic-----
wd0 <- getwd()
if(!dir.exists(cisTargetPath)){
  dir.create(file.path('.', cisTargetPath), showWarnings = TRUE,recursive = TRUE)
  setwd(cisTargetPath)
  cat('Downloading cisTarget DBs...(this might take a while)\n')
  cisTarget_url <- c('https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather',
                     'https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather')
  for(i in seq_along(database_url)){
    aux <- paste('wget',database_url[i])
    system(aux)
  }
  setwd(wd0)
}

# (0.2) scenic directory structure ----
if(dir.exists(scenicPath)){
  setwd(file.path('.', scenicPath))
  scenicOptions<-readRDS('int/scenicOptions.Rds')
  if(!is.null(resumeFrom)){
    scenicOptions@status$current <- resumeFrom
  }
}else{
  scenicOptions <- initializeScenic(org="mgi", dbDir=paste0(getwd(),'/',cisTargetPath), nCores=numcores)
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  scenicOptions@inputDatasetInfo$colVars  <- "int/colVars.Rds"
  
  dir.create(file.path('.', scenicPath,'int'), showWarnings = TRUE,recursive = TRUE)
  setwd(file.path('.', scenicPath))
  saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
}
setwd(wd0)

# (0.3) Auxiliary functions ----
#function para calcular la moda de una variable numerica continua
estimate_mode <- function(x) {
  d <- density(x)
  out <- d$x[which.max(d$y)]
  return(out)
}

# filter ssce object by week and pseudotime value
condition_by_week_pseudotime <- function(ssce, week, pstime){
  
  condition <- rep(NA, ncol(ssce))
  
  for (i in 1:length(week)){
    name_cond <- paste(paste(week[[i]], collapse = ""), paste(round(pstime[[i]], 2), collapse = ":"),sep = "_")
    tmp <- ssce$week %in% week[[i]] & ssce$pseudoUMAP3D > pstime[[i]][1] & ssce$pseudoUMAP3D < pstime[[i]][2]
    condition[tmp] <- name_cond
  }
  
  return(condition)  
}



# (1.0) Load DS2 ----
(load('data/05_DS2_seuratLT.Rdata'))
(load('data/07_DS2_pseudot_umap.Rdata'))

ptime0        <- rep(NA,ncol(ssce))
names(ptime0) <- colnames(ssce)
ptime0[names(ptime)] <- ptime
colData(ssce)[names(ptime0),'pseudoUMAP3D'] <- ptime0

#Remove cells without pseudotime
table(!is.na(ssce$pseudoUMAP3D))
ssce <- ssce[,!is.na(ssce$pseudoUMAP3D)]

cellsOfInterest <- names(table(colData(ssce)[,label])[table(colData(ssce)[,label]) > 0])
ssce$clusP2      <- factor(ssce$clusP2, levels = cellsOfInterest)

exprMat           <- counts(ssce)
rownames(exprMat) <- rowData(ssce)$mgi_symbol

cellInfo <- colData(ssce)

# (2) scenic ----
setwd(scenicPath)
saveRDS(cellInfo,file=getDatasetInfo(scenicOptions, "cellInfo"))


# Color to assign to the variables (same format as for NMF::aheatmap)
if(scenicOptions@status$current==0){
  ctypes <- levels(cellInfo$cohort_pstime)
  colVars <- list(CellType=rainbow(length(ctypes))) 
  names(colVars[['CellType']]) <- ctypes       
  saveRDS(colVars, file="int/colVars.Rds")
  #plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))
  
  
  # (2.1) gene filtering ----
  # (Adjust minimum values according to your dataset)
  genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions=scenicOptions,
                             minCountsPerGene=3*.01*ncol(exprMat),
                             minSamples=ncol(exprMat)*.01)
  
  # (2.2) TF set
  library(dorothea)
  ddf <- data.frame(dorothea_mm)
  tfs <- unique(ddf$tf)
  
  genesOK   <- intersect(genesKept,rowData(ssce)$mgi_symbol[rowData(ssce)$isHVG])
  genesKept <- union(genesOK,tfs[tfs%in%genesKept])
  saveRDS(genesKept,file='int/1.1_genesKept.Rds')
  
}else{
  genesKept <- readRDS('int/1.1_genesKept.Rds')
}

exprMat_filtered <- exprMat[genesKept,]

# (2.2) Co-expression network ----
if(scenicOptions@status$current<1){
  runCorrelation(as.matrix(exprMat_filtered), scenicOptions)
}

exprMat_filtered_log <- log2(exprMat_filtered+1) 
exprMat_log          <- log2(exprMat+1)

# (2.3) GENIE3 in R ----
if(scenicOptions@status$current < 1){
  scenicOptions@settings$seed <- 123457  
  runGenie3(as.matrix(exprMat_filtered_log), scenicOptions)
  
  # (2.3.1) Create coexpression modules from regression modules ----
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  saveRDS(scenicOptions,'int/scenicOptions.Rds')
}

# (2.4) Regulons  ----
# (2.4.1) Create regulons (motif enrichment & prunning) ----
if(scenicOptions@status$current < 2){
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,onlyPositiveCorr = bOnlyPosCorr) 
  saveRDS(scenicOptions,'int/scenicOptions.Rds')
}
# (2.5) AUCell ----
if(scenicOptions@status$current < 3){
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, as.matrix(exprMat_log))
  scenicOptions@status$current <- 3  
  saveRDS(scenicOptions,'int/scenicOptions.Rds')
}

# (3) RSS ----

# (3.1) Determine local maximums from pseudotime density distribution ----
#1st peak
(p1 <- estimate_mode(ssce$pseudoUMAP3D[ssce$week %in% c("2w", "3w", "4w") & ssce$pseudoUMAP3D > 0 & ssce$pseudoUMAP3D < 5]))
#2nd peak
(p2 <- estimate_mode(ssce$pseudoUMAP3D[ssce$week %in% c("2w", "3w", "4w") & ssce$pseudoUMAP3D > 5 & ssce$pseudoUMAP3D < 9]))
#3rd peak
(p3 <- estimate_mode(ssce$pseudoUMAP3D[ssce$week %in% c("3w", "4w", "5w") & ssce$pseudoUMAP3D > 10 & ssce$pseudoUMAP3D < 25]))
#4th peak
(p4 <- estimate_mode(ssce$pseudoUMAP3D[ssce$week %in% c("4w", "5w", "8w") & ssce$pseudoUMAP3D > 25 & ssce$pseudoUMAP3D < 35]))

w = 2.5 #width around the peaks

ssce$picos <- NA
colData(ssce)[ssce$pseudoUMAP3D > (p1 - w) & ssce$pseudoUMAP3D < (p1 + w), "picos"] <- "p1"
colData(ssce)[ssce$pseudoUMAP3D > (p2 - w) & ssce$pseudoUMAP3D < (p2 + w), "picos"] <- "p2"
colData(ssce)[ssce$pseudoUMAP3D > (p3 - w) & ssce$pseudoUMAP3D < (p3 + w), "picos"] <- "p3"
colData(ssce)[ssce$pseudoUMAP3D > (p4 - w) & ssce$pseudoUMAP3D < (p4 + w), "picos"] <- "p4"

#cells by peak
table(colData(ssce)[, "picos"],useNA = 'ifany')

# (3.2) Calculate RSS ----
cellInfo$picos <- factor(ssce$picos, levels = c("p1", "p2", "p3", "p4"))
#cellInfo_s     <- cellInfo[!is.na(cellInfo$picos),]
cellInfo_s <- cellInfo[cellInfo$picos %in% c("p3", "p4"),] #pairwise comparison

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

gauc       <- getAUC(regulonAUC)
gauc       <- gauc[,rownames(cellInfo_s)]
rss        <- calcRSS(AUC=gauc, cellAnnotation=cellInfo_s[, "picos"])

# (3.3) Bubble plot RRS ----
rss_non_extended <- rss[onlyNonDuplicatedExtended(rownames(rss)), ]
#rss_non_extended <- rss_non_extended[, c("p1", "p2", "p3", "p4")] #sort peaks
rss_non_extended <- rss_non_extended[, c("p3", "p4")] #sort peaks
rssPlot_peaks   <- plotRSS(rss_non_extended, cluster_columns = FALSE)

if(bpdf){
  fname <- paste0(wd0,'/data_figs/10_',sname,'.pdf')
  rssPlot_peaks$plot
  ggsave(fname,width=4,height=16)
}else{
  rssPlot_peaks$plot
}


