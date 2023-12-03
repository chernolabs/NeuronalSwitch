
# TF expression profile cascade
#
# input: 
#   data/05_DS2_seuratLT.Rdata    : DS2
#   data/07_DS2_pseudot_umap.Rdata: pseudo time
# output:
#  data_figs/09_DS2_tfsCascade.pdf


require(SingleCellExperiment)
require(scran)
require(dorothea)
library(pheatmap)
library(ComplexHeatmap)
library(dendsort)
options(bitmapType = 'cairo')


# (0) Load data and initialize ----
(load('data/05_DS2_seuratLT.Rdata'))
(load('data/07_DS2_pseudot_umap.Rdata'))

# add pseudotime to colData
ssce$pt_umap                          <- rep(NA,ncol(ssce))
colData(ssce)[names(ptime),'pt_umap'] <- ptime

# (0.1) initialize things ----
bpdf      <- TRUE    #export pdf figure

bdendsort <- FALSE             #use dendsort to reorder dendrograms
sname     <- paste0('DS2_tfsCascade') # figure id
label     <- 'clusP2'          # lets consider this clusterization label

ssce      <- ssce[,!colData(ssce)[,label]%in%'NB2'] # Discard 7 nuclei dubiously mapped to NB2 by seuratLT
clusoi    <- c('RGL','NPC','NB','GCimm1','GCimm2','GCyoung','GCmat1') #clusters of interests [-c(1:2)]
saux                         <- as.character(colData(ssce)[,label])
saux[saux%in%c('NB1','NB2')] <- 'NB'
colData(ssce)[,label]        <- factor(saux,levels=clusoi)

# (0.2) keep clusterized and pseudotime-assigned nuclei  ----
iin <- !(is.na(colData(ssce)[,label]) | is.na(ssce$pt_umap))
table(iin,useNA='ifany')
# FALSE  TRUE 
#  2173 24536
if(sum(iin)<length(iin)) ssce <- ssce[,iin]

# (0.3) dorothea TFs  ----
require(dorothea)
tfs  <- unique(dorothea_mm$tf)
tfs  <- tfs[tfs%in%rownames(ssce)]
length(tfs)
# [1] 700

# (1) Sort nuclei according to estimated pseudotime ----
ssce  <- ssce[,order(colData(ssce)[,'pt_umap'],decreasing=FALSE)]


# (2) difftest ----

# Parallel implementation of TSCAN difftest
# Honestly, do not remember if this hack was mine. If you know the author
# please let us know.
pdifftest<-function(data, TSCANorder, df = 3,numcores=15){ 
  library(future.apply)
  plan(multicore,workers=numcores)
  
  ptime <- 1:length(TSCANorder)
  pval <- future_apply(data[, TSCANorder], 1, function(x) {
    if (sum(x) == 0) {
      1
    }
    else {
      model <- mgcv::gam(x ~ s(ptime, k = 3))
      pchisq(model$null.deviance - model$deviance, model$df.null - 
               model$df.residual, lower.tail = F)
    }
  })
  qval <- p.adjust(pval, method = "fdr")
  data.frame(pval = pval, qval = qval)
}

# This step could take a while....
bcalcul <- !file.exists(paste0('data_aux/09_',sname,'_diffTest','.Rdata'))
if(bcalcul){
  cat('Calculating difftest:\n')
  dtest  <- pdifftest(as.matrix(assay(ssce,"logcounts")[tfs,]),1:ncol(ssce),numcores=40)
  # * save ----
  save(dtest,file=paste0('data/09_',sname,'_diffTest','.Rdata'))
}else{
  (load(paste0('data_aux/09_',sname,'_diffTest','.Rdata')))
} 

# (2.1) keep significatives ----
qval <- 0.05
goi  <- rownames(dtest)[dtest$qval < qval]  # 478 genes of interest


# (3) Smoothing transcriptional profiles ----
# ____nfilt parameter ----
nfilt <- 601 
f21   <- rep(1/nfilt,nfilt)
x1    <- as.matrix(assay(ssce,'logcounts')[goi,])
xx1   <- t(apply(x1,1,function(y){stats::filter(y,f21,sides=2)}))
xx1[is.na(xx1)]<-0  
colnames(xx1)  <- colnames(ssce)

# (4) Inter profile distances ----
dcor <- 1-(1+cor(t(xx1)))*0.5

# (5) Cascade ----
# (5.0.1) disregard border time points (because of smoothing) ----
nout   <- (nfilt-1)/2
icells <- (nout+1):(ncol(xx1)-(nout+1))

# (5.0.2) Standarize ----
aux           <- aggregate(t(xx1[,icells,drop=FALSE]),
                           list(colData(ssce)[colnames(xx1)[icells],label]),
                           mean)
rownames(aux) <- aux[,1]
meanXpression <- t(aux[,-1])
zz <- xx1[,icells,drop=FALSE] * 1/apply(meanXpression[,,drop=FALSE],1,max)

# (5.1) check for ON/OFF switching behavior----

# ____minCorValue       : minimum correlation value wrt switch on/off template
# ____minExpressionLevel: minimum cluster expression level
# ____minChangeLevel    : minimum change expected between quartile80 - quartile20
p <- list()
p$nfilt              <- nfilt
p$qval               <- qval
p$minCorValue        <- 0.65
p$minExpressionLevel <- NA     # criterio gok1
p$minChangeLevel     <- NA     # criterio gok2 
p$addCorThreshold    <- 0.9
p$bAddMaxCor         <- FALSE  # add genes presenting template correlations above p$addCorThreshold value
p$bslope             <- TRUE   # criterio gok4: diff(range(y)) > (0.5*quantile(y,0.1))
p$brankcor           <- TRUE   # criterio gok5: de spearman + expression >0.5 +   diff(range(y)) > (0.5*quantile(y,0.1)) 
p$bFilterPre         <- TRUE   # deal with TFs presenting more than one switch-off pattern along pseudotime (keep only first transition)
p$bFilterPost        <- TRUE   # deal with TFs presenting more than one switch-off pattern along pseudotime (keep only last transition)

lres <- list()
for(groupOnOff in c('ON','OFF')){  
  xx11 <- xx1[,icells]
  
  
  geneSlopeMat <- c()
  lonoff <- list()
  for(iClusA in 1:length(clusoi)){
    gok   <- c()
    if(groupOnOff=='ON'){  # if looking for switch-ON, start from clusterA and consider 'subsequent' clusters
      rangeClusB <- (iClusA):length(clusoi)
    }else{                # if looking for switch-OFF, start from clusterA and consider 'previous' clusters
      rangeClusB <- (iClusA):1
    }
    maxCorReceived <- c()
    dfaux          <- c()
    dfMaxCor       <- c()
    for(iClusB in rangeClusB){
      if(iClusB > max(rangeClusB)) next
      cluster <- clusoi[c(iClusA:iClusB)] 
      
      template <- rep(0,length=length(icells))
      cclus <- as.character(colData(ssce)[colnames(xx11),label])[icells]
      cclus[cclus%in%c('NB1','NB2')] <- 'NB'
      
      template[cclus%in%cluster] <- 1
      
      a     <- cor(template,t(xx11))
      
      dfaux <- cbind(dfaux,a[1,])
      
      
      gok   <- unique(c(gok,colnames(a)[a>p$minCorValue]))
      
      if(length(maxCorReceived)==0){
        maxCorReceived <- a[1,]
      }else{
        maxCorReceived <- pmax(maxCorReceived,a[1,])
      }
    }
    
    colnames(dfaux) <- clusoi[rangeClusB]
    a<-t(apply(dfaux,1,function(x){
      ii <- which.max(x)
      c(clusoi[iClusA],colnames(dfaux)[ii],x[ii])
    }))  
    dfMaxCor <- data.frame(a[,1],a[,2],as.numeric(a[,3]))
    
    gok <- gok[!is.na(gok)]
    
    genes<-gok1<-gok2<-gok3<-gok4<-gok5<-filtersOK <- c()
    if(length(gok)>0){ # depuracion de seleccion
      
      # (5.1.1) me cercioro que efectivamente haya expresion en el luster iClusA ----
      if(!is.na(p$minExpressionLevel)){
        # (podria ser, si el cluster es muy chico, que domine la senial global)
        
        zz1  <- apply(zz[gok,cclus%in%clusoi[iClusA],drop=FALSE],1,mean)
        gok1 <- names(zz1)[zz1>p$minExpressionLevel]
        filtersOK <- c(filtersOK,'gok1')
      }
      
      # (5.1.2) alternativamente me fijo si se produjo un encendido relevante en iClusA ----
      # (cambio en x' de mas de p$minChangeLevel)
      if(!is.na(p$minChangeLevel)){   
        zz2 <- apply(zz[gok,cclus%in%clusoi[iClusA],drop=FALSE],1,function(x){
          diff(quantile(x,probs=c(0.2,0.8))) 
        })
        gok2 <- names(zz2)[zz2>p$minChangeLevel]
        filtersOK <- c(filtersOK,'gok2')
      }
      
      
      if(p$bslope){
        if(length(gok)>0){
          auxClus <- clusoi[iClusA]
          #if(auxClus%in%'NB') auxClus <- c('NB1','NB2')
          bcell <- colData(ssce)[colnames(zz),label]%in%auxClus
          aux   <- apply(zz[gok,bcell,drop=FALSE],1,
                         function(y){
                           x <- seq_along(y)
                           return(diff(range(y)) > (0.5*quantile(y,0.1)) )
                         })
          gok4 <- gok[aux]
          filtersOK <- c(filtersOK,'gok4')
        }
      }
      
      # (5.1.3) genes Spearman > 0.9 +diff(range(y)) > (0.5*quantile(y,0.1))   ----
      # para ver monoticidad
      if(p$brankcor){  
        auxClus <- clusoi[iClusA]
        #if(auxClus%in%'NB') auxClus <- c('NB1','NB2')
        bcell <- colData(ssce)[colnames(zz),label]%in%auxClus
        if(FALSE){
          bcell <- rep(TRUE,ncol(zz))
          if(groupOnOff=='ON'){
            auxSpearman <- aggregate(t(zz[gok,bcell,drop=FALSE]),
                                     list(colData(ssce)[colnames(zz)[bcell],label]),
                                     function(y){
                                       x <- seq_along(y)
                                       return(cor(x,y,method='spearman') > 0.9 &
                                                diff(range(y)) > (0.5*quantile(y,0.1)) )
                                     })
          }else{
            auxSpearman <- aggregate(t(zz[gok,bcell,drop=FALSE]),
                                     list(colData(ssce)[colnames(zz)[bcell],label]),
                                     function(y){
                                       x <- seq_along(y)
                                       return(-cor(x,y,method='spearman') > 0.9 &
                                                diff(range(y)) > (0.5*quantile(y,0.1)) )
                                     })
          }
          
          rownames(auxSpearman)<- auxSpearman[,1]
          auxSpearman          <- auxSpearman[,-1,drop=FALSE]
          # considero el NB<-max(NB1,NB2)
          if('NB2'%in%clusoi){
            auxSpearman['NB1',]   <- apply(auxSpearman[c('NB1','NB2'),],2,function(x){max(x)})
            auxSpearman <- auxSpearman[-which(clusoi%in%'NB2'),,drop=FALSE]
          }
          rownames(auxSpearman)[which(clusoi%in%'NB1')] <- 'NB'
          
          #expresion media en cluster de interes
          zz1  <- apply(zz[gok,bcell,drop=FALSE],1,mean)
          #gok5 <- gok[auxSpearman>0 & zz1[gok]>0.5]
          gok5 <- gok[auxSpearman[clusoi[iClusA],]>0 & zz1[gok]>0.5]
          
        }else{
          auxClus <- clusoi[iClusA]
          #if(auxClus%in%'NB') auxClus <- c('NB1','NB2')
          bcell <- colData(ssce)[colnames(zz),label]%in%auxClus
          
          if(groupOnOff=='ON'){
            auxSpearman <- apply(zz[gok,bcell,drop=FALSE],1,function(y){
              x <- seq_along(y)
              return(cor(x,y,method='spearman') > 0.9 & diff(range(y)) > (0.5*quantile(y,0.1)))
            })
          }else{
            auxSpearman <- apply(zz[gok,bcell,drop=FALSE],1,function(y){
              x <- seq_along(y)
              return(-cor(x,y,method='spearman') > 0.9 & diff(range(y)) > (0.5*quantile(y,0.1)))
            })
          }
        }
        
        
        #expresion media en cluster de interes
        zz1  <- apply(zz[gok,bcell,drop=FALSE],1,mean)
        #gok5 <- gok[auxSpearman>0 & zz1[gok]>0.5]
        gok5 <- gok[auxSpearman[names(zz1)]>0 & zz1>0.2]
        
        if(p$brankcor>0){
          #genes <- unique(c(genes,gok5))
          filtersOK <- c(filtersOK,'gok5')
        }
      }
      
      
    }
    
    
    
    if(FALSE){
      ii<-14
      plot(template,main=paste(genes[ii],signif(a[1,genes[ii]],2)))
      points(xx11[genes[ii],]/max(xx11[genes[ii],]))
    }
    
    genes <- c()
    for(i in seq_along(filtersOK)){
      if(i==1){
        genes <- get(filtersOK[i])
      }else{
        genes <- intersect(genes,get(filtersOK[i]))
      }
    }

    aa <- names(maxCorReceived)[maxCorReceived>p$addCorThreshold]
    aa <- aa[!aa%in%genes]
    
    if(p$bAddMaxCor) genes <- unique(c(names(maxCorReceived)[maxCorReceived>p$addCorThreshold],genes))
    
    lonoff[[clusoi[iClusA]]] <- genes
  }
  
  
  # (5.2) depurado de listas ----
  if(groupOnOff=='ON'){
    if(p$bFilterPost){
      if(FALSE){
        # saco genes detectados en un estadio
        # si los mismos aparecen prendidos en estadios inmediatamente anteriores
        for(ilon in length(lonoff):2){
          iout <- which(lonoff[[ilon]]%in%lonoff[[ilon-1]])
          if(length(iout)==0) next
          lonoff[[ilon]] <- lonoff[[ilon]][-iout]
        }
      }else{
        # si hay mas de una detecttion me quedo con el primero
        laux   <- reverseSplit(lonoff)
        laux   <- lapply(laux,function(x){x[which.min(match(x,clusoi))]})
        lonoff <- reverseSplit(laux)[clusoi]
      }
    }
    lres[['ON']] <- lonoff
  }else{
    if(p$bFilterPost){
     if(FALSE){
       # saco genes detectados en un estadio
       # si los mismos aparecen prendidos en estadios inmediatamente posteriores
       for(ilon in 1:(length(lonoff)-1)){
        iout <- which(lonoff[[ilon]]%in%lonoff[[ilon+1]])
        if(length(iout)==0) next
        lonoff[[ilon]] <- lonoff[[ilon]][-iout]
      }
     }else{
       # si hay mas de una detecttion me quedo con el ultimo
       laux   <- reverseSplit(lonoff)
       laux   <- lapply(laux,function(x){x[which.max(match(x,clusoi))]})
       lonoff <- reverseSplit(laux)[clusoi]
     }
    }
    lres[['OFF']] <- lonoff
  }
}

# (6) Heatmap ----
on  <- factor(unlist(reverseSplit(lres$ON)) ,levels=clusoi)
off <- factor(unlist(reverseSplit(lres$OFF)),levels=clusoi)
u   <- unique(c(unlist(lres$OFF),unlist(lres$ON)))

# (6.1)  define temporal classes to organize heatmap rows based on on-off patterns ----
df  <- matrix(NA,ncol=2,nrow=length(u))
rownames(df) <- u
df[names(on),1]  <- on
df[names(off),2] <- off
df <- data.frame(df)

dfs     <- df[do.call(order,df),]
pattern <- apply(dfs,1,paste0,collapse='_')
table(pattern)
pattern
# 1_NA  2_3  2_5  2_6  2_7 2_NA 3_NA  4_6  4_7 4_NA  5_6 5_NA 6_NA 7_NA NA_1 NA_2 NA_3 NA_5 NA_6 NA_7 
#    3    2    4    5    1   32    2    1    2    4    1    8   12    1   12   53   16   10    2    2 
lclass <- list()
lclass[[1]]<-names(pattern)[pattern%in%c('NA_1','NA_2','1_NA')]
lclass[[2]]<-names(pattern)[pattern%in%c('NA_3','2_3','2_4','1_5')]
lclass[[3]]<-names(pattern)[pattern%in%c('2_5','2_6','2_7','2_NA','3_NA','NA_5','NA_6','4_6','4_7')]
lclass[[4]]<-names(pattern)[pattern%in%c('4_NA','5_NA','NA_7')]
lclass[[5]]<-names(pattern)[pattern%in%c('6_NA','7_NA')]
names(lclass)<-as.character(1:length(lclass))

# (6.2) intra-class hierarchical clustering ----
ggenes <- c()
for(i in 1:length(lclass)){
  genes <- lclass[[i]]
  if(!exists('dcor')) dcor <- 1-(1+cor(t(xx1[genes,])))*0.5
  
  if(length(genes)>1){
    hc   <- hclust(as.dist(dcor[genes,genes]),method='ward.D2')
    
    dd <- dendsort(as.dendrogram(hc),type='average')
    hc  <- as.hclust(dd)
    
    ggenes <- c(ggenes,hc$labels[hc$order])
  }else{
    ggenes <- c(ggenes,genes)
  }
}
genes     <- ggenes
geneClass <- rep(names(lclass),unlist(lapply(lclass,length))) 

ncolouts <- ((nfilt-1)/2+1)
icols    <-  ncolouts:(ncol(xx1)-ncolouts)

# (6.3) nuclei subsampling for visualization ----
set.seed(12337)
icols <- sort(sample(icols,length(icols)*0.2),decreasing=FALSE)
if(FALSE) plot(colData(ssce)[,'pt_umap'][icols])

# (6.4) standarization ----
gxprofile <- t(apply(xx1[genes,icols],1,function(x){(x-min(x))/(max(x)-min(x))}))


# (6.5) heatmap annotation ----
annotation_row = data.frame(MaxExp = apply(xx1[rownames(gxprofile),],1,max))
annotation_col = data.frame(ptime   = rank(colData(ssce)[colnames(xx1)[icols],'pt_umap']),
                            cluster = colData(ssce)[colnames(xx1)[icols],label])
newCols       <- c("F5BA6C","EB7604","C6ADCF","A75528","DC2A85","2E9737","DEA709")
newCols       <- paste0("#",newCols)
names(newCols)<- c("RGL","NPC","NB","GCimm1","GCimm2","GCyoung","GCmat1")
annoCol        <- list(cluster = newCols)

ph <- pheatmap(gxprofile,cluster_rows=FALSE,cluster_cols = FALSE,
                 show_colnames = FALSE,
                 use_raster    = TRUE,
                 row_gap      = unit(0.5,'mm'),
                 column_split = factor(colData(ssce)[colnames(gxprofile),label],levels=clusoi),
                 column_gap = unit(0.5,'mm'),
                 annotation_row    = annotation_row,
                 annotation_col    = annotation_col, 
                 annotation_colors = annoCol)



if(bpdf){
  fname <- paste0('data_figs/09_',sname,'.pdf')
  pdf(fname,width=10,height=max(4,0.2*nrow(gxprofile)))
  draw(ph)
  dev.off()
}else{
  draw(ph)
}





