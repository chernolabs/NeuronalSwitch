#From limma
strsplit2 <- function (x, split, ...) 
{
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i]) 
      out[i, 1:nc[i]] <- s[[i]]
  }
  out
}


# Funcion para construir matriz de distancias.
# Es posible trabajar con un subset de las columnas de sce  [nodes1] (se devuelve dist:nodes1 x nodes1)
# o especificar un segundo grupo de nodos (nodes2) y estimar la dist: nodes1 x nodes2
#
# Input:
#   sce: swingleCellExperiment object
#   nodes1: vector de nombres (colnames(ssce)) para especificar nodos query
#   nodes2: vector de nombres (colnames(ssce)) para especificar nodos de refs

buildDistMatrix <- function(sce,nodes1=NULL,nodes2=NULL,
                            use.reducedDim=c(NA,"PCA")[2],num.PCA=NA,mode=c("pearson","euclidean")[1],
                            bverbose=FALSE){
  
  # Estimo matriz de distancias
  t1<-Sys.time()
  
  if(!is.null(nodes2) & is.null(nodes1)){
    warning('nodes2 cant be specified without nodes1\n')
    return(NA)
  }
  
  berror<-FALSE
  if(!is.null(nodes1)){
    if(any(!nodes1%in%colnames(sce))) berror<-TRUE
    if(!is.null(nodes2) & any(!nodes2%in%colnames(sce))) berror<-TRUE
  }
  if(berror) warning('No nodes1 or nodes2 found in colnames(sce)')
  
  if(is.na(use.reducedDim)){
    x <- t(logcounts(sce))
  }else{
    if(!use.reducedDim%in%reducedDimNames(sce)){
      stop("Non recognized reducedDim.\n")
    }else{
      x <- reducedDim(sce,use.reducedDim)
      if(!is.na(num.PCA)){
        x <- x[,1:(min(num.PCA,ncol(x)))]
      }
    }  
  }
  
  #me quedo con sets de interes
  nodes <- rownames(x)
  if(!is.null(nodes1)) nodes <- nodes1
  if(!is.null(nodes2)) nodes <- c(nodes,nodes2)
  x <- x[nodes,]
  
  if(mode=="pearson"){
    if(class(x)[1]=="dgCMatrix"){
      dist.mat <- 1- (sparse.cor4(t(x))$cor+1)/2
    }else{
      dist.mat <- 1- (cor(t(x))+1)/2
    }
  }else{
    #plot(pca$sdev);abline(v=numPCA)
    dist.mat <- as.matrix(dist(x))
  }
  if(is.null(rownames(dist.mat))){
    colnames(dist.mat)<-rownames(dist.mat)<-colnames(sce)
  }
  
  #me quedo con el subset de correlaciones necesario
  if(!is.null(nodes1)){
    dist.mat <- dist.mat[nodes1,]
    if(!is.null(nodes2)){
      dist.mat <- dist.mat[,nodes2]
    }else{
      dist.mat <- dist.mat[,nodes1]
    }
  }
  t2<-Sys.time()
  
  if(bverbose) cat('Build distance Matrix:', t2-t1, '\n')  
  return(dist.mat)
}

# Funcion para construir redes/listas de k primeros vecinos mutuos.
# Es posible trabajar con un subset de las columnas de sce  [nodes1]
# o identificar enlaces mutuos entre grupos disjuntos de nodos (nodes1 y nodes2)
# En este ultimo caso se estiman la lista de los primeros vecinos mutuos (del set nodes2)
# para los nodos del set nodes1.
#
# Input:
#   sce: swingleCellExperiment object
#   dist.mat: matriz de distancias pre-calculada. Puede ser no-cuadrada. Si se especifica se toman
#             rownames como nodes1 y colnames como nodes2
#   nodes1: vector de nombres (colnames(ssce)) para especificar nodos query
#   nodes2: vector de nombres (colnames(ssce)) para especificar nodos de refs
#   mutualK: numero de prmeros vecinos
#   bvalue : si TRUE la lista de primero svecinos se construye como lista de distancias nombrada
#            si FALSE se cosntruye solo con los nombres de los primeros vecinos
#   bOutGraph: si FALSE decvuelve lista de primeros vecinos, 
#              si TRUE devueve grafo de primeros vecinos mutuos
#   bWeighted: devuelve el grafo pesado por similaridad en el espacio reducido utilizado
#   bCalculateNodeSims: si TRUE, calcula similaridad de Jaccard, Dice e invLogWeight
#   
buildMKNN     <-function(sce=NULL,mutualK=30,
                         dist.mat=NULL, nodes1=NULL,nodes2=NULL,bvalue=FALSE,
                         use.reducedDim=c(NA,"PCA")[2],num.PCA=NA,mode=c("pearson","euclidean")[1],
                         bOutGaph=TRUE,bWeighted=FALSE,bCalculateNodeSims=FALSE,bverbose=FALSE){
  
  if(is.null(sce) & is.null(dist.mat)){
    warning('Either sce or dist.mat should be specified.\n')
  }
  
  # Estimo matriz de distancias
  if(is.null(dist.mat)){
    dist.mat <- buildDistMatrix(sce,nodes1,nodes2,use.reducedDim,num.PCA,mode,bverbose)
    if(FALSE){
      t1<-Sys.time()
      
      if(!is.null(nodes2) & is.null(nodes1)){
        warning('nodes2 cant be specified without nodes1\n')
        return(NA)
      }
      
      berror<-FALSE
      if(!is.null(nodes1)){
        if(any(!nodes1%in%colnames(sce))) berror<-TRUE
        if(!is.null(nodes2) & any(!nodes2%in%colnames(sce))) berror<-TRUE
      }
      if(berror) warning('No nodes1 or nodes2 found in colnames(sce)')
      
      if(is.na(use.reducedDim)){
        x <- t(logcounts(sce))
      }else{
        if(!use.reducedDim%in%reducedDimNames(sce)){
          stop("Non recognized reducedDim.\n")
        }else{
          x <- reducedDim(sce,use.reducedDim)
          if(!is.na(num.PCA)){
            x <- x[,1:(min(num.PCA,ncol(x)))]
          }
        }  
      }
      
      #me quedo con sets de interes
      nodes <- rownames(x)
      if(!is.null(nodes1)) nodes <- nodes1
      if(!is.null(nodes2)) nodes <- c(nodes,nodes2)
      x <- x[nodes,]
      
      if(mode=="pearson"){
        if(class(x)[1]=="dgCMatrix"){
          dist.mat <- 1- (sparse.cor4(t(x))$cor+1)/2
        }else{
          dist.mat <- 1- (cor(t(x))+1)/2
        }
      }else{
        #plot(pca$sdev);abline(v=numPCA)
        dist.mat <- as.matrix(dist(x))
      }
      if(is.null(rownames(dist.mat))){
        colnames(dist.mat)<-rownames(dist.mat)<-colnames(sce)
      }
      
      #me quedo con el subset de correlaciones necesario
      if(!is.null(nodes1)){
        dist.mat <- dist.mat[nodes1,]
        if(!is.null(nodes2)){
          dist.mat <- dist.mat[,nodes2]
        }else{
          dist.mat <- dist.mat[,nodes1]
        }
      }
      
      if(bverbose) cat('Build distance Matrix:', difftime(Sys.time(),t1,units = 'secs')[[1]],' secs.\n')  
    }
  }else{
    nodes1 <- rownames(dist.mat)
    nodes2 <- colnames(dist.mat)
  }
  
  # KNN
  t1<-t0<-Sys.time()
  iknn <- 2:(mutualK+1)
  if(is.null(dist.mat) & !is.null(nodes2)) iknn <- iknn-1
  knn <- apply(dist.mat,1,function(x){
    #saux <- names(sort(rank(x),decreasing=FALSE))[2:(mutualK+1)]
    saux <- names(sort(x,decreasing=FALSE)[iknn])
    return(saux) #devuelve el nombre de k_primeros_vecinos
  })
  
  if(bverbose) cat('...Build nodes1 knn matrix',difftime(Sys.time(),t1,units = 'secs')[[1]],' secs.\n')
  
  #mutualKNN
  t1<-Sys.time()
  lmknn<-list()
  if(is.null(nodes2)){
    if(bvalue){
      for(i in 1:ncol(knn)){
        vecinos <- knn[,i]
        vecinos <- vecinos[unlist(apply(knn[,vecinos],2,function(x){colnames(knn)[i]%in%x}))]
        xx        <- dist.mat[colnames(knn)[i],vecinos]
        lmknn[[colnames(knn)[i]]] <- xx
      }
    }else{
      for(i in 1:ncol(knn)){
        vecinos <- knn[,i]
        vecinos <- vecinos[unlist(apply(knn[,vecinos],2,function(x){colnames(knn)[i]%in%x}))]
        #xx        <- dist.mat[colnames(knn)[i],vecinos]
        lmknn[[colnames(knn)[i]]] <- vecinos #xx
      }
    }
  }else{
    knn2 <- apply(dist.mat,2,function(x){
      #saux <- names(sort(rank(x),decreasing=FALSE))[2:(mutualK+1)]
      saux <- names(sort(x,decreasing=FALSE)[iknn])
      return(saux) #devuelve el nombre de k_primeros_vecinos
    })
    t11<-Sys.time()
    
    #Manera vieja de calcularlo
    if(FALSE){
      t1<-Sys.time()
      for(i in 1:ncol(knn)){
        vecinos <- knn[,i]
        vecinos   <- vecinos[unlist(apply(knn2[,vecinos],2,function(x){colnames(knn)[i]%in%x}))]
        #xx        <- dist.mat[colnames(knn)[i],vecinos]
        #names(xx) <- vecinos
        #lmknn[[colnames(knn)[i]]] <- xx 
        lmknn[[colnames(knn)[i]]] <- vecinos
      }
      
      if(bverbose) cat(difftime(Sys.time(),t1,units = 'secs')[[1]],' secs.\n')
    }else{  #Manera nueva
      lknn       <- as.list(data.frame(knn2))
      names(lknn)<- sub('.',':',colnames(knn2),fixed=TRUE)
      lNeighOf2 <- Biobase::reverseSplit(lknn)
      lmknn <- list()
      if(bvalue){
        for(i in seq_along(lNeighOf2)){
          vecinos <- intersect(lNeighOf2[[i]],knn[,names(lNeighOf2)[i]])
          #lmknn[[names(lNeighOf2)[i]]] <- vecinos
          xx        <- dist.mat[names(lNeighOf2)[i],vecinos]
          names(xx) <- vecinos
          lmknn[[names(lNeighOf2)[i]]] <- xx 
        } 
      }else{
        for(i in seq_along(lNeighOf2)){
          vecinos <- intersect(lNeighOf2[[i]],knn[,names(lNeighOf2)[i]])
          lmknn[[names(lNeighOf2)[i]]] <- vecinos
        }
      }
    }
  }
  if(bverbose){
    if(is.null(nodes2)){
      cat('...Build list KMNN:', difftime(Sys.time(),t1,units = 'secs')[[1]], ' secs.\n')
    }else{
      cat('...Build nodes2 knn matrix:', difftime(t11,t1,units = 'secs')[[1]], ' secs.\n')
      cat('...Build list KMNN:', difftime(Sys.time(),t11,units = 'secs')[[1]], ' secs.\n')
    }
    cat('Total:',difftime(Sys.time(),t0,units = 'secs')[[1]],' secs.\n')
  }
  
  
  if(FALSE){
    table(unlist(lapply(lmknn,length)))
    hist(unlist(lapply(lmknn,length)),breaks=39)
    a<-unlist(lapply(names(lmknn),function(x){
      xx<-lmknn[[x]]
      if(length(xx)>0){
        return(dist.mat[x,xx[length(xx)]])
      }else{
        return(-1)
      }
    }))
  }
  
  if(!bOutGaph){
    return(lmknn)
  }
  
  #Armo edgelist
  t1<-Sys.time()
  el <- c()
  for(i in seq_along(lmknn)){
    if(bvalue){
      a<-names(lmknn[[i]])
    }else{
      a<-lmknn[[i]]
    }
    if(length(a)>0)
      el<-rbind(el,cbind(rep(names(lmknn)[i],length(a)),a))
  }
  a<-t(apply(el,1,sort))
  el<-a[!duplicated(a),]
  
  g <- igraph::graph_from_edgelist(el,directed=FALSE)
  g <- igraph::simplify(g)
  
  if(bCalculateNodeSims){
    t11<-Sys.time()
    simCellJ <- igraph::similarity(g,method="jaccard")
    simCellD <- igraph::similarity(g,method="dice")
    simCellW <- igraph::similarity(g,method = "invlogweighted")
    colnames(simCellJ)<-colnames(simCellD)<-colnames(simCellW)<-igraph::V(g)$name
    rownames(simCellJ)<-rownames(simCellD)<-rownames(simCellW)<-igraph::V(g)$name
    
    #para asignar las ismilitudes a los enlaces
    eattr<-t(apply(el,1,function(x){
      c(x,simCellJ[x[1],x[2]],simCellD[x[1],x[2]],simCellW[x[1],x[2]])
    }))    
    colnames(eattr) <- c("Source","Target","simJ","simD","simW")
    rownames(eattr)  <- apply(eattr[,1:2],1,function(x){paste0(sort(x),collapse=":")})
    
    #pongo las similutdes calculadas como atrubuto de edges del grafo
    saux <- apply(igraph::ends(g,igraph::E(g)),1,function(x){paste0(sort(x),collapse=":")})
    igraph::E(g)$simJ <- as.numeric(eattr[saux,"simJ"])
    igraph::E(g)$simD <- as.numeric(eattr[saux,"simD"])
    igraph::E(g)$simW <- as.numeric(eattr[saux,"simW"])
    t22<-Sys.time()
    if(bverbose) cat('    Node similarities:', t22-t11, '\n')
  }
  if(bWeighted){
    t11<-Sys.time()
    a <- apply(igraph::get.edgelist(g),1,function(x){
      return(dist.mat[x[1],x[2]])
    })
    igraph::E(g)$weight <-  1/a
    t22<-Sys.time()
    if(bverbose) cat('    Edge weights assigned:', t22-t11, '\n')
  }
  if(is.na(use.reducedDim)){
    saux <- paste0("K",mutualK,"_Dim",ncol(dist.mat))
  }else{
    saux <- paste0("K",mutualK,"_reduced",use.reducedDim,"_Dim",ncol(dist.mat),"_mode",mode)
  }
  g <- igraph::set.graph.attribute(g,"params",saux)
  
  if(bverbose) cat('Build graph:', difftime(Sys.time(),t1,units = 'secs')[[1]], ' secs.\n')
  if(bverbose) cat('TOTAL TIME ELAPSED:', difftime(Sys.time(),t0,units = 'secs')[[1]], ' secs.\n')
  return(g) 
}

# Calculo de clusterizacion sobre el grafo g. El resultado queda almacenado en colData(sce)
# Inputs:
#   sce: single cell object (usado para calcular previamente el grafo g)
#   g: grafo igraph
#   lAlgorithm: es una lista de funciones igraph de deteccion de comunidades
#       lAlgorithms <- list(louvain=igraph::cluster_louvain, 
#                           infomap=igraph::cluster_infomap,
#                           labProp=igraph::cluster_label_prop,
#                           fgreedy=igraph::cluster_fast_greedy,
#                           walktrap=igraph::cluster_walktrap)
#   prefix: prefijo para armar e nombre de las columnas de colData(sce)
# Outpus:
#   singleCellExperiment object con resultados de clusterizacion en colData(sce)
doClustering<-function(sce,g,lAlgorithms,prefix="label_"){
  lclus<-list()
  #layout(matrix(seq_along(lAlgorithms),2,2))
  for(i in seq_along(lAlgorithms)){
    lclus[[i]]<-lAlgorithms[[i]](g)
    #    (tt<-table(colData(sce)[,"cell_type"],membership(lclus[[i]])[rownames(colData(sce))]))
    #    (inf<-infoPartition(tt))
    #    plot(inf,as.numeric(table(lclus[[i]]$membership)),xlim=c(0,1),xlab="coherence",ylab="size",log="y",main=names(lAlgorithms)[i])
  }
  names(lclus)<-names(lAlgorithms)
  
  laux <- lapply(lclus,igraph::membership)
  aux <- matrix(unlist(laux),ncol=length(lclus))
  colnames(aux)<-names(lclus)
  rownames(aux)<-names(laux[[1]])
  
  labels          <- matrix(0,ncol=length(lclus),nrow=ncol(sce))
  rownames(labels)<- colnames(sce)
  labels[rownames(aux),]<-aux
  colnames(labels)<-paste0(prefix,colnames(aux))
  labels<-apply(labels,2,function(x){   #pongo los ceros como clusters de un elemento
    i0     <- which(x==0)
    newVal <- max(x)+seq_along(i0)
    x[i0]  <- newVal
    return(x)
  })  
  
  colData(sce) <- cbind(colData(sce),labels)
  cat("Clustering done\n")
  return(sce)
}

#Exporta el grafo gg a formato Gephi. colData de sce se exporta como atributos de nodos 
# Inputs:
#   sce: singleCellExperiment
#   gg : igraph network to be exported
#   file.name: nombre para el file a exportar
exportGraphToGephi<-function(sce,gg,file.name="foo"){
  aux  <- edge_attr(gg)
  maux <- matrix(unlist(aux),ncol=length(aux))
  colnames(maux)<-names(aux)
  eattr    <- data.frame(ends(gg,E(gg)),maux)
  colnames(eattr)<-c("Source","Target",names(aux))
  
  #edges 
  write.table(eattr,sep=",",row.names=FALSE,col.names=TRUE,
              file=paste0(file.name,".adjList.csv"),quote=FALSE)
  
  el <- get.edgelist(gg)
  el <- t(apply(el,1,sort))
  
  #node attributes
  maux<-colData(sce)
  linkedNodes <- unique(c(el[,1],el[,2]))
  a<-maux[,1]%in%linkedNodes
  if(all(!a)) a<-rownames(maux)%in%linkedNodes
  maux<-maux[a,]
  #colnames(maux)[1]<-"Id"
  maux <- cbind(Id=rownames(maux),maux)
  
  aux <- grep("percent_top",colnames(maux))
  if(length(aux)>0) maux <- maux[,-aux]
  
  #pongo los atributos de nodo en el grafo
  for(i in 2:(ncol(maux))){
    if(class(maux[1,i])=="character"){
      gg <- set.vertex.attribute(gg,colnames(maux)[i],V(gg),maux[,i])
    }else{
      gg <- set.vertex.attribute(gg,colnames(maux)[i],V(gg),as.numeric(maux[,i]))
    }
  }
  
  write.table(maux,sep=",",row.names=FALSE,col.names=TRUE,file=paste0(file.name,".nodes.csv"),quote=FALSE)
  
  cat("Graph exported.\n")
}

