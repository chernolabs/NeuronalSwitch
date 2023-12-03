# Cluster and  density plots over KMNN graphs 
#
# input: 
#   data/06_DS1gdf_xy.txt: Layout computed using the ForceLayout2 algorithm in Gephi
#
# output:
#  data/06_DS1_MKNN_partition.pdf
#  data_figs/06_DS1_MKNN_cohort_*_.pdf

options(bitmapType = 'cairo')



# (0) Load data ----
data <- read.csv('data_aux/06_DS1gdf_xy.txt',header=TRUE,sep="\t")

yfactor <- 0.5
bpdf    <- TRUE

# (1) Clusters ----
# (1.1) Define custom colors
natCol <- c("F5BA6C","EB7604","C6ADCF","663E8E","A75528","DC2A85","2E9737","DEA709","1B77B5","AACC83","8FCCC4","A2C8DC","F8F4B5")
natCol <- paste0("#",natCol)
names(natCol)<- c("RGL","NPC","NB1","NB2","GCimm1","GCimm2","GCyoung","GCmat1","GCmat2","Astro","OPC","Oligo","Peri")

# (1.2) Select main path  
ipoints <- data$cluster %in% c("NB2","NB1","GCimm1","GCmat1","RGL","NPC","Astro","GCyoung","GCimm2","GCmat2")
df      <- cbind(x=data$x,y=data$y*yfactor,cluster=data$cluster)[ipoints,]

ccol  <- natCol[df[,'cluster']]
ccol1 <- paste0(ccol,"0D")
ccol2 <- ccol
ccol3 <- paste0(ccol,'0D')

if(bpdf) pdf(file='data_figs/06_DS1_MKNN_partition.pdf',width=7,height=7*748/996,compress=TRUE)
plot(df[,1],df[,2],pch=21,cex=.4,
     xlab="",ylab="",axes=FALSE,
     bg=ccol2,col=ccol3,lwd=.2)
if(bpdf)dev.off()


# (2) Cohorts   Kernel density estimators ----
library(ks)
week <- c("1w","2w","4w","8w")

w <- unlist(lapply(strsplit(rownames(data),"_"),function(x){x[2]}))
data$week <- w

layout(matrix(1:4,2,2,byrow=TRUE))
for(i in seq_along(week)){
  fname   <- paste0('data_figs/06_DS1_MKNN_cohort_',week[i],'.pdf')
  
  if(bpdf)pdf(file=fname,width=7,height=7*748/996)#650/1200,compress=TRUE)
  
  ccol    <- rep("#999999",nrow(data))
  ipoints <- which(data$week%in%week[i])
  df <- cbind(x=data$x,y=data$y)[ipoints,]
  
  H    <- Hpi.diag(x=df)
  fhat <- kde(x=df, H=H)
  
  par(bty='n')
  par(mar=c(1,1,1,1))
  plot(data$x,data$y,pch=20,cex=.4,col="#DDDDDD",xlab="",ylab="",axes=FALSE)
  plot(fhat, display="filled.contour",add=TRUE,alpha=1,cont=c(5,10,20,40,60,80,90,95),frame=FALSE)
  if(bpdf)dev.off()
  
}



