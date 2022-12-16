
p.cores=0.35
library(dplyr)
#doParallel::registerDoParallel(round(parallel::detectCores()*p.cores,0))
doParallel::registerDoParallel(34)
foreach::getDoParWorkers()
library(foreach)

args=commandArgs(trailingOnly = TRUE)

t1<-read.table(args[1],sep=",",header=TRUE)
dim(t1)
names(t1)<-c("Genome_1","Genome_2","ANI")
t1$Genome_1<-gsub(".gb.fasta","",t1$Genome_1)
t1$Genome_2<-gsub(".gb.fasta","",t1$Genome_2)

length(unique(t1$Genome_1))
length(unique(t1$Genome_2))
length(unique(union(t1$Genome_1,t1$Genome_2)))


print("Now multiplying by shared orthogroups and then normalizing")
if(!file.exists("temp.rda")){
  print("c1")
  m1<-matrix(nrow=length(unique(union(t1$Genome_1,t1$Genome_2)))
             ,ncol=length(unique(union(t1$Genome_1,t1$Genome_2))))
  print("matrix dimensions")
  
  print(dim(m1))
  m1[,]<-0
  
  rownames(m1)<-unique(union(t1$Genome_1,t1$Genome_2))
  colnames(m1)<-unique(union(t1$Genome_1,t1$Genome_2))
  m1[t1$Genome_1[3],t1$Genome_2[4]]
  a=0
  #foreach(i=1:round(length(t1$Genome_1)/1000,0),.verbose=TRUE) %dopar%  {
  for(i in seq_along(t1$Genome_1)){
    print("i")
    print(i)
    if(t1$ANI[i]>0) try(m1[t1$Genome_1[i],t1$Genome_2[i]]<-t1$ANI[i])
    if(is.na(m1[t1$Genome_1[i],t1$Genome_2[i]])) m1[t1$Genome_1[i],t1$Genome_2[i]]<-0
    a=a+1
  }
  print("matrix max")
  
  print(max(m1))
  m1=(m1/max(m1))
}

print("c2")
max(m1)
median(m1)

(sum(m1>0.1)/(ncol(m1)*ncol(m1)))

print(m1[1:10,1:10])

mean(m1[,])
if(!file.exists(paste0(args[2],"hclust.rda"))){
  mean(m1[,])
  d<-dist(m1)
  print("check")
  hc <- hclust(d, method = "average")
  pdf(paste0(args[2],"hclust.pdf"),height=10,width=165)
  plot(hc)
  dev.off()
  
  save(hc,file="hclust.rda")
}
stop()

load(paste0(args[2],"hclust.rda"))
print("c3")

dmat_nr<-m1
genome=unique(union(t1$Genome_1,t1$Genome_2))
genome_clusters <- data.frame(genome)

if (!file.exists("./cache/df_cluster_strict.rda")) {
  print("Optimising number of clusters, using the non-redundant data set.")
  
  ##CHANGES: all tails_small->tails_nr, all dmat->dmat_nr
  
  h <- hc #average instead of the default complete
  
  homogeneity <- function(mat, threshold = 0.7) {
    mat <- as.matrix(mat)
    if (nrow(mat) == 1) {
      ho <- NA
    }
    if (nrow(mat) >1) {
      hit <- sum(mat >= threshold)-nrow(mat)
      all <- nrow(mat)^2-nrow(mat)
      ho <- round(hit/all, 3)
    }
    return(ho)
  }
  
  completeness <- function(dmat_nr, index, threshold = 0.7) {
    if (length(index) == 0){
      stop("cluster must contain at least one element. Please provide indices.")
    }
    imat <- dmat_nr[index, index]
    emat <- dmat_nr[index, -index]
    amat <- dmat_nr[index, ]
    hit <- sum(emat >= threshold)
    all <- sum(amat >= threshold)-length(index)-(sum(imat >= threshold)-length(index))/2
    co <- round(1-hit/all, 3)
    return(co)
  }
  
  df_cluster <- data.frame(
    nclust = seq(from = 1, to = nrow(genome_clusters), by = 1)
  )
  
  df_cluster$shoco <- parallel::mclapply(df_cluster$nclust, function(x){
    ho = vector()
    co = vector()
    clusters <-  cutree(tree = h, k = x)
    genome_clusters$cluster <- paste0("cluster_",clusters)
    for (j in unique(genome_clusters$cluster)) {
      index <- which(genome_clusters$cluster == j)
      smat <- dmat_nr[index, index]
      hoexpr=length(index)*homogeneity(smat, threshold = 0.9)
      #print(length(index))
      #if (hoexpr < 0.01& length(index)>1) (hoexpr= -1)
      #if (hoexpr > 0.9 & length(index)>1) (hoexpr= length(index)**2)
      
      ho <- c(ho, hoexpr)
      co <- c(co, length(index)*completeness(dmat_nr, index, threshold = 0.7))
    }
    #print(ho)
    shoco <- sum(ho, na.rm = TRUE)
    return(shoco)
  }, mc.cores = 34 )
  
  #
  save(df_cluster, file = "./cache/df_cluster_strict.rda")
}

load("./cache/df_cluster_strict.rda")

df_cluster$shoco=unlist(df_cluster$shoco)
#with(df_cluster, plot(nclust, shoco))
best_cut <- which(df_cluster$shoco == max(df_cluster$shoco))

clusters <-  cutree(tree = hc, k = best_cut)

genome_clusters$cluster <- paste0("cluster_",sprintf("%03d", clusters))
load("
table(genome_clusters$cluster)
save(genome_clusters,file="cache/genome_clusters.rda")
#file.remove("cache/df_cluster_strict.rda")
