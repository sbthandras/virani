
library(tidyverse)
library(foreach)
library(doParallel)
library(grid)

args = commandArgs(trailingOnly=TRUE)


genomesdir="genomes/"
outputdir="results/"
threads=args[1]
genomesdir=args[2]
outputdir=args[3]
if(!dir.exists(genomesdir)) dir.create(genomesdir)
if(!dir.exists(outputdir)) dir.create(outputdir)


print(threads)
print(genomesdir)
print(outputdir)

cl <- makeCluster(as.numeric(threads))
print(cl)
registerDoParallel(cl)


fils<-list.files(genomesdir,pattern = "fasta")
print(fils)
anidf <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("seq1", "seq2", "ANI")
colnames(anidf) <- x
count=0

anidf = foreach(i=1:length(fils), .combine=rbind) %dopar% {
  totaldf <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("seq1", "seq2", "ANI")
  colnames(anidf) <- x
  print(fils[i])

  count=count+1
  print(count)
  #if(count==2) break()
  for(o in i:length(fils)){
    print(fils[o])
    if(file.info(paste0(genomesdir,fils[o]))$size>1000000) next
    
    system(paste0("rm ",outputdir,"tmptbl",i),ignore.stderr = TRUE)
    system(paste0('blastn  -subject ',genomesdir,fils[o],' -query ',genomesdir,fils[i],'    -out ',outputdir,'tmptbl',i,' -max_target_seqs 1    -outfmt "6 qseqid sseqid qstart qend sstart send qseq sseq evalue bitscore pident qlen slen" ')
           ,ignore.stderr = TRUE)
    
    if(file.exists(paste0(outputdir,"tmptbl",i))){
      if(file.info(paste0(outputdir,"tmptbl",i))$size==0) next
      
      tmptbl<-read.table(paste0(outputdir,"tmptbl",i),sep="\t",header=F)
      names(tmptbl)<-c("qseqid" ,"sseqid" ,"qstart", "qend", "sstart", 
"send","qseq","sseq", "evalue", "bitscore" ,"pident" ,"qlen", "slen")
      
      if(tmptbl$qlen[1]>tmptbl$slen[1]){
        tmptbl$alignlen<-0
        for(a in 1:length(tmptbl$qstart)){
          tmptbl$alignlen[a]<-length(tmptbl$qstart[a]:tmptbl$qend[a])
          if(tmptbl$pident[a]>100) tmptbl$pident[a]=100
        }
        
        tmptbl$biggerlen<-tmptbl$qlen
        
tmptbl$weightedpid<-tmptbl$pident*abs(tmptbl$alignlen/tmptbl$biggerlen)
        
        for(ii in seq_along(tmptbl$qstart)){
          
          for(xx in ii:length(tmptbl$qstart)){
            
print(length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$qstart[xx]:tmptbl$qend[xx])))
            print("alignlen")
            print(tmptbl$qlen[ii])
            print(tmptbl$qlen[xx])
            
            
if(length(intersect(tmptbl$qstart[ii]:tmptbl$qend[ii],tmptbl$qstart[xx]:tmptbl$qend[xx]))>=1){
              if(tmptbl$weightedpid[ii]>tmptbl$weightedpid[xx]) 
tmptbl$weightedpid[xx]<-tmptbl$weightedpid[xx]*
                  (
                    
(length(tmptbl$qstart[xx]:tmptbl$qend[xx])-length(intersect(
                      
tmptbl$qstart[ii]:tmptbl$qend[ii],tmptbl$qstart[xx]:tmptbl$qend[xx]))
                    )/
                      length(tmptbl$qstart[xx]:tmptbl$qend[xx]))
              if(tmptbl$weightedpid[ii]<tmptbl$weightedpid[xx]) 
tmptbl$weightedpid[ii]<-tmptbl$weightedpid[ii]*
                  
((length(tmptbl$qstart[ii]:tmptbl$qend[ii])-length(intersect(tmptbl$qstart[ii]:tmptbl$qend[ii],tmptbl$qstart[xx]:tmptbl$qend[xx])))/
                     length(tmptbl$qstart[ii]:tmptbl$qend[ii]))
            }
            
          }
        }
      } else{
        tmptbl$alignlen<-0
        for(a in 1:length(tmptbl$qstart)){
          tmptbl$alignlen[a]<-length(tmptbl$sstart[a]:tmptbl$send[a])
          if(tmptbl$pident[a]>100) tmptbl$pident[a]=100
          
        }
        tmptbl$biggerlen<-tmptbl$slen
        
tmptbl$weightedpid<-tmptbl$pident*abs(tmptbl$alignlen/tmptbl$biggerlen)
        for(ii in seq_along(tmptbl$sstart)){
          for(xx in ii:length(tmptbl$sstart)){
            
print(length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$sstart[xx]:tmptbl$send[xx])))
            print("alignlen")
            print(tmptbl$slen[ii])
            print(tmptbl$slen[xx])
            
            
if(length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$sstart[xx]:tmptbl$send[xx]))>=1){
              if(tmptbl$weightedpid[ii]>tmptbl$weightedpid[xx]) 
tmptbl$weightedpid[xx]<-tmptbl$weightedpid[xx]*(
                
(length(tmptbl$sstart[xx]:tmptbl$send[xx])-length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$sstart[xx]:tmptbl$send[xx])))/
                  length(tmptbl$sstart[xx]:tmptbl$send[xx])
              )
              if(tmptbl$weightedpid[ii]<tmptbl$weightedpid[xx]) 
tmptbl$weightedpid[ii]<-tmptbl$weightedpid[ii]*(
                
(length(tmptbl$sstart[ii]:tmptbl$send[ii])-length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$sstart[xx]:tmptbl$send[xx])))/
                  length(tmptbl$sstart[ii]:tmptbl$send[ii])
              )
            }
            
          }
        }
        
      }
      
      
      

      tmpdf <- tmptbl[1,c("qseqid","sseqid")]
      tmpdf$ANI<-sum(tmptbl$weightedpid)
      tmpdf$qseqid<-fils[o]
      tmpdf$sseqid<-fils[i]
      totaldf<-rbind(totaldf,tmpdf)
      
      
      
    }
  
  }
  totaldf
}


write.table(anidf,file=paste0(outputdir,"tempres.tsv"),quote=F,row.names=F,sep=",")

anim<-matrix(nrow=length(unique(anidf$qseqid)),ncol=length(unique(anidf$sseqid)))

row.names(anim)=unique(anidf$qseqid)
colnames(anim)=unique(anidf$qseqid)

for(i in seq_along(anidf$qseqid)){
  anim[anidf$qseqid[i],anidf$sseqid[i]]<-anidf$ANI[i]
  anim[anidf$sseqid[i],anidf$qseqid[i]]<-anidf$ANI[i]
  
  
}

#Personal cases
row.names(anim)<-gsub("vB_Aba._","",gsub(".fasta","",row.names(anim))) 
colnames(anim)<-gsub("vB_Aba._","",gsub(".fasta","",colnames(anim))) 

row.names(anim)<-gsub("\\_Hun","",row.names(anim)) 
colnames(anim)<-gsub("\\_Hun","",row.names(anim)) 

row.names(anim)<-gsub("\\_"," ",row.names(anim)) 
colnames(anim)<-gsub("\\_"," ",row.names(anim)) 


row.names(anim)<-gsub("phiAbaA1 Rocket","Rocket",row.names(anim)) 
colnames(anim)<-gsub("phiAbaA1 Rocket","Rocket",row.names(anim)) 

row.names(anim)<-gsub("Nayy","Navy",row.names(anim)) 
colnames(anim)<-gsub("Navy","Navy",row.names(anim)) 

row.names(anim)<-gsub("Storma","Storm",row.names(anim)) 
colnames(anim)<-gsub("Storma","Storm",row.names(anim)) 


<<<<<<< HEAD
p$width <- 1400
p$height <- 1400
=======
>>>>>>> 210083d (Changes regarding newlines in systemcalls)



anim[is.na(anim)] <- 0

pdf(paste0(outputdir,"tempheat30.pdf"), height =55, width = 55)
pheatmap::pheatmap(anim,
                   treeheight_row=0,
                   treeheight_col=0,
                   cellwidth = 20, 
                   cellheight = 20,
                   cluster_row = T,
                   cluster_col = T,
                   show_rownames = T,
                   show_colnames = T,
                   cex = 2.5,
                   legend = T,
                   fontsize =2,
                   color = colorRampPalette(
                     RColorBrewer::brewer.pal(n = 7, name = 
"Purples"))(100))
dev.off()


pdf(paste0(outputdir,"tempheat60.pdf"), height =90, width = 90)
e=pheatmap::pheatmap(anim,
                   treeheight_row=0,
                   treeheight_col=0,
                   cellwidth = 30, 
                   cellheight = 30,
                   cluster_row = T,
                   cluster_col = T,
                   show_rownames = T,
                   show_colnames = T,
                   cex = 2.5,
                   legend = T,
                   fontsize =4,
                   color = colorRampPalette(
                     RColorBrewer::brewer.pal(n = 7, name = "PuRd"))(10))


dev.off()

pdf(paste0(outputdir,"tempheat60.pdf"), height =90, width = 90)
pheatmap::pheatmap(anim,
                   treeheight_row=0,
                   treeheight_col=0,
                   cellwidth = 20, 
                   cellheight = 20,
                   cluster_row = T,
                   cluster_col = T,
                   show_rownames = T,
                   show_colnames = T,
                   cex = 2.5,
                   legend = T,
                   fontsize =5,
                   color = colorRampPalette(
                     RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100))
dev.off()

