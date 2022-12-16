#system("rm genomes/*;cp ../Metaphage_outputs_whole/cd-hit/vOTUs_consensus/* genomes/")
#TODO: do blastp mode :)

library(foreach)
library(doParallel)
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

#system(paste0("rm ",genomesdir,"*; ls ../Metaphage_outputs_whole/cd-hit/vOTUs_consensus/* | while read f; do cp $f ",genomesdir,"/; done"))


# file_filtd=read.table("../iglike_git/iglike2022/results/count_meta_filtd.tsv",sep="\t",header=T)
# fils<-list.files("genomes/")
# fils
# file_filtd$ViralOTU<-paste0(file_filtd$ViralOTU,".fasta")
# file_filtd$ViralOTU
# for(i in fils){
#   if(i %in% file_filtd$ViralOTU){
#     print(i)
#   }else{
#     system(paste0("rm genomes/",i))
#   }
# }

fils<-list.files(genomesdir,pattern = "fasta")
system('mv collected_ANI collected_ANI.legacy')

anidf <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("seq1", "seq2", "ANI")
colnames(anidf) <- x
count=0

### csak egyet csinaljon meg 
anidf = foreach(i=1:length(fils), .combine=rbind) %dopar% {
  totaldf <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("seq1", "seq2", "ANI")
  colnames(anidf) <- x
  print(fils[i])
  if(file.info(paste0(genomesdir,fils[i]))$size>1000000) next
  
  count=count+1
  print(count)
  #if(count==2) break()
  for(o in i:length(fils)){
    print(fils[o])
    if(file.info(paste0(genomesdir,fils[o]))$size>1000000) next
    
    #system(paste0("~/fastANI -q genomes/",fils[i]," -r ../genomes/",fils[o]," -o tmp"))
    system(paste0("rm ",outputdir,"tmptbl",i),ignore.stderr = TRUE)
    system(paste0('blastn  -subject ',genomesdir,fils[o],' -query ',genomesdir
                  ,fils[i],'    -out ',outputdir,'tmptbl',i,' -max_target_seqs 1    -outfmt "6 qseqid sseqid qstart qend sstart send qseq sseq evalue bitscore pident qlen slen" ')
           ,ignore.stderr = TRUE)
    #system(paste0("sed -i '1s/^/",fils[i],"\t",fils[o],"\t","/' tmp"))
    
    if(file.exists(paste0(outputdir,"tmptbl",i))){
      if(file.info(paste0(outputdir,"tmptbl",i))$size==0) next
      
      tmptbl<-read.table(paste0(outputdir,"tmptbl",i),sep="\t",header=F)
      names(tmptbl)<-c("qseqid" ,"sseqid" ,"qstart", "qend", "sstart", "send","qseq","sseq", "evalue", "bitscore" ,"pident" ,"qlen", "slen")
      
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
              if(tmptbl$weightedpid[ii]>tmptbl$weightedpid[xx]) tmptbl$weightedpid[xx]<-tmptbl$weightedpid[xx]*
                  (
                    (length(tmptbl$qstart[xx]:tmptbl$qend[xx])-length(intersect(
                      tmptbl$qstart[ii]:tmptbl$qend[ii],tmptbl$qstart[xx]:tmptbl$qend[xx]))
                    )/
                      length(tmptbl$qstart[xx]:tmptbl$qend[xx]))
              if(tmptbl$weightedpid[ii]<tmptbl$weightedpid[xx]) tmptbl$weightedpid[ii]<-tmptbl$weightedpid[ii]*
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
              if(tmptbl$weightedpid[ii]>tmptbl$weightedpid[xx]) tmptbl$weightedpid[xx]<-tmptbl$weightedpid[xx]*(
                (length(tmptbl$sstart[xx]:tmptbl$send[xx])-length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$sstart[xx]:tmptbl$send[xx])))/
                  length(tmptbl$sstart[xx]:tmptbl$send[xx])
              )
              if(tmptbl$weightedpid[ii]<tmptbl$weightedpid[xx]) tmptbl$weightedpid[ii]<-tmptbl$weightedpid[ii]*(
                (length(tmptbl$sstart[ii]:tmptbl$send[ii])-length(intersect(tmptbl$sstart[ii]:tmptbl$send[ii],tmptbl$sstart[xx]:tmptbl$send[xx])))/
                  length(tmptbl$sstart[ii]:tmptbl$send[ii])
              )
            }
            
          }
        }
        
      }
      
      
      
      #print(paste0("slen is ",tmptbl$slen[1],"qlen is",tmptbl$qlen[1],"weightedpid is ",sum(tmptbl$weightedpid),"biggerlen is ",tmptbl$biggerlen[1]))  
      #print(sum(tmptbl$weightedpid))
      #print(coverage)
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
slackr::slackr_bot("Finished with ANIs")
