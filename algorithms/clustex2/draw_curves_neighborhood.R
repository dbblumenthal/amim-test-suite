################
#build time
#function: draw several plots
#(1) one plot contain 20 subplot(can it be that many?)
#(2) one subplot contains 4 curves, including
#     a) clustered genes
#     b) clustered seed genes
#     c) largest module size
#     d) seed genes in the largest module
#     perhaps along with e) second largest module size and seed genes in it
#(3) input file
# 1.neighborhood
# 2.clustered genes
# 3.clustered seed genes
# 4.largest module size
# 5.seed genes in largest
# 6.2nd largest module size  
# 7.seed genes in 2nd largest
################

#setwd("C:/Users/zding/workspace/projects/bio_summer/DiMo/clustex2_64_win")
#setwd("C:/Users/zding/workspace/projects/bio_summer/DiMo/clustex2_32_win/")
#setwd("C:/Users/zding/workspace/projects/bio_summer/DiMo/dimo_lapack_eigen/dimo_lapack_eigen")


args <- commandArgs(trailingOnly=TRUE)

data = as.matrix( read.table(args[1],header=T,sep="\t",quote="",row.names=NULL) )

#data = as.matrix( read.table("TNF_neighborhood.txt",header=T,sep="\t",quote="",row.names=NULL) )

simi = unique(data[,1])
simi = sort(simi, decreasing=T)
levels <- length(simi)

vectors <- c()
if( levels == 6 ){
  vectors <- seq( 0.05, 0.3, by = 0.05)
}
if( levels == 15 ){
  vectors <- seq( 0.02, 0.3, by = 0.02)
}
if( levels == 30 ){
  vectors <- seq( 0.01, 0.3, by =0.01)
}
if( levels == 99 ){
  vectors <- seq(0.01, 0.99, by =0.01)
}



pdf(args[2],width=14,height=7)
#pdf("fix_neighborhood.pdf")


for( i in 1:levels ){
  
  #par(cex=0.85)
  vec = (as.numeric(data[,1])==simi[i])
  temp_mat = data[vec,]
  
  #start to plot
  temp_str_2 = paste( "=",vectors[i],sep="")
  temp_str_1 = paste("Module expansion process: neighborhood",temp_str_2,sep="")
  #win.graph(height=7,width=14)
  split.screen(c(1,2))
  #split.screen(matrix(c(0,.7,  .7,1,  0,0,  1,1), ncol=4))
  if( is.matrix(temp_mat) ){
    

    if(max(as.numeric(temp_mat[,3])) > max(as.numeric(temp_mat[,4])) )
    {
      screen(1)
      plot(temp_mat[,2],temp_mat[,3],col="green","l",ylab="number of genes", xlab="number of clustered genes",main=temp_str_1,lab=c(20,25,5))#green: seed genes
      lines(temp_mat[,2],temp_mat[,4],col="red","l")# black: largest module
      lines(temp_mat[,2],temp_mat[,5],col="red","l",lty=2)# yellow: largest module seeds
      lines(temp_mat[,2],temp_mat[,6],col="blue","l")# blue: 2nd largest module
      lines(temp_mat[,2],temp_mat[,7],col="blue","l",lty=2)# gray: 2nd largest module seeds
      grid()
      legend("topleft", legend=c("clusterd seed genes", "largest module", "largest module seeds","2nd largest module","2nd largest module seeds"), col=c("green","red","red","blue","blue"),lty=c(1,1,2,1,2) )
      screen(2)
      plot( as.numeric(temp_mat[,5])/as.numeric(temp_mat[,4]),temp_mat[,4],ylab="number of genes in the largest module", xlab="seed gene fraction  in largest module","l",main=temp_str_1,lab=c(20,25,10))
      grid()
      #plot(temp_mat[,5]/temp_mat[,4], temp_mat[,2],  ylab="number of clustered genes", xlab="seed gene fraction  in largest module",col="red","l",main=temp_str_1)
      
    }
    if( max(as.numeric(temp_mat[,3])) <= max(as.numeric(temp_mat[,4])) )
    {
      screen(1)
      plot(temp_mat[,2],temp_mat[,4],col="red","l",ylab="number of genes", xlab="number of clustered genes",main=temp_str_1,lab=c(20,25,5))
      lines(temp_mat[,2],temp_mat[,3],col="green","l")
      lines(temp_mat[,2],temp_mat[,5],col="red","l",lty=2)
      lines(temp_mat[,2],temp_mat[,6],col="blue","l")
      lines(temp_mat[,2],temp_mat[,7],col="blue","l",lty=2)
      #legend
      grid()
      legend("topleft", legend=c("clusterd seed genes", "largest module", "largest module seeds","2nd largest module","2nd largest module seeds"), col=c("green","red","red","blue","blue"),lty=c(1,1,2,1,2) )
      screen(2)
      plot( as.numeric(temp_mat[,5])/as.numeric(temp_mat[,4]),temp_mat[,4],ylab="number of genes in the largest module", xlab="seed gene fraction  in largest module","l",main=temp_str_1,lab=c(20,25,10))
      grid()
      #plot(temp_mat[,5]/temp_mat[,4], temp_mat[,2],  ylab="number of clustered genes", xlab="seed gene fraction in largest module","l",main=temp_str_1)
      
    }
    
  }
  if( !is.matrix(temp_mat) ){
    if(max(as.numeric(temp_mat[3])) > max(as.numeric(temp_mat[4])) )
    {
      screen(1)
      plot(temp_mat[2],temp_mat[3],col="green","l",xlab="number of genes", ylab="number of clustered genes",main=temp_str_1,lab=c(20,25,5))#green: seed genes
      lines(temp_mat[2],temp_mat[4],col="red","l")# black: largest module
      lines(temp_mat[2],temp_mat[5],col="red","l",lty=2)# yellow: largest module seeds
      lines(temp_mat[2],temp_mat[6],col="blue","l")# blue: 2nd largest module
      lines(temp_mat[2],temp_mat[7],col="blue","l",lty=2)# gray: 2nd largest module seeds
      grid()
      legend("topleft", legend=c("clusterd seed genes", "largest module", "largest module seeds","2nd largest module","2nd largest module seeds"), col=c("green","red","red","blue","blue"),lty=c(1,1,2,1,2) )
      screen(2)
      plot( as.numeric(temp_mat[5])/as.numeric(temp_mat[4]),temp_mat[4],ylab="largest gene module size", xlab="seed gene fraction  in largest module","l",main=temp_str_1,lab=c(20,25,10))
      grid()
      #plot(temp_mat[,5]/temp_mat[,4], temp_mat[,2],  ylab="number of clustered genes", xlab="seed gene fraction  in largest module",col="red","l",main=temp_str_1)
      
    }
    if( max(as.numeric(temp_mat[3])) <= max(as.numeric(temp_mat[4])) )
    {
      screen(1)
      plot(temp_mat[2],temp_mat[4],col="red","l",xlab="number of genes", ylab="number of clustered genes",main=temp_str_1,lab=c(20,25,5))
      lines(temp_mat[2],temp_mat[3],col="green","l")
      lines(temp_mat[2],temp_mat[5],col="red","l",lty=2)
      lines(temp_mat[2],temp_mat[6],col="blue","l")
      lines(temp_mat[2],temp_mat[7],col="blue","l",lty=2)
      #legend
      grid()
      legend("bottomright", legend=c("clusterd seed genes", "largest module", "largest module seeds","2nd largest module","2nd largest module seeds"), col=c("green","red","red","blue","blue"),lty=c(1,1,2,1,2) )
      screen(2)
      plot( as.numeric(temp_mat[5])/as.numeric(temp_mat[4]),temp_mat[4],ylab="largest gene module size", xlab="seed gene fraction  in largest module","l",main=temp_str_1,lab=c(20,25,10))
      grid()
      #plot(temp_mat[,5]/temp_mat[,4], temp_mat[,2],  ylab="number of clustered genes", xlab="seed gene fraction in largest module","l",main=temp_str_1)
      
    }
  }
  close.screen(all.screens=TRUE)
} 



dev.off()


# pdf("temp.pdf")
# for( i in 1:levels ){
#   
#   
#   vec = (data[,1]==simi[i])
#   temp_mat = data[vec,]
#   
#   #start to plot
#   temp_str_2 = paste( vectors[i],"%",sep="")
#   temp_str_1 = paste(temp_str_2,"neighborhood",sep=" ")
#   
#   plot(as.numeric(temp_mat[,7])/as.numeric(temp_mat[,6]), temp_mat[,2],  ylab="number of clustered genes", xlab="seed gene fraction in largest module","l",main=temp_str_1)
#   
#   
# }
# dev.off()



