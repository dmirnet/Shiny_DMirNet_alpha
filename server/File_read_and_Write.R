## Read dataset from csv file
Read_Scale<-function(dataset,root){
  data=NULL
  if(is.null(dataset)){
    file=paste0(root,"//sample_datasets/Sample_31miRNAs_1151mRNAs_expr.csv")
    data<-read.csv(file, header=TRUE, sep=",")
  }
  else{
    data=read.csv(dataset, header=TRUE, sep=",")
  }
  data=scale(data)
  data=round(data,2)
  data=as.matrix(data)
  return(data)
}

write_file<-function (data,filename){
  write.table(data, filename, sep=",", row.names=FALSE, col.names = TRUE)
}

read_file<-function(file,scale_data=TRUE){
  data<-read.csv(file, header=TRUE, sep=",",encoding = "UTF-8")
  if(scale_data==TRUE){
    data<-scale(data)
  }
  data=as.matrix(data)
  return(data)
}

write_description<-function(directory,fname,content){
  temp=getwd()
  setwd(directory)
  fileConn<-file(fname)
  writeLines(content, fileConn,useBytes=T)
  close(fileConn)
  setwd(temp)
}
write_description_append<-function(directory,fname,content){
  temp=getwd()
  setwd(directory)
  write(content,file="readme.txt",append=TRUE)
  setwd(temp)
}
