#Ensembling the selected files
ensemble_files<-function(data_list,type,irm_topk,filename,colname){
  Result=NULL
  if(type=="ensemble_irm_mean"){
    irm_result <- IRM(do.call(cbind,data_list),"mean")
    Result <- uptri2mat(colname, irm_result)
  }
  else if(type=="ensemble_irm_median"){
    irm_result <- IRM(do.call(cbind,data_list),"median")
    Result <- uptri2mat(colname, irm_result)
  }
  else if(type=="ensemble_topk_mean"){
    irm_result <- IRM_topk(do.call(cbind,data_list),irm_topk,"mean")
    Result <- uptri2mat(colname, irm_result)
  }else{
    irm_result <- IRM_topk(do.call(cbind,data_list),irm_topk,"median")
    Result <- uptri2mat(colname, irm_result)
  }
  write_file(Result,paste0(dir_ensemble,filename))
}

result_analysis<-function(data_analysis, type, m,mirna,fname){

  if(type=="analysis_mpairs"){
    result=Extopk_mimR(data_analysis,m,mirna)
    write_file(result,paste0(dir_analysis,"(top_m)Analysis_of_",fname))
    #write description
    content=c(paste0("(top_m)Analysis_of_",fname,":- Analysis type=Select top m pairs",", top m=",m,",number of miRNAs=",mirna))
    write_description_append(dir_analysis,"readme.txt",content)
  }else{
    result=Extopk_mimR_miRk(data_analysis,m,mirna)
    write_file(result,paste0(dir_analysis,"(top_m_mRNA)Analysis_of_",fname))
    #write description
    content=c(paste0("(top_m)Analysis_of_",fname,":- Analysis type=Select top m mRNAs for each miRNA",", top m= ",m,",number of miRNAs=",mirna))
    write_description_append(dir_analysis,"readme.txt",content)
  }
}

result_validation<-function(selected_data, valdata,fname){
  result=Validation(selected_data,valdata)
  write.table(result, paste0(dir_validation,"Validation_of_",fname),sep = ",",row.names = FALSE,col.names = TRUE)
}