# Setup the working directory for the output files 

#set working directory for first run of the experment 
setWd_first_run<-function(dir,run_iter){
  setwd(dir)
  dir.create("DMirNet_Data")
  setwd(paste0(getwd(),"/DMirNet_Data/"))
  dir.create(paste0("Experment_Results_",run_iter))
  setwd(paste0(getwd(),"/Experment_Results_1/"))
  create_directory()
}

#set working directory for first run of the experment 
setWd_for_reRun<-function(dataiter){
  dirnew=paste0("Experment_Results_",dataiter)
  dir.create(dirnew)
  dirnew=paste0("/",dirnew)
  setwd(paste0(getwd(),dirnew))
  create_directory()
}

#create directory
create_directory<-function(){
  dir_experment<<-getwd()
  dir.create("Direct_Correlation_and_Bootstrapping")
  dir_direct_bootstrap<<-paste0(getwd(),"/Direct_Correlation_and_Bootstrapping/")
  write_description(dir_direct_bootstrap,"readme.txt","------Experment Description----")
  dir.create("Ensemble_Results")
  dir_ensemble<<-paste0(getwd(),"/Ensemble_Results/")
  write_description(dir_ensemble,"readme.txt","------Experment Description----")
  dir.create("Result_Analysis")
  dir_analysis<<-paste0(getwd(),"/Result_Analysis/")
  write_description(dir_analysis,"readme.txt","------Experment Description----")
  dir.create("Validation_Results")
  dir_validation<<-paste0(getwd(),"/Validation_Results/")
  write_description(dir_validation,"readme.txt","------Experment Description----")
  setwd(dir_direct_bootstrap)
  dir.create("uppertri")
  dir_direct_bootstrap_uppertri<<-paste0(getwd(),"/uppertri/")
  setwd(dir_experment)
}

setup_dirs_only<-function(){
  dir_experment<<-getwd()
  dir_direct_bootstrap<<-paste0(getwd(),"/Direct_Correlation_and_Bootstrapping/")
  dir_ensemble<<-paste0(getwd(),"/Ensemble_Results/")
  dir_analysis<<-paste0(getwd(),"/Result_Analysis/")
  dir_validation<<-paste0(getwd(),"/Validation_Results/")
  setwd(dir_direct_bootstrap)
  dir_direct_bootstrap_uppertri<<-paste0(getwd(),"/uppertri/")
  setwd(dir_experment)
}
