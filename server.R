
## Define server logic required
shinyServer(function(input, output,session) {
  #Observe Run Button
  observeEvent(input$run,{
    tryCatch({
      #run progress
      withProgress(message = 'Computing DMirNet...', style="notification", value = 0.01, {
        #Get the values of input parameters
        if(input$parallel_method=="core"){
          clusters<<-as.numeric(input$cores)
        }else{
          i=1
          hosts=list()
          while(i<=(count-1)){
            val=paste("addmore",i,sep = "")
            hosts=c(hosts,as.character(input[[val]]))
            i=i+1
          }
          clusters<<-hosts
        }
        
        iterations<<-as.numeric(input$iteration)
        sample.percentage<<-as.numeric(input$rate)/100
        miRNAs<<-as.numeric(input$miRNAs)
        dataset<<-input$dataset$datapath
        seed_num<<-as.numeric(input$seed)
        dir=as.character(readDirectoryInput(session, 'directory'))
        dir=normalizePath(dir)
        working_dir<<-paste0(dir,"/DMirNet_Data/")
        
        #run the experment 
        disable("page")
        #check the run type and set the working directory
        if(rerun==FALSE){
          setWd_first_run(dir,rerun_iter)
          rerun_iter<<-rerun_iter+1
        }
        if(input$wd_list=="New Experiment"){
          setWd_for_reRun(rerun_iter)
          rerun_iter<<-rerun_iter+1
          ensemble_iter<<-1
        }
        incProgress(0.3)
        progress_iteration=0.1
        if(input$corpcor&&input$space&&input$corND&&input$ida){
          progress_iteration=0.13
        }else if((input$corpcor&&input$space&&input$corND) || (input$space&&input$corND&&input$ida) || (input$corpcor&&input$corND&&input$ida) ){
          progress_iteration=0.17
        }else if((input$corpcor&&input$space)|| (input$corND&&input$ida)|| (input$corpcor&& input$corND) || (input$space&&input$ida)||(input$corpcor&&input$ida)||(input$space&&input$corND)){
          progress_iteration=0.25
        }else {
          progress_iteration=0.5
        }
        #Read data and scale 
        data<<-Read_Scale(dataset,root_dir)
        ### Step 1. Perform experiment Direct corelation and bootstrapping ###
        #corpcor 
        if(input$corpcor){
          if(input$bootstrap_action=="Disable"){
            setProgress(value = NULL, message =NULL , detail = "Performing Corpcor Direct corelation",session = session)
            buildCorpcor(data,input$lambda)
            incProgress(progress_iteration)
          }else{
            setProgress(value = NULL, message =NULL , detail = "Bootstrapping with Corpcor",session = session)
            paramsc=c(as.numeric(input$lambda))
            corpcor_bootstrap(data,input$bootstrap_method,iterations,sample.percentage,as.numeric(input$bootstrap_topk),paramsc)
            incProgress(progress_iteration)
          }
          #write description
          content=c(paste0("Corpcor parameter setting:- lambda=",as.numeric(input$lambda)))
          write_description_append(dir_direct_bootstrap,"readme.txt",content)
        }
        #Space
        if(input$space){
          if(input$bootstrap_action=="Disable"){
            setProgress(value = NULL, message =NULL , detail = "Performing Space Direct corelation",session = session)
            buildSpace(data,input$lam1,input$lam2,input$w,input$iter)
            incProgress(progress_iteration)
          }else{
            setProgress(value = NULL, message =NULL , detail = "Bootstrapping with Space",session = session)
            paramss=c(input$lam1,input$lam2,input$w,input$iter)
            space_bootstrap(data,input$bootstrap_method,iterations,sample.percentage,as.numeric(input$bootstrap_topk),paramss)
            incProgress(progress_iteration)
          }
          #write description
          content=c(paste0("Space parameter setting:- lam1=",input$lam1,",lam2=",input$lam2,", weight=",input$w,",iteration=",input$iter))
          write_description_append(dir_direct_bootstrap,"readme.txt",content)
        }
        #CorND
        if(input$corND){
          if(input$bootstrap_action=="Disable"){
            setProgress(value = NULL, message =NULL , detail = "Performing CorND Direct corelation",session = session)
            buildCorND(data,input$corND_alpha,input$corND_beta)
            incProgress(progress_iteration)
          }else{
            setProgress(value = NULL, message =NULL , detail = "Bootstrapping with corND",session = session)
            paramsco=c(input$corND_alpha,input$corND_beta)
            corND_bootstrap(data,input$bootstrap_method,iterations,sample.percentage,as.numeric(input$bootstrap_topk),paramsco)
            incProgress(progress_iteration)
          }
          #write description
          content=c(paste0("CorND parameter setting:- alpha=",input$corND_alpha,", beta=",input$corND_beta))
          write_description_append(dir_direct_bootstrap,"readme.txt",content)
        }
        #IDA
        if(input$ida){
          if(input$bootstrap_action=="Disable"){
            setProgress(value = NULL, message =NULL , detail = "Performing IDA Direct corelation",session = session)
            buildIDA(data,input$alpha_idea)
          }else{ 
            setProgress(value = NULL, message =NULL , detail = "Bootstrapping with IDA",session = session)
            paramsi=c(input$alpha_idea)
            IDA_bootstrap(data,input$bootstrap_method,iterations,sample.percentage,as.numeric(input$bootstrap_topk),paramsi)
            incProgress(progress_iteration)
          }
          #write description
          content=c(paste0("IDA parameter setting:- alpha=",input$alpha_idea))
          write_description_append(dir_direct_bootstrap,"readme.txt",content)
          incProgress(progress_iteration)
        }
        #write description
        content=NULL
        if(input$bootstrap_action!="Disable"){
          if(input$bootstrap_method=="bootstrap_irm_median"||input$bootstrap_method=="bootstrap_irm_mean"){
            content=c(paste0("Bootstrapping parameter setting:- Bootstrapping method=",input$bootstrap_method,", Iteration=",iterations,", Sample Rate(%)=",input$rate))
          }else{
            content=c(paste0("Bootstrapping parameter setting:- Bootstrapping method=",input$bootstrap_method,", Iteration=",iterations,", Sample Rate(%)=",input$rate,", topk=",as.numeric(input$bootstrap_topk)))
          }
        }
        content=c(content,"-------------------------------------------------------")
        write_description_append(dir_direct_bootstrap,"readme.txt",content)
        #Result show
        setProgress(value = NULL, message =NULL , detail = "Preparing Result",session = session)
        result_show("tab1")
        incProgress(0.1)
        Sys.sleep(2)
        updateActionButton(session,'run',"Re-Run")  
        enable("page")
        setwd(working_dir)
        #update list of experment 
        wd_dirs=list.files(recursive = FALSE)
        rerun_iter<<-length(wd_dirs)+1
        shinyjs::show("wd_list")
        i=length(wd_dirs)
        j=i+1
        wd_dirs[j]="New Experiment"
        if(input$wd_list=="New Experiment" || rerun==FALSE){
          updateRadioButtons(session,'wd_list',choices = wd_dirs,selected = wd_dirs[i])
        }else{
          updateRadioButtons(session,'wd_list',choices = wd_dirs,selected = input$wd_list)
        }
        rerun<<-TRUE
        incProgress(0.1)
        Sys.sleep(1)
      })#end of progress 
     },error=function(cond){
       showModal(modalDialog(
         title = "ERROR!",
         paste0("Unable to perform experment.Please try again. Error detail=",cond),
         easyClose=TRUE,
         fade=TRUE))
     },finally = {
       enable("page")
     })#end of tryCatch
  })#end of run observer
  
  #Observer for tab selection
  observeEvent(input$tabpanel,{
    if(!is.null(working_dir)){  
      if(input$tabpanel=="Ensemble"){
        if(input$wd_list=="New Experiment" || input$wd_list=="None"){  
          updateRadioButtons(session,'ensemble_result_list',NULL,choices =list("None"),selected = FALSE)
          updateTextAreaInput(session,"ensemble_readme","Expernment Description",value = "-")
          updateTextInput(session,"ensemble_prev_file","File",value = "-")
          updateCheckboxGroupInput(session,'ag_list',NULL,choices =list("NULL")) 
          output$table_ensemble_result <- DT::renderDataTable(DT::datatable({NULL}))
        }else{
          temp=getwd()
          setwd(dir_direct_bootstrap_uppertri)
          result_files=list.files(pattern ="\\.csv$",recursive = FALSE)
          i=1
          file_name=list()
          while(i<=length(result_files)){
            file_name[[i]] <-result_files[i]
            i <- i + 1
          }
          if(i==1){
            file_name[[1]]<-"None"
          }
          updateCheckboxGroupInput(session,'ag_list',label=NULL,choices = file_name,selected = c(file_name[1],file_name[2]))
          setwd(temp)
          result_show("tab2")
        }
      }else if(input$tabpanel=="Result Analysis"){
        if(input$wd_list=="New Experiment" || input$wd_list=="None"){
          updateRadioButtons(session,'analysis_result_list',NULL,choices =list("None"),selected = FALSE)
          updateTextAreaInput(session,"analysis_readme","Expernment Description",value = "-")
          updateTextInput(session,"analysis_prev_file","File",value = "-")
          updateRadioButtons(session,'analysis_file',NULL,choices =list("NULL"))
          output$table_analysis_result <- DT::renderDataTable(DT::datatable({NULL}))
        }else{
          temp=getwd()
          setwd(dir_direct_bootstrap)
          result_files=list.files(pattern ="\\.csv$",recursive = FALSE)
          i=1
          file_name=list()
          while(i<=length(result_files)){
            file_name[[i]] <-result_files[i]
            i <- i + 1
          }
          setwd(temp)
          setwd(dir_ensemble)
          result_files2=list.files(pattern ="\\.csv$",recursive = FALSE)
          j=length(result_files2)+length(result_files)
          z=1
          while(i<=j){
            file_name[[i]] <-result_files2[z]
            z<-z+1
            i <- i + 1
          }
          if(i==1){
            file_name[[1]]<-"None"
          }
          updateRadioButtons(session,'analysis_file',label=NULL,choices = file_name,selected = file_name[1])
          setwd(temp)
          result_show("tab3")
        }
      }else if(input$tabpanel=="Validation"){
        if(input$wd_list=="New Experiment" || input$wd_list=="None"){
          updateRadioButtons(session,'validation_result_list',NULL,choices =list("None"),selected = FALSE)
          updateTextAreaInput(session,"validation_readme","Expernment Description",value = "-")
          updateTextInput(session,"validation_prev_file","File",value = "-")
          updateRadioButtons(session,'validation_file',NULL,choices =list("NULL"))
          output$table_validation_result <- DT::renderDataTable(DT::datatable({NULL}))
        }else{
          temp=getwd()
          setwd(dir_analysis)
          result_files=list.files(pattern ="\\.csv$",recursive = FALSE)
          i=1
          file_name=list()
          while(i<=length(result_files)){
            file_name[[i]] <-result_files[i]
            i <- i + 1
          }
          if(i==1){
            file_name[[1]]<-"None"
          }
          updateRadioButtons(session,'validation_file',label=NULL,choices = file_name,selected = file_name[1])
          setwd(temp)
          result_show("tab4")
        }
      }
    }else{
      update_display()
    }
  })
  #Observer for running Ensemble 
  observeEvent(input$run_ensemble,{
    tryCatch({
      withProgress(message = 'Computing DMirNet...', style="notification", value = 0.01, {
        setProgress(value = NULL, message =NULL , detail = "IRM Ensemble for selected files",session = session)
        disable("run_ensemble")
        temp=getwd()
        setwd(dir_direct_bootstrap_uppertri)
        data_list=list()
        filename=paste0("Ensemble_file_result_",ensemble_iter,".csv")
        readme_add=paste0(filename,":- ensembled files=")
        i=1
        while(i<=length(input$ag_list)){
          data_list[[i]]=read_file(input$ag_list[i])
          readme_add=paste0(readme_add,input$ag_list[i],",")
          i=i+1
        }
        setwd(temp)
        setwd(dir_direct_bootstrap)
        colname=colnames(read_file(input$ag_list[1]))
        setwd(temp)
        incProgress(0.5)
        ensemble_files(data_list,input$ensemble,as.numeric(input$irm_topk),filename,colname)
        ensemble_iter<<-ensemble_iter+1
        setwd(dir_ensemble)
        if(input$ensemble=="ensemble_irm_mean"||input$ensemble=="ensemble_irm_median"){
          readme_add=paste0(readme_add,"; Ensemble Method=",input$ensemble)
        }else{
          readme_add=paste0(readme_add,"; Ensemble Method=",input$ensemble,",topk=",input$irm_topk)
        }
        write(readme_add,file="readme.txt",append=TRUE)
        setwd(temp)
        setProgress(value = NULL, message =NULL , detail = "Preparing Result",session = session)
        result_show("tab2")
        incProgress(1)
        Sys.sleep(2)
      })
    },error=function(cond){
      showModal(modalDialog(
        title = "ERROR!",
        paste0("Unable to perform ensemble.Please try again. Error detail=",cond),
        easyClose=TRUE,
        fade=TRUE))
    },finally={
      enable("run_ensemble")
    })
  }) 
  
  #Observer for analysis run
  observeEvent(input$run_analysis,{
    withProgress(message = 'Computing DMirNet...', style="notification", value = 0.01, {
      setProgress(value = NULL, message =NULL , detail = "Performing analysis for selected file",session = session)
      tryCatch({
        temp=getwd()
        disable("run_analysis")
        data_analysis=NULL
        setwd(dir_direct_bootstrap)
        if(file.exists(input$analysis_file)){
          data_analysis=read_file(input$analysis_file)
        }else{
          setwd(dir_ensemble)
          data_analysis=read_file(input$analysis_file)
        }
        incProgress(0.5)
        result_analysis(data_analysis,input$analysis_criteria,as.numeric(input$analysis_topk),as.numeric(input$miRNAs),input$analysis_file)
        setwd(temp)
        setProgress(value = NULL, message =NULL , detail = "Preparing Result",session = session)
        result_show("tab3")
        enable("run_analysis")
        incProgress(1)
        Sys.sleep(2)
        
      },error=function(cond){
        showModal(modalDialog(
          title = "ERROR!",
          paste0("Unable to perform analysis. Please try again. Error detail=",cond),
          easyClose=TRUE,
          fade=TRUE))
      },finally={
        enable("run_analysis")
      })
    })
  }) 
  #Observer for validation run
  observeEvent(input$run_validation,{
    withProgress(message = 'Computing DMirNet...', style="notification", value = 0.01, {
      setProgress(value = NULL, message =NULL , detail = "Performing validation for selected file",session = session)
      tryCatch({
        temp=getwd()
        disable("run_validation")
        setwd(dir_analysis)
        val_dataset=input$validationdata$datapath
        if (is.null(val_dataset)){
          val_dataset=paste0(root_dir,"/sample_datasets/Sample_ValidationDataset.csv")
          val_dataset_name="Sample_ValidationDataset.csv"
        }else{
          val_dataset_name=basename(val_dataset)
        }
        dataset_validation=read.csv( val_dataset, header=TRUE, sep=",")
        data_validation=read_file(input$validation_file,FALSE)
        incProgress(0.5)
        result_validation(data_validation,dataset_validation,input$validation_file)
        #write description
        content=c(paste0("Validation_of_",input$validation_file,":-Data used for Result Validation=",val_dataset_name))
        write_description_append(dir_validation,"readme.txt",content)
        setwd(temp)
        setProgress(value = NULL, message =NULL , detail = "Preparing Result",session = session)
        result_show("tab4")
        enable("run_validation")
        incProgress(1)
        Sys.sleep(2)
      },error=function(cond){
        showModal(modalDialog(
          title = "ERROR!",
          paste0("No matching record found",cond),
          easyClose=TRUE,
          fade=TRUE))
      },finally = {
        enable("run_validation")
      })
    })
  }) 
  
  #output result show function
  result_show<-function(tab){
    #shows result for expermental tab result
    if(tab=="tab1"){
      temp=getwd()
      setwd(dir_direct_bootstrap)
      result_files=list.files(pattern ="\\.csv$",recursive = FALSE)
      i=1
      file_name=list()
      while(i<=length(result_files)){
        file_name[[i]] <-result_files[i]
        i <- i + 1
      }
      if(i==1){
        file_name[[1]]<-"None"
      }
      updateRadioButtons(session,'dir_boot_list',label=NULL,choices = file_name,selected = file_name[1])
      readme_view1<<-readLines("readme.txt", file.info("readme.txt")$size)
      updateTextAreaInput(session,"exp_readme",value =readme_view1 )
      updateTextInput(session,"ex_prev_file",value=file_name[1])
      disable("ex_prev_file")
      setwd(temp)
    }else if(tab=="tab2"){
      #shows result for ensembling result
      temp=getwd()
      setwd(dir_ensemble)
      result_files=list.files(pattern ="\\.csv$",recursive = TRUE)
      i=1
      file_name=list()
      while(i<=length(result_files)){
        file_name[[i]] <-result_files[i]
        i <- i + 1
      }
      if(i==1){
        file_name[[1]]<-"None"
      }
      updateRadioButtons(session,'ensemble_result_list',label=NULL,choices = file_name,selected = file_name[1])
      readme_view2<<-readLines("readme.txt", file.info("readme.txt")$size)
      updateTextAreaInput(session,"ensemble_readme",value =readme_view2 )
      updateTextInput(session,"ensemble_prev_file",value=file_name[1])
      disable("ensemble_prev_file")
      setwd(temp)
    }else if(tab=="tab3"){
      #shows result for result analysis
      temp=getwd()
      setwd(dir_analysis)
      result_files=list.files(pattern ="\\.csv$",recursive = TRUE)
      i=1
      file_name=list()
      while(i<=length(result_files)){
        file_name[[i]] <-result_files[i]
        i <- i + 1
      }
      if(i==1){
        file_name[[1]]<-"None"
      }
      updateRadioButtons(session,'analysis_result_list',label=NULL,choices = file_name,selected = file_name[1])
      readme_view3<<-readLines("readme.txt", file.info("readme.txt")$size)
      updateTextAreaInput(session,"analysis_readme",value =readme_view3 )
      updateTextInput(session,"analysis_prev_file",value=file_name[1])
      disable("analysis_prev_file")
      setwd(temp)
    }else if(tab=="tab4"){
      #shows result for result validation
      temp=getwd()
      setwd(dir_validation)
      result_files=list.files(pattern ="\\.csv$",recursive = TRUE)
      i=1
      file_name=list()
      while(i<=length(result_files)){
        file_name[[i]] <-result_files[i]
        i <- i + 1
      }
      if(i==1){
        file_name[[1]]<-"None"
      }
      updateRadioButtons(session,'validation_result_list',label=NULL,choices = file_name,selected = file_name[1])
      readme_view4<<-readLines("readme.txt", file.info("readme.txt")$size)
      updateTextAreaInput(session,"validation_readme",value =readme_view4 )
      updateTextInput(session,"validation_prev_file",value=file_name[1])
      disable("validation_prev_file")
      setwd(temp)
    }
  }
  #file preview selection observer
  observeEvent(input$dir_boot_list,{
    if(input$dir_boot_list!="None"){
      updateTextInput(session,"ex_prev_file",value=input$dir_boot_list)
      disable("ex_prev_file")
      output$table_experment_result <- DT::renderDataTable(DT::datatable({NULL}))
    }
  })
  observeEvent(input$ensemble_result_list,{
    if(input$ensemble_result_list!="None"){
      updateTextInput(session,"ensemble_prev_file",value=input$ensemble_result_list)
      disable("ensemble_prev_file")
      output$table_ensemble_result <- DT::renderDataTable(DT::datatable({NULL}))
    }
  })
  observeEvent(input$analysis_result_list,{
    if(input$analysis_result_list!="None"){
      updateTextInput(session,"analysis_prev_file",value=input$analysis_result_list)
      disable("analysis_prev_file")
      output$table_analysis_result <- DT::renderDataTable(DT::datatable({NULL}))
    }
  })
  observeEvent(input$validation_result_list,{
    if(input$validation_result_list!="None"){
      updateTextInput(session,"validation_prev_file",value=input$validation_result_list)
      disable("validation_prev_file")
      output$table_validation_result <- DT::renderDataTable(DT::datatable({NULL}))
    }
  })
  #Function for Display result in table  
  display_table<-function(directory,file_name,tab=""){
    tryCatch({
      if(tab=="tab1"&&file_name!="None"){
        data=read.csv(paste0(directory,file_name), header=TRUE, sep=",")
        output$table_experment_result <- DT::renderDataTable(DT::datatable({
          data
        }))
      }else if(tab=="tab2"&&file_name!="None"){
        data=read.csv(paste0(directory,file_name), header=TRUE, sep=",")
        output$table_ensemble_result <- DT::renderDataTable(DT::datatable({
          data
        }))
      }else if(tab=="tab3"&&file_name!="None"){
        data=read.csv(paste0(directory,file_name), header=TRUE, sep=",")
        output$table_analysis_result <- DT::renderDataTable(DT::datatable({
          data
        }))
      }else if(tab=="tab4"&&file_name!="None"){
        data=read.csv(paste0(directory,file_name), header=TRUE, sep=",")
        output$table_validation_result <- DT::renderDataTable(DT::datatable({
          data
        }))
      }
    },error=function(cond){
      showModal(modalDialog(
        title = "ERROR!",
        paste0("Cannot display output file. Error type=",cond),
        easyClose=TRUE,
        fade=TRUE))
    })
  }
  
  #---------------------------------------- UI observers----------------------
  observe({
    path=as.character(readDirectoryInput(session, 'directory'))
    if (file.exists(paste0(path,"/DMirNet_Data"))){
      showModal(modalDialog(
        title = "WARNING!",
        "The working directory has beed used before. You can see the list of output from previouse run on the preview area"
      ))
      working_dir<<-paste0(path,"/DMirNet_Data/")
      setwd(working_dir)
      wd_dirs=list.files(recursive = FALSE)
      rerun_iter<<-length(wd_dirs)+1
      updateActionButton(session,'run',"Re-Run") 
      rerun<<-TRUE
      shinyjs::show("wd_list")
      j=length(wd_dirs)+1
      wd_dirs[j]="New Experiment"
      updateRadioButtons(session,'wd_list',choices = wd_dirs,selected = wd_dirs[1])
      wd_list_show(wd_dirs[1])
    }
  })
  #Observer for the running methods
  observeEvent(input$parallel_method,{
    if(input$parallel_method=="core"){
      enable('run')
      shinyjs::show("no_cores")
      shinyjs::hide("no_machines")
      shinyjs::hide('host_input')
    }else{
      if(run_disable==TRUE){
        disable('run')
      }else{
        enable('run')
      }
      shinyjs::show('host_input')
      shinyjs::hide("no_cores")
      shinyjs::show("no_machines")
    } 
  })
  #Observer for displaying Bootstrapping settings
  observeEvent(input$bootstrap_action,{
    if(input$bootstrap_action=="Enable"){
      shinyjs::show("boot_div")
    }else{
      shinyjs::hide("boot_div")
    }
  })
  #Observer for Selection of Editing parameters for Direct Association Methods 
  observeEvent(input$dir_corl_edit,{
    if(input$dir_corl_edit){
      shinyjs::show("edit_div")
    }else{
      shinyjs::hide("edit_div")
    }
  })
  #Observer for Direct Correlation Methods selections
  observeEvent(c(input$corpcor,input$space,input$corND,input$ida,input$bootstrap_method,input$ensemble,input$ag_list,input$analysis_file,input$validation_file,input$validationdata),{
    #Upon selection show input parameters of Direct Association Methods 
    if(input$corpcor){
      shinyjs::show("divcorp")
    }else{
      shinyjs::hide("divcorp")
    }
    if(input$space){
      shinyjs::show("divspace")
      dat=NULL
      if(!is.null(input$dataset$datapath)){
        dat<-read.csv(input$dataset$datapath, header=TRUE, sep=",")
      }else{
        dat<-read.csv(paste0(root_dir,"/sample_datasets/Sample_31miRNAs_1151mRNAs_expr.csv"), header=TRUE, sep=",")  
      }
      l1=(1/sqrt(nrow(dat))*qnorm(1-1/(2*ncol(dat)^2)))*nrow(dat)*.161
      updateNumericInput(session,'lam1',label="lam1",value=l1,min=1,max=10)
      col=ncol(dat)
      row=nrow(dat)
      updateNumericInput(session,'weight',label="weight for each data point",value=0,min=0,max=row)
    }else{
      shinyjs::hide("divspace")
    }
    if(input$corND){
      shinyjs::show("divcorND")
    }else{
      shinyjs::hide("divcorND")
    }
    if(input$ida){
      shinyjs::show("divida")
    }else{
      shinyjs::hide("divida")
    }
    if(!input$space&&!input$corpcor&&!input$ida&&!input$corND){
      disable("run")
    }else{
      enable("run")
    }
    if(input$bootstrap_method=="bootstrap_topk_mean"||input$bootstrap_method=="bootstrap_topk_median"){
      shinyjs::show("bootstrap_topk")
    }else{
      shinyjs::hide("bootstrap_topk")
    }
    if(input$ensemble=="ensemble_topk_mean"||input$ensemble=="ensemble_topk_median"){
      shinyjs::show("irm_topk")
    }else{
      shinyjs::hide("irm_topk")
    }
    if(length(input$ag_list)<2){
      disable("run_ensemble")
    }else{
      enable("run_ensemble")
    }
    if(input$analysis_file=="NULL"){
      disable("run_analysis")
    }else{
      enable("run_analysis")
    }
  })
  
  #Observer for Directory input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        # condition prevents handler execution on initial app launch
        # launch the directory selection dialog with initial path read from the widget
        path = choose.dir(default = readDirectoryInput(session, 'directory'))
        # update the widget value
        updateDirectoryInput(session, 'directory', value = path)
        #check for directory existance
        if (file.exists(paste0(path,"/DMirNet_Data"))){
          showModal(modalDialog(
            title = "WARNING!",
            "The working directory has beed used before. You can see the list of output from previouse run on the preview area"
          ))
          working_dir<<-paste0(path,"/DMirNet_Data/")
          setwd(working_dir)
          wd_dirs=list.files(recursive = FALSE)
          rerun_iter<<-length(wd_dirs)+1
          updateActionButton(session,'run',"Re-Run") 
          rerun<<-TRUE
          shinyjs::show("wd_list")
          j=length(wd_dirs)+1
          wd_dirs[j]="New Experiment"
          updateRadioButtons(session,'wd_list',choices = wd_dirs,selected = wd_dirs[1])
          wd_list_show(wd_dirs[1])
        }else{
          working_dir<<-NULL
          rerun_iter<<-1
          rerun<<-FALSE
          shinyjs::hide("wd_list")
          updateActionButton(session,'run',"Run")
          updateRadioButtons(session,'wd_list',choices = list("None"),selected = "None")
          update_display()
        }
        
      }
    }
  )
  #Observer for working directory lists
  observeEvent(input$wd_list,{
    output$table_experment_result <- DT::renderDataTable(DT::datatable({NULL}))
    output$table_validation_result <- DT::renderDataTable(DT::datatable({NULL}))
    output$table_ensemble_result <- DT::renderDataTable(DT::datatable({NULL}))
    output$table_analysis_result <- DT::renderDataTable(DT::datatable({NULL}))
    if(input$wd_list!="None" && input$wd_list!="New Experiment"){
      wd_list_show(input$wd_list)
    }else{
      update_display()
      dir=as.character(readDirectoryInput(session, 'directory'))
      dir=normalizePath(dir)
      working_dir<<-paste0(dir,"/DMirNet_Data/")
    }
  })
  
  wd_list_show<-function(wd){
    tryCatch({
      setwd(working_dir)
      setwd(paste0(getwd(),"/",wd,"/"))
      setup_dirs_only()
      setwd(dir_ensemble)
      en_files=list.files(recursive = FALSE,pattern ="\\.csv$")
      ensemble_iter<<-length(en_files)+1
      setwd(working_dir)
      result_show("tab1")
      },error=function(cond){
        showModal(modalDialog(
          title = "ERROR!",
          paste0("Unable to change directory.Please try again. Error detail=",cond),
          easyClose=TRUE,
          fade=TRUE))
      })#end of 
  }
  #Observer for input datasets validity 
  observeEvent(input$dataset,{
    ext<-tools::file_ext(input$dataset$datapath)
    #validate file type
    if(ext!="csv"){
      showModal(modalDialog(
        title = "ERROR!",
        "The file format of the input dataset must be in .csv format. Please select the correct dataset"
      ))
    }else{
      d=read.csv(input$dataset$datapath)
      e=i=1
      if(any(is.na(d))){
        showModal(modalDialog(
          title = "ERROR!",
          "The input file contains 'NA' values. Please select the correct dataset"
        ))
      }
      coln=colnames(d)
      while(i<=length(coln)){
        if(coln[i]==" "||is.na(coln[i])||coln[i]=="X"){
          e=2
          break
        }
        i=i+1
      }
      if(e==2){
        showModal(modalDialog(
          title = "ERROR!",
          "The input file contains empty column names. Please select the correct dataset"
        ))
      }
      dat=NULL
      if(!is.null(input$dataset$datapath)){
        dat<-read.csv(input$dataset$datapath, header=TRUE, sep=",")
      }else{
        dat<-read.csv(paste0(root_dir,"/sample_datasets/Sample_31miRNAs_1151mRNAs_expr.csv"), header=TRUE, sep=",")  
      }
      l1=(1/sqrt(nrow(dat))*qnorm(1-1/(2*ncol(dat)^2)))*nrow(dat)*.161
      updateNumericInput(session,'lam1',label="lam1",value=l1,min=1,max=10)
      col=ncol(dat)
      row=nrow(dat)
      updateNumericInput(session,'weight',label="weight for each data point",value=0,min=0,max=row)
    }
  })
  observeEvent(input$validationdata,{
    ext<-tools::file_ext(input$validationdata$datapath)
    #validate file type
    if(ext!="csv"){
      showModal(modalDialog(
        title = "ERROR!",
        "The file format of the input dataset must be in .csv format. Please select the correct dataset"
      ))
    }else{
      d=read.csv(input$validationdata$datapath)
      e=i=1
      if(any(is.na(d))){
        showModal(modalDialog(
          title = "ERROR!",
          "The input file contains 'NA' values. Please select the correct dataset"
        ))
      }
      coln=colnames(d)
      while(i<=length(coln)){
        if(coln[i]==" "||is.na(coln[i])||coln[i]=="X"){
          e=2
          break
        }
        i=i+1
      }
      if(e==2){
        showModal(modalDialog(
          title = "ERROR!",
          "The input file contains empty column names. Please select the correct dataset"
        ))
      }
    }
  })
  
  #observer for display output file action
  observeEvent(c(input$dis_ana,input$dis_ens,input$dis_exp,input$dis_val),{
    if(input$dis_ana){
      display_table(dir_analysis,input$analysis_prev_file,"tab3")
    }
    if(input$dis_ens){
      display_table(dir_ensemble,input$ensemble_prev_file,"tab2")
    }
    if(input$dis_exp){
      display_table(dir_direct_bootstrap,input$ex_prev_file,"tab1")
    }
    if(input$dis_val){
      display_table(dir_validation,input$validation_prev_file,"tab4")
    }
  })
  update_display<-function(){
    updateRadioButtons(session, 'dir_boot_list',NULL,choices =list("None"),selected = FALSE)
    updateTextInput(session,"ex_prev_file","File",value = "-")
    updateTextAreaInput(session,"exp_readme","Expernment Description",value = "-")
    output$table_experment_result <- DT::renderDataTable(DT::datatable({NULL}))
  }
  
  observeEvent(input$addbtn,{
    count<<-1
    html_txt=NULL
    updateTextAreaInput(session,"host_input","List of added computers",value = "Warning: Please add all the computer's hostname of the computers")
    disable('run')
    run_disable<<-TRUE
    while (count<=input$addmore) {
      html_txt=paste(html_txt,'<input type="text" style="width=60%!important;" name="addmore',sep ="")
      html_txt=paste(html_txt,count,sep="")
      html_txt=paste(html_txt,'"class="form-control" placeholder="Enter the hostname of the computer"/><br>',sep = "")
      count<<-count+1
    }
    showModal( modalDialog(
      title = "Add the hostname of computers",
      HTML(html_txt),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok", "ADD", class="btn btn-success")
      )
       ))
  })
  
  observeEvent(input$ok,{
    removeModal()
    i=1
    hosts=NULL
    while(i<=(count-1)){
      val=paste("addmore",i,sep = "")
      if((input[[val]])==""){
        hosts=NULL
        break
      }else{
        hosts=paste(hosts,i,sep = "  ")
        hosts=paste(hosts,") ",sep = "")
        hosts=paste(hosts,as.character(input[[val]]),sep="")
      }
      i=i+1
    }
    if(is.null(hosts)){
      updateTextAreaInput(session,"host_input","List of added computers",value = "Warning: Please add all the computer's hostname of the computers")
      disable('run')
      run_disable<<-TRUE
    }else{
      enable('run')
      run_disable<<-FALSE
      updateTextAreaInput(session,"host_input","List of added computers",value = hosts)
    }
    shinyjs::show('host_input')
  })
  
})
