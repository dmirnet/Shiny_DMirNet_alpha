
#Attempt to detect the number of CPU cores on the current host
noCores<-detectCores(all.tests = FALSE, logical = TRUE)-1 # to avoid the computer to stand still until R task finishes 
shinyUI( 
  fluidPage(
    shinyjs::useShinyjs(),
    id="page",
    theme=shinytheme("flatly"),
    navbarPage("DMirNet"),
    fluidRow(
      div(
        tabsetPanel(type = "tabs",
                    id="tabpanel",
                    tabPanel("Experiments", 
                             column(6,
                                    
                                    h4("Experiment"),
                                    directoryInput('directory', label = 'select a working directory', value = '~'),
                                    hidden(
                                      h4("List of Experments"),
                                      radioButtons('wd_list',"Choose Experiment",choices = c("None"))
                                    ),
                                    fileInput("dataset",label = "*Select Dataset",placeholder = "Sample_31miRNAs_1151mRNAs_expr.csv",
                                              accept=c('.xls','.xlsx','.csv'),width = '65%'
                                    ),
                                    radioButtons("parallel_method","Choose a parallelism method for running the experiment",choiceNames = c("Cluster of Cores","Cluster of Computers"),choiceValues = c("core","server"), selected = "core"),
                                    hidden(
                                      div(id="no_cores",
                                      selectInput("cores","Number of cores",choices =(1:noCores),selected = "1",multiple = FALSE,
                                                  selectize = TRUE, width = '60%')
                                      )
                                    ),
                                    hidden(
                                      div(id="no_machines",
                                        tags$div( 
                                        class="input-group control-group after-add-more",
                                        tags$p("Enter the number of Computers"),
                                        tags$input(type="number",required="required", style="width:70%!important;",value="1", name= "addmore",class="form-control"),
                                        br(),
                                        tags$div(  
                                         actionLink("addbtn","Add Computers",class="btn btn-default")
                                          ),
                                        br()
                                        )
                                       ),
                                        disabled(
                                          textAreaInput("host_input","List of added computers",value = "Warning: Please add all the computers hostname",width ="150%", rows =4)
                                        )
                                      
                                    )
                             ),
                             column(3,
                                    height="30%",
                                    h4("Direct Correlation Methods"),
                                    checkboxInput('space', 'Space',TRUE),
                                    checkboxInput('corND', 'CorND'),
                                    checkboxInput('corpcor', 'Corpcor'),
                                    checkboxInput('ida', 'IDA'),
                                    br(),
                                    checkboxInput('dir_corl_edit', 'Edit the parameters for the selected direct correlation methods')
                             ),
                             column(3,
                                    h4("Bootstrapping"),
                                    radioButtons("bootstrap_action","Choose bootstrapping action",choices = c("Enable","Disable"),selected = "Disable"),
                                    hidden(
                                      div( id="boot_div",
                                           numericInput("iteration","Iteration","50",min =2, max=100,000),
                                           numericInput("seed","Seed",value = 123),
                                           sliderInput("rate","Sampling rate (%)","98",min=1,max=100),
                                           radioButtons('bootstrap_method','Select Bootstrapping Method',choiceNames = c("Simple aggregation-Mean","Simple aggregation-Median","Top k aggregation-Mean","Top k aggregation-Median"), choiceValues = c("bootstrap_irm_mean","bootstrap_irm_median","bootstrap_topk_mean","bootstrap_topk_median"),selected = "bootstrap_irm_mean"),
                                           numericInput("bootstrap_topk","Number of Top K","100",min =1, max=1189)
                                      )
                                    )
                             ),
                             hidden(
                             div(
                            id="edit_div",
                             column(12,br(),
                                    hr(),
                             column(3, 
                                    hidden(
                                      div(
                                        id="divspace",
                                        h4("Parameters selection for Space (Optional)"),
                                        numericInput("lam1","lam1",value="0",min=0.0,max=10),
                                        numericInput("lam2","lam2",value="0",min=0.0,max=10),
                                        numericInput("w","Weight",value="2",min=0,max=10),
                                        numericInput("iter","Iter",value="3",min=2,max=100)
                                      )
                                    )
                             ),
                             column(3, 
                                    hidden(
                                      div(
                                        id="divcorND",
                                        h4("Parameters selection for CorND (Optional)"),
                                        numericInput("corND_alpha","Alpha","1",min=0,max=10),
                                        sliderInput("corND_beta","Beta (0-1)","0.9",min=0.0,max=1)
                                        
                                      )
                                    )
                             ),
                             column(3, 
                                    hidden(
                                      div(
                                        id="divcorp",
                                        h4("Parameter selection for Corpcor (Optional)"),
                                        numericInput("lambda","Lambda (0-1)(optional)",value=0,min=0.0,max=1)
                                      )
                                    )
                             ),
                             column(3,
                                    hidden(
                                      div(
                                        id="divida",
                                        h4("Parameters selection for IDA (Optional)"),
                                        numericInput("alpha_idea","Alpha","0.01",min=0,max=10)
                                      )
                                    )
                             )
                             )
                             )
                             ),
                             column(12, style="margin-left:50%",
                                actionButton('run',"Run", class="btn btn-success", icon = NULL, width = NULL)
                             ),
                             div(
                               column(12,hr(), h4("Output Preview"),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 300px",
                                               h4("Direct Correlation and Bootstrapping Result"),
                                               radioButtons('dir_boot_list',NULL,choices =list("None"),selected = FALSE)
                                             )
                                      ),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 300px",
                                               textAreaInput("exp_readme","Expernment Description",cols =15, rows = 10)
                                             )
                                      ),
                                      column(6,
                                             
                                             h4("Preview"),
                                             disabled(
                                               textInput("ex_prev_file","File",value = "-",width =1000)
                                             ),
                                             actionButton("dis_exp","Display File in table"),  
                                             br(),br(),
                                             div(
                                               style = "overflow-y:scroll; max-height: 900px",
                                               DT::dataTableOutput("table_experment_result")
                                             ),
                                             br(),br()
                                      )
                               )
                             )
                             
                    ),
                    tabPanel("Ensemble",
                             id="ensemble_tab",
                             column(6,
                                    style = "overflow-y:scroll; max-height: 350px",
                                    div(
                                      h4("Select File"),
                                      checkboxGroupInput('ag_list',NULL,choices =list("NULL")) 
                                    )
                             ),
                             column(6,
                                    id="ensemble_2",
                                    radioButtons('ensemble','Ensemble Method ',choiceNames = c("Simple aggregation – Mean","Simple aggregation – Median","Top k aggregation – Mean","Top k aggregation – Median"), choiceValues= c("ensemble_irm_mean","ensemble_irm_median","ensemble_topk_mean","ensemble_topk_median"),selected = "ensemble_irm_mean"),
                                    numericInput("irm_topk","Number of Top K","10000",min =1),
                                    actionButton('run_ensemble',"Ensemble", class="btn btn-success")
                                   
                             ),
                             div(
                               column(12,hr(), h4("Output Preview"),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 320px",
                                               h4("Ensemble Result"),
                                               radioButtons('ensemble_result_list',NULL,choices =list("None"),selected = FALSE)
                                             )
                                      ),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 300px",
                                               textAreaInput("ensemble_readme","Expernment Description",cols =15, rows = 10)
                                             )
                                      ),
                                      column(6,
                                             
                                             h4("Preview"),
                                             disabled(
                                               textInput("ensemble_prev_file","File",value = "-",width =1000)
                                             ),
                                             actionButton("dis_ens","Display File in table"), 
                                             br(),br(),
                                             div(
                                               style = "overflow-y:scroll; max-height: 900px",
                                               DT::dataTableOutput("table_ensemble_result")
                                             ),
                                             br(),br()
                                      )
                               )
                             )
                             
                    ),
                    tabPanel("Result Analysis",
                             column(6,
                                    style = "overflow-y:scroll; max-height: 320px",
                                    div(
                                      id="result_analysis",
                                      h4("Select file for Analysis"),
                                      radioButtons('analysis_file',NULL,choices =list("NULL"))
                                    )
                             ),
                             column(6,
                                    height="30%",
                                    numericInput("miRNAs","Number of miRNAs","10",min =1, max=1189),
                                    radioButtons("analysis_criteria","Criterion for Selecting Highly-Ranked Pairs",choiceNames = c("Select Top m pairs:","Select Top m mRNAs for each miRNA"),choiceValues=c("analysis_mpairs","analysis_mmpairs"),selected="analysis_mpairs"),
                                    numericInput("analysis_topk","Top m:","1000",min =1, max=1189),
                                    actionButton('run_analysis',"Run", class="btn btn-success")
                             ),
                             div(
                               column(12,hr(), h4("Output Preview"),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 320px",
                                               h4("Analysis Result"),
                                               radioButtons('analysis_result_list',NULL,choices =list("None"),selected = FALSE)
                                             )
                                      ),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 300px",
                                               textAreaInput("analysis_readme","Expernment Description",cols =15, rows = 10)
                                             )
                                      ),
                                      column(6,
                                             
                                             h4("Preview"),
                                             disabled(
                                               textInput("analysis_prev_file","File",value = "-",width =1000)
                                             ),
                                             actionButton("dis_ana","Display File in table"), 
                                             br(),br(),
                                             div(
                                               style = "overflow-y:scroll; max-height: 900px",
                                               DT::dataTableOutput("table_analysis_result")
                                             ),
                                             br(),br()
                                      )
                               )
                             )
                             
                    ),
                    tabPanel("Validation",
                             column(6,
                                    style = "overflow-y:scroll; max-height: 320px",
                                    div(
                                      id="result_validation",
                                      h4("Select which output file to validate"),
                                      radioButtons('validation_file',NULL,choices =list("NULL"))
                                    )
                             ),
                             column(6,
                                    height="30%",
                                    fileInput("validationdata",label = "Select a Dataset for Validating the Results",placeholder = "Sample_ValidationDataset.csv",
                                              accept=c('.csv','.xls','.xlsx')
                                              
                                    ),
                                    actionButton('run_validation',"Run", class="btn btn-success")
                             ),
                             div(
                               column(12,hr(), h4("Output Preview"),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 300px",
                                               h4("Validation Result"),
                                               radioButtons('validation_result_list',NULL,choices =list("None"),selected = FALSE)
                                             )
                                      ),
                                      column(3,
                                             div(
                                               style = "overflow-y:scroll; max-height: 300px",
                                               textAreaInput("validation_readme","Expernment Description",cols =15, rows = 10)
                                             )
                                      ),
                                      column(6,
                                             
                                             h4("Preview"),
                                             disabled(
                                               textInput("validation_prev_file","File",value = "-",width =1000)
                                             ),
                                             actionButton("dis_val","Display File in table"), 
                                             br(),br(),
                                             div(
                                               style = "overflow-y:scroll; max-height: 900px",
                                               DT::dataTableOutput("table_validation_result")
                                             ),
                                             br(),br()
                                      )
                               )
                             )
                    ),
                    selected = "Experiments"
        )
      )
    )
  )
)
