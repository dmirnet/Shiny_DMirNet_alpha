##### DMirNet Common Functions #####

# Matrix normalization based on absolutle values (range: -1 ~ 1)
norm_mat <- function(mat){
  
  abs_mat <- abs(mat)
  min_mat <- min(abs_mat[row(abs_mat)!=col(abs_mat)])
  max_mat <- max(abs_mat[row(abs_mat)!=col(abs_mat)])
  scale_abs_mat <- (abs_mat-min_mat)/(max_mat-min_mat)
  result <- sign(mat)*scale_abs_mat
  diag(result) <- 1
  
  return(result)
}

### Ensemble aggregation method - IRM ###
# Inversed Rank Aggregation with Mean, Median
IRM <- function(adjacencyList, method = c("mean", "median")){
  
  # Change the correlation values by the ranking
  abs_rank_list <- apply(abs(adjacencyList), 2, function(x){rank(x,ties.method = "min")})
  cor_rank_list <- abs_rank_list*sign(adjacencyList) 
  
  # Aggregate with user specified method (default: mean)
  if(method == "median"){
    agg_result <- apply(cor_rank_list, 1, median)}
  else {
    agg_result <- apply(cor_rank_list, 1, mean)}
  
  # Linearly mapping the Ensemble aggregation result to be between -1 and 1
  abs_agg_result <- abs(agg_result)
  scale_agg_result <- (abs_agg_result-min(abs_agg_result))/(max(abs_agg_result)-min(abs_agg_result))
  result <- sign(agg_result)*scale_agg_result
  
  return(result)
  
}


## Inversed Rank Aggregation with Mean, Median (top k version)
IRM_topk <- function(adjacencyList, topk, method = c("mean", "median")){
  
  cutoff <- nrow(adjacencyList)-topk
  # Change the correlation values by the ranking
  abs_rank_list <- apply(abs(adjacencyList), 2, function(x){ifelse(rank(x,ties.method = "max")<cutoff, 0, rank(x,ties.method = "min"))})
  cor_rank_list <- abs_rank_list*sign(adjacencyList) 
  
  # Aggregate with user specified method (default: mean)
  if(method == "median"){
    agg_result <- apply(cor_rank_list, 1, median)}
  else {
    agg_result <- apply(cor_rank_list, 1, mean)}
  
  # Linearly mapping the Ensemble aggregation result to be between -1 and 1
  abs_agg_result <- abs(agg_result)
  scale_agg_result <- (abs_agg_result-min(abs_agg_result))/(max(abs_agg_result)-min(abs_agg_result))
  result <- sign(agg_result)*scale_agg_result
  
  return(result)
  
}


### Two methods for selecting Top K interactions ###

## Extract top k associations of a miRNA
## It returns k*(number of miRNA) miRNA-mRNA pairs based on the absolute value of direct correlation coefficients
Extopk_mimR_miRk <- function(cormat,topk,num_miR){
  mirname<-colnames(cormat)
  row.names(cormat)=mirname
  mirname <- mirname[1:num_miR]
  Result<-matrix(-1,ncol=5,nrow=num_miR*topk)
  n<-0
  for (i in 1:num_miR){ # for each miRNA
    corlist=cormat[(num_miR+1):nrow(cormat),i] # extract the correlation values
    corlist=corlist[order(-abs(corlist))]
    corlist=corlist[1:topk]
    rn=names(corlist)
    for (j in 1:topk){
      n <- (n+1)
      Result[n, 1] <- mirname[i] # miR
      Result[n, 2] <- rn[j] # mR
      Result[n, 3] <- corlist[j] # correlation
      Result[n, 4] <- abs(as.numeric(Result[n, 3])) # Calculate absolute value
      Result[n, 5] <- sign(corlist[j])
    }
  }
  rownames(Result)=NULL
  colnames(Result)=c("miRNA","mRNA","Correlation","abs_corr","Sign")
  Sorted_Result <- Result[sort.list(as.numeric(Result[,4]),decreasing=TRUE),] #Rank the whole interations by decreasing the absolute value
  
  return(Sorted_Result) # Return top k miRNA-mRNA interactions for each miRNA
}      

## Extract top k miR-mR associations
## It returns k miRNA-mRNA pairs based on the absolute value of direct correlation coefficients
Extopk_mimR <- function(cormat,topk,num_miR){
  mir<-ncol(cormat)
  gene<-nrow(cormat)
  mirname<-colnames(cormat)
  genename<-mirname
  
  Result<-matrix(-1,ncol=5,nrow=num_miR*(mir-num_miR))
  n<-0
  for(i in 1:num_miR){
    for(j in (num_miR+1):mir){
      n<-(n+1)
      Result[n, 1] <- mirname[i] # miR
      Result[n, 2] <- genename[j] # mR
      Result[n, 3] <- cormat[j,i] # correlation
      Result[n, 4] <- abs(as.numeric(Result[n, 3])) # Calculate absolute value
      Result[n, 5] <- sign(cormat[j,i])
    }
  } 
  
  rownames(Result)=NULL
  colnames(Result)=c("miRNA","mRNA","Correlation","abs_corr","Sign")
  SortedResult <- Result[sort.list(as.numeric(Result[,4]),decreasing=TRUE),] #Rank the whole interations by decreasing the absolute value
  
  return(SortedResult[1:topk,]) # Return top k miRNA-mRNA interactions
  
}     


### Validation ### 

## Validate the targets of a miRNA

Validation=function(topkList, valdata_csv){
  
  # Preprocessing the list of topk results 
  topkList = as.matrix(topkList, ncol=3);
  if(nrow(topkList)<1) stop("No data in the input")
  
  # Read the validation data from file 
  dt=valdata_csv
  
  # Get the validation lists from the data
  dt = paste(dt[, 1], dt[, 2], sep=" ");
  tmp= paste(topkList[, 1], topkList[, 2], sep=" ");  
  tmp<-gsub("\\.", "-", tmp)
  result=topkList[which(tmp %in% dt), ] 
  NoConfimred=NULL
  if(is.matrix(result)){
    result=data.frame(result)
    result[,3]=as.numeric(as.character(result[,3]))
    NoConfimred=nrow(result)
  } else {
    result=as.matrix(result)
    result=t(result)
    result=data.frame(result)
    result[,3]=as.numeric(as.character(result[,3]))
    NoConfimred=nrow(result)
  }
  return (result)
  #return(list(result, NoConfimred))
}


### Matrix conversion function ###

## Convert from an upper triangular vector to a full matrix.
uptri2mat <- function(genes, tri, diag=1){
  mat <- matrix(0, nrow=length(genes), ncol=length(genes))
  colnames(mat) <- rownames(mat) <- genes
  diag(mat) <- diag
  
  #extract the indices of all upper-triangular elements
  uti <- matrix(c(row(mat)[upper.tri(mat)], col(mat)[upper.tri(mat)]), ncol = 2)
  mat[uti] <- tri
  
  # make a matrix symmetric
  mat[matrix(c(col(mat)[upper.tri(mat)], row(mat)[upper.tri(mat)]), ncol = 2)] <- mat[matrix(c(row(mat)[upper.tri(mat)], col(mat)[upper.tri(mat)]), ncol = 2)]
  
  return(mat)
}
### Bootstrapping experiment function ###

# Bootstrap the reconstruction of a network
bootstrap <- function(data, direct_fun, ensemble_fun, sample.percentage, iterations,topk,params){
  if (typeof(direct_fun) != "character"){
    stop("You must provide the character name of the function you want to bootstrap. For instance, fun=\"buildSpace\"")
  }
  fun <- get(direct_fun)
  if (typeof(ensemble_fun) != "character"){
    stop("You must provide the character name of the ensemble method you want to ensemble. For instance, fun=\"IRM\"")
  }
  rand.seed=seed_num
  # Bootstrapping function
  funWrapper <- function(rand.seed, fun, data, sample.percentage,params, ...){
    set.seed(rand.seed)
    sampledData <- data[sample(1:nrow(data),round(sample.percentage * nrow(data))),]
    net=NULL
    if(direct_fun=="buildCorpcor"){
      net=fun(sampledData,params[1],0)
    }else if(direct_fun=="buildSpace"){
      net=fun(sampledData,params[1],params[2],params[3],params[4],0)
    }else if(direct_fun=="buildCorND"){
      net=fun(sampledData,params[1],params[2],0)
    }else if(direct_fun=="buildIDA"){
      net=fun(sampledData,params[1],0)
    }
    return(net[upper.tri(net)])
  }
  cl=makePSOCKcluster(clusters,setup_timeout=120)
  setDefaultCluster(cl)
  vars=list(direct_fun,"norm_mat","write_file","dir_direct_bootstrap","dir_direct_bootstrap_uppertri","pcor.shrink","space.joint","pc_stable","gaussCItest","getNextSet","udag2pdagRelaxed","idaFast")
  clusterExport(cl, vars, envir = .GlobalEnv)
  result <- parLapply(cl, 1:iterations, funWrapper, fun, data, sample.percentage,params)	
  stopCluster(cl)
  # Run bootstrapping experiment
  #result <- lapply(1:iterations, funWrapper, fun, data, sample.percentage,params)		
  result <- as.data.frame(result)
  boot_ens_result=NULL
  # Ensemble aggregation
  if(ensemble_fun=="bootstrap_irm_mean"){
    ensfun <- get("IRM")
    boot_ens_result <- ensfun(result,"mean")
  }else if(ensemble_fun=="bootstrap_irm_median"){
    ensfun <- get("IRM")
    boot_ens_result <- ensfun(result,"median")
  }else if(ensemble_fun=="bootstrap_topk_mean"){
    ensfun <- get("IRM_topk")
    boot_ens_result <- ensfun(result,topk,"mean")
  }else{
    ensfun <- get("IRM_topk")
    boot_ens_result <- ensfun(result,topk,"median")
  }
  return(boot_ens_result)
}
