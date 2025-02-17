#################################################################################################################
### Title: Project 2 - Sim 1 V3
### Current Analyst: Alexander Alsup
### Last Updated: (11/12/2024) 
### Notes: 
#################################################################################################################

# Note: The function below clears the global environment 
rm(list=ls())
##################################################################################################
# PACKAGES -----
##################################################################################################
library(ggplot2)
library(compositions)
library(tidyr)
library(ggfortify)
library(openxlsx)
library(readxl)
library(stringr)
library(dplyr)
library(writexl)
  # ANCOM BC2 packages
library(ANCOMBC)
library(phyloseq)
library(DirichletReg)
##################################################################################################
# GLOBAL ARGUMENTS ----
##################################################################################################
## PATHS----
# Controls whether to save simulation results locally as an excel document
Outputs=TRUE
# Path for output document
Outputs_path = ""
##################################################################################################
# READ ME ------
##################################################################################################
## Simulation Parameters -----
# data.cells = Simulated Cell Proportions
# group = A vector
# ncell = The number of cell types used in simulation
# Nsims = Number of simulations to run per configuration
# lambda1 = The vector of Poisson coefficients used to generate count data for group 1
# lambda2 = The vector of Poisson coefficients used to generate count data for group 2
# sample.size = A vector of sample sizes to use in simulation
##################################################################################################
# FUNCTIONS ----
##################################################################################################
## Data Generation: ----
## Generates Cell Counts, Cell Proportions, and Patient Group Based on inputs
data.generation <- function(ncell=ncell,lambda1=lambda1,lambda2=lambda2,N=N){
  # Vector of "Cell Type" Names for naming
  cell.type <- paste("Cell_",seq(from=1,to=ncell),sep="")
  # Generate Count data
  pois1 <- sapply(1:ncell, function(x) rpois(n=N/2, lambda=lambda1[x]))+0.05
  pois2 <- sapply(1:ncell, function(x) rpois(n=N/2, lambda=lambda2[x]))+0.05
  data.cells <- rbind(pois1,pois2); colnames(data.cells) <- cell.type
  # Convert to Proportion data
  data.proportions <- data.cells/rowSums(data.cells)
  return(list(cell.type,data.cells,data.proportions))
}
## Effect Size Calculations -----
Poisson_effect <- function(lambda11=lambda11,lambda22=lambda22){
  Cohens_D <- (lambda11-lambda22)/(sqrt(lambda11+lambda22))
  return(Cohens_D)
}
## TEST: ANCOM Test Basic: Simple T-test equivalent -----
ANCOM.test <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- ncell
  model_test <- NULL
  model_pvals <- NULL
  model_pvals_corrected <- NULL
  for(i in 1:ncell){
    y1 <- log(data.proportions)
    y = apply(y1, 2, function(x) x - y1[,i])[,-i]  
    x <- as.factor(group) 
    lm <-summary(lm(formula = y ~ x, data = data.frame(X = x)))
    model_pvals <- c(unlist(lapply(coef(lm), function(x)x[2,4])))
    model_pvals_corrected <- p.adjust(model_pvals,"BH")
    W <- sum(model_pvals_corrected <= 0.05)
    model_test[i] <- ifelse(W>= W.star,1,0)   
  } 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="ANCOM","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## TEST: ANCOM BC2: Using the ANCOMBC Package ----
ANCOMBC2.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N,power_cell=power_cell,fdr_cell=fdr_cell){
  # Results objects
  model_test <- NULL
  model_pvals <- NULL
  model_pvals_corrected <- NULL
  # Phyloseq object
  data.proportions.ancombc2 <- data.proportions
  rownames(data.proportions.ancombc2) <- paste0("patient_",1:nrow(data.proportions.ancombc2))
  sample_data <- data.frame(sample=rownames(data.proportions.ancombc2),group=group)
  rownames(sample_data) <- rownames(data.proportions.ancombc2)
  otu_table <- data.proportions.ancombc2
  ps <- phyloseq(otu_table(otu_table,taxa_are_rows=FALSE),sample_data(sample_data))
  # ANCOM BC2
  output = ancombc2(data = ps, assay_name = "counts",
                    fix_formula = "group",
                    p_adj_method = "holm", pseudo_sens = FALSE,
                    prv_cut = 0.10, s0_perc = 0.05,
                    group = "group",
                    alpha = 0.05,
                    iter_control = list(tol = 1e-5, max_iter = 10,verbose = FALSE),
                    em_control = list(tol = 1e-5, max_iter = 10),
                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 50))
  res_prim = output$res; model_pvals=res_prim$q_groupB
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="ANCOMBC2","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## TEST: Median Test: New Test -----
  # f even
f <- Vectorize(function(m, n) {
  if (isTRUE(m==1)) { # (Optional)
    y <- n^2 * choose(2*n, n) * (1/2)^(2*n-1) / (2*n-1)
    return(y)
  }
  g <- 2 * log(n) + lchoose(2*n, n) # A common factor, potentially very large
  j <- seq_len(n) - 1               # The sum indices from 0 through n-1
  x <- (n+j) * log(m / 2) - log(n+j) + lchoose(n-1, j) + (n-1-j) * log(1-m)
  sum(exp(g + x))
}, "m")
  # Median p-value distribution for even-sized sample
f.M <- function(m, n, theta=0) {
  x <- pmax(0, pmin(1, m - theta)) # Focus on the interval [0,1]
  x <- ifelse(x > 1/2, 1-x, x)     # Apply symmetry
  f(2*x, n)*2                      # Rescale by a factor of 2
}
  # P-value for even sample size
p_value_even <- function(m,n){
  if(is.na(m)){
    NA
  } else{
    integrate(f.M,lower=0,upper=m,n=n)$value
  }
}
  # Median Test
median_test <- function(median_length=median_length,pval_median=pval_median){
  if(median_length %% 2 == 1){
    shape <- (median_length-1)/2+1
    mediantest_pvalues <- pbeta(pval_median,shape,shape,lower.tail=TRUE)
  } else {
    shape <- median_length/2
    mediantest_pvalues <- as.numeric(lapply(pval_median,function(x) p_value_even(x,shape)))
  }
  return(mediantest_pvalues)
}
# Omission function
omit_median <- function(pval_median=pval_median,alr_pvals_test=alr_pvals_test){
  omit_index <- which(pval_median==min(pval_median,na.rm=TRUE))
  alr_pvals_test[omit_index,] <- NA ; alr_pvals_test[,omit_index] <- NA 
  return(alr_pvals_test)
}
  # Full Function
Median_Test_Full <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N,cell.type=cell.type,power_cell=power_cell,fdr_cell=fdr_cell){
  # Output Objects
  model_test <- NULL
  model_pvals <- NULL
  model_pvals_corrected <- NULL
  # MEDIAN TEST
  # ALR P-values
  alr_pvals<- data.frame(matrix(NA,nrow=ncol(data.proportions),ncol=ncol(data.proportions)))
  dimnames(alr_pvals) <- list(cell.type,cell.type)
  for(i in 1:ncell){
    y1 <- log(data.proportions)
    y = apply(y1, 2, function(x) x - y1[,i])  
    x <- as.factor(group) 
    alr_lm <- summary(lm(formula = y ~ x, data = data.frame(X = x)))
    alr_pvals[i,] <- c(unlist(lapply(coef(alr_lm), function(x)x[2,4])))
  }
  diag(alr_pvals) <- NA
  # ALR P-value Correction
  alr_pvals <- apply(alr_pvals, 1,function(x) p.adjust(x,method="BH")) #Correction on the Sub-hypothesis test p-values
  # Save ALR P-value Object for test
  alr_pvals_test <- alr_pvals
  ## Median Subcomposition Test ----
  sequence_matrix <- matrix(NA,nrow=ncell,ncol=ncell)
  dimnames(sequence_matrix) <- list(paste("Step",seq(1,ncell,1)),cell.type)
  median_test_pvals <- c(0,0,0);i <- 1
  # Loop
  while(min(median_test_pvals,na.rm=TRUE)<=0.05){
    pval_median <- apply(alr_pvals_test, 1,function(x) median(x,na.rm=TRUE))
    median_length <- length(na.omit(pval_median))
    median_test_pvals <- median_test(median_length=median_length,pval_median=pval_median)
    median_test_pvals <- p.adjust(median_test_pvals,method="bonferroni") #Correction on the Median test p-values 
    sequence_matrix[i,] <- median_test_pvals
    alr_pvals_test <- omit_median(pval_median=pval_median,alr_pvals_test=alr_pvals_test)
    i <- i+1
  } 
  #sequence_matrix <- as.data.frame(sequence_matrix[rowSums(sequence_matrix,na.rm=TRUE)!=0,])
  model_test <- apply(sequence_matrix,2,function(x) min(x,na.rm=TRUE))
  model_test <- ifelse(model_test<= 0.05,1,0)
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Median","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## TEST: Current Standard Test: T-test with Correction for multiple testing and log-transformed proportions -----
Standard.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N,power_cell=power_cell,fdr_cell=fdr_cell){
  model_pvals <- NULL
  x <- as.factor(group)
  y1 <- log(data.proportions)
  for(i in 1:ncell){
    # Trad model
    lm <- summary(lm(formula = y1[,i] ~ x, data = data.frame(X = x)))
    model_pvals[i] <- coef(lm)[2,4]
  }   
  #model_pvals_corrected <- p.adjust(model_pvals,"BH")
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Standard","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## TEST: Poisson Count Test: A test using Cell Counts and a Poisson GLM ----
Poisson.test <- function(data.cells=data.cells,group=group,ncell=ncell,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- ncell
  model_pvals <- NULL
  x <- as.factor(group)
  y1 <- data.cells
  for(i in 1:ncell){
    lm <- summary(glm(formula = data.cells[,i] ~ x, data = data.frame(X = x), family = poisson(link = "log")))
    model_pvals[i] <- coef(lm)[2,4]
  }
  #model_pvals_corrected <- p.adjust(model_pvals,"BH")
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Poisson on Counts","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## TEST: Mann-Whitney U Test ----- 
  # A non-parametric test on Cell Proportions 
MWU.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- ncell
  model_pvals <- NULL
  x <- as.factor(group)
  y1 <- data.proportions
  for(i in 1:ncell){
    lm <- wilcox.test(data.proportions[,i] ~ x, exact = FALSE)
    model_pvals[i] <-lm$p.value
  }
  #model_pvals_corrected <- p.adjust(model_pvals,"BH")
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Mann-Whitney U","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## TEST: Dirichlet Regression on Proportions -----
# Description Pending 
Dirichlet.Regression <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N,power_cell=power_cell,fdr_cell=fdr_cell){
  model_pvals <- NULL
  result <- NULL
  group <- c(rep("A", N/2), rep("B", N/2))
  x <- as.factor(group)
  data <- data.frame(x=x,data.proportions)
  data$y <- DR_data(data[,2:(ncell+1)])
  result <-tryCatch({
    lm <- summary(DirichReg(y~x,data))
    coef_matrix <- lm$coef.mat; model_pvals <- coef_matrix[seq(2,nrow(coef_matrix),2),4]
    model_pvals_corrected <- p.adjust(model_pvals,"BH")
    model_test <- ifelse(model_pvals<= 0.05,1,0) 
    model_power <- sum(model_test[power_cell])
    model_fdr <- sum(model_test[(fdr_cell)])
    data.frame("Sample_Size"=N,"Method"="Dirichlet Regression","Power"=model_power,"FDR"=model_fdr)  
  },error=function(e){
    data.frame("Sample_Size"=N,"Method"="Dirichlet Regression","Power"=NA,"FDR"=NA) 
  })
  return(result)
}
## TEST: Fisher Sub-composition Scanning ------
# Fisher Statistic Calculation 
fisher_statistic_calc <- function(p_values){
  if(length(na.omit(p_values))==0){
    fisher_statistic=NA
  } else{
    p_values=p_values[!is.na(p_values)]
    fisher_statistic=-2*sum(log(p_values)) 
  }
  return(fisher_statistic)
}
# Fisher Test Full
Fisher_Test_Full <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N,cell.type=cell.type,power_cell=power_cell,fdr_cell=fdr_cell,FDR_control=FDR_control,FDR_control_type=FDR_control_type){
  # Output Objects
  model_test <- NULL
  model_pvals <- NULL
  model_pvals_corrected <- NULL
  # SUB HYPOTHESIS P-VALUES
  alr_pvals<- data.frame(matrix(NA,nrow=ncol(data.proportions),ncol=ncol(data.proportions)))
  dimnames(alr_pvals) <- list(cell.type,cell.type)
  for(i in 1:ncell){
    y1 <- log(data.proportions)
    y = apply(y1, 2, function(x) x - y1[,i])  
    x <- as.factor(group) 
    alr_lm <- summary(lm(formula = y ~ x, data = data.frame(X = x)))
    alr_pvals[i,] <- c(unlist(lapply(coef(alr_lm), function(x)x[2,4])))
  }
  diag(alr_pvals) <- NA
  alr_pvals_test <- alr_pvals # Save ALR P-value Object for test
  if(FDR_control_type == "Sub Hypothesis" & FDR_control!="None"){
    alr_pvals_test <- apply(alr_pvals, 1,function(x) p.adjust(x,method=FDR_control))
  }
  ## Subcomposition Test ----
  sequence_matrix <- matrix(NA,nrow=ncell,ncol=ncell)
  dimnames(sequence_matrix) <- list(paste("Step",seq(1,ncell,1)),cell.type)
  test_pvals <- c(0,0,0);i <- 1
  # Loop
  while(min(test_pvals,na.rm=TRUE)<=0.05){
    fisher_statistic_row <- apply(alr_pvals_test, 1,function(x) fisher_statistic_calc(x))
    df <- 2*length(na.omit(fisher_statistic_row))
    test_pvals <- pchisq(fisher_statistic_row,df,lower.tail=FALSE)
    if(FDR_control_type == "Fisher Test" & FDR_control!="None"){
      test_pvals <- p.adjust(test_pvals,method=FDR_control)
    }
    sequence_matrix[i,] <- test_pvals
    index_omit <- which(fisher_statistic_row==max(fisher_statistic_row,na.rm=TRUE))
    alr_pvals_test[index_omit,] <- NA;alr_pvals_test[,index_omit] <- NA
    i <- i+1
  } 
  #sequence_matrix <- as.data.frame(sequence_matrix[rowSums(sequence_matrix,na.rm=TRUE)!=0,])
  model_test <- apply(sequence_matrix,2,function(x) max(x,na.rm=TRUE))
  model_test <- ifelse(model_test<= 0.05,1,0)
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Fisher","Power"=model_power,"FDR"=model_fdr)
  return(result) 
}
## RESULTS: Packaging -----
results_package <- function(Abundance=NA){
  Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
  #
  sim_full=simulation[[2]]%>%
    mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2)
  Full_Results<-bind_rows(Full_Results,sim_full)
  #
  sim=simulation[[1]]%>%
    mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2)
  Results<-bind_rows(Results,sim)
}
## SIMULATION: Adjust for tests to be included -----
Simulation <- function(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell){
  time_0 <- Sys.time()
  ticks <- seq(0,Nsims,by=Nsims/100)
  K_compare <- data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA)
  Full_Results <- data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA)
  test_index= c(power_cell,fdr_cell)
  # K Level - Cycling through Sample Size
  for(k in 1:length(sample.size)){
    N <- sample.size[k]
    group <- c(rep("A", N/2), rep("B", N/2))
    J_results <- data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA)
    cat(paste0("Sample Size: ",N," \n "))
  # J Level - Simulations
    for(j in 1:Nsims){
      data <- data.generation(ncell=ncell,lambda1=lambda1,lambda2=lambda2,N=N)
      cell.type <- data[[1]];data.cells <- data[[2]];data.proportions <- data[[3]]
      J_results <- bind_rows(J_results,
                             ANCOMBC2.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N,power_cell=power_cell,fdr_cell=fdr_cell),
                             Standard.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N,power_cell=power_cell,fdr_cell=fdr_cell),
                             #Median_Test_Full(ncell=ncell,data.proportions=data.proportions,group=group,N=N,cell.type=cell.type,power_cell=power_cell,fdr_cell=fdr_cell),
                             Dirichlet.Regression(data.proportions=data.proportions,group=group,ncell=ncell,N=N,power_cell=power_cell,fdr_cell=fdr_cell),
                             Fisher_Test_Full(ncell=ncell,data.proportions=data.proportions,group=group,N=N,cell.type=cell.type,power_cell=power_cell,fdr_cell=fdr_cell,FDR_control="BH",FDR_control_type="Sub Hypothesis")
                             )
      if(j %in% ticks){
        cat(paste0("|"))
      }
    }
  # Saving Full Results
    Full_Results <- bind_rows(Full_Results,J_results)
  # Calculating Mean Power and FDR across J Simulations
    J_results <- J_results %>%
      group_by(Method)%>%
      filter(!is.na(Power))%>%
      mutate(Power=mean(Power,na.rm=TRUE),
             FDR=mean(FDR,na.rm=TRUE))%>%
      ungroup()%>% distinct()
  # Appending results to K-Level Dataframe
    K_compare <- bind_rows(K_compare,J_results)
  # Posting Progress
  #cat(paste0(round(k/length(sample.size)*100,1),"%","||"))
    cat(paste0("\n"))
  }
  time_1 <- Sys.time()
  cat(paste0("Simulation Complete. Duration ",round(time_1-time_0,2)," minutes"))
  return(list(K_compare,Full_Results))
}
#######################################################
# SIMULATIONS: Real-World Example ----
#######################################################
## Simulation: ONE DIFFERENTIALLY ABUNDANT CELL TYPE/ 5 CELL TYPES -----
### Global Simulation Parameters 
set.seed(8895)
Nsims=500
ncell=5
sample.size=c(20,40,60,80,100,120,140,160,180,200)
W.star=ceiling(ncell*0.8)
power_cell=1
fdr_cell=2
Full_Results=data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA,"Abundance"=NA,"Effect"=NA,"Effect_2"=NA)
Results=data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA,"Abundance"=NA,"Effect"=NA,"Effect_2"=NA)
### Set 1: Minor Effect Size
lambda1 <- c(2.8, 2.4, 0.44, 0.32, 0.04)
lambda2 <- c(lambda1[1]-0.5, 2.4, 0.44, 0.32, 0.04);Effect=Poisson_effect(lambda1[1],lambda2[1])
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)
  #
Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)
# Set 2: Moderate Effect Size
lambda1 <- c(2.8, 2.4, 0.44, 0.32, 0.04)
lambda2 <- c(lambda1[1]-1.0, 2.4, 0.44, 0.32, 0.04);Effect=Poisson_effect(lambda1[1],lambda2[1])
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)
#
Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)
# Set 3: Large Effect Size
lambda1 <- c(2.8, 2.4, 0.44, 0.32, 0.04)
lambda2 <- c(lambda1[1]-1.5, 2.4, 0.44, 0.32, 0.04);Effect=Poisson_effect(lambda1[1],lambda2[1])
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)
#
Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

## Final Results ----
Results <- Results %>%
  filter(!is.na(Power))%>%
  arrange(Effect,Sample_Size,Method)
Full_Results <- Full_Results %>%
  filter(!is.na(Power))%>%
  arrange(Effect,Sample_Size,Method)

## Save Final Results -----
if(Outputs==TRUE){
  data_list <- list(Results=Results,Results_Full=Full_Results)
  filename <- "Sim1_RealWorld_1Diff_5Types_V2.xlsx"
  write_xlsx(data_list,paste0(Outputs_path,filename))
}

#######################################################
# SIMULATIONS: General Form ----
#######################################################
## Simulation: ONE DIFFERENTIALLY ABUNDANT CELL TYPE/ 15 CELL TYPES -----
### Global Simulation Parameters 
set.seed(8895)
Nsims=500
ncell=15
sample.size=c(20,40,60,80,100,120,140,160,180,200)
power_cell=1
fdr_cell=14
Full_Results=data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA,"Abundance"=NA,"Effect"=NA,"Effect_2"=NA)
Results=data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA,"Abundance"=NA,"Effect"=NA,"Effect_2"=NA)

# Set 1: Low Abundance/ Small Effect Size
lambda1 = c(250          ,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(lambda1[1]-25,250,250,250,500,250,500,500,500,500,1000,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="Low",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="Low",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

# Set 2: Low Abundance/ Moderate Effect Size
lambda1 = c(250          ,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(lambda1[1]-50,250,250,250,500,250,500,500,500,500,1000,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="Low",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="Low",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

# Set 3: Low Abundance/ Moderate Effect Size
lambda1 = c(250          ,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(lambda1[1]-75,250,250,250,500,250,500,500,500,500,1000,1000,1000,1000,1000)
Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="Low",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="Low",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

### Global Simulation Parameters 
power_cell=11
fdr_cell=14

# Set 4: High Abundance/ Small Effect Size
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000          ,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,500,lambda1[11]-25,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

# Set 5: High Abundance/ Moderate Effect Size
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000           ,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,500,lambda1[11]-50,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

# Set 6: High Abundance/ High Effect Size
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000           ,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,500,lambda1[11]-75,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

## Final Results ----
Results <- Results %>%
  filter(!is.na(Power))%>%
  arrange(Effect,Sample_Size,Method)
Full_Results <- Full_Results %>%
  filter(!is.na(Power))%>%
  arrange(Effect,Sample_Size,Method)

## Save Final Results -----
if(Outputs==TRUE){
  data_list <- list(Results=Results,Results_Full=Full_Results)
  filename <- "Sim1_General_1Diff_15Types_V2.xlsx"
  write_xlsx(data_list,paste0(Outputs_path,filename))
}

#######################################################
# SIMULATIONS: General Form - Multiple DA ----
#######################################################
## Simulation: TWO DIFFERENTIALLY ABUNDANT CELL TYPE/ 15 CELL TYPES -----
### Global Simulation Parameters 
set.seed(8895)
Nsims=500
ncell=15
sample.size=c(20,40,60,80,100,120,140,160,180,200)
power_cell=11
fdr_cell=1
Full_Results=data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA,"Abundance"=NA,"Effect"=NA,"Effect_2"=NA)
Results=data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA,"Abundance"=NA,"Effect"=NA,"Effect_2"=NA)
### Global Simulation Parameters 
# Set 1: High Abundance/ Small Effect Size
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000          ,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,450,lambda1[11]-50,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)
# Set 2: High Abundance/ Moderate Effect Size
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000           ,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,450,lambda1[11]-75,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)
# Set 3: High Abundance/ High Effect Size
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000           ,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,450,lambda1[11]-100,1000,1000,1000,1000)
sim <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,power_cell=power_cell,fdr_cell=fdr_cell)

Effect=Poisson_effect(lambda1[power_cell],lambda2[power_cell]);Effect_2=paste0(lambda1[power_cell],"/",lambda2[power_cell])
sim_full=sim[[2]];sim_full<-sim_full%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Full_Results<-bind_rows(Full_Results,sim_full)
sim=sim[[1]];sim<-sim%>%mutate(Abundance="High",Effect=Effect,Effect_2=Effect_2);Results<-bind_rows(Results,sim)

## Final Results ----
Results <- Results %>%
  filter(!is.na(Power))%>%
  arrange(Effect,Sample_Size,Method)
Full_Results <- Full_Results %>%
  filter(!is.na(Power))%>%
  arrange(Effect,Sample_Size,Method)

## Save Final Results ----
if(Outputs==TRUE){
  data_list <- list(Results=Results,Results_Full=Full_Results)
  filename <- "Sim1_General_2Diff_15Types_V2.xlsx"
  write_xlsx(data_list,paste0(Outputs_path,filename))
}









