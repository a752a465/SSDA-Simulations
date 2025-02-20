#################################################################################################################
### Title: SubComposition Scanning for Differential Abundance testing - Functions & Demonstration
### Current Analyst: Alexander Alsup
### Last Updated: (02/17/2025) 
### Notes: 
#################################################################################################################

# Note: The function below clears the global environment 
rm(list=ls())
##################################################################################################
# PACKAGES -----
##################################################################################################
library(ggplot2)
library(tidyr)
library(dplyr)
library(lmerTest)
##################################################################################################
# READ ME ------
##################################################################################################
## ON PACKAGES -----
# Packages are included for ease of use, but are not required for the central SSDA system
# The system runs on base R when a matrix of P-values is provided

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
## Data Generation: Simple ----
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
## TEST: SSDA Function ------
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
# SSDA 
SSDA <- function(test_matrix=test_matrix){
  # Checks
  if(ncol(test_matrix) != nrow(test_matrix)){
    stop("Error: Test matrix has differing number of rows and columns")
  }
  if(isSymmetric(test_matrix)==FALSE){
    stop("Error: Test matrix is not symmetric")
  }
  # Multiple Testing Correction
  test_matrix1 <- apply(test_matrix, 1,function(x) p.adjust(x,method="BH"))
  test_matrix=test_matrix1
  # Output Preparation
  SSDA_Sequence_Matrix <- matrix(NA,nrow=ncol(test_matrix),ncol=ncol(test_matrix))
  dimnames(SSDA_Sequence_Matrix) <- list(paste("Step",seq(1,ncol(test_matrix),1)),colnames(test_matrix))
  test_pvals <- c(0,0,0);i <- 1
  # Subcomposition Scanning Loop
  while(min(test_pvals,na.rm=TRUE)<=0.05){
    # Fisher Statistic
    fisher_statistic_row <- apply(test_matrix, 1,function(x) fisher_statistic_calc(x))
    # Degrees of Freedom
    df <- 2*length(na.omit(fisher_statistic_row))
    # Fisher's Method Test
    test_pvals <- pchisq(fisher_statistic_row,df,lower.tail=FALSE)
    # Outputting results to Sequence Matrix
    SSDA_Sequence_Matrix[i,] <- test_pvals
    # Index of Ommision Cell Type
    index_omit <- which(fisher_statistic_row==max(fisher_statistic_row,na.rm=TRUE))
    # Cell type ommision
    test_matrix[index_omit,] <- NA;test_matrix[,index_omit] <- NA
    i <- i+1
  }
  SSDA_Sequence_Matrix <- as.data.frame(SSDA_Sequence_Matrix[rowSums(SSDA_Sequence_Matrix,na.rm=TRUE)!=0,])
  # Establishing DA and Non-DA Elements
  DA_indices <- which(is.na(SSDA_Sequence_Matrix[nrow(SSDA_Sequence_Matrix),]))
  non_DA_indices <- which(!is.na(SSDA_Sequence_Matrix[nrow(SSDA_Sequence_Matrix),]))
  # Combined Probality with Non-DA Matrix
  ALR_Test_Matrix_DA <- data.frame(test_matrix1);ALR_Test_Matrix_DA[,DA_indices] <- NA
  df <- 2*rowSums(!is.na(ALR_Test_Matrix_DA))
  fisher <- apply(ALR_Test_Matrix_DA, 1,function(x) fisher_statistic_calc(x))
  ALR_Test_Matrix_DA$F_statistic=fisher;ALR_Test_Matrix_DA$df=df
  ALR_Test_Matrix_DA$`Pvalue` <- pchisq(ALR_Test_Matrix_DA$F_statistic,ALR_Test_Matrix_DA$df,lower.tail=FALSE)
  #ALR_Test_Matrix_DA=ALR_Test_Matrix_DA[,-DA_indices]
  # Packaging Results
  results <- list("Sequence Matrix"=SSDA_Sequence_Matrix,"Fisher's Method With non-DA Sub-composition"=ALR_Test_Matrix_DA)
  return(results)
}

##################################################################################################
# EXAMPLES ----
##################################################################################################
## OVERVIEW ------

## SIMPLE TWO GROUP COMPARISON -----


## MULTIPLE REGRESSION MODEL -----
### Notes -----

### Data Generation -----

# Define the sample
set.seed(15155)
n <- 200
K <- 5
celltype_names <- paste0("CellType_",1:K)
Age <- rnorm(n,50,10)
Treatment <- rbinom(n,1,0.5)

# Define baseline log-abundance 
log_baseline <- c(6, 6.5, 5.5, 6.2, 6)

# Define effect sizes for Treatment (binary) and Age (continuous)
beta_treatment <- c(0,0,0.1,-0.6,0)
beta_age <- c(0.05,-0.02,0.03,-0.04,0.01)

# Generate log cell counts based on model
log_counts <- matrix(NA,nrow=n, ncol=K)
for(k in 1:K){
  log_counts[,k] <- log_baseline[k] +
    beta_treatment[k] * Treatment +
    beta_age[k] * Age +
    rnorm(n,0,0.3)
}
# Generate count data
counts <- matrix(rpois(n*K,lambda=exp(log_counts)),nrow=n,ncol=K)
colnames(counts) <- celltype_names

# Convert Counts to proportions
data.proportions <- counts/rowSums(counts)

### ALR Hypothesis Testing -----

#### Manual Hypothesis Testing ------

# Let's start by defining a simple function which will give us any ALR we want
# In the function below, if we set k=1 and k=2, we will be provided with the ALR
# alr_12 = log(y_1/y_2)
alr_kk <- function(k=k,k_prime=k_prime,data=data){
  data <- log(data)
  y = apply(data, 2, function(x) x - data[,k])[,k_prime]
  return(y)
}
# Example 
  #alr_kk(k=1,k_prime=2,data=data.proportions)

# Now we can build a matrix of alr hypothesis test p_values
ALR_Test_Matrix <- matrix(data=NA,nrow=K,ncol=K)

# Let's Start by Looking at the hypothesis test for the alr_12
model_12 <- glm(alr_kk(k=1,k_prime=2,data=data.proportions)~Treatment+Age)
  #summary(model_12)

# For this seed, the p-value for the treatment coefficient should be 0.5659
# We'll place that p-value in our ALR Test Matrix
# Note, the p-value for alr_12 is the same as alr_21
ALR_Test_Matrix[1,2] <- 0.5659 ; ALR_Test_Matrix[2,1] <- 0.5659 

# However, we can notice that there is a typical location for our p-value of interest in these models
summary(model_12)$coefficients[2,4]

#### Automatic Hypothesis Testing ------
# Instead of Manually extracting our p-values, we can use a for-loop to calculate and store them
ALR_Test_Matrix2 <- matrix(data=NA,nrow=K,ncol=K)
colnames(ALR_Test_Matrix2)=celltype_names;rownames(ALR_Test_Matrix2)=celltype_names
for(k in 1:K){
  for(k_prime in 1:K){
    model <- glm(alr_kk(k=k,k_prime=k_prime,data=data.proportions)~Treatment+Age)
    ALR_Test_Matrix2[k,k_prime] <- summary(model)$coefficients[2,4]
  }
}
diag(ALR_Test_Matrix2) <- NA

### SSDA -----
# Now we can apply SSDA. 
# In our simulation, we specified that cell type 3 and cell type 4 would be differentially abundant based on Treatment. 
# Now we can test that using the SSDA() function.  

SSDA(test_matrix=ALR_Test_Matrix2)

# Here, we can observe in the sequence matrix ($`Sequence Matrix`) that cell type 3 and cell type 4 were omitted in steps 2 and 3 respectively
# Comparing the ALRs of all cell types against the non-DA matrix (in $`Fisher's Method With non-DA Sub-composition`),
# We can observethat only cell types 3 and 4 have a statistically significant Fisher's Method p-value as well

## MIXED EFFECTS MODEL ------
### Notes -----

### Data Generation -----

# Define the sample
set.seed(15155)
n_patients <- 80
n_samples_per_pt <- 4
K <- 5
celltype_names <- paste0("CellType_",1:K)

# Patient data
patient_data <- data.frame(
  Patient=rep(1:n_patients,each=n_samples_per_pt),
  Observation=rep(1:n_samples_per_pt,n_patients),
  Age=rep(rnorm(n_patients,mean=50,sd=10),each=n_samples_per_pt),
  Treatment=rep(rbinom(n_patients,1,0.5),each=n_samples_per_pt)
)

# Define baseline log-abundance 
log_baseline <- c(6, 6.5, 5.5, 6.2, 6)

# Define effect sizes for Treatment (binary) and Age (continuous)
beta_treatment <- c(0,0,0.2,-0.3,0)
beta_age <- c(0.05,-0.02,0.03,0.02,0.01)

# Random intercept per cell type per patient
intercepts <- rep(rnorm(n_patients*K, mean=0,sd=0.2),each=n_samples_per_pt)
random_intercepts <- matrix(intercepts,nrow=n_patients*n_samples_per_pt,ncol=K)

# Generate log cell counts based on model
log_counts <- matrix(NA,nrow=n_patients*n_samples_per_pt, ncol=K)
for(k in 1:K){
  log_counts[,k] <-log_baseline[k] +
    beta_treatment[k] * patient_data$Treatment +
    beta_age[k] * patient_data$Age +
    random_intercepts[,k]+
    rnorm(n_patients*n_samples_per_pt,0,0.02)
}

# Generate count data
counts <- matrix(rpois(n_patients*n_samples_per_pt*K,lambda=exp(log_counts)),nrow=n_patients*n_samples_per_pt,ncol=K)
for(i in 1:n_patients*n_samples_per_pt){
  counts[i,]  <- rpois(K,lambda=exp(log_counts[i,]))
}

# Convert Counts to proportions
data.proportions <- counts/rowSums(counts)

### ALR Hypothesis Testing (Automatic) -----
ALR_Test_Matrix2 <- matrix(data=NA,nrow=K,ncol=K)
colnames(ALR_Test_Matrix2)=celltype_names;rownames(ALR_Test_Matrix2)=celltype_names
for(k in 1:K){
  for(k_prime in 1:K){
    if(k!=k_prime){
      df_analysis <- data.frame(patient_data,"ALR"=alr_kk(k=k,k_prime=k_prime,data=data.proportions))
      ALR_Test_Matrix2[k,k_prime] <- summary(lmer(ALR~Treatment+Age+(1| Patient),data=df_analysis))$coefficients[2,5]
    }
  }
}
diag(ALR_Test_Matrix2) <- NA

### SSDA -----
# Now we can apply SSDA. 
# In our simulation, we specified that cell type 3 and cell type 4 would be differentially abundant based on Treatment. 
# Now we can test that using the SSDA() function.  

SSDA(test_matrix=ALR_Test_Matrix2)








