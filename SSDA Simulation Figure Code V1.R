#################################################################################################################
### Title: SSDA Simulation Figure Code V1
### Current Analyst: Alexander Alsup
### Last Updated: (11/30/2024) 
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
# Path for output documents
Outputs_path = ""
# Effect Size calculation function 
Poisson_effect <- function(lambda11=lambda11,lambda22=lambda22){
  Cohens_D <- (lambda11-lambda22)/(sqrt(lambda11+lambda22))
  return(Cohens_D)
}
#######################################################
# PLOTS ----
#######################################################
## PLOT: SIMULATION 1 REAL WORLD -----
# Data Prep
Results <- read_excel(paste0(Outputs_path,"Sim1_RealWorld_1Diff_5Types_V3.xlsx"),sheet="Results")
#
sample.size <- unique(Results$Sample_Size)
df_graph <- Results %>%
  pivot_longer(cols=c("Power","FDR"),names_to="Metric")%>%arrange(Method,Metric,Sample_Size)%>%
  #filter(Method!="Dirichlet Regression")%>%
  #mutate(Effect_2=factor(Effect_2,levels=Effect_2_levels))%>%
  mutate(CohenD=paste0("d = ",round(Effect,2)),Effect=factor(Effect))%>%
  separate(Effect_2,into=c("Cell","Effect"),sep="/")%>%
  mutate(Cell=as.numeric(Cell),Effect=as.numeric(Effect),
         Effect3=Cell-Effect)%>%
  mutate(Effect=paste0("mu = ",Effect3/50))%>%
  mutate(Method=case_when(Method=="Fisher"~"SSDA",Method=="Standard"~"Proportion Regression",TRUE~Method))
### Plot: With Mu
plot1 <- ggplot(df_graph,aes(x=Sample_Size,y=value,color=Method,linetype=Method,shape=Method))+
  geom_line(size=1.5,alpha=0.5)+
  geom_point(size=2.5,alpha=0.9)+
  labs(title="Simulation 1 Real-World Example",
       x="Sample Size",
       y="Metric Value")+
  scale_color_manual(values=c("ANCOMBC2"="#4075d6","Proportion Regression"="#631307","SSDA"="#15bf3a","Dirichlet Regression"="purple"))+
  scale_shape_manual(values=c("ANCOMBC2"=17,"Proportion Regression"=16,"SSDA"=15,"Dirichlet Regression"=4))+
  facet_grid(Effect~Metric,scales="fixed",labeller=labeller(Metric=c(Power="Power",FDR="FDR")))+
  geom_hline(yintercept=0.05,color="black",alpha=0.2)+
  scale_y_continuous(breaks=c(0,0.05,0.25,0.5,0.75,1.0),limits=c(0,1),minor_breaks=NULL)+
  scale_x_continuous(breaks=sample.size,minor_breaks=NULL)+
  #theme_prism(palette = "winter_bright", base_size = 16,base_family = "serif")+
  theme(legend.position="right",strip.text=element_text(size=12))
plot1
# Saving Plot
if(Outputs==TRUE){
  name <-"Figures/Sim1_V3.png"
  file_name<-paste0(Outputs_path,name)
  ggsave(file_name,plot=plot1, width = 7, height = 8, units="in",dpi=300)
}

## PLOT: SIMULATION 2 GENERAL FORM 1 DA -----
### PLOT: SIM 2 DIFFERENCE & ABUNDANCE (High Abundance) -----
# Data Prep
Results <- read_excel(paste0(Outputs_path,"Sim1_General_1Diff_15Types_V2.xlsx"),sheet="Results")
sample.size <- unique(Results$Sample_Size)
#
df_graph <- Results %>%
  pivot_longer(cols=c("Power","FDR"),names_to="Metric")%>%arrange(Method,Metric,Sample_Size)%>%
  #filter(Method!="Dirichlet Regression")%>%
  #mutate(Effect_2=factor(Effect_2,levels=Effect_2_levels))%>%
  mutate(CohenD=paste0("d = ",round(Effect,2)),Effect=factor(Effect))%>%
  separate(Effect_2,into=c("Cell","Effect"),sep="/")%>%
  mutate(Cell=as.numeric(Cell),Effect=as.numeric(Effect),
         Effect=case_when(Cell-Effect==25~"mu = 25",
                          Cell-Effect==50~"mu = 50",
                          Cell-Effect==75~"mu = 75",TRUE~NA),
         Abundance=case_when(Abundance=="High"~"High Abun",
                             Abundance=="Low"~"Low Abun",TRUE~NA))%>%
  mutate(Method=case_when(Method=="Fisher"~"SSDA",Method=="Standard"~"Proportion Regression",TRUE~Method))
df_graph_high <- df_graph %>% filter(Abundance=="High Abun")
df_graph_low <- df_graph %>% filter(Abundance=="Low Abun")

max <- df_graph %>%
  filter(Method=="Dirichlet Regression")

# Plot: High Abundance
caption = "Simulation 2 Results: High Abundance \n Power and False Discovery Rate"
plot1 <- ggplot(df_graph_high,aes(x=Sample_Size,y=value,color=Method,linetype=Method,shape=Method))+
  geom_line(size=1.5,alpha=0.5)+
  geom_point(size=2.5,alpha=0.9)+
  labs(title=caption,x="Sample Size",y="Metric Value")+
  scale_color_manual(values=c("ANCOMBC2"="#4075d6","Proportion Regression"="#631307","SSDA"="#15bf3a","Dirichlet Regression"="purple"))+
  scale_shape_manual(values=c("ANCOMBC2"=17,"Proportion Regression"=16,"SSDA"=15,"Dirichlet Regression"=4))+
  facet_grid(Effect~Metric,scales="fixed",labeller=labeller(Metric=c(Power="Power",FDR="FDR")))+
  geom_hline(yintercept=0.05,color="black",alpha=0.3)+
  scale_y_continuous(breaks=c(0,0.05,0.25,0.5,0.75,1.0),limits=c(0,1),minor_breaks=NULL)+
  scale_x_continuous(breaks=sample.size,minor_breaks=NULL)+
  #theme_prism(palette = "winter_bright", base_size = 16,base_family = "serif")+
  theme(legend.position="right",strip.text=element_text(size=12))
plot1

# Saving Plot
if(Outputs==TRUE){
  name <-"Figures/Sim2_High_V3.png"
  file_name<-paste0(Outputs_path,name)
  ggsave(file_name,plot=plot1, width = 7, height = 8, units="in",dpi=300)
}

# Plot: Low Abundance
caption = "Simulation 2 Results: Low Abundance \n Power and False Discovery Rate"
plot1 <- ggplot(df_graph_low,aes(x=Sample_Size,y=value,color=Method,linetype=Method,shape=Method))+
  geom_line(size=1.5,alpha=0.5)+
  geom_point(size=2.5,alpha=0.9)+
  labs(title=caption,x="Sample Size",y="Metric Value")+
  scale_color_manual(values=c("ANCOMBC2"="#4075d6","Proportion Regression"="#631307","SSDA"="#15bf3a","Dirichlet Regression"="purple"))+
  scale_shape_manual(values=c("ANCOMBC2"=17,"Proportion Regression"=16,"SSDA"=15,"Dirichlet Regression"=4))+
  facet_grid(Effect~Metric,scales="fixed",labeller=labeller(Metric=c(Power="Power",FDR="FDR")))+
  geom_hline(yintercept=0.05,color="black",alpha=0.2)+
  scale_y_continuous(breaks=c(0.05,0.25,0.5,0.75,1.0),limits=c(0,1),minor_breaks=NULL)+
  scale_x_continuous(breaks=sample.size,minor_breaks=NULL)+
  #theme_prism(palette = "winter_bright", base_size = 16,base_family = "serif")+
  theme(legend.position="right",strip.text=element_text(size=12))
plot1

# Saving Plot
if(Outputs==TRUE){
  name <-"Figures/Sim2_Low_V3.png"
  file_name<-paste0(Outputs_path,name)
  ggsave(file_name,plot=plot1, width = 7, height = 8, units="in",dpi=300)
}

### PLOT: SIM 2 EFFECT SIZE ----
# Data Prep
Results <- read_excel(paste0(Outputs_path,"Sim1_General_1Diff_15Types_V2.xlsx"),sheet="Results")
sample.size <- unique(Results$Sample_Size)
#
df_graph <- Results %>%
  pivot_longer(cols=c("Power","FDR"),names_to="Metric")%>%arrange(Method,Metric,Sample_Size)%>%
  #filter(Method!="Dirichlet Regression")%>%
  #mutate(Effect_2=factor(Effect_2,levels=Effect_2_levels))%>%
  mutate(Effect=paste0("d = ",round(Effect,2)),Effect=factor(Effect))
# Plot
caption = "Simulation 2 Results: Power and False Discovery Rate (Cohen's D)"
plot1 <- ggplot(df_graph,aes(x=Sample_Size,y=value,color=Method,linetype=Method,shape=Method))+
  geom_line(size=1.5,alpha=0.5)+
  geom_point(size=2.5,alpha=0.9)+
  labs(title=caption,x="Sample Size",y="Metric Value")+
  scale_color_manual(values=c("ANCOMBC2"="#4075d6","Standard"="#631307","Fisher"="#15bf3a","Dirichlet Regression"="purple"))+
  scale_shape_manual(values=c("ANCOMBC2"=17,"Standard"=16,"Fisher"=15,"Dirichlet Regression"=4))+
  facet_grid(Effect~Metric,scales="fixed",labeller=labeller(Metric=c(Power="Power",FDR="FDR")))+
  geom_hline(yintercept=0.05,color="black",alpha=0.2)+
  scale_y_continuous(breaks=c(0,0.05,0.25,0.5,0.75,1.0),limits=c(0,1),minor_breaks=NULL)+
  scale_x_continuous(breaks=sample.size,minor_breaks=NULL)+
  #theme_prism(palette = "winter_bright", base_size = 16,base_family = "serif")+
  theme(legend.position="right",strip.text=element_text(size=12))
plot1
# Saving Plot
if(Outputs==TRUE){
  name <-"Figures/Sim1_General_1Diff_15Types_EffectMain_V2.png"
  file_name<-paste0(Outputs_path,name)
  ggsave(file_name,plot=plot1, width = 7, height = 8, units="in",dpi=300)
}


## PLOT: SIMULATION 3 GENERAL FORM 2 DA -----
# Data Prep
Results <- read_excel(paste0(Outputs_path,"Sim1_General_2Diff_15Types_V2.xlsx"),sheet="Results")
sample.size <- unique(Results$Sample_Size)
#
df_graph <- Results %>%
  pivot_longer(cols=c("Power","FDR"),names_to="Metric")%>%arrange(Method,Metric,Sample_Size)%>%
  #filter(Method!="Dirichlet Regression")%>%
  #mutate(Effect_2=factor(Effect_2,levels=Effect_2_levels))%>%
  mutate(CohenD=paste0("d = ",round(Effect,2)),Effect=factor(Effect))%>%
  separate(Effect_2,into=c("Cell","Effect"),sep="/")%>%
  mutate(Cell=as.numeric(Cell),Effect=as.numeric(Effect),
         Effect=case_when(Cell-Effect==50~"mu = 50",
                          Cell-Effect==75~"mu = 75",
                          Cell-Effect==100~"mu = 100",TRUE~NA),
         Effect=factor(Effect,levels=c("mu = 50","mu = 75","mu = 100")),
         Abundance=case_when(Abundance=="High"~"High Abun",
                             Abundance=="Low"~"Low Abun",TRUE~NA))%>%
  mutate(Method=case_when(Method=="Fisher"~"SSDA",Method=="Standard"~"Proportion Regression",TRUE~Method))
# Plot
caption = "Simulation 3 Results: Power and False Discovery Rate"
plot1 <- ggplot(df_graph,aes(x=Sample_Size,y=value,color=Method,linetype=Method,shape=Method))+
  geom_line(size=1.5,alpha=0.5)+
  geom_point(size=2.5,alpha=0.9)+
  labs(title=caption,x="Sample Size",y="Metric Value")+
  scale_color_manual(values=c("ANCOMBC2"="#4075d6","Proportion Regression"="#631307","SSDA"="#15bf3a","Dirichlet Regression"="purple"))+
  scale_shape_manual(values=c("ANCOMBC2"=17,"Proportion Regression"=16,"SSDA"=15,"Dirichlet Regression"=4))+
  facet_grid(Effect~Metric,scales="fixed",labeller=labeller(Metric=c(Power="Power",FDR="FDR")))+
  geom_hline(yintercept=0.05,color="black",alpha=0.2)+
  scale_y_continuous(breaks=c(0,0.05,0.25,0.5,0.75,1.0),limits=c(0,1),minor_breaks=NULL)+
  scale_x_continuous(breaks=sample.size,minor_breaks=NULL)+
  #theme_prism(palette = "winter_bright", base_size = 16,base_family = "serif")+
  theme(legend.position="right",strip.text=element_text(size=12))
plot1
# Saving Plot
if(Outputs==TRUE){
  name <-"Figures/Sim1_General_2Diff_15Types_EffectMain_V2.png"
  file_name<-paste0(Outputs_path,name)
  ggsave(file_name,plot=plot1, width = 7, height = 8, units="in",dpi=300)
}

### Plot: Base Abundance and Difference
df_graph %>%
  mutate(Difference=Effect_addntl)

plot1 <- ggplot(df_graph,aes(x=Sample_Size,y=value,color=Method,linetype=Method,shape=Method))+
  geom_line(size=1.5,alpha=0.5)+
  geom_point(size=2.5,alpha=0.9)+
  labs(title="Simulation 1 Real-World Example V1",
       x="Sample Size",
       y="Metric Value")+
  scale_color_manual(values=c("ANCOMBC2"="#4075d6","Standard"="#631307","Fisher"="#15bf3a","Dirichlet Regression"="#592353"))+
  scale_shape_manual(values=c("ANCOMBC2"=17,"Standard"=16,"Fisher"=15,"Dirichlet Regression"=18))+
  facet_grid(Abundance_base+Difference~Metric,scales="fixed",labeller=labeller(Metric=c(Power="Power",FDR="FDR")))+
  geom_hline(yintercept=0.05,color="black",alpha=0.2)+
  scale_y_continuous(breaks=c(0,0.05,0.25,0.5,0.75,1.0),limits=c(0,1),minor_breaks=NULL)+
  scale_x_continuous(breaks=sample.size,minor_breaks=NULL)+
  #theme_prism(palette = "winter_bright", base_size = 16,base_family = "serif")+
  theme(legend.position="right",strip.text=element_text(size=12))
plot1


