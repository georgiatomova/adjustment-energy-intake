rm(list=ls())
#setwd("~/University of Leeds/Peter Tennant - Leeds causal inference group/Research/Projects and papers/Nutrition Adjustment")
library(progress); library(devtools); library(dagitty)
#install_github("jtextor/dagitty/r") #update dagitty to latest version if required
options(scipen=999) 

# DAG does not change in parts 1 and 2, and is the same as in the original study
# DAG.i differs between parts 1 and 2 in the following way:
  # Part 1: ALL dietary variables have equal total causal effects on the outcome (5.0kg/100kcal)
  # Part 2: The total causal effect on the outcome is 5.0kg/100kcal for the exposure and 3.0kg/100kcal for each of the other nutrients

# where .i is added as a suffix, it indicates identical causal effects among nutrients

#---PART-1------------------

############################
# Different causal effects #
############################

DAG <- dagitty("dag {
                
                NMES  -> WT   [beta = 0.25]
                CRB   -> WT   [beta = 0.33]
                FBR   -> WT   [beta = -0.02]
                UF    -> WT   [beta = 0.24]
                PRO   -> WT   [beta = 0.15]
                SF    -> WT   [beta = 0.175]
                ALC   -> WT   [beta = 0.09]
                
                U     -> NMES [beta = 0]
                U     -> CRB  [beta = 0]
                U     -> FBR  [beta = 0]
                U     -> UF   [beta = 0]
                U     -> PRO  [beta = 0]
                U     -> SF   [beta = 0]
                U     -> ALC  [beta = 0]
                
                }") 


############################
# Identical causal effects #
############################

# All path coefficients correspond to a total causal effect of 5.0kg/100kcal, once rescaled

DAG.i <- dagitty("dag {
                
                NMES  -> WT   [beta = 0.25] 
                CRB   -> WT   [beta = 0.5]  
                FBR   -> WT   [beta = 0.1]  
                UF    -> WT   [beta = 0.4]  
                PRO   -> WT   [beta = 0.3]  
                SF    -> WT   [beta = 0.25]
                ALC   -> WT   [beta = 0.1]
                
                U     -> NMES [beta = 0]
                U     -> CRB  [beta = 0]
                U     -> FBR  [beta = 0]
                U     -> UF   [beta = 0]
                U     -> PRO  [beta = 0]
                U     -> SF   [beta = 0]
                U     -> ALC  [beta = 0]
                
                }") 

set.seed(96)

# Set up dataframes to store means 
SimulationMeans <- data.frame(mod_relative_allcomp=numeric(0), mod_relative_allcomp.i=numeric(0), mod_relative_standard=numeric(0), mod_relative_standard.i=numeric(0),
                              mod_total_allcomp=numeric(0), mod_total_allcomp.i=numeric(0), mod_total_partition=numeric(0), mod_total_partition.i=numeric(0))

# Conduct simulations 
Nsims <- 100000
Nobs  <- 1000
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

for (i in 1:Nsims) {
  
  SimData   <- simulateSEM(DAG,  N=Nobs)
  SimData.i <- simulateSEM(DAG.i, N=Nobs)
  
# Rescale variables in all datasets based on plausible values
SimNames    <- c("SimData", "SimData.i")
SimDatasets <- list(SimData, SimData.i)
  
Rescaled <- lapply(SimDatasets, function(x) {
    x$NMES <- x$NMES*125+250
    x$CRB  <- x$CRB*250+500   
    x$FBR  <- x$FBR*50+100  
    x$UF   <- x$UF*200+400    
    x$PRO  <- x$PRO*150+300   
    x$SF   <- x$SF*125+275    
    x$ALC  <- x$ALC*50+175   
    x$WT   <- x$WT*25+80;
    return(x)
  })
names(Rescaled) <- paste0(SimNames)
list2env(Rescaled, envir = .GlobalEnv)

# Calculate compositional variables
SimDatasets <- list(SimData, SimData.i)
CompVar <- lapply(SimDatasets, function(y) {
  y$TotalEnergy    <- y$NMES + y$CRB + y$FBR + y$UF + y$PRO +  y$SF + y$ALC
  y$ResidualEnergy <- y$CRB + y$FBR + y$UF + y$PRO + y$SF + y$ALC;
  return(y)
})
names(CompVar) <- paste0(SimNames)
list2env(CompVar, envir = .GlobalEnv)

###########################################
### The average relative causal effect  ###
###########################################

# 1. The all-components model
mod_relative_allcomp    <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
mod_relative_allcomp.i  <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData.i)

#Calculate weights:
wi <- c(mean(SimData$CRB)/mean(SimData$ResidualEnergy),
        mean(SimData$FBR)/mean(SimData$ResidualEnergy),
        mean(SimData$UF)/mean(SimData$ResidualEnergy),
        mean(SimData$PRO)/mean(SimData$ResidualEnergy),
        mean(SimData$SF)/mean(SimData$ResidualEnergy),
        mean(SimData$ALC)/mean(SimData$ResidualEnergy))

wi.i <- c(mean(SimData.i$CRB)/mean(SimData.i$ResidualEnergy),
        mean(SimData.i$FBR)/mean(SimData.i$ResidualEnergy),
        mean(SimData.i$UF)/mean(SimData.i$ResidualEnergy),
        mean(SimData.i$PRO)/mean(SimData.i$ResidualEnergy),
        mean(SimData.i$SF)/mean(SimData.i$ResidualEnergy),
        mean(SimData.i$ALC)/mean(SimData.i$ResidualEnergy))


mean_relative_allcomp  <- mod_relative_allcomp$coefficients[2]*100-
  (wi[1]*mod_relative_allcomp$coefficients[3]*100+
     wi[2]*mod_relative_allcomp$coefficients[4]*100+
     wi[3]*mod_relative_allcomp$coefficients[5]*100+
     wi[4]*mod_relative_allcomp$coefficients[6]*100+
     wi[5]*mod_relative_allcomp$coefficients[7]*100+
     wi[6]*mod_relative_allcomp$coefficients[8]*100)

mean_relative_allcomp.i  <- mod_relative_allcomp.i$coefficients[2]*100-
  (wi.i[1]*mod_relative_allcomp.i$coefficients[3]*100+
     wi.i[2]*mod_relative_allcomp.i$coefficients[4]*100+
     wi.i[3]*mod_relative_allcomp.i$coefficients[5]*100+
     wi.i[4]*mod_relative_allcomp.i$coefficients[6]*100+
     wi.i[5]*mod_relative_allcomp.i$coefficients[7]*100+
     wi.i[6]*mod_relative_allcomp.i$coefficients[8]*100)

# 2. The standard model
mod_relative_standard    <- lm(WT ~ NMES + TotalEnergy, data=SimData)
mod_relative_standard.i  <- lm(WT ~ NMES + TotalEnergy, data=SimData.i)

mean_relative_standard   <- mod_relative_standard$coefficients[2]*100 
mean_relative_standard.i <- mod_relative_standard.i$coefficients[2]*100

################################
### The total causal effect  ###
################################

# 1. The all-components model
mod_total_allcomp    <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
mod_total_allcomp.i  <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData.i)

mean_total_allcomp   <- mod_total_allcomp$coefficients[2]*100 
mean_total_allcomp.i <- mod_total_allcomp.i$coefficients[2]*100

# 2. The energy partition model

mod_total_partition   <- lm(WT ~ NMES + ResidualEnergy, data=SimData)
mod_total_partition.i <- lm(WT ~ NMES + ResidualEnergy, data=SimData.i)

mean_total_partition   <- mod_total_partition$coefficients[2]*100 
mean_total_partition.i <- mod_total_partition.i$coefficients[2]*100

# Save local vectors of means into dataframe:
SimulationMeans[nrow(SimulationMeans)+1,]  <- c(mean_relative_allcomp, mean_relative_allcomp.i, mean_relative_standard, mean_relative_standard.i,
                                                mean_total_allcomp, mean_total_allcomp.i, mean_total_partition, mean_total_partition.i)


# Display simulation progress
pb$tick()

}

### Save median and centiles for key coefficients ####

SummaryMeans <- data.frame(model=character(0), 
                           lower=numeric(0), 
                           point=numeric(0), 
                           upper=numeric(0))

modelname <- list("Relative Allcomp", "Relative Allcomp.i", "Relative Standard", "Relative Standard.i",
                  "Total Allcomp", "Total Allcomp.i", "Total Partition", "Total Partition.i")

for (j in 1:8) {
  
  centiles                             <- c(round(quantile(SimulationMeans[,j], 0.025), digits=2), round(quantile(SimulationMeans[,j], 0.5), digits=2), round(quantile(SimulationMeans[,j], 0.975),digits=2))
  SummaryMeans[nrow(SummaryMeans)+1,]  <- c(modelname[j], unname(as.list(centiles)))
  rm(centiles)
}


View(SummaryMeans)
#write.csv(SummaryMeans, "InfoLossMeans.csv")



#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())


#---PART-2------------------

############################
# Different causal effects #
############################

DAG <- dagitty("dag {
                
                NMES  -> WT   [beta = 0.25]
                CRB   -> WT   [beta = 0.33]
                FBR   -> WT   [beta = -0.02]
                UF    -> WT   [beta = 0.24]
                PRO   -> WT   [beta = 0.15]
                SF    -> WT   [beta = 0.175]
                ALC   -> WT   [beta = 0.09]
                
                U     -> NMES [beta = 0]
                U     -> CRB  [beta = 0]
                U     -> FBR  [beta = 0]
                U     -> UF   [beta = 0]
                U     -> PRO  [beta = 0]
                U     -> SF   [beta = 0]
                U     -> ALC  [beta = 0]
                
                }") 


############################
# Identical causal effects #
############################

# All path coefficients, except the exposure, correspond to a total causal effect of 3.0kg/100kcal, once rescaled

DAG.i <- dagitty("dag {
                
                NMES  -> WT   [beta = 0.25] 
                CRB   -> WT   [beta = 0.3]  
                FBR   -> WT   [beta = 0.06]  
                UF    -> WT   [beta = 0.24]  
                PRO   -> WT   [beta = 0.18]  
                SF    -> WT   [beta = 0.15]
                ALC   -> WT   [beta = 0.06]
                
                U     -> NMES [beta = 0]
                U     -> CRB  [beta = 0]
                U     -> FBR  [beta = 0]
                U     -> UF   [beta = 0]
                U     -> PRO  [beta = 0]
                U     -> SF   [beta = 0]
                U     -> ALC  [beta = 0]
                
                }") 

set.seed(96)

# Set up dataframes to store means 
SimulationMeans <- data.frame(mod_relative_allcomp=numeric(0), mod_relative_allcomp.i=numeric(0), mod_relative_standard=numeric(0), mod_relative_standard.i=numeric(0),
                              mod_total_allcomp=numeric(0), mod_total_allcomp.i=numeric(0), mod_total_partition=numeric(0), mod_total_partition.i=numeric(0))

# Conduct simulations 
Nsims <- 100000
Nobs  <- 1000
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

for (i in 1:Nsims) {
  
  SimData   <- simulateSEM(DAG,  N=Nobs)
  SimData.i <- simulateSEM(DAG.i, N=Nobs)
  
  # Rescale variables in all datasets based on plausible values
  SimNames    <- c("SimData", "SimData.i")
  SimDatasets <- list(SimData, SimData.i)
  
  Rescaled <- lapply(SimDatasets, function(x) {
    x$NMES <- x$NMES*125+250
    x$CRB  <- x$CRB*250+500   
    x$FBR  <- x$FBR*50+100  
    x$UF   <- x$UF*200+400    
    x$PRO  <- x$PRO*150+300   
    x$SF   <- x$SF*125+275    
    x$ALC  <- x$ALC*50+175   
    x$WT   <- x$WT*25+80;
    return(x)
  })
  names(Rescaled) <- paste0(SimNames)
  list2env(Rescaled, envir = .GlobalEnv)
  
  # Calculate compositional variables
  SimDatasets <- list(SimData, SimData.i)
  CompVar <- lapply(SimDatasets, function(y) {
    y$TotalEnergy    <- y$NMES + y$CRB + y$FBR + y$UF + y$PRO +  y$SF + y$ALC
    y$ResidualEnergy <- y$CRB + y$FBR + y$UF + y$PRO + y$SF + y$ALC;
    return(y)
  })
  names(CompVar) <- paste0(SimNames)
  list2env(CompVar, envir = .GlobalEnv)
  
  ###########################################
  ### The average relative causal effect  ###
  ###########################################
  
  # 1. The all-components model
  mod_relative_allcomp    <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
  mod_relative_allcomp.i  <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData.i)
  
  #Calculate weights:
  wi <- c(mean(SimData$CRB)/mean(SimData$ResidualEnergy),
          mean(SimData$FBR)/mean(SimData$ResidualEnergy),
          mean(SimData$UF)/mean(SimData$ResidualEnergy),
          mean(SimData$PRO)/mean(SimData$ResidualEnergy),
          mean(SimData$SF)/mean(SimData$ResidualEnergy),
          mean(SimData$ALC)/mean(SimData$ResidualEnergy))
  
  wi.i <- c(mean(SimData.i$CRB)/mean(SimData.i$ResidualEnergy),
            mean(SimData.i$FBR)/mean(SimData.i$ResidualEnergy),
            mean(SimData.i$UF)/mean(SimData.i$ResidualEnergy),
            mean(SimData.i$PRO)/mean(SimData.i$ResidualEnergy),
            mean(SimData.i$SF)/mean(SimData.i$ResidualEnergy),
            mean(SimData.i$ALC)/mean(SimData.i$ResidualEnergy))
  
  
  mean_relative_allcomp  <- mod_relative_allcomp$coefficients[2]*100-
    (wi[1]*mod_relative_allcomp$coefficients[3]*100+
       wi[2]*mod_relative_allcomp$coefficients[4]*100+
       wi[3]*mod_relative_allcomp$coefficients[5]*100+
       wi[4]*mod_relative_allcomp$coefficients[6]*100+
       wi[5]*mod_relative_allcomp$coefficients[7]*100+
       wi[6]*mod_relative_allcomp$coefficients[8]*100)
  
  mean_relative_allcomp.i  <- mod_relative_allcomp.i$coefficients[2]*100-
    (wi.i[1]*mod_relative_allcomp.i$coefficients[3]*100+
       wi.i[2]*mod_relative_allcomp.i$coefficients[4]*100+
       wi.i[3]*mod_relative_allcomp.i$coefficients[5]*100+
       wi.i[4]*mod_relative_allcomp.i$coefficients[6]*100+
       wi.i[5]*mod_relative_allcomp.i$coefficients[7]*100+
       wi.i[6]*mod_relative_allcomp.i$coefficients[8]*100)
  
  # 2. The standard model
  mod_relative_standard    <- lm(WT ~ NMES + TotalEnergy, data=SimData)
  mod_relative_standard.i  <- lm(WT ~ NMES + TotalEnergy, data=SimData.i)
  
  mean_relative_standard   <- mod_relative_standard$coefficients[2]*100 
  mean_relative_standard.i <- mod_relative_standard.i$coefficients[2]*100
  
  ################################
  ### The total causal effect  ###
  ################################
  
  # 1. The all-components model
  mod_total_allcomp    <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
  mod_total_allcomp.i  <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData.i)
  
  mean_total_allcomp   <- mod_total_allcomp$coefficients[2]*100 
  mean_total_allcomp.i <- mod_total_allcomp.i$coefficients[2]*100
  
  # 2. The energy partition model
  
  mod_total_partition   <- lm(WT ~ NMES + ResidualEnergy, data=SimData)
  mod_total_partition.i <- lm(WT ~ NMES + ResidualEnergy, data=SimData.i)
  
  mean_total_partition   <- mod_total_partition$coefficients[2]*100 
  mean_total_partition.i <- mod_total_partition.i$coefficients[2]*100
  
  # Save local vectors of means into dataframe:
  SimulationMeans[nrow(SimulationMeans)+1,]  <- c(mean_relative_allcomp, mean_relative_allcomp.i, mean_relative_standard, mean_relative_standard.i,
                                                  mean_total_allcomp, mean_total_allcomp.i, mean_total_partition, mean_total_partition.i)
  
  
  # Display simulation progress
  pb$tick()
  
}

### Save median and centiles for key coefficients ####

SummaryMeans2 <- data.frame(model=character(0), 
                           lower=numeric(0), 
                           point=numeric(0), 
                           upper=numeric(0))

modelname <- list("Relative Allcomp", "Relative Allcomp.i", "Relative Standard", "Relative Standard.i",
                  "Total Allcomp", "Total Allcomp.i", "Total Partition", "Total Partition.i")

for (j in 1:8) {
  
  centiles                             <- c(round(quantile(SimulationMeans[,j], 0.025), digits=2), round(quantile(SimulationMeans[,j], 0.5), digits=2), round(quantile(SimulationMeans[,j], 0.975),digits=2))
  SummaryMeans2[nrow(SummaryMeans2)+1,]  <- c(modelname[j], unname(as.list(centiles)))
  rm(centiles)
}


View(SummaryMeans2)
#write.csv(SummaryMeans2, "InfoLossMeans2.csv")
