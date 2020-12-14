###------------------------------------------------------------------------------------------------------
# In the original study simulations, an apparent 'information loss' was observed, i.e. when total energy intake was adjusted for,
#the estimate was biased even in the absence of any confounding.

# We hypothesised that this is a result of combining variables, that all have distinct effects on the outcome, into one.
# Therefore, if we simulate all nutrients to have the same effects on the outcome, then we expect the estimates to be unbiased in the absence of confounding.

# The following simulations explore scenarios in which either:
#    (1) all nutrients, including the exposure, have the same effects on the outcome;
# or (2) all nutrients, excluding the exposure, have the same effects on the outcome.
###------------------------------------------------------------------------------------------------------

# Load required packages
library(progress); library(devtools); library(dagitty)

# Update dagitty directly from GitHub to the latest version if required
install_github("jtextor/dagitty/r")

###------------------------------------------------------------------------------------------------------
# DATA STRUCTURE

# Data are sumulated based on two different DAGs:
# (1) DAG: 
#     - identical to the one used in the original simulations
#     - does not change in parts I and II 
#     - each nutrient has a different causal effect on the outcome
# (2) DAG.i:
#     - in Part I, all dietary variable have equal total causal effects on the outcome (5.0kg/100kcal)
#     - in Part II, the total causal effect of the exposure NMES on the outcome is 5.00kg/100kcal, while
# the total causal effects of all other nutrients on the outcome are 3.00kg/100kcal


# All scenarios are simulated without any confounding, to make demonstrating the issues at hand easier
###------------------------------------------------------------------------------------------------------


###------------------------------------------------------------------------------------------------------
# PART I
###------------------------------------------------------------------------------------------------------

# Define data structures for simulation and set up dataframe to store means 

# in this DAG, all nutrients have different causal effect on the outcome
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


# in this DAG, all nutrients have the same causal effect on the outcome - 5.00kg/100kcal after re-scaling
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

# Set seed to ensure reproducible results
set.seed(96)

# Set up empty data frames to store means later in the simulation
simulation_means <- data.frame(mod_relative_allcomp=numeric(0), mod_relative_allcomp.i=numeric(0), mod_relative_standard=numeric(0), mod_relative_standard.i=numeric(0),
                              mod_total_allcomp=numeric(0), mod_total_allcomp.i=numeric(0), mod_total_partition=numeric(0), mod_total_partition.i=numeric(0))

###----------------------------------------------------------------------------------------------
# DATA SIMULATION

# The simulation involves the following steps:
# 1. Re-scale of all variables based on plausible mean and standard deviation values
# 2. Calculate total and residual energy intake from the simulated nutrients
# 3. Estimate the average relative causal effect of NMES on WT using:
#   3.1. The all-components model
#   3.2. The standard model
# 4. Estimate the total causal effect of NMES on WT using:
#   4.1. The all-components model
#   4.2. The energy partition model
# 5. Save exposure coefficient estimates from each model and store in a vector
###----------------------------------------------------------------------------------------------

# Run 100,000 simulations with 1000 observations of data each
Nsims <- 100000
Nobs  <- 1000

# Create simulation progress bar
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

# Run simulation
for (i in 1:Nsims) {
  
  SimData   <- simulateSEM(DAG,  N=Nobs)
  SimData.i <- simulateSEM(DAG.i, N=Nobs)
  
# Rescale variables based on plausible values (the same values used in the original simulation)
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

# Calculate total and residual energy intake
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
simulation_means[nrow(simulation_means)+1,]  <- c(mean_relative_allcomp, mean_relative_allcomp.i, mean_relative_standard, mean_relative_standard.i,
                                                mean_total_allcomp, mean_total_allcomp.i, mean_total_partition, mean_total_partition.i)


# Display simulation progress
pb$tick()

}

# Create empty dataframes to store point coefficients and simulation intervals for each model
summary_means <- data.frame(model=character(0), 
                           lower=numeric(0), 
                           point=numeric(0), 
                           upper=numeric(0))

# Create a list of model names
modelname <- list("Relative Allcomp", "Relative Allcomp.i", "Relative Standard", "Relative Standard.i",
                  "Total Allcomp", "Total Allcomp.i", "Total Partition", "Total Partition.i")

# For each model, extract the point estimate and the 95% simulation intervals and store in a vector
for (j in 1:8) {
  
  centiles                             <- c(round(quantile(simulation_means[,j], 0.025), digits=2), round(quantile(simulation_means[,j], 0.5), digits=2), round(quantile(simulation_means[,j], 0.975),digits=2))
  summary_means[nrow(summary_means)+1,]  <- c(modelname[j], unname(as.list(centiles)))
  rm(centiles)
}


# View results
View(summary_means)

# Save results (optional)
#write.csv(summary_means, "InfoLossMeans.csv")



###------------------------------------------------------------------------------------------------------
# PART II
###------------------------------------------------------------------------------------------------------

# empty the environment from Part I
rm(list=ls())

# Define data structures for simulation and set up dataframe to store means 

# in this DAG, all nutrients have different causal effect on the outcome
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


# in this DAG, all nutrients have the same causal effect on the outcome - 3.00kg/100kcal after re-scaling, except the exposure (5.00kg/100kcal)
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

# Set seed to ensure reproducible results
set.seed(96)

# Set up empty data frames to store means later in the simulation
simulation_means <- data.frame(mod_relative_allcomp=numeric(0), mod_relative_allcomp.i=numeric(0), mod_relative_standard=numeric(0), mod_relative_standard.i=numeric(0),
                              mod_total_allcomp=numeric(0), mod_total_allcomp.i=numeric(0), mod_total_partition=numeric(0), mod_total_partition.i=numeric(0))

# Run 100,000 simulations with 1000 observations of data each
Nsims <- 100000
Nobs  <- 1000

# Create simulation progress bar
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

# Run simulation
for (i in 1:Nsims) {
  
  SimData   <- simulateSEM(DAG,  N=Nobs)
  SimData.i <- simulateSEM(DAG.i, N=Nobs)
  
  # Rescale variables based on plausible values (the same values used in the original simulation)
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
  
  # Calculate total and residual energy intake
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
  simulation_means[nrow(simulation_means)+1,]  <- c(mean_relative_allcomp, mean_relative_allcomp.i, mean_relative_standard, mean_relative_standard.i,
                                                  mean_total_allcomp, mean_total_allcomp.i, mean_total_partition, mean_total_partition.i)
  
  
  # Display simulation progress
  pb$tick()
  
}

# Create empty dataframes to store point coefficients and simulation intervals for each model
summary_means_2 <- data.frame(model=character(0), 
                           lower=numeric(0), 
                           point=numeric(0), 
                           upper=numeric(0))

# Create a list of model names
modelname <- list("Relative Allcomp", "Relative Allcomp.i", "Relative Standard", "Relative Standard.i",
                  "Total Allcomp", "Total Allcomp.i", "Total Partition", "Total Partition.i")

# For each model, extract the point estimate and the 95% simulation intervals and store in a vector
for (j in 1:8) {
  
  centiles                             <- c(round(quantile(simulation_means[,j], 0.025), digits=2), round(quantile(simulation_means[,j], 0.5), digits=2), round(quantile(simulation_means[,j], 0.975),digits=2))
  summary_means_2[nrow(summary_means_2)+1,]  <- c(modelname[j], unname(as.list(centiles)))
  rm(centiles)
}

# View results
View(summary_means_2)

# Save results (optional)
#write.csv(summary_means_2, "InfoLossMeans2.csv")
