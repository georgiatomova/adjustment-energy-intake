# Load required packages
library(progress); library(devtools); library(dagitty)

# Update dagitty directly from GitHub to the latest version if required
install_github("jtextor/dagitty/r")

###----------------------------------------------------------------------------------------------
# DATA STRUCTURE

# Two datasets are simulated - without and with proxy confounding by common causes of diet,
# based on DAG and DAGU, respectively. 
###----------------------------------------------------------------------------------------------

# Define data structures for simulation and set up dataframe to store means 

# No confounding
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
                U     -> SF   [beta = 0]
                U     -> UF   [beta = 0]
                U     -> PRO  [beta = 0]
                U     -> ALC  [beta = 0]
                
                }") 

# Confounding
DAGU <- dagitty("dag {
                
                NMES  -> WT   [beta = 0.25]
                CRB   -> WT   [beta = 0.33]
                FBR   -> WT   [beta = -0.02]
                UF    -> WT   [beta = 0.24]
                PRO   -> WT   [beta = 0.15]
                SF    -> WT   [beta = 0.175]
                ALC   -> WT   [beta = 0.09]
                
                U     -> NMES [beta = 0.5]
                U     -> CRB  [beta = 0.25]
                U     -> FBR  [beta = -0.5]
                U     -> SF   [beta = 0.5]
                U     -> UF   [beta = 0.25]
                U     -> PRO  [beta = 0.25]
                U     -> ALC  [beta = 0.5]
                
                }") 

# Set seed to ensure reproducible results
set.seed(96)

# Set up empty data frames to store means later in the simulation
SimulationMeans <- data.frame(model0=numeric(0),
                              model0U=numeric(0),
                              model1=numeric(0),
                              model1U=numeric(0), 
                              model2=numeric(0),
                              model2U=numeric(0),
                              model3a=numeric(0),
                              model3aU=numeric(0),
                              model3b=numeric(0),
                              model3bU=numeric(0),
                              model4=numeric(0),
                              model4U=numeric(0),
                              model5a=numeric(0),
                              model5aU=numeric(0),
                              model5b=numeric(0),
                              model5bU=numeric(0))


###----------------------------------------------------------------------------------------------
# DATA SIMULATION

# The simulation involves the following steps:
# 1. Re-scale of all variables based on plausible mean and standard deviation values
# 2. Calculate total and residual energy intake from the simulated nutrients
# 3. Run the different models that adjust for energy intake: 
#   3.0. The unadjusted model
#   3.1. The energy partition model
#   3.2. The standard model
#   3.3. The nutrient density model, and the multivariable nutrient density model
#   3.4. The residual model
#   3.5. The alternative 'all-components' model
# 4. Save exposure coefficient estimates from each model and store in a vector
###----------------------------------------------------------------------------------------------

# Run 100,000 simulations with 1000 observations of data each
Nsims <- 100000
Nobs  <- 1000

# Create simulation progress bar
pb <- progress_bar$new(total = Nsims, format = ":bar :percent eta: :eta")

# Run simulation
for (i in 1:Nsims) {
  
  SimData  <- simulateSEM(DAG,  N=Nobs)
  SimUData <- simulateSEM(DAGU, N=Nobs)
  
  # Rescale variables based on plausible values, making sure that each target causal effect corresponds to the correct standardised path coefficient before re-scaling:
  
  # 5.0Kg per 100 calories = 5/25*(125/100) = 0.25 for NMES
  # 3.3Kg per 100 calories = 3/25*(250/100) = 0.33 for CRB
  # -1.0Kg per 100 calories = -1/25*(50/100) = -0.02 for FBR
  # 3.0kg per 100 calories = 3/25*(200/100) = 0.24 for UF
  # 2.5Kg per 100 calories = 2.5/25*(150/100) = 0.15 for PRO
  # 3.5Kg per 100 calories = 3.5/25*(125/100) = 0.175 for SF
  # 4.5Kg per 100 calories = (4.5/25)*(50/100) = 0.09 for ALC
  
  SimData$NMES   <- SimData$NMES*125+250  
  
  SimData$CRB    <- SimData$CRB*250+500   
  SimData$FBR    <- SimData$FBR*50+100    
  SimData$UF     <- SimData$UF*200+400    
  SimData$PRO    <- SimData$PRO*150+300   
  SimData$SF     <- SimData$SF*125+275    
  SimData$ALC    <- SimData$ALC*50+175    
  SimData$WT     <- SimData$WT*25+80  
  
  SimUData$NMES  <- SimUData$NMES*125+250
  
  SimUData$CRB   <- SimUData$CRB*250+500  
  SimUData$FBR   <- SimUData$FBR*50+100  
  SimUData$UF    <- SimUData$UF*200+400  
  SimUData$PRO   <- SimUData$PRO*150+300 
  SimUData$SF    <- SimUData$SF*125+275
  SimUData$ALC   <- SimUData$ALC*50+175
  SimUData$WT    <- SimUData$WT*25+80  
  
  # Calculate compositional variables
  
  SimData$TotalEnergy      <- SimData$NMES + SimData$CRB + SimData$FBR + SimData$UF + SimData$PRO +  SimData$SF + SimData$ALC
  SimData$ResidualEnergy   <- SimData$CRB + SimData$FBR + SimData$UF + SimData$PRO + SimData$SF + SimData$ALC
  SimUData$TotalEnergy     <- SimUData$NMES + SimUData$CRB + SimUData$FBR + SimUData$UF + SimUData$PRO + SimUData$SF + SimUData$ALC
  SimUData$ResidualEnergy  <- SimUData$CRB + SimUData$FBR + SimUData$UF + SimUData$PRO + SimUData$SF + SimUData$ALC
  
  # Run each model for energy intake adjustment  
  
### 0. The unadjusted model
  
  mod0   <- lm(WT ~ NMES, data=SimData)
  mod0U  <- lm(WT ~ NMES, data=SimUData)
  
  mean0  <- mod0$coefficients[2]*100
  mean0U <- mod0U$coefficients[2]*100
  
  
### 1. The energy partition model 
  
  mod1   <- lm(WT ~ NMES + ResidualEnergy, data=SimData)
  mod1U  <- lm(WT ~ NMES + ResidualEnergy, data=SimUData)
  
  mean1  <- mod1$coefficients[2]*100 
  mean1U <- mod1U$coefficients[2]*100
  
### 2. Тhe standard model
  mod2   <- lm(WT ~ NMES + TotalEnergy, data=SimData)
  mod2U  <- lm(WT ~ NMES + TotalEnergy, data=SimUData)
  
  mean2  <- mod2$coefficients[2]*100 
  mean2U <- mod2U$coefficients[2]*100

  
### 3. The nutrient density model
  
  # Create nutrient density variable for the exposure (as a percentage)
  SimData$NMESdens  <- SimData$NMES/SimData$TotalEnergy*100
  SimUData$NMESdens <- SimUData$NMES/SimUData$TotalEnergy*100
  
  mod3a  <- lm(WT ~ NMESdens, data=SimData);  summary(mod3a);  confint(mod3a)
  mod3aU <- lm(WT ~ NMESdens, data=SimUData); summary(mod3aU); confint(mod3aU)
  
  mean3a  <- mod3a$coefficients[2]
  mean3aU <- mod3aU$coefficients[2]
  
  # The multivariable nutrient density model
  
  mod3b  <- lm(WT ~ NMESdens + TotalEnergy, data=SimData);  summary(mod3b);  confint(mod3b)
  mod3bU <- lm(WT ~ NMESdens + TotalEnergy, data=SimUData); summary(mod3bU); confint(mod3bU)

  mean3b  <- mod3b$coefficients[2]
  mean3bU <- mod3bU$coefficients[2]
  
### 4. The residual model 
  
  Residual  <- lm(NMES ~ TotalEnergy, data=SimData); summary(Residual)
  ResidualU <- lm(NMES ~ TotalEnergy, data=SimUData); summary(ResidualU)
  
  mod4  <- lm(SimData$WT ~ Residual$residuals)  
  mod4U <- lm(SimUData$WT ~ ResidualU$residuals)
  
  mean4  <- mod4$coefficients[2]*100
  mean4U <- mod4U$coefficients[2]*100
  
### 5. The all-components model
  
  # total causal effect
  mod5a  <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
  mod5aU <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimUData)
  mean5a  <- mod5a$coefficients[2]*100
  mean5aU <- mod5a$coefficients[2]*100
  
  # average relative causal effect
  mod5b  <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
  mod5bU <- lm(WT ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimUData)
  
  #Calculate weights:
  wi <- c(mean(SimData$CRB)/mean(SimData$ResidualEnergy),
          mean(SimData$FBR)/mean(SimData$ResidualEnergy),
          mean(SimData$UF)/mean(SimData$ResidualEnergy),
          mean(SimData$PRO)/mean(SimData$ResidualEnergy),
          mean(SimData$SF)/mean(SimData$ResidualEnergy),
          mean(SimData$ALC)/mean(SimData$ResidualEnergy))
  
  mean5b  <- mod5b$coefficients[2]*100-
    (wi[1]*mod5b$coefficients[3]*100+
       wi[2]*mod5b$coefficients[4]*100+
       wi[3]*mod5b$coefficients[5]*100+
       wi[4]*mod5b$coefficients[6]*100+
       wi[5]*mod5b$coefficients[7]*100+
       wi[6]*mod5b$coefficients[8]*100)
  
  mean5bU  <- mod5bU$coefficients[2]*100-
    (wi[1]*mod5bU$coefficients[3]*100+
       wi[2]*mod5bU$coefficients[4]*100+
       wi[3]*mod5bU$coefficients[5]*100+
       wi[4]*mod5bU$coefficients[6]*100+
       wi[5]*mod5bU$coefficients[7]*100+
       wi[6]*mod5bU$coefficients[8]*100)
  
  
  
  # Save local vectors of means into dataframe:
  SimulationMeans[nrow(SimulationMeans)+1,]  <- c(mean0, mean0U, mean1, mean1U, mean2, mean2U,
                                                  mean3a, mean3aU, mean3b, mean3bU, mean4, mean4U,
                                                  mean5a, mean5aU, mean5b, mean5bU)

  rm(mod0, mod0U, mod1, mod1U, mod2, mod2U, mod3a, mod3aU, mod3b, mod3bU, mod4, mod4U, mod5a, mod5aU, mod5b, mod5bU,
     mean0, mean0U, mean1, mean1U, mean2, mean2U, mean3a, mean3aU, mean3b, mean3bU, mean4, mean4U, mean5a, mean5aU, mean5b, mean5bU, wi)

  # Display simulation progress
  pb$tick()
  
}


# Create empty dataframes to store point coefficients and simulation intervals for each model
SummaryMeans <- data.frame(model=character(0), 
                             lower=numeric(0), 
                             point=numeric(0), 
                             upper=numeric(0))

# Create a list of model names
modelname <- list("Model0", "Model0U", "Model1", "Model1U", "Model2", "Model2U", "Model3a", "Model3aU",
                  "Model3b", "Model3bU", "Model4", "Model4U", "Model5a", "Model5aU", "Model5b", "Model5bU")


# For each model, extract the point estimate and the 95% simulation intervals and store in a vector
for (j in 1:16) {
  
  centiles                             <- c(round(quantile(SimulationMeans[,j], 0.025), digits=2), round(quantile(SimulationMeans[,j], 0.5), digits=2), round(quantile(SimulationMeans[,j], 0.975),digits=2))
  SummaryMeans[nrow(SummaryMeans)+1,]  <- c(modelname[j], unname(as.list(centiles)))
  rm(centiles)
}

# View results
View(SummaryMeans)

# Save results (optional)
#write.csv(SummaryMeans, "SummaryMeans.csv")


