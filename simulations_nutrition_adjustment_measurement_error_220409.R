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
                
                NMES_T  -> GLUC [beta = 0.25]
                CRB_T   -> GLUC [beta = 0.33]
                FBR_T   -> GLUC [beta = -0.02]
                UF_T    -> GLUC [beta = 0.24]
                PRO_T   -> GLUC [beta = 0.15]
                SF_T    -> GLUC [beta = 0.175]
                ALC_T   -> GLUC [beta = 0.09]
              
                error   ->  NMES_e [beta = 0.5] 
                error   ->  CRB_e  [beta = 0.5]
                error   ->  FBR_e  [beta = 0.5]
                error   ->  UF_e   [beta = 0.5]
                error   ->  PRO_e  [beta = 0.5]
                error   ->  SF_e   [beta = 0.5]
                error   ->  ALC_e  [beta = 0.5]
                
                U     -> NMES_T [beta = 0]
                U     -> CRB_T  [beta = 0]
                U     -> FBR_T  [beta = 0]
                U     -> SF_T   [beta = 0]
                U     -> UF_T   [beta = 0]
                U     -> PRO_T  [beta = 0]
                U     -> ALC_T  [beta = 0]

                }") 

# Confounding
DAGU <- dagitty("dag {
                
                NMES_T  -> GLUC [beta = 0.25]
                CRB_T   -> GLUC [beta = 0.33]
                FBR_T   -> GLUC [beta = -0.02]
                UF_T    -> GLUC [beta = 0.24]
                PRO_T   -> GLUC [beta = 0.15]
                SF_T    -> GLUC [beta = 0.175]
                ALC_T   -> GLUC [beta = 0.09]
                
                error   ->  NMES_e [beta = 0.5] 
                error   ->  CRB_e  [beta = 0.5]
                error   ->  FBR_e  [beta = 0.5]
                error   ->  UF_e   [beta = 0.5]
                error   ->  PRO_e  [beta = 0.5]
                error   ->  SF_e   [beta = 0.5]
                error   ->  ALC_e  [beta = 0.5]
                
                U     -> NMES_T [beta = 0.5]
                U     -> CRB_T  [beta = 0.25]
                U     -> FBR_T  [beta = -0.5]
                U     -> SF_T   [beta = 0.5]
                U     -> UF_T   [beta = 0.25]
                U     -> PRO_T  [beta = 0.25]
                U     -> ALC_T  [beta = 0.5]
                
                U     -> error  [beta = 0.5] 
                
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
# 2. Calculate total and remaining energy intake from the simulated nutrients
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
  
  # Construct 'observed' variables by adding 66% of the variance of the true value and 33% of the variance of the error term, to produce the same final mean and SDs as the original simulation:
  
  SimData$NMES   <- (SimData$NMES_T*102.01 + SimData$NMES_e*72.13)+250   
  SimData$CRB    <- (SimData$CRB_T*204.02  + SimData$CRB_e*144.26)+500   
  SimData$FBR    <- (SimData$FBR_T*40.80   + SimData$FBR_e*28.85)+100 
  SimData$UF     <- (SimData$UF_T*163.22  + SimData$UF_e*115.41)+400   
  SimData$PRO    <- (SimData$PRO_T*122.41 + SimData$PRO_e*86.56)+300  
  SimData$SF     <- (SimData$SF_T*102.01  + SimData$SF_e*72.13)+275 
  SimData$ALC    <- (SimData$ALC_T*40.80  + SimData$ALC_e*28.85)+175  
  SimData$GLUC   <- SimData$GLUC*25+80
  
  SimUData$NMES   <- (SimUData$NMES_T*102.01 + SimUData$NMES_e*72.13)+250   
  SimUData$CRB    <- (SimUData$CRB_T*204.02  + SimUData$CRB_e*144.26)+500   
  SimUData$FBR    <- (SimUData$FBR_T*40.80   + SimUData$FBR_e*28.85)+100 
  SimUData$UF     <- (SimUData$UF_T*163.22  + SimUData$UF_e*115.41)+400   
  SimUData$PRO    <- (SimUData$PRO_T*122.41 + SimUData$PRO_e*86.56)+300  
  SimUData$SF     <- (SimUData$SF_T*102.01  + SimUData$SF_e*72.13)+275 
  SimUData$ALC    <- (SimUData$ALC_T*40.80  + SimUData$ALC_e*28.85)+175 
  SimUData$GLUC  <- SimUData$GLUC*25+80
  
  # Calculate compositional variables
  
  SimData$TotalEnergy      <- SimData$NMES + SimData$CRB + SimData$FBR + SimData$UF + SimData$PRO +  SimData$SF + SimData$ALC
  SimData$RemainingEnergy   <- SimData$CRB + SimData$FBR + SimData$UF + SimData$PRO + SimData$SF + SimData$ALC
  SimUData$TotalEnergy     <- SimUData$NMES + SimUData$CRB + SimUData$FBR + SimUData$UF + SimUData$PRO + SimUData$SF + SimUData$ALC
  SimUData$RemainingEnergy  <- SimUData$CRB + SimUData$FBR + SimUData$UF + SimUData$PRO + SimUData$SF + SimUData$ALC
  
  # Run each model for energy intake adjustment  
  
### 0. The unadjusted model
  
  mod0   <- lm(GLUC ~ NMES, data=SimData)
  mod0U  <- lm(GLUC ~ NMES, data=SimUData)
  
  mean0  <- mod0$coefficients[2]*100
  mean0U <- mod0U$coefficients[2]*100
  
  
### 1. The energy partition model 
  
  mod1   <- lm(GLUC ~ NMES + RemainingEnergy, data=SimData)
  mod1U  <- lm(GLUC ~ NMES + RemainingEnergy, data=SimUData)
  
  mean1  <- mod1$coefficients[2]*100 
  mean1U <- mod1U$coefficients[2]*100
  
### 2. Тhe standard model
  mod2   <- lm(GLUC ~ NMES + TotalEnergy, data=SimData)
  mod2U  <- lm(GLUC ~ NMES + TotalEnergy, data=SimUData)
  
  mean2  <- mod2$coefficients[2]*100 
  mean2U <- mod2U$coefficients[2]*100

  
### 3. The nutrient density model
  
  # Create nutrient density variable for the exposure (as a percentage)
  SimData$NMESdens  <- SimData$NMES/SimData$TotalEnergy*100
  SimUData$NMESdens <- SimUData$NMES/SimUData$TotalEnergy*100
  
  mod3a  <- lm(GLUC ~ NMESdens, data=SimData);  summary(mod3a);  confint(mod3a)
  mod3aU <- lm(GLUC ~ NMESdens, data=SimUData); summary(mod3aU); confint(mod3aU)
  
  mean3a  <- mod3a$coefficients[2]
  mean3aU <- mod3aU$coefficients[2]
  
  # The multivariable nutrient density model
  
  mod3b  <- lm(GLUC ~ NMESdens + TotalEnergy, data=SimData);  summary(mod3b);  confint(mod3b)
  mod3bU <- lm(GLUC ~ NMESdens + TotalEnergy, data=SimUData); summary(mod3bU); confint(mod3bU)

  mean3b  <- mod3b$coefficients[2]
  mean3bU <- mod3bU$coefficients[2]
  
### 4. The residual model 
  
  Residual  <- lm(NMES ~ TotalEnergy, data=SimData); summary(Residual)
  ResidualU <- lm(NMES ~ TotalEnergy, data=SimUData); summary(ResidualU)
  
  mod4  <- lm(SimData$GLUC ~ Residual$residuals)  
  mod4U <- lm(SimUData$GLUC ~ ResidualU$residuals)
  
  mean4  <- mod4$coefficients[2]*100
  mean4U <- mod4U$coefficients[2]*100
  
### 5. The all-components model
  
  # total causal effect
  mod5a  <- lm(GLUC ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
  mod5aU <- lm(GLUC ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimUData)
  mean5a  <- mod5a$coefficients[2]*100
  mean5aU <- mod5aU$coefficients[2]*100
  
  # average relative causal effect
  mod5b  <- lm(GLUC ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimData)
  mod5bU <- lm(GLUC ~ NMES + CRB + FBR + UF + PRO + SF + ALC, data=SimUData)
  
  #Calculate weights:
  wi <- c(mean(SimData$CRB)/mean(SimData$RemainingEnergy),
          mean(SimData$FBR)/mean(SimData$RemainingEnergy),
          mean(SimData$UF)/mean(SimData$RemainingEnergy),
          mean(SimData$PRO)/mean(SimData$RemainingEnergy),
          mean(SimData$SF)/mean(SimData$RemainingEnergy),
          mean(SimData$ALC)/mean(SimData$RemainingEnergy))
  
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


