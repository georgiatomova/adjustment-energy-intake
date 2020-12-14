# Adjustment for energy intake in nutritional research: a causal inference perspective  

## Table of Contents  
1. Introduction  
2. Data and setup
3. Simulation
4. Sensitivity analyses 

## 1. Introduction

Energy intake is routinely adjusted for in nutritional research. There are four strategies for adjustment that are most commonly used: the standard model, the energy partition model, the nutrient density model, and the residual model.
There has been considerable debate as to which is most appropriate, but there does not seem to be enough clarity in practice. It is also underappreciated that each of these models corresponds to a distinct causal effect estimand, meaning that these models, although 
ostensibly similar in purpose, might actually estimate very different effects.  
  
This study uses simulations to demonstrate the causal effect estimands that are associated with each model, and explores the performance of these models in reducing confounding by common causes of diet. We also explore an alternative model, that we call 'the all-components model', 
and compare its performance to the performance of the approaches that are currently most commonly used in practice.   

## 2. Data and setup  
All data used in this study are simulated, and the code for data simulation is included. Once simulated, all variables are re-scaled so that their mean and standard deviation values are plausible, to aid illustration. Values used for re-scaling have been informed by the [NDNS Data Year 7-8](https://www.gov.uk/government/statistics/ndns-results-from-years-7-and-8-combined), 
a summary of which is available in the NDNS_yr_7_to_8_statistics.xlsx file.  

All analyses were conducted using R version 4.0.3.  

Please have the following packages installed before running the simulations:  
```  
install.packages("devtools")
install_github("jtextor/dagitty/r")
install.packages("progress")
```  

You might also prefer to remove scientific notation:  
```
options(scipen=999)  
```  

## 3. Simulation  

Two datasets are simulated - with and without the presence of confounding by common causes of diet (i.e. determinants of dietary intake), but otherwise identical. The simulated variables are:  
* non-milk extrinsic sugars (NMES) - **the exposure**  
* carbohydrates (CRB)  
* fibre (FBR)  
* unsaturated fat (UF)  
* saturated fat (SF)  
* protein (PRO)  
* alcohol (ALC)  
* body weight (WT) - **the outcome**  

The following variables are calculated, as opposed to simulated:  
* total energy intake - the sum of energy from all nutrients, *including* the exposure  
* residual energy intake - the sum of energy from all nutrients, *excluding* the exposure  

The following models for energy intake adjustment are then run:  
  
  0. The unadjusted model  
  1. The energy partition model  
  2. The standard model  
  3. The nutrient density model, and the multivariable nutrient density model  
  4. The residual model  
  5. The alternative 'all-components' model  

Exposure coefficient estimates from each model are stored in a vector, the point estimate and the corresponding 95% simulation intervals of each model are extracted and saved as the primary study results.  

## 4. Sensitivity analyses  

In the original study simulations, an apparent 'information loss' was observed, i.e. when total energy intake was adjusted for, the estimate was biased even in the absence of any confounding.  

We hypothesised that this is a result of combining variables, that all have distinct effects on the outcome, into one. Therefore, if we simulate all nutrients to have the same effects on the outcome, then we expect the estimates to be unbiased in the absence of confounding.

The following simulations explore scenarios in which either:  
(1) all nutrients, including the exposure, have the same effects on the outcome; or  
(2) all nutrients, excluding the exposure, have the same effects on the outcome.  
