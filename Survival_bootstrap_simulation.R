######################################################################
######################################################################
# Simulate survival data, and estimate the effect 
# of walking pace -> outcome that is mediated through BMI,
# with standard error derived using bootstrapping.
######################################################################
######################################################################

library(MASS)
library(readr)
library(survival)
library(simsurv)
library(flexsurv)


######################################################################
######################################################################
# Step 1: Simulate walking pace, 
#         and two confounders (ldl cholesterol and smoking).
######################################################################
######################################################################

# assume centred variables.
# define mean and covariance.

mean.vec <- c(0,0,0)
cov.matrix <-  matrix( c( 1, -0.3, -0.2,
			    -0.3, 1, 0.1,
 			    -0.2, 0.1, 1),		
			   nrow = 3, 
			   dimnames = list(c("wp","ldl","smoke"),
					     c("wp","ldl","smoke")))


# set sample size for simulation.

sample_size <- 10000

# simulate walking pace and confounders.

sample_data <- as.data.frame(mvrnorm(n = sample_size,
                               mu = mean.vec, 
                               Sigma = cov.matrix  ))



######################################################################
######################################################################
# Step 2: Simulate bmi as a mediator (walking pace -> bmi -> outcome.).
######################################################################
######################################################################

# bmi is negatively correlated with walking pace.

sample_data$bmi <- -sqrt(0.5)*sample_data$wp + sqrt(0.8)*rnorm(sample_size)

# add id column

sample_data$id <- 1:sample_size
sample_data <- sample_data[,c(5,1:4)]


######################################################################
######################################################################
# Step 3: Simulate survival outcome using simsurv.
######################################################################
######################################################################

# define effect sizes of variables on outcome.

betas = c(wp = 0.4 , bmi = -0.3, ldl = -0.1, smoke = -0.3)

simdat <- simsurv(dist = "weibull", lambdas = 0.05, gammas = 1, betas = betas ,
                  x = sample_data, maxt = 12)

# combine with covariate data.

surv_data <- cbind(sample_data, simdat[,-1])

######################################################################
######################################################################
# Step 4: Define function to estimate the proportion of the total effect
#         of walking pace -> outcome that is mediated through bmi.
######################################################################
######################################################################


mediation.model <- function(surv_data){


# stage 1 identify direct effect of walking pace -> outcome, independent of bmi:

cox.mod.direct <- coxph(Surv(eventtime, status) ~ wp+bmi+ldl+smoke,  data = surv_data) 
direct_effect <- cox.mod.direct$coefficients[1]


# stage 2 identify total effect of walking pace -> outcome,
# and identify proportion mediated.

cox.mod.total <- coxph(Surv(eventtime, status) ~ wp+ldl+smoke,  data = surv_data) 
total_effect <- cox.mod.total$coefficients[1]
proportion_mediated <- (total_effect-direct_effect)/total_effect 


proportion_mediated


}


######################################################################
######################################################################
# Step 5: Evaluate function: estimate proportion mediated.
######################################################################
######################################################################


proportion.mediated <- mediation.model(surv_data)
print(proportion.mediated)

# proportion mediated = 34%


######################################################################
######################################################################
# Step 6: Derive confidence interval for proportion mediated
#         using bootstrap simulation on the log scale.
######################################################################
######################################################################


# create function to resample rows of data.frame

bootstrap_df <- function(surv_data){
		
		boot_data <- surv_data[sample(1:nrow(surv_data), size = nrow(surv_data), replace = TRUE),]
		boot_data$id <- 1:nrow(surv_data)
		boot_data
}



# create data frame to store bootstrap simulation results.
# 100 simulations.


B <- 100

results <- data.frame("iteration" = 1:B,"proportion.mediated" = rep(0,B) )

mediation.model(bootstrap_df(surv_data))

for(i in 1:B){

results$proportion.mediated[i] <- mediation.model(bootstrap_df(surv_data))

}

se.log.prop.mediated <- sd(log(results$proportion.mediated))

prop.mediated.lower.95 <-  exp(log(proportion.mediated)-1.96*se.log.prop.mediated)
prop.mediated.upper.95 <-  exp(log(proportion.mediated)+1.96*se.log.prop.mediated)

# print 95% confidence interval:

print(c(proportion.mediated,prop.mediated.lower.95, prop.mediated.upper.95 ))


# proportion mediated = 34%, 95% CI (30% - 38%)


