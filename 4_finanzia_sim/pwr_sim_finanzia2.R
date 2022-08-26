if (!require("pacman")) install.packages("pacman")
pacman::p_load(purrr, dplyr, discreteRV, estimatr, doParallel, foreach, tidyverse)

# SET WORKING DIRECTORY
 setwd('C:/Users/isaac/Downloads')

# Hyperparameters for simulation
# Number of replications
m <<- 250
# Significance level
significancia <<- 0.0500

# Calibrated parameters
n_obs <<-143916 
  #
mu_baseline <<- 983.8833
sd_baseline <<- 19.97
trend1 <<- 37.76
sd_trend1 <<- 3.46
trend2 <<- -.325
sd_trend2 <<- .124


# Effect size grid 
MP_effects <<- c(-10, -20, -30, -40)
sd_mp <<- 50

R_effects <<- c(-10, -20, -30, -40)
sd_tasa <<- 50


################################################################################

mnthRV <<- RV(outcomes = c(1:27), probs = c(0, 0,	0,	0.0009,	0.0229,	0.0191,	0.0284,	0.0297,	0.0293,
                                            0.0241,	0.0263,	0.0367,	0.0445,	0.0559,	0.0346,	0.0343,	0.0327,
                                            0.0673,	0.0485,	0.0426,	0.0423,	0.047,	0.0383,	0.092,	0.0646,
                                            0.0518,	0.0862))
mnthRV_c <<- RV(outcomes = c(1:27), probs = c(0.0887,	0.0397,	0.0414,	0.0387,	0.0486,	0.0407,	0.0431,	0.0421,	
                                              0.0315,	0.0342,	0.0414,	0.0281,	0.0421,	0.0339,	0.0291,	0.0305,
                                              0.0397,	0.0377,	0.0397,	0.0404,	0.0233,	0.0346,	0.0332,	0.0284,
                                              0.0315,	0.0164,	0.0212))


pwrMP <- matrix(,nrow=27,ncol=length(MP_effects))
pwrMPh <- matrix(,nrow=27,ncol=length(MP_effects))
pwrMPhh <- matrix(,nrow=27,ncol=length(MP_effects))
pwrMPl <- matrix(,nrow=27,ncol=length(MP_effects))
pwrMPll <- matrix(,nrow=27,ncol=length(MP_effects))

pwrR <- matrix(,nrow=27,ncol=length(R_effects))
pwrRh <- matrix(,nrow=27,ncol=length(R_effects))
pwrRhh <- matrix(,nrow=27,ncol=length(R_effects))
pwrRl <- matrix(,nrow=27,ncol=length(R_effects))
pwrRll <- matrix(,nrow=27,ncol=length(R_effects))


pwr_mppvals <- function(mp_effect) {
  
  id <- c(1:n_obs)
  #Exogenous randomization
  MP <- rbinom(n_obs,1,0.5)
  tasa <- rdunif(n_obs,3,0)/3
  
  #Simulate weights
  probab <- sample(c(0.6,0.65,1.56,1.69,3.75,9.75,9.84,10.66,61.5),n_obs,replace=TRUE)
  
  #Baseline dep var (default) probability (Bernoulli)
  y_ <- rnorm(n_obs,.1800389,.017701) + rnorm(n_obs,0.010,0.010/4)*MP +
    rnorm(n_obs,0.02, 0.020/4)*tasa
  y_ <- as.integer((runif(n_obs)<y_))
  #Time of default if default
  mnth <- as.integer(rsim(mnthRV,n_obs))*y_
  #Baseline dep var (cancel) probability (Bernoulli)
  c_ <- rnorm(n_obs,.1624131,.00245255) + rnorm(n_obs,0.018,0.018/4)*MP +
    rnorm(n_obs,0.028,0.028/4)*tasa
  c_ <- as.integer((runif(n_obs)<c_))
  #Time of cancel if cancel & not default
  mnth_c <- as.integer(rsim(mnthRV_c,n_obs))*c_*(1-y_)
  #Baseline dep var (sdo_prom) (Normal)
  s_ <- rnorm(n_obs,mu_baseline, sd_baseline)
  
  #Expand
  data <- data.frame(id, MP, tasa, probab, mnth, mnth_c,s_) %>% slice(rep(1:n(), each = 27)) %>% 
    group_by(id) %>% mutate(month = row_number())
  remove(id, MP, tasa, probab, mnth, mnth_c, y_, c_, s_)
  
  #Simulate dep var	
  data <- data %>% mutate(s = s_ + rnorm(n(),trend1, sd_trend1)*month 
                           + rnorm(n(),trend2, sd_trend2)*month^2
                           + rnorm(n(),250,15)*MP 
                           + rnorm(n(),mp_effect,sd_mp)*MP*month 
                           + rnorm(n(),-100, 20)*tasa
                           + rnorm(n(),-10, sd_tasa)*tasa*month) %>% 
    mutate(s = replace(s, (month >= mnth & mnth!=0) |  (month >= mnth_c & mnth_c!=0), NA)) %>%
    mutate(s = replace(s, month < 3 + mnth_c & month >= mnth_c & mnth_c!=0, 0)) 
  
  
  #REG
  model <- lm_robust(s ~ as.factor(month):MP + as.factor(month):tasa, data=data, 
                     weights = probab, fixed_effects = ~ as.factor(probab):month,
                     se_type = "HC1", try_cholesky = TRUE)
  pvalues <- model$p.value[1:27]
  remove(data, model)
  gc()
  return(as.numeric(pvalues))
}
pwr_rpvals <- function(tasa_effect) {
  
  id <- c(1:n_obs)
  #Exogenous randomization
  MP <- rbinom(n_obs,1,0.5)
  tasa <- rdunif(n_obs,3,0)/3
  
  #Simulate weights
  probab <- sample(c(0.6,0.65,1.56,1.69,3.75,9.75,9.84,10.66,61.5),n_obs,replace=TRUE)
  
  #Baseline dep var (default) probability (Bernoulli)
  y_ <- rnorm(n_obs,.1800389,.017701) + rnorm(n_obs,0.010,0.010/4)*MP +
    rnorm(n_obs,0.02, 0.020/4)*tasa
  y_ <- as.integer((runif(n_obs)<y_))
  #Time of default if default
  mnth <- as.integer(rsim(mnthRV,n_obs))*y_
  #Baseline dep var (cancel) probability (Bernoulli)
  c_ <- rnorm(n_obs,.1624131,.00245255) + rnorm(n_obs,0.018,0.018/4)*MP +
    rnorm(n_obs,0.028,0.028/4)*tasa
  c_ <- as.integer((runif(n_obs)<c_))
  #Time of cancel if cancel & not default
  mnth_c <- as.integer(rsim(mnthRV_c,n_obs))*c_*(1-y_)
  #Baseline dep var (sdo_prom) (Normal)
  s_ <- rnorm(n_obs,mu_baseline, sd_baseline)
  
  #Expand
  data <- data.frame(id, MP, tasa, probab, mnth, mnth_c,s_) %>% slice(rep(1:n(), each = 27)) %>% 
    group_by(id) %>% mutate(month = row_number())
  remove(id, MP, tasa, probab, mnth, mnth_c, y_, c_, s_)
  
  #Simulate dep var	
  data <- data %>% mutate(s = s_ + rnorm(n(),trend1, sd_trend1)*month 
                          + rnorm(n(),trend2, sd_trend2)*month^2
                          + rnorm(n(),250,15)*MP 
                          + rnorm(n(),-30, sd_mp)*MP*month 
                          + rnorm(n(),-100, 20)*tasa
                          + rnorm(n(),tasa_effect, sd_tasa)*tasa*month) %>% 
    mutate(s = replace(s, (month >= mnth & mnth!=0) |  (month >= mnth_c & mnth_c!=0), NA)) %>%
    mutate(s = replace(s, month < 3 + mnth_c & month >= mnth_c & mnth_c!=0, 0)) 
  
  
  #REG
  model <- lm_robust(s ~ as.factor(month):MP + as.factor(month):tasa, data=data, 
                     weights = probab, fixed_effects = ~ as.factor(probab):month,
                     se_type = "HC1", try_cholesky = TRUE)
  
  pvalues <- model$p.value[28:54]
  remove(data, model)
  gc()
  return(as.numeric(pvalues))
}

# Use the detectCores() function to find the number of cores in system
no_cores <- detectCores() - 1
registerDoParallel(makeCluster(no_cores))
#registerDoParallel(makeCluster(no_cores, type = "FORK"))

i = 0
for (mp_effect in MP_effects) {
  i = i + 1
  mppvals <- foreach(t = 1:m, .combine =rbind, 
                   .packages=c("dplyr", "estimatr", "purrr", "discreteRV"))  %dopar%  {
                     pwr_mppvals(mp_effect)
                     }
  gc()
  pwrMP[1:27,i] <- colMeans((mppvals<=significancia))  
  pwrMPh[1:27,i] <- colMeans((mppvals<=significancia+0.025))  
  pwrMPhh[1:27,i] <- colMeans((mppvals<=significancia+0.0125))  
  pwrMPl[1:27,i] <- colMeans((mppvals<=significancia-0.025))  
  pwrMPll[1:27,i] <- colMeans((mppvals<=significancia-0.0125))  
}

df_pwrMP <- data.frame(pwrMP, pwrMPh, pwrMPhh, pwrMPl, pwrMPll)
write_csv(df_pwrMP, "./pwr_simMPsdo.csv")

i = 0
for (tasa_effect in R_effects) {
  i = i + 1
  rpvals <- foreach(t = 1:m, .combine =rbind, 
                     .packages=c("dplyr", "estimatr", "purrr", "discreteRV"))  %dopar%  {
                       pwr_rpvals(tasa_effect)
                     }
  gc()
  pwrR[1:27,i] <- colMeans((rpvals<=significancia))  
  pwrRh[1:27,i] <- colMeans((rpvals<=significancia+0.025))  
  pwrRhh[1:27,i] <- colMeans((rpvals<=significancia+0.0125))  
  pwrRl[1:27,i] <- colMeans((rpvals<=significancia-0.025))  
  pwrRll[1:27,i] <- colMeans((rpvals<=significancia-0.0125))  
}

df_pwrR <- data.frame(pwrR, pwrRh, pwrRhh, pwrRl, pwrRll)
write_csv(df_pwrR, "./pwr_simRsdo.csv")

stopImplicitCluster()
gc()


		
  