# Instrumental Variables Estimation

We compute power for different sample size, and effect size, and stores them in a matrix whose columns are the effect size, and the rows the sample size. Each cell in each matrix stores the statistical power.

A 2SLS (IV) model is considered
			
$$X = \alpha_{FS} + \beta_{FS} Z + \epsilon_{FS} \ldots\ldots\ldots (1)$$

$$Y = \alpha_{2SLS} + \beta_{2SLS} \hat{X} + \epsilon_{2SLS} \ldots\ldots\ldots (2)$$

where $\hat{X}$ are the predicted values from (1), also known as First Stage.

As well as a 2SLS (CF) model
			
$$P(X = 1 | Z) = \Phi(z\beta_{FS}) \ldots\ldots\ldots (1)$$

$$Y = \alpha_{CF} + \beta_{CF} X + \delta_{CF}(X\lambda(z\beta_{FS})-(1-X)\lambda(-z\beta_{FS})) + u \ldots\ldots\ldots	(2)$$
				
where the first equation denotes a probit model and \lambda(.) is the well-known inverse Mills ratio - X\lambda(z\beta_{FS})-(1-X)\lambda(-z\beta_{FS}) is usually known as a 'generalized error'. As noted by Wooldrige (2015), bootstraped s.e. are computed. We need to do this because we are running these as separate	regressions and, in effect, the residual from the first regression is an estimated value. 
				
The variables denote the following
			Z - exogenous randomization (treatment variable)
			X - take-up variable
			Y - dummy dependent variable/continuous dependent variable
			
The interpretation is as follows:

Receiving treatment (Z) (with probability 'p_treat') increases the probability of take-up (X) by $\beta_{FS}$ in average from a baseline probability of 'p_take_up' which in turn increases the probability of 

  1. Being controlled (dummy)
      We suppose that the r.v. follows a Bernoulli distribution with mean 'p_baseline' = 0.40 (and a s.d. of 'sd' = 0.24)

  2. hba1c (continuous)
      We suppose that the r.v. follows a Log-normal distribution with mean 'p_baseline' = 9.74 and a s.d. of 'sd' = 2.70 

  3. BMI (continuous)
      We suppose that the r.v. follows a Normal distribution with mean 'p_baseline' = 30.76 and a s.d. of 'sd' = 6.36. By the CLT, we can argue that such assumption is reasonable and valid. 
      
by $\beta_{2SLS}$ ($\beta_{CF}$) in average from a baseline probability of 'p_baseline'. Note that in the case of the continuous variable the effect is as percentage of the mean, while it is pp for the binay outcome.
	
The effect size (of both FS and IV) is drawn from a multi-Normal distribution	with the first stage s.d. equal to 'sd_fs' = 0.30; the iv s.d. equal to 'sd_iv' = 0.25; and the correlation equal to 'cor' = 0.20.

Within each cell 'm' replications are run to estimate the power.
