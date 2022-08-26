# Linear power simulation


Computes power for different sample size and effect size, and stores them in a matrix whose columns are the effect size and the rows the sample size. Each cell in each matrix stores the statistical power. A first-stage model (also can be interpreted as a reduced-form model) is considered
			
$$X = \alpha_{FS} + \beta_{FS} Z + \epsilon_{FS}$$
			
The variables denote the following
			Z - exogenous randomization (treatment variable)
			X - dummy dependent variable
			
The interpretation is as follows:

Receiving treatment (Z) (with probability 'p_treat') increases the probability of the dependent variable (X) by $\beta_{FS}$ in average from a baseline probability of 'p_baseline'. Within each cell 'm' replications are run to estimate the power.
