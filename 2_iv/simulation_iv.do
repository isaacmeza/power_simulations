********************
version 17.0
********************
 
/*******************************************************************************
* Name of file:	simulation_iv
* Author:	Isaac M
* Machine:	Isaac M 											
* Date of creation:	February. 2, 2022
* Last date of modification: 
* Modifications: 
* Files used:     
		- 
* Files created:  
* Purpose: 


		Power calculations with moving parts.

		We compute power for different sample size, and effect size, and stores them in a
		matrix whose columns are the effect size, and the rows the sample size.
		Each cell in each matrix stores the statistical power.

		A 2SLS (IV) model is considered
			
			X = \alpha_{FS} + \beta_{FS} Z + \epsilon_{FS}				    (1)
			Y = \alpha_{2SLS} + \beta_{2SLS} \hat{X} + \epsilon_{2SLS}		(2)

		where \hat{X} are the predicted values from (1), also known as First Stage.

		As well as a 2SLS (CF) model
			
			P(X = 1 | Z) = \Phi(z\beta_{FS})  										(1)
			Y = \alpha_{CF} + \beta_{CF} X + 
				\delta_{CF}(X\lambda(z\beta_{FS})-(1-X)\lambda(-z\beta_{FS})) + u	(2)
				
		where the first equation denotes a probit model and \lambda(.) is the well-known
		inverse Mills ratio - X\lambda(z\beta_{FS})-(1-X)\lambda(-z\beta_{FS}) is usually
		known as a 'generalized error'. As noted by Wooldrige (2015), bootstraped s.e. 
		are computed. We need to do this because we are running these as separate
		regressions and, in effect, the residual from the first regression is an 
		estimated value. 
				

		The variables denote the following

			Z - exogenous randomization (treatment variable)
			X - take-up variable
			Y - dummy dependent variable/continuous dependent variable
			
		The interpretation is as follows:

		Receiving treatment (Z) (with probability `p_treat') increases the probability 
		of take-up (X) by \beta_{FS} in average from a baseline
		probability of `p_take_up' which in turn increases the probability of 

			1) Being controlled (dummy)
				We suppose that the r.v. follows a Bernoulli distribution with 
				mean `p_baseline' = 0.40 (and a s.d. of `sd' = 0.24)
			2) hba1c (continuous)
				We suppose that the r.v. follows a Log-normal distribution with 
				mean `p_baseline' = 9.74 and a s.d. of `sd' = 2.70 
			3) BMI (continuous)
				We suppose that the r.v. follows a Normal distribution with 
				mean `p_baseline' = 30.76 and a s.d. of `sd' = 6.36
				By the CLT, we can argue that such assumption is reasonable and valid.

		 by \beta_{2SLS} (\beta_{CF})
		in average from a baseline probability of `p_baseline'. Note that in the case
		of the continuous variable the effect is as percentage of the mean, while 
		it is pp for the binay outcome.

		The effect size (of both FS and IV) is drawn from a multi-Normal distribution
		with the first stage s.d. equal to `sd_fs' = 0.30; the iv s.d. equal to
		`sd_iv' = 0.25; and the correlation equal to `cor' = 0.20

		Within each cell `m' replications are run to estimate the power.

*/


clear all
set more off
	
timer on 1

*Sample size
local nn "500(500)5000"
local length_nn=10
*Number of replications
local m=900
*Number of bootstrap replications
local b_reps=500
*Baseline dependent variable (Y)
	*calibrated from data
local p_baseline=9.74
*Standard deviation of dep variable (Y)
	*calibrated from data
local sd=2.70
*Baseline for first stage
local p_take_up=0.20
*SD for first stage (heterogeneity)
local sd_take_up=0.2
*Probability of treatment
local p_treat=0.5
*Significance level
local signif=0.050

*Grid for 2SLS effect (ATT)
/*
local depvar="binary"
local beta_iv "0.05 0.1 0.15 0.2"
local length_iv=4
local sd_iv=0.25
*/

local depvar="continuous"
local beta_iv "0.05 0.10 0.15 0.20 0.25 0.30"
local length_iv=6
local sd_iv=0.25

*Grid for beta_fs (FS)
local beta_fs "0.20 .30 .40 .50"
local length_fs=4
local sd_fs=0.35
*Correlation between fs and 2sls
*If exclution restriction is met, then `cor'=0
local cor=0.20


********************************************************************************
********************************************************************************

local np=`length_nn'*`length_iv'*`length_fs'*`m'
di "--------"
di "Número de procesos:"
di `np'
di "--------"
		
		
local k=0
local pb=0	
qui {	
foreach mu_x in `beta_fs' {
	local k=`k'+1
	*Define matrices in which power is stored 
	matrix pwr_iv_`k'= J(`length_nn',`length_iv',0)
	*matrix pwr_cf_`k'= J(`length_nn',`length_iv',0)
	
	local j=0
	foreach mu_y in `beta_iv' {
		local j=`j'+1

		local i=0
		forvalues n = `nn' {
			local i=`i'+1

			forvalues t=1/`m' {
				local pb=`pb'+1
				
				clear 
				set obs `n'
				
				*Simulate random noise
				gen epsilon_fs=rnormal(0,`sd_take_up')
				gen epsilon_y=rnormal(0,`sd')
				 
				*Simulate a multivariate normal from where FS & ATE effects
				* are going to be drawn
				* If exclution restriction is met, then `cor'=0
				matrix C = (1 , `cor'  \ `cor' , 1 )
				matrix m = (`mu_x',`mu_y')
				matrix sd = (`sd_fs',`sd_iv')
				drawnorm xx yy , n(`n') means(m) sds(sd) corr(C)

				*Exogenous randomization
				gen zz=(uniform()<=`p_treat')

				*Simulate fs
				gen x=.
					*Baseline probability of take-up 
				replace x=(uniform()<=`p_take_up'+epsilon_fs) if zz==0
					*Take.up (encouragement design)
					*Note that for Z = 1 probability of having X = 1 increases by xx
					*(following a marginal normal distribution).					
				replace x=(uniform()<=`p_take_up'+epsilon_fs+xx) if zz==1


				*Simulate dep var
				gen y=.
				
				if "`depvar'"=="binary" {	
						*Baseline dep var probability
					replace y=(uniform()<=`p_baseline'+epsilon_y) if x==0
						*ATT
					replace y=(uniform()<=`p_baseline'+epsilon_y+yy) if x==1
					}
				
				else {	
						*Baseline dep var probability
							*Note that  
							* X\sim log\mathcal{N} (\mu, \sigma) \leftrightarrow
							*	Y=\ln(X)\sim \mathcal{N} (\mu, \sigma)
					replace y=rnormal(`p_baseline'+epsilon_y,`sd') if x==0
						*ATT
					replace y=rnormal((`p_baseline'+epsilon_y)*(1+yy),`sd') if x==1
					}
								
				*2SLS
					*IV
				ivregress 2sls y (x=z), r  
				test _b[x]=0
				matrix pwr_iv_`k'[`i',`j']=pwr_iv_`k'[`i',`j']+(`r(p)'<=`signif')
					*CF
				probit x z, r	
				predict xb, xb
				*Generalized residuals
				gen gen_resid = cond(x == 1, normalden(xb)/normal(xb), -normalden(xb)/(1-normal(xb)))	
				reg y x gen_resid,  vce(bootstrap, reps(`b_reps'))   
				*p-value
				local pval=2*ttail(e(df_r), abs(_b[_x]/_se[_x]))
				matrix pwr_cf_`k'[`i',`j']=pwr_cf_`k'[`i',`j']+(`pval'<=`signif')
				
				
								*Progress bar
				noi {
				if `pb'==1 {
					di "Progress"
					di "--------"
					}
				if `pb'==floor(`np'/10) {
					di "10%"
					}
				if `pb'==floor(`np'*2/10) {
					di "20%"
					}
				if `pb'==floor(`np'*3/10) {
					di "30%"
					}
				if `pb'==floor(`np'*4/10) {
					di "40%"
					}
				if `pb'==floor(`np'*5/10) {
					di "50%"
					}
				if `pb'==floor(`np'*6/10) {
					di "60%"
					}
				if `pb'==floor(`np'*7/10) {
					di "70%"
					}
				if `pb'==floor(`np'*8/10) {
					di "80%"
					}
				if `pb'==floor(`np'*9/10) {
					di "90%"
					}
				if `pb'==floor(`np'-10) {
					di "100%"
					di "--------"
					di "        "
					}
		
				}
				}
			matrix pwr_iv_`k'[`i',`j']=pwr_iv_`k'[`i',`j']/`m'
			matrix pwr_cf_`k'[`i',`j']=pwr_cf_`k'[`i',`j']/`m'
			}
		
	
		}
	clear
	 svmat pwr_iv_`k'
	export delimited using "$directorio\Simulations\c1pwr_iv_`k'.csv", replace novarnames

	clear
	 svmat pwr_cf_`k'
	export delimited using "$directorio\Simulations\c1pwr_cf_`k'.csv", replace novarnames
	
	}
	}

timer off 1
timer list



*POWER CALCULATION GRAPHS
forvalues kk=1/`k' {
	
	import delimited "$directorio\Simulations\c1pwr_cf_`kk'.csv", clear
	*Power calculation graph
	gen sample_size=(_n-1)*500+2500

	twoway (line v1 sample_size, lwidth(medthick)) ///
		(line v2 sample_size, lwidth(medthick)) ///
		(line v3 sample_size, lwidth(medthick)) ///
		(line v4 sample_size, lwidth(medthick)) ///
			, graphregion(color(white)) xtitle("Sample size") ytitle("Statistical power") ///
			 /*legend(order(1 "5%" 2 "10%" 3 "15%" 4 "20%") rows(2))*/ ///
			 legend(order(1 "15%" 2 "20%" 3 "25%" 4 "30%") rows(2))  ///			 
			 xlabel(2500(500)7000)
	graph export "$directorio\Simulations\c1stat_pwr_cf_`kk'.pdf", replace
	
	*---------------------------------------------------------------------------
	
	import delimited "$directorio\Simulations\c1pwr_iv_`kk'.csv", clear
	*Power calculation graph
	gen sample_size=(_n-1)*500+500


	twoway (line v1 sample_size, lwidth(medthick)) ///
		(line v2 sample_size, lwidth(medthick)) ///
		(line v3 sample_size, lwidth(medthick)) ///
		(line v4 sample_size, lwidth(medthick)) ///
		(line v5 sample_size, lwidth(medthick)) ///
		(line v6 sample_size, lwidth(medthick)) ///		
			, graphregion(color(white)) xtitle("Sample size") ytitle("Statistical power") ///
			  /*legend(order(1 "5%" 2 "10%" 3 "15%" 4 "20%") rows(2))*/  ///
			 legend(order(1 "5%" 2 "10%" 3 "15%" 4 "20%" 5 "25%" 6 "30%") rows(2))  ///			 
			 xlabel(500(500)5000)
	graph export "$directorio\Simulations\c1stat_pwr_iv_`kk'.pdf", replace
	}

