********************
version 17.0
********************
 
/*******************************************************************************
* Name of file:	
* Author:	Isaac M
* Machine:	Isaac M 											
* Date of creation:	February. 2, 2022
* Last date of modification: 
* Modifications: 
* Files used:     
		- 
* Files created:  

* Purpose: 

	Power calculations for success in notification for the 'Notifiers experiment'.

	Data generation process is as follows: 
	We have r Poisson(`mu_r') number of cases per each day-region to be loaded for notification. Each casefile itself has a Poisson(`mu_d') number of defendants. Hence, the number of 'diligencias' follows a compound Poisson distribution per day (\sum^{\sum_{r} Poisson(`mu_r')} Poisson(`mu_d')). 

	Casefiles are assigned randomly (with probability `p_treat') within regions to the treatment arm.

	Baseline probability ( `p_baseline'=a/(a+b) , `var_p_baseline'=ab/({a+b}^2(a+b+1)) ) of succesful notification follows a Beta(a,b) distribution.

	Each casefile has associated a region r, which has a differential fixed effect of \alpha_r , and a notifier n which has a differential fixed effect of \gamma_n.

	Assignment to treatment (ATT) has a (random) treatment effect of N(\beta,\sigma). 

	In sum, the model is 
		Y_i = 1( runiform_i < `p_baseline'_i + alfa_r + gamma_n + N(\beta,\sigma)*1(rotador) )
		


	We compute power for different sample size and effect size, and store them in a
	matrix whose columns are the effect size and the rows the sample size.
	Each cell in each matrix stores the statistical power.

	The following model is considered
		
		Y_{i} = \alpha_{r} + \beta 1(Rotador)_i + \epsilon_{i}    	(1)
		
	where i indexes casefile, r - region, and n - rotating notifier
		\alpha_r     - Region FE
		1(Rotador)_i - Binary variable indicating when casefile i is rotator vs fixed
		
and we cluster standard errors at the region level.		
		
	The interpretation is as follows:
	Casefile i in time t, assigned to Rotating arm (with probability `p_treat') increases
	the probability of succesful notification Y by \beta in average.

	Within each cell `m' replications are run to estimate the power.
********************************************************************************
*/


	
clear all
timer on 1

*Hyperparameters for simulation
	*Number of replications
	local m = 1500
	*Significance level
	local signif = 0.0501

*Calibrated parameters
	
	*Mean # of casefiles per regions
	local mu_1 = 7
	local mu_2 = 2
	local mu_3 = 5
	local mu_4 = 7
	local mu_5 = 7
	local mu_6 = 12
	local mu_7 = 3
	local mu_8 = 4
	local mu_9 = 5
	local mu_10 = 5
	local mu_11 = 7
	local mu_12 = 9
	local mu_13 = 3
	local mu_14 = 5
	local mu_15 = 6
	local mu_16 = 6
	local mu_17 = 5
	local mu_18 = 9
	local mu_19 = 8
	local mu_20 = 6
	local mu_21 = 6
	local mu_22 = 7
	local mu_23 = 4
	local mu_24 = 4
	local mu_25 = 3
	
	*Region fixed effects
	local alfa_1 = 0.00
	local alfa_2 = 0.0373512
	local alfa_3 = 0.0905405
	local alfa_4 = 0.0491221
	local alfa_5 = 0.006141
	local alfa_6 = 0.0262904
	local alfa_7 = 0.119692
	local alfa_8 = 0.0564404
	local alfa_9 = 0.0731232
	local alfa_10 = 0.096254
	local alfa_11 = -0.168103
	local alfa_12 = -0.2185041
	local alfa_13 = -0.1551029
	local alfa_14 = -0.131526
	local alfa_15 = -0.1454873
	local alfa_16 = -0.1528162
	local alfa_17 = -0.2019781
	local alfa_18 = -0.1557296
	local alfa_19 = -0.155013
	local alfa_20 = -0.1934067
	local alfa_21 = -0.1919277
	local alfa_22 = -0.1320546
	local alfa_23 = -0.1197123
	local alfa_24 = -0.0816916
	local alfa_25 = -0.1858485

	*Mean # of defendants per casefile
	local mu_d_1 = 4
	local mu_d_2 = 6
	local mu_d_3 = 7
	local mu_d_4 = 4
	local mu_d_5 = 4
	local mu_d_6 = 4
	local mu_d_7 = 4
	local mu_d_8 = 3
	local mu_d_9 = 7
	local mu_d_10 = 3
	local mu_d_11 = 4
	local mu_d_12 = 13
	local mu_d_13 = 6
	local mu_d_14 = 4
	local mu_d_15 = 3
	local mu_d_16 = 3
	local mu_d_17 = 4
	local mu_d_18 = 4
	local mu_d_19 = 3
	local mu_d_20 = 4
	local mu_d_21 = 4
	local mu_d_22 = 4
	local mu_d_23 = 3
	local mu_d_24 = 7
	local mu_d_25 = 10

	*Notifier fixed effects
	forvalues r = 1/50 {
		local gama_`r' = rnormal(.0220244, 0.1449842^2)
	}
	*Probability of treatment
	local p_treat=0.5
	*Baseline probability of succesful notification
	local p_baseline = 0.2753778
	local var_p_baseline = (0.075)^2
	local a = (`p_baseline'^2)*((1-`p_baseline')/`var_p_baseline'-1/`p_baseline')
	local b = `a'*(1/`p_baseline'-1)
				
*Moving arms
	*Number of days until rotations
	local dr = 5
	*Sample size grid (number of working days - increased by number of days until rotation)
	local numdays "5(`dr')100"
	local length_nn = 0
	forvalues nn = `numdays' {
		local length_nn = `length_nn'+1
	}

	*Effect size grid (as % of the mean) - \beta = beta_eff*`p_baseline'
	local beta_eff  = ".025 .05 .075 .10"
	local length_eff = 4
		*We compute the variance of beta_eff such that the effect has +-5% of error
	local var_eff = (0.05/1.96)^2


********************************************************************************
********************************************************************************
		
		
local pb=0
				
*Define matrices in which power is stored 
matrix pwr = J(`length_nn',`length_eff',0)
matrix pwrh = J(`length_nn',`length_eff',0)
matrix pwrhh = J(`length_nn',`length_eff',0)
matrix pwrl = J(`length_nn',`length_eff',0)
matrix pwrll = J(`length_nn',`length_eff',0)
matrix obs = J(`length_nn',`length_eff',0)

qui {
	
local j=0
	foreach beta in `beta_eff' {
		local j=`j'+1
		di "`beta'"
	
		local i=0
		forvalues nn = `numdays' {
			local i=`i'+1
			
			forvalues t = 1/`m' {
				local pb = `pb'+1

				clear
				*Number of rotations in the `nn' days
				local nr = round(`nn'/`dr')
				set obs 25

				gen num_c = .
				forvalues r = 1/25 {
					*Simulate number of casefiles per day-region (Poisson distributed)
					replace num_c = rpoisson(`mu_`r''*`nn') + 1 in `r'
				}
				
				*Identification of notifier per rotation period
				gen not_1 = _n 
				forvalues d = 2/`nr' {
					cap drop u 
					gen u = runiform()
					sort u
					gen not_`d' = _n
				}
				sort not_1
				
				*Total number of observations (Sum of compound Poisson distributions)				
				gen tot_num_obs = sum(num_c)
				local num_obs = tot_num_obs[_N]
				set obs `num_obs'
				matrix obs[`i',`j'] = obs[`i',`j'] + `num_obs'		
				
				gen region = 1 if _n<=tot_num_obs[1]
				gen alfa = `alfa_1' if _n<=tot_num_obs[1]
				forvalues r = 2/25 {
					*Identification of region FE
					replace alfa = `alfa_`r'' if inrange(_n,tot_num_obs[`r'-1]+1,tot_num_obs[`r'])
					*Identification of region
					replace region = `r' if inrange(_n,tot_num_obs[`r'-1]+1,tot_num_obs[`r'])		
				}
				
				
				*Identification of notifier				
				bysort region : gen rotator_period = 1 if (_n<=_N/`nr')
				forvalues d = 2/`nr' {
					bysort region : replace rotator_period = `d' if inrange(_n,(_N/`nr')*(`d'-1),(_N/`nr')*(`d'))
				}
				
				gen not_r = .
				forvalues d = 1/`nr' {
					forvalues r = 1/25 {
						replace not_r = not_`d'[`r'] if region==`r' & rotator_period==`d'
					}
				}
				

				*Identification of casefiles				
				gen casefile = _n 
				
				*Exogenous randomization (clustered)
				gen rotador = .
				bysort region : replace rotador = (runiform()<=`p_treat')
				
				*Identification of notifier FE
					*Fixed notifiers
				replace not_r = 25 + region if rotador==0
				gen gama = .
				forvalues r = 1/50 {
					replace gama = `gama_`r'' if not_r==`r'
				}

				
				*Expansion of dataset to include different defendants
				gen num_d = .
				forvalues r = 1/25 {
					replace num_d = rpoisson(`mu_d_`r'') if region==`r'
				}
				*Truncate
				replace num_d = 15 if num_d>=15 & !missing(num_d)
				expand num_d

				*Simulate a normal dist from where treatment effect is going to be drawn
				local bta = `p_baseline'*`beta'
				local sd_bta = sqrt(`var_p_baseline'*`var_eff'+`var_p_baseline'*`beta'^2+`var_eff'*`p_baseline'^2)
				gen beta = rnormal(`bta',`sd_bta')
			
				*Simulate dep var (success)
				gen y = .
				gen pbaseline = rbeta(`a', `b')
				gen unif = runiform()

					*Control mean
				replace y = (unif<=pbaseline + alfa + gama) if rotador==0
					*ATT
				replace y = (unif<=pbaseline + alfa + gama + beta) if rotador==1

				*REG
				reghdfe y rotador, absorb(region) cluster(casefile)
				local sign_rot = sign(_b[rotador])
				*p-value one-sided test - H_0 : \beta<=0
				local pvalue = ttail(`e(df_r)',`sign_rot'*sqrt(`e(F)'))
				di `pvalue'
				matrix pwr[`i',`j'] = pwr[`i',`j'] + (`pvalue'<=`signif')
				matrix pwrh[`i',`j'] = pwrh[`i',`j'] + (`pvalue'<=`signif'+0.025)
				matrix pwrhh[`i',`j'] = pwrhh[`i',`j'] + (`pvalue'<=`signif'+0.0125)
				matrix pwrl[`i',`j'] = pwrl[`i',`j'] + (`pvalue'<=`signif'-0.025)
				matrix pwrll[`i',`j'] = pwrll[`i',`j'] + (`pvalue'<=`signif'-0.0125)
				
					*Progress bar

				if `pb'==1 {
					local np=`length_nn'*`length_eff'*`m'
					noi di "--------"
					noi di "Número de procesos:" `np'
					noi di "--------"
					noi di "Progress"
					noi di "--------"
				}
				if `pb'==floor(`np'/10) {
					noi di "10%"
				}
				if `pb'==floor(`np'*2/10) {
					noi di "20%"
				}
				if `pb'==floor(`np'*3/10) {
					noi di "30%"
				}
				if `pb'==floor(`np'*4/10) {
					noi di "40%"
				}
				if `pb'==floor(`np'*5/10) {
					noi di "50%"
				}
				if `pb'==floor(`np'*6/10) {
					noi di "60%"
				}
				if `pb'==floor(`np'*7/10) {
					noi di "70%"
				}
				if `pb'==floor(`np'*8/10) {
					noi di "80%"
				}
				if `pb'==floor(`np'*9/10) {
					noi di "90%"
				}
				if `pb'==floor(`np'-1) {
					noi di "100%"
					noi di "--------"
					noi di "        "
				}
			
			}
			
			*Power
			matrix pwr[`i',`j']=pwr[`i',`j']/`m'
			*Upper estimate
			matrix pwrh[`i',`j']=pwrh[`i',`j']/`m'
			matrix pwrhh[`i',`j']=pwrhh[`i',`j']/`m'
			*Lower estimate
			matrix pwrl[`i',`j']=pwrl[`i',`j']/`m'		
			matrix pwrll[`i',`j']=pwrll[`i',`j']/`m'
			*Obs casefiles
			matrix obs[`i',`j'] = obs[`i',`j']/`m'
		}
	}
}				 
				 
				 
timer off 1
timer list				 

********************************************************************************
********************************************************************************

svmat pwr 
svmat pwrh 
svmat pwrhh
svmat pwrl
svmat pwrll
svmat obs
 
gen wk_days = _n*5 if !missing(pwr1)
egen obs = rowmean(obs*)
replace obs = round(obs)
save "$directorio/DB/pwr_sim_2.dta", replace



* Power calculation plot (# of working days)
twoway (rarea pwrh1 pwrl1 wk_days, color(navy%15)) (rarea pwrhh1 pwrll1 wk_days, color(navy%25)) (line pwr1 wk_days, xlabel(5(10)100) color(navy) lwidth(thick)) ///
	(rarea pwrh2 pwrl2 wk_days, color(maroon%15)) (rarea pwrhh2 pwrll2 wk_days, color(maroon%25)) (line pwr2 wk_days, xlabel(5(10)100) color(maroon) lwidth(thick)) ///
	(rarea pwrh3 pwrl3 wk_days, color(dkgreen%15)) (rarea pwrhh3 pwrll3 wk_days, color(dkgreen%25)) (line pwr3 wk_days, xlabel(5(10)100) color(dkgreen) lwidth(thick)) ///
	(rarea pwrh4 pwrl4 wk_days, color(gs5%15)) (rarea pwrhh4 pwrll4 wk_days, color(gs5%25)) (line pwr4 wk_days, xlabel(5(10)100) color(gs5) lwidth(medthick)) ///
	, xtitle("# of working days") ytitle("Statistical power") graphregion(color(white)) legend(order (3 "2.5%" 6 "5%" 9 "7.5%" 12 "10%") rows(1) title("Effect", size(small) color(black)))
graph export "$directorio/Figuras/pwr_sim_graph_1.pdf", replace
	
	
* Power calculation plot (# of casefiles)
twoway (rarea pwrh1 pwrl1 obs, color(navy%15)) (rarea pwrhh1 pwrll1 obs, color(navy%25)) (line pwr1 obs, xlabel(1000(2000)15000) color(navy) lwidth(thick)) ///
	(rarea pwrh2 pwrl2 obs, color(maroon%15)) (rarea pwrhh2 pwrll2 obs, color(maroon%25)) (line pwr2 obs, xlabel(1000(2000)15000) color(maroon) lwidth(thick)) ///
	(rarea pwrh3 pwrl3 obs, color(dkgreen%15)) (rarea pwrhh3 pwrll3 obs, color(dkgreen%25)) (line pwr3 obs, xlabel(1000(2000)15000) color(dkgreen) lwidth(thick)) ///
	(rarea pwrh4 pwrl4 obs, color(gs5%15)) (rarea pwrhh4 pwrll4 obs, color(gs5%25)) (line pwr4 obs, xlabel(1000(2000)15000) color(gs5) lwidth(medthick)) ///
	, xtitle("# of casefiles") ytitle("Statistical power") graphregion(color(white)) legend(order (3 "2.5%" 6 "5%" 9 "7.5%" 12 "10%") rows(1) title("Effect", size(small) color(black)))
graph export "$directorio/Figuras/pwr_sim_graph_2.pdf", replace
	
