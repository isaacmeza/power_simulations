********************
version 17.0
********************
 
/*******************************************************************************
* Name of file:	
* Author:	Isaac M
* Machine:	Isaac M 											
* Date of creation:	May. 27, 2022
* Last date of modification: 
* Modifications: 
* Files used:     
		- 
* Files created:  
* Purpose: 
	Power calculations for MDE in Finanzia (RESTUD) paper
	Data generation process is as follows:
	We have 143916 individuals followed during 27 months in the experimental phase. Each individual has one of 9 payment profile and tenure with the bank, as well as a $N(p_baseline,sd_baseline)$ probability of default. 

	The experiment has two arms in MP and 4 arms for different interests rates. Each arm has a $N(\cdot, \cdot)$ effect in the probability of default. The simuation is such that the MP and the interest rate have separable and additive effects on the propensity of default. This is the MDE.
		
	We compute power for different effect size, and store them in a
	matrix whose columns are the effect size and the rows the month.
	Each cell in each matrix stores the statistical power.
	The following model is considered
		
		Y_{it} = \alpha_{t} + \beta_t 1(MP_i = 10%) + \gamma_t(45%-r_i)/30% + \epsilon_{it}    	(1)
		
	where i indexes individual, t - month
		
robust standard errors are computed.		
		
********************************************************************************
*/

clear all
timer on 1

*Hyperparameters for simulation
*Number of replications
local m = 250
*Significance level
local signif = 0.0500

*Calibrated parameters
local n = 143916
local p_baseline = .1800389
local sd_baseline = .017701

local p_0 = 0
local p_1 = 0
local p_2 = 0
local p_3 = 0
local p_4 = 0.0009
local p_5 = 0.0238
local p_6 = 0.0429
local p_7 = 0.0713
local p_8 = 0.101
local p_9 = 0.1303
local p_10 = 0.1544
local p_11 = 0.1807
local p_12 = 0.2174
local p_13 = 0.2619
local p_14 = 0.3178
local p_15 = 0.3524
local p_16 = 0.3867
local p_17 = 0.4194
local p_18 = 0.4867
local p_19 = 0.5352
local p_20 = 0.5778
local p_21 = 0.6201
local p_22 = 0.6671
local p_23 = 0.7054
local p_24 = 0.7974
local p_25 = 0.862
local p_26 = 0.9138
local p_27 = 1

*Effect size grid 
local MP_effect = ".0025 .005 .0075 .010 .0125 .015 .0175 .02"
local sd_mp = 0.010/4
local length_MP = 8

local R_effect = "-.005 -.01 -.015 -.02 -.025 -0.03 -0.035 -.04"
local sd_tasa = 0.020/4
local length_R = 8

********************************************************************************

local np=`length_MP'*`m'
di "--------"
di "Número de procesos:"
di `np'
di "--------"		
		
*Define matrices in which power is stored 
matrix pwrMP= J(27,`length_MP',0)
matrix pwrR= J(28,`length_R',0)

matrix pwrMPh= J(27,`length_MP',0)
matrix pwrRh= J(28,`length_R',0)
matrix pwrMPhh= J(27,`length_MP',0) 
matrix pwrRhh= J(28,`length_R',0)

matrix pwrMPl= J(27,`length_MP',0)
matrix pwrRl= J(28,`length_R',0)
matrix pwrMPll= J(27,`length_MP',0)
matrix pwrRll= J(28,`length_R',0)
	
		
	local pb=0	
	local i=0
	foreach mp_effect in `MP_effect' {
		local i = `i' + 1
		
		qui {
		forvalues t = 1/`m' {
			local pb = `pb' + 1
				
			clear 
			*Number of individuals*month
			set obs `n'
			gen id = _n
			
			*Exogenous randomization
			gen MP = (uniform()<=0.5)
			gen unif = uniform()
			gen tasa = 0 if inrange(unif,0,0.25)
			replace tasa = 1/3 if inrange(unif,0.25,0.5)
			replace tasa = 2/3 if inrange(unif,0.5,0.75)
			replace tasa = 1 if inrange(unif,0.75,1)
			drop unif
			
			*Simulate weights
			gen unif = uniform()
			gen probab = 0.6 if inrange(unif,0,1/9)
			replace probab = 0.65 if inrange(unif,1/9,2/9)
			replace probab = 1.56 if inrange(unif,2/9,3/9)
			replace probab = 1.69 if inrange(unif,3/9,4/9)
			replace probab = 3.75 if inrange(unif,4/9,5/9)
			replace probab = 9.75 if inrange(unif,5/9,6/9)
			replace probab = 9.84 if inrange(unif,6/9,7/9)
			replace probab = 10.66 if inrange(unif,7/9,8/9)
			replace probab = 61.5 if inrange(unif,8/9,1)				
			drop unif
			
			egen strata = group(probab)
			
			*Expand
			expand 27
			bysort id : gen month = _n 
			
			*Simulate dep var	
				*Baseline dep var probability (Bernoulli)
			bysort id : gen propensity = rnormal(`p_baseline',`sd_baseline') if _n==1
				*ATT
			bysort id : replace propensity = propensity + rnormal(`mp_effect',`sd_mp') if MP==1 & _n==1
			bysort id : replace propensity = propensity + rnormal(0.02,`sd_tasa')*tasa if tasa!=0 & _n==1
				
			bysort id : gen y_ = (uniform()<propensity) if _n==1
			
			bysort id : gen unif = uniform() if _n==1
			gen mnth_ = .
			forvalues k = 1/27 {
				replace mnth = `k' if inrange(unif,`p_`=`k'-1'',`p_`k'') & y_==1
			}
				
			gen y = 0
			bysort id : replace y = 1 if y_[1]==1 & month>=mnth[1]
			
				
			*REG
			keep y month MP tasa probab strata month
			reghdfe y ibn.month#1.MP ibn.month#c.tasa [pw = probab], absorb(strata#month) nocons
				
			forvalues k = 5/27 {
				test `k'.month#1.MP = 0
				matrix pwrMP[`k',`i']=pwrMP[`k',`i']+(`r(p)'<=`signif')
				matrix pwrMPh[`k',`i'] = pwrMPh[`k',`i'] + (`r(p)'<=`signif'+0.025)
				matrix pwrMPhh[`k',`i'] = pwrMPhh[`k',`i'] + (`r(p)'<=`signif'+0.0125)
				matrix pwrMPl[`k',`i'] = pwrMPl[`k',`i'] + (`r(p)'<=`signif'-0.025)
				matrix pwrMPll[`k',`i'] = pwrMPll[`k',`i'] + (`r(p)'<=`signif'-0.0125)
				}

				
							*Progress bar
			if `pb'==1 {
			noi di " "
			noi _dots 0, title(Loop through replications) reps(`np')
			noi di " "
			}
			noi _dots `pb' 0
	
			}
			
		
		forvalues k = 5/27 {
			*Power
			matrix pwrMP[`k',`i']=pwrMP[`k',`i']/`m'
			*Upper estimate
			matrix pwrMPh[`k',`i']=pwrMPh[`k',`i']/`m'		
			matrix pwrMPhh[`k',`i']=pwrMPhh[`k',`i']/`m'		
			*Lower estimate
			matrix pwrMPl[`k',`i']=pwrMPl[`k',`i']/`m'	
			matrix pwrMPll[`k',`i']=pwrMPll[`k',`i']/`m'	
			}
		}
		}

timer off 1
timer list				 

preserve
clear
svmat pwrMP 
svmat pwrMPh 
svmat pwrMPhh
svmat pwrMPl
svmat pwrMPll

save "C:\Users\isaac\Downloads\pwr_simMP.dta", replace
********************************************************************************

timer on 1
	local pb=0	
	local i=0
	foreach tasa_effect in `R_effect' {
		local i = `i' + 1

		qui {
		forvalues t = 1/`m' {
			local pb = `pb' + 1
				
			clear 
			*Number of individuals*month
			set obs `n'
			gen id = _n
			
			*Exogenous randomization
			gen MP = (uniform()<=0.5)
			gen unif = uniform()
			gen tasa = 0 if inrange(unif,0,0.25)
			replace tasa = 1/3 if inrange(unif,0.25,0.5)
			replace tasa = 2/3 if inrange(unif,0.5,0.75)
			replace tasa = 1 if inrange(unif,0.75,1)
			drop unif
			
			*Simulate weights
			gen unif = uniform()
			gen probab = 0.6 if inrange(unif,0,1/9)
			replace probab = 0.65 if inrange(unif,1/9,2/9)
			replace probab = 1.56 if inrange(unif,2/9,3/9)
			replace probab = 1.69 if inrange(unif,3/9,4/9)
			replace probab = 3.75 if inrange(unif,4/9,5/9)
			replace probab = 9.75 if inrange(unif,5/9,6/9)
			replace probab = 9.84 if inrange(unif,6/9,7/9)
			replace probab = 10.66 if inrange(unif,7/9,8/9)
			replace probab = 61.5 if inrange(unif,8/9,1)				
			drop unif
			
			egen strata = group(probab)
			
			*Expand
			expand 27
			bysort id : gen month = _n 
			
			*Simulate dep var	
				*Baseline dep var probability (Bernoulli)
			bysort id : gen propensity = rnormal(`p_baseline',`sd_baseline') if _n==1
				*ATT
			bysort id : replace propensity = propensity + rnormal(0.01,`sd_mp') if MP==1 & _n==1
			bysort id : replace propensity = propensity + rnormal(`tasa_effect',`sd_tasa')*tasa if tasa!=0 & _n==1
				
			bysort id : gen y_ = (uniform()<propensity) if _n==1
			
			bysort id : gen unif = uniform() if _n==1
			gen mnth_ = .
			forvalues k = 1/27 {
				replace mnth = `k' if inrange(unif,`p_`=`k'-1'',`p_`k'') & y_==1
			}
				
			gen y = 0
			bysort id : replace y = 1 if y_[1]==1 & month>=mnth[1]
			
				
			*REG
			keep y month MP tasa probab strata month
			reghdfe y ibn.month#1.MP ibn.month#c.tasa [pw = probab], absorb(strata#month) nocons
				
			forvalues k = 5/27 {
				test `k'.month#c.tasa = 0
				matrix pwrR[`k',`i']=pwrR[`k',`i']+(`r(p)'<=`signif')
				matrix pwrRh[`k',`i'] = pwrRh[`k',`i'] + (`r(p)'<=`signif'+0.025)
				matrix pwrRhh[`k',`i'] = pwrRhh[`k',`i'] + (`r(p)'<=`signif'+0.0125)
				matrix pwrRl[`k',`i'] = pwrRl[`k',`i'] + (`r(p)'<=`signif'-0.025)
				matrix pwrRll[`k',`i'] = pwrRll[`k',`i'] + (`r(p)'<=`signif'-0.0125)
				}

			reghdfe y i.tasa [pw = probab] if month==27 & MP==0 & inlist(tasa,0,1), absorb(strata)	
			test 1.tasa = 0
				matrix pwrR[28,`i']=pwrR[28,`i']+(`r(p)'<=`signif')
				matrix pwrRh[28,`i'] = pwrRh[28,`i'] + (`r(p)'<=`signif'+0.025)
				matrix pwrRhh[28,`i'] = pwrRhh[28,`i'] + (`r(p)'<=`signif'+0.0125)
				matrix pwrRl[28,`i'] = pwrRl[28,`i'] + (`r(p)'<=`signif'-0.025)
				matrix pwrRll[28,`i'] = pwrRll[28,`i'] + (`r(p)'<=`signif'-0.0125)
							*Progress bar

										*Progress bar
			if `pb'==1 {
			noi di " "
			noi _dots 0, title(Loop through replications) reps(`np')
			noi di " "
			}
			noi _dots `pb' 0
	
			}
		
		forvalues k = 5/27 {
			*Power
			matrix pwrR[`k',`i']=pwrR[`k',`i']/`m'
			*Upper estimate
			matrix pwrRh[`k',`i']=pwrRh[`k',`i']/`m'		
			matrix pwrRhh[`k',`i']=pwrRhh[`k',`i']/`m'		
			*Lower estimate
			matrix pwrRl[`k',`i']=pwrRl[`k',`i']/`m'	
			matrix pwrRll[`k',`i']=pwrRll[`k',`i']/`m'	
			}
		
		}
		}

timer off 1
timer list			
********************************************************************************
********************************************************************************

clear
svmat pwrR 
svmat pwrRh 
svmat pwrRhh
svmat pwrRl
svmat pwrRll

save "C:\Users\isaac\Downloads\pwr_simR.dta", replace
********************************************************************************

use "C:\Users\isaac\Downloads\pwr_simMP.dta", clear

gen month = _n if _n>=4
twoway (rarea pwrMPh2 pwrMPl2 month, color(navy%10)) ///
	(rarea pwrMPhh2 pwrMPll2 month, color(navy%30)) ///
	(line pwrMP2 month, lwidth(medthick) color(navy)) ///
	(rarea pwrMPh4 pwrMPl4 month, color(maroon%10)) ///
	(rarea pwrMPhh4 pwrMPll4 month, color(maroon%30)) ///
	(line pwrMP4 month, lwidth(medthick) color(maroon)) ///
	(rarea pwrMPh6 pwrMPl6 month, color(dkgreen%10)) ///
	(rarea pwrMPhh6 pwrMPll6 month, color(dkgreen%30)) ///
	(line pwrMP6 month, lwidth(medthick) color(dkgreen)) ///
	(rarea pwrMPh8 pwrMPl8 month, color(purple%10)) ///
	(rarea pwrMPhh8 pwrMPll8 month, color(purple%30)) ///
	(line pwrMP8 month, lwidth(medthick) color(purple)) ///
	,graphregion(color(white)) legend(order(3 "{&beta}{subscript:MP} = 0.005" ///
	6 "{&beta}{subscript:MP} = 0.010" ///
	9 "{&beta}{subscript:MP} = 0.015" ///
	12 "{&beta}{subscript:MP} = 0.020") rows(1) size(small)) xtitle("Month since treatment") ///
	ytitle("Statistical power")
graph export "C:\Users\isaac\Downloads\pwr_simMP.pdf", replace

********************************************************************************

use "C:\Users\isaac\Downloads\pwr_simR.dta", clear

gen month = _n if inrange(_n,4,27)
twoway (rarea pwrRh2 pwrRl2 month, color(navy%10)) ///
	(rarea pwrRhh2 pwrRll2 month, color(navy%30)) ///
	(line pwrR2 month, lwidth(medthick) color(navy)) ///
	(rarea pwrRh4 pwrRl4 month, color(maroon%10)) ///
	(rarea pwrRhh4 pwrRll4 month, color(maroon%30)) ///
	(line pwrR4 month, lwidth(medthick) color(maroon)) ///
	(rarea pwrRh6 pwrRl6 month, color(dkgreen%10)) ///
	(rarea pwrRhh6 pwrRll6 month, color(dkgreen%30)) ///
	(line pwrR6 month, lwidth(medthick) color(dkgreen)) ///
	(rarea pwrRh8 pwrRl8 month, color(purple%10)) ///
	(rarea pwrRhh8 pwrRll8 month, color(purple%30)) ///
	(line pwrR8 month, lwidth(medthick) color(purple)) ///
	,graphregion(color(white)) legend(order(3 "{&gamma}{subscript:R} = -.01" ///
	6 "{&gamma}{subscript:R} = -.02" ///
	9 "{&gamma}{subscript:R} = -.03" ///
	12 "{&gamma}{subscript:R} = -.04") rows(1) size(small)) xtitle("Month since treatment") ///
	ytitle("Statistical power")
graph export "C:\Users\isaac\Downloads\pwr_simR.pdf", replace

********************************************************************************

keep if _n==28
gen id = _n
reshape long pwrR pwrRhh pwrRll pwrRh pwrRl,  i(id) j(j)

replace j = -.005*_n 

twoway (rcap pwrRh pwrRl j, color(navy%20)) ///
	(rcap pwrRhh pwrRll j, color(navy%50)) ///
	(scatter pwrR j, msymbol(O) msize(medthick) color(navy)) ///
	(line pwrR j, color(navy%30)) ///
	,graphregion(color(white)) legend(off) xtitle("Effect size in R: {&gamma}{subscript:R}") ///
	ytitle("Statistical power")
graph export "C:\Users\isaac\Downloads\pwr_simRpooled.pdf", replace