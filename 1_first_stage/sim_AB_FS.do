********************
version 17.0
********************
 
/*******************************************************************************
* Name of file:	sim_AB_FS
* Author:	Isaac M
* Machine:	Isaac M 											
* Date of creation:	
* Last date of modification: 
* Modifications: 
* Files used:     
		- 
* Files created:  
* Purpose: 

		Power calculations with a moving part.

		Computes power for different sample size and effect size, and stores them in a
		matrix whose columns are the effect size and the rows the sample size.
		Each cell in each matrix stores the statistical power.

		A first-stage model (also can be interpreted as a reduced-form model) is considered
			
			X = \alpha_{FS} + \beta_{FS} Z + \epsilon_{FS}				(1)

			
		The variables denote the following

			Z - exogenous randomization (treatment variable)
			X - dummy dependent variable
			
		The interpretation is as follows:

		Receiving treatment (Z) (with probability `p_treat') increases the probability 
		of the dependent variable (X) by \beta_{FS} in average from a baseline 
		probability of `p_baseline'.

		Within each cell `m' replications are run to estimate the power.

*/


clear all
timer on 1

*Sample size
local nn "200(50)1000"
local length_nn=17
*Number of replications
local m=1000
*Baseline dependent variable
local p_baseline=0.64
*Probability of treatment
local p_treat=0.5
*Probability of attrition
local p_attrition=0.32
local diff_attrition=-0.0
*Significance level
local signif=0.051
*Grid for beta_fs (FS)
local beta_fs "0(-0.02)-0.26"
local length_fs=14
local sd_fs=0.074

********************************************************************************
********************************************************************************


local np=`length_nn'*`length_fs'*`m'
di "--------"
di "Número de procesos:"
di `np'
di "--------"
		
		
local k=0
local pb=0


	*Define matrices in which power is stored 
	matrix pwr= J(`length_nn',`length_fs',0)
	
	local j=0
	forvalues mu_x = `beta_fs' {
		local j=`j'+1

		local i=0
		forvalues n = `nn' {
			local i=`i'+1

			forvalues t=1/`m' {
				local pb=`pb'+1
				
				clear 
				qui set obs `n'
				
				*Simulate a normal dist from where FS is going to be drawn
				qui gen yy=rnormal(`mu_x',`sd_fs')
				
				*Exogenous randomization
				qui gen zz=(uniform()<=`p_treat')

				*Simulate dep var
				qui gen y=.
					*Baseline dep var probability
				qui replace y=(uniform()<=`p_baseline') if zz==0
					*ATT
				qui replace y=(uniform()<=`p_baseline'+yy) if zz==1

				*Drop observations with attrition
				qui replace y=. if  (uniform()<=`p_attrition') & zz==0
					*Differential attrition
				qui replace y=. if  (uniform()<=`p_attrition'+`diff_attrition') & zz==1
				
				*IVREG
				qui reg y zz, r  
				qui test _b[zz]=0
				matrix pwr[`i',`j']=pwr[`i',`j']+(`r(p)'<=`signif')
				
				
								*Progress bar

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
				
			matrix pwr[`i',`j']=pwr[`i',`j']/`m'
			}
		}
	
	
	clear
	qui svmat pwr
	
	


timer off 1
timer list

