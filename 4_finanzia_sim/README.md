# Finanzia power simulation 

We run power simulations for the main specification in [Expanding Financial Access Via Credit Cards:Evidence from Mexico](www.diegojimenezh.com/assets/pdf/creditcards_main.pdf)


--


For different effect sizes $\beta_t$, and $\gamma_t$, we will be producing synthetic data and estimating the main specification on this simulated data to approximate its power.

$$
Y_{it} = \alpha_{t} + \beta_t 1(MP_i = 10) + \gamma_t(45-r_i)/30 + \epsilon_{it}
$$
	
We will be assuming that the effects of minimum payments and interest rates are separable and that the effect of interest rate changes is linear.


--

## DGP

Data generation process is as follows:	We have 143916 individuals followed during 27 months in the experimental phase. Each of this individuals will have one of nine equally probable possible payment profiles. This profiles will serve as both strata and weights when estimating equation (\ref{main_spec}). Moreover, every individual will be assigned to a minimum payment arm (5\%/10\%) and an interest rate arm (15\%/25\%/35\%/45\%) - each of this drawn from a discrete uniform distribution.


## $Y_{it}$ - Cumulative default


When $Y_{it}$ is defined to be the cumulative default, the model for the DGP can be described as follows:


1. Default ($Y_i$) is drawn from a $\operatorname{Ber}(p_i)$ distribution for every individual, where the propensity to default $p_i$ is described as

$$p_i = \mathcal{N}(\mu_{\text{baseline}},\sigma_{\text{baseline}})+ \mathcal{N}(\beta_{MP},\sigma_{MP})\times 1(MP_i = 10)+ \mathcal{N}(\gamma_{r},\sigma_{r})\times (45-r_i)/30$$

    
2. Conditional on default, the month ($\bar{t}\in[1,27]$) where individual $i$ defaults is simulated from a categorical distribution. We assume that the time of default is independent (or at least has no correlation) of treatment status.

    
3. Finally, $Y_{it}$ is defined as

$$Y_{it} = Y_i\times 1(t\geq\bar{t})$$

i.e., it represents the cumulative default status.


## $Y_{it}$ - Average balance


When $Y_{it}$ is defined to be the average balance, the model for the DGP gets more complicated. We first have to take into account the cumulative default, as well as the cumulative cancellation status, since this is essentially attrition. Moreover, cumulative default and cumulative cancellation are outcomes affected by treatment.


1. Cumulative default ($D_{it}$) is simulated as in the previous section:

$$D_{it} = \operatorname{Ber}(p_i)\times 1(t\geq\bar{t})$$

where 

$$p_i = \mathcal{N}(0.18,0.017)+ \mathcal{N}(0.01,\frac{0.01}{4})\times 1(MP_i = 10)+ \mathcal{N}(0.02,\frac{0.02}{4})\times (45-r_i)/30$$

    
2. Cumulative cancellation ($C_{it}$) is simulated as

$$C_{it} = \operatorname{Ber}(q_i)\times 1(t\geq\tilde{t})\times(1-D_{it})$$

where 

$$q_i = \mathcal{N}(0.016,0.0024)+ \mathcal{N}(0.018,\frac{0.018}{4})\times 1(MP_i = 10)+ \mathcal{N}(0.028,\frac{0.028}{4})\times (45-r_i)/30$$

and $\tilde{t}$ indicates the month of cancellation.
    
3. The average balance $S_{it}$ is modeled as

$$
\begin{align}
    S_{it} &=1(t\neq\tilde{t})\times\left\lbrace\mathcal{N}(\mu_{\text{baseline}},\sigma_{\text{baseline}})+ t\mathcal{N}(\mu_{t},\sigma_{t})+t^2\mathcal{N}(\mu_{t^2},\sigma_{t^2}) \right.\\
    & \quad\quad + \mathcal{N}(250,15)\times 1(MP_i = 10)+ t\mathcal{N}(\beta_{MP},\sigma_{MP})\times 1(MP_i = 10) \\
   & \left.\quad\quad + \mathcal{N}(-100,20)\times \frac{45-r_i}{30}+ t\mathcal{N}(\gamma_{r},\sigma_{r})\times \frac{45-r_i}{30} \right\rbrace \;\quad\text{ if }\;C_{it}=0, D_{it}=0
\end{align}
$$
    
note that the average balance is only defined when there is no attrition ($C_{it}=0, D_{it}=0$), and when individual cancels its debt, then the average is zero. Note also that we are assuming a quadratic trend in time for the average balance.




Because we know the real effect size $\beta_{MP},\gamma_{r}$, we can calculate power for each effect size as the fraction of estimated effects that are statistically different from zero at the 95% confidence level from a total of 250 replications.


    
Finally, we also pool all treatment arms in a single regression and at the individual level

$$Y_{i} = \beta_{mp,r} 1(MP_i = mp)1(r_i = r)$$

and we compute power for the null

$$H_0 : \beta_{5,15} =  \beta_{5,45}$$


--

This folder includes implementation in STATA and also in R. R makes uses of parallelization packages.
