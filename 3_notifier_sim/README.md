# Notifiers experiment power simulation

We have $r=1,\ldots, 25$, $\operatorname{Poisson}(\mu_r)$, number of cases per each day-region to be loaded for notification. Each casefile itself has a $\operatorname{Poisson}(\mu_d)$ number of defendants. Hence, the number of 'diligencias' per working day follows a compound Poisson distribution : $\sum^{\sum_{r} \operatorname{Poisson}(\mu_r)} \operatorname{Poisson}(\mu_d)$. 

Casefiles are assigned randomly (with probability $p_{\text{treat}}$) within regions to the treatment arm. 

Baseline probability, $p_{baseline}$, of successful notification follows a $\text{Beta}(a,b)$ distribution, with parameters calibrated such that 

$$\mathbb{E}[p_{\text{baseline}}]=\frac{a}{a+b}$$ 

, and 

$$\text{V}[p_{\text{baseline}}]=\frac{ab}{(a+b)^2(a+b+1)}$$. 

Moreover, each casefile has associated a region r, which has a differential fixed effect of $\bar{\alpha_r}$, and a notifier n which has a differential fixed effect of $\bar{\gamma_n}$. 

Assignment to treatment (ATT) has a (random) treatment effect that is normally distributed $N(\mu_\beta,\sigma^2_{\beta})$. 

In sum, the model for the DGP is 

$$Y_i = 1( U[0,1]_i < \text{Beta}(a,b)_i + \bar{\alpha_r} + \bar{\gamma_n} + N(\mu_\beta,\sigma^2_{\beta})_i 1(i \text{ is Rotator}) )$$

and we estimate it using the following specification

$$Y_{i} = \alpha_{r}  + \beta 1(i \text{ is Rotator}) + \epsilon_{i}$$

clustering standard errors at the region level.
