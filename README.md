# Power simulation 

Statistical power $(1-\beta)$ refers to the likelihood/probability of detecting an effect where a non-zero effect is present, i,e. given a null $H_0 : \theta = 0$, and its alternative $H_1 : \theta \neq 0$, the power is defined as

$$1-\beta = \Pr(\text{reject } H_0 | H_1  \text{ is true})$$

To estimate the power of our main specification, we use Monte Carlo simulation by replicating the data generating process of our experiment:

$$(X,Y)\sim F$$

where $F$ is a very general distribution and can be as complex as desired.

Given this synthetic data, together with an identification method we will estimate the parameters of interest of the model and compute its power.
For instance, we might use a linear conditional mean model:

$$\mathbb{E}[Y|X]=\theta^{\mathsf{T}}X$$


![Type error](https://raw.githubusercontent.com/isaacmeza/power_simulations/main/table_type_error.png)


In general, we will describe (approximate) power as a function of sample size and effect size, and it will be a non-decreasing function in this two arguments.

In this repository you can find several examples of power simulations of varying complexity.
