# Power simulation 

Statistical power $(1-\beta)$ refers to the likelihood/probability of detecting an effect where a non-zero effect is present, i,e. given a null $H_0 : \theta = 0$, and its alternative $H_1 : \theta \neq 0$, the power is defined as

$$1-\beta = \Pr(\text{reject } H_0 \;|\;H_1  \text{ is true})$$

To estimate the power of our main specification, we use Monte Carlo simulation by replicating the data generating process of our experiment:

$$(X,Y)\sim F$$

Given this synthetic data, together with an identification method we will estimate the parameters of interest of the model and compute its power.
For instance, we might use a linear conditional mean model:

$$\mathbb{E}[Y|X]=\theta^{\mathsf{T}}X$$


\begin{table}[H]
  \centering
    \begin{tabular}{l|cc}
    \toprule
    \multicolumn{1}{r}{Conclusion} & \multicolumn{2}{c}{Reality} \\
    \midrule
    \midrule
     & $H_0$ TRUE & $H_0$ FALSE \\
     \midrule
    \multicolumn{1}{c|}{$H_0$ TRUE} & $(1-\alpha)$ & Type I : $\beta$ \\
    $H_0$ FALSE & Type II: $\alpha$ & $(1-\beta)$ \\
    \bottomrule
    \end{tabular}%
\end{table}%
