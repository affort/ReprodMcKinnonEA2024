# Independent project: Reproducing McKinnonEA2024

This independent project aims to verify the power of the statistical test proposed in the following paper:
>K.A. McKinnon, I.R. Simpson, A.P. Williams, The pace of change of summertime temperature extremes, *Proc. Natl. Acad. Sci. U.S.A.* 121 (42) e2406143121, https://doi.org/10.1073/pnas.2406143121 (2024).

The goal is reproduce figures S2(d) and S2(e) of the paper. 
I then introduce two new cases to check whether the test retains its power under alternative specifications of the data.

The code produced and shared by the authors, `synthetic_data_tests.py`, 
was used as a reference. My reproduced code is in `Reprod_Synthetic_McKinnonEA2024.r`. An Rmarkdown version is in `Reprod_Synthetic_McKinnonEA2024.rmd`, and the output is visible in `Reprod_Synthetic_McKinnonEA2024.pdf`.

Please email me (Teng Wei Yeo) at yeotengwei@gmail.com if you have any questions.

## Summary
McKinnon et al. (2024, Fig. S2) propose a non-parametric test for detecting trends in the magnitude of the difference between extremes and the median over time. 
The null hypothesis is that there is no change in the variance. 
For each year from 1959 to 2023 (for a total of 65 years), they find the annual maximum deviation between the highest recorded temperature of that summer, and the median temperature of that summer.
They then rank these deviations from 1 to 65 in order of magnitude, then do the following linear regression: 

$E([\textnormal{rank}]_i) = \beta_0 + \beta_1[\textnormal{year}]_i, i = 1,...,65$.

They check whether the trend in these ranks are statistically significant or not by checking for the significance of $\hat{\beta_1}$. 
They do so by finding the 95% confidence interval based on $N$ simulations from the null hypothesis: simulate data from a normal distribution that has no change in variance over time, calculate the order of the ranks for each of the $N$ simulations, calculate $\hat{\beta_1}^{(\textnormal{sim})}$ for each simulation, then find the 2.5 and 97.5 percentile of the $\hat{\beta_1}^{(\textnormal{sim})}$ values. 
Reject $\hat{\beta_1}$ if $\hat{\beta_1} \leq \hat{\beta_1}^{(\textnormal{sim, 2.5pct})}$ or $\hat{\beta_1} \geq \hat{\beta_1}^{(\textnormal{sim, 97.5pct})}$.

McKinnon et al. (2024) show that the test's power is strong, in that it is sensitive enough to reject the null hypothesis when the data has linearly-increasing variance over time. 
I.e., if the variance is $x$% higher in 2023 than it was in 1959, and if the variance increased linearly from $v$ in 1959 to $v\left(1+\frac{x}{100}\right)$ in 2023, the test will reject the null almost 100% of the time.
I reproduced their simulations and got the same result.
I then introduced two new test cases --- I check whether the test's power remains strong when detecting a (i) quadratically increasing variance over time, and (ii) exponentially increasing variance over time, for various values of $x$.

I show that, for a given rejection region with 5% Type I error rate, the test’s power is strong if the data’s variance increases linearly or quadratically over time, but weakens if the increase is exponential over time. 
This is visualized in `Main_Plot.pdf`, or in page 9 of `Reprod_Synthetic_McKinnonEA2024.pdf`.
