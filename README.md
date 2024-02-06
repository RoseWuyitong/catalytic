
<!-- README.md is generated from README.Rmd. Please edit that file -->

# catalytic

<!-- badges: start -->
<!-- badges: end -->

The goal of catalytic is to build framework for catalytic prior
distributions, as proposed in the paper “Catalytic prior distributions
with application to generalized linear models” by Huang et al. The goal
is to stabilize high-dimensional working models by shrinking them
towards simplified models. This is achieved by supplementing observed
data with weighted synthetic data generated from a predictive
distribution under the simpler model.

## Installation

You can install the development version of catalytic from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RoseWuyitong/catalytic")
```

## Example

This is a basic example which shows you how to solve a logistic
regression problem:

``` r
library(catalytic)
#> 
#> Attaching package: 'catalytic'
#> The following object is masked _by_ '.GlobalEnv':
#> 
#>     initialization
set.seed(1)

## Generate dummy data sets for a logistic regression model
N <- 25 # Number of observations data
p <- 10 # Number of predictors for the data
beta_true <- c(1, 2, -3, 4, -5, 6, -7, 8, -9, 10, 11) # true coefficients

obs_x <- matrix(rnorm(N * p), N, p)
obs_y <- rbinom(N, 1, stats::plogis(obs_x %*% beta_true[-1] + beta_true[1]))

## Initialization
data_ <- initialization(obs_x, obs_y, synthetic_response_generation_family = "binomial", regression_type = "logistic")
#> generating synthetic covariate
#> Column 1 is processing
#> Column 2 is processing
#> Column 3 is processing
#> Column 4 is processing
#> Column 5 is processing
#> Column 6 is processing
#> Column 7 is processing
#> Column 8 is processing
#> Column 9 is processing
#> Column 10 is processing
#> generating synthetic response
```

The synthetic covariate and response will be generated after
initialization:

``` r
head(data_$syn_x)
#>           V1         V2          V3         V4         V5          V6
#> 1 -0.6212406  1.3586796  1.43302370 -0.5425200  1.4322822 -0.30097613
#> 2  1.1249309  0.5566632 -0.13505460  0.2914462  1.3430388  0.53149619
#> 3  0.9438362 -0.6887557 -0.61202639  1.0630998  1.7672873  0.01739562
#> 4 -0.8204684 -0.3942900  1.43302370  1.0630998 -0.4616447 -0.52827990
#> 5  0.1836433 -1.4707524  0.47550953 -0.4734006  0.3841854 -1.11592011
#> 6  0.9189774 -0.4781501  0.02800216  1.1604026 -0.6203667 -0.75081900
#>            V7         V8          V9         V10
#> 1 -0.14439960  0.2199248 -0.63030033  0.19719344
#> 2 -0.92936215  0.5210227 -1.15657236 -0.05972328
#> 3 -0.31806837 -0.1587546 -0.01339952  0.56082073
#> 4  0.45699881 -1.0479844 -1.37020788  1.09677704
#> 5 -0.07715294 -0.7660820  0.64224131  1.16296456
#> 6  2.20610246 -1.2080828 -0.63030033 -0.25567071
```

``` r
head(data_$syn_y)
#> [1] 1 1 1 0 1 0
```

And all the information:

``` r
names(data_)
#>  [1] "obs_x"           "obs_y"           "comb_x"          "comb_y"         
#>  [5] "x_deg"           "syn_n"           "reg_type"        "recd_x"         
#>  [9] "obs_n"           "fit_fam"         "syn_y_gen_model" "x_col_n"        
#> [13] "syn_x"           "syn_y"
```

Then you can choose different way to get the point estimate:

``` r
## Single tau
data_ <- frequentist_prior(data_, taus = 10)
#> Will only get the point estimate for single tau value for faster process.

## Frequentist
data_ <- frequentist_prior(data_,
  taus = 1:10,
  risk_estimate_method = "steinian_estimate",
  discrepancy_function = "logistic_deviance"
)

data_ <- frequentist_prior(data_,
  taus = 1:10,
  risk_estimate_method = "cross_validation",
  discrepancy_function = "logistic_deviance"
)

## Bayesian
data_ <- bayesian_hyperpriors(data_)
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.34 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.491 seconds (Warm-up)
#> Chain 1:                0.393 seconds (Sampling)
#> Chain 1:                0.884 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.7e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.17 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.469 seconds (Warm-up)
#> Chain 2:                0.435 seconds (Sampling)
#> Chain 2:                0.904 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.8e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.473 seconds (Warm-up)
#> Chain 3:                0.309 seconds (Sampling)
#> Chain 3:                0.782 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1.7e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.17 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.492 seconds (Warm-up)
#> Chain 4:                0.447 seconds (Sampling)
#> Chain 4:                0.939 seconds (Total)
#> Chain 4:

data_ <- bayesian_hyperpriors(data_,
  joint = TRUE,
  chains = 3,
  iter = 500,
  warmup = 100
)
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.36 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 1: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1: Iteration: 101 / 500 [ 20%]  (Sampling)
#> Chain 1: Iteration: 150 / 500 [ 30%]  (Sampling)
#> Chain 1: Iteration: 200 / 500 [ 40%]  (Sampling)
#> Chain 1: Iteration: 250 / 500 [ 50%]  (Sampling)
#> Chain 1: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 1: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.017 seconds (Warm-up)
#> Chain 1:                0.056 seconds (Sampling)
#> Chain 1:                0.073 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.6e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: WARNING: There aren't enough warmup iterations to fit the
#> Chain 2:          three stages of adaptation as currently configured.
#> Chain 2:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 2:          the given number of warmup iterations:
#> Chain 2:            init_buffer = 15
#> Chain 2:            adapt_window = 75
#> Chain 2:            term_buffer = 10
#> Chain 2: 
#> Chain 2: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 2: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 2: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 2: Iteration: 101 / 500 [ 20%]  (Sampling)
#> Chain 2: Iteration: 150 / 500 [ 30%]  (Sampling)
#> Chain 2: Iteration: 200 / 500 [ 40%]  (Sampling)
#> Chain 2: Iteration: 250 / 500 [ 50%]  (Sampling)
#> Chain 2: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 2: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 2: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 2: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 2: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 2:                0.057 seconds (Sampling)
#> Chain 2:                0.078 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.6e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: WARNING: There aren't enough warmup iterations to fit the
#> Chain 3:          three stages of adaptation as currently configured.
#> Chain 3:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 3:          the given number of warmup iterations:
#> Chain 3:            init_buffer = 15
#> Chain 3:            adapt_window = 75
#> Chain 3:            term_buffer = 10
#> Chain 3: 
#> Chain 3: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 3: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 3: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 3: Iteration: 101 / 500 [ 20%]  (Sampling)
#> Chain 3: Iteration: 150 / 500 [ 30%]  (Sampling)
#> Chain 3: Iteration: 200 / 500 [ 40%]  (Sampling)
#> Chain 3: Iteration: 250 / 500 [ 50%]  (Sampling)
#> Chain 3: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 3: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 3: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 3: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 3: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 3:                0.051 seconds (Sampling)
#> Chain 3:                0.072 seconds (Total)
#> Chain 3:

## Results
data_$est_beta
#> [[1]]
#> [[1]]$method
#> [1] "single point estimate: tau =  10"
#> 
#> [[1]]$beta
#> (Intercept)          V1          V2          V3          V4          V5 
#> -0.12898872 -0.09277015  0.46968137  0.39755787 -0.45472511  0.53436596 
#>          V6          V7          V8          V9         V10 
#> -0.27720667  0.60650005 -0.89820432  0.80758469  0.87946610 
#> 
#> [[1]]$taus
#> [1] 10
#> 
#> [[1]]$risk
#> [1] NA
#> 
#> 
#> [[2]]
#> [[2]]$method
#> [1] "frequentist: risk estimate method = steinian_estimate ; discrepancy function = logistic_deviance ; bootstrap iteration times = 100 ; tau 0 =  1 ; cv fold num =  5"
#> 
#> [[2]]$optimal_tau
#> [1] 2
#> 
#> [[2]]$beta
#> (Intercept)          V1          V2          V3          V4          V5 
#>  -0.1179750   0.3057027   1.1470119   0.6804794  -1.2325746   0.9597947 
#>          V6          V7          V8          V9         V10 
#>  -0.6676486   1.9824179  -1.9746820   2.1169091   2.3036149 
#> 
#> [[2]]$taus
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[2]]$risk
#>  [1] 5.210595 4.988431 5.194170 5.528119 5.898986 6.270088 6.625710 6.959717
#>  [9] 7.270623 7.559158
#> 
#> 
#> [[3]]
#> [[3]]$method
#> [1] "frequentist: risk estimate method = cross_validation ; discrepancy function = logistic_deviance ; bootstrap iteration times = 100 ; tau 0 =  1 ; cv fold num =  5"
#> 
#> [[3]]$optimal_tau
#> [1] 1
#> 
#> [[3]]$beta
#> (Intercept)          V1          V2          V3          V4          V5 
#>  -0.1194295   0.5074289   1.5211205   0.8208972  -1.5967941   1.1866903 
#>          V6          V7          V8          V9         V10 
#>  -0.8663213   2.7051720  -2.5754138   2.8441552   3.0738984 
#> 
#> [[3]]$taus
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[3]]$risk
#>  [1] 1.732437 1.762122 2.274482 2.477139 2.305302 3.270654 2.220932 1.842939
#>  [9] 2.404525 2.334275
#> 
#> 
#> [[4]]
#> [[4]]$method
#> [1] "bayesian: joint = FALSE ; tau =  2.5 ; chains =  4 ; iter =  2000 ; warmup =  1000"
#> 
#> [[4]]$beta_median
#>  [1]  0.4017043  2.6774610  1.5247548 -2.2613643  2.5360893 -1.8316026
#>  [7]  4.1066374 -4.0810769  4.6537178  5.0941590
#> 
#> [[4]]$beta_lower_iqr
#>  [1] -2.3292451 -1.9558797 -1.3093614 -6.5748463 -1.1065238 -5.9512497
#>  [7]  0.8337115 -9.2513612  1.7071399  1.9488505
#> 
#> [[4]]$beta_higher_iqr
#>  [1]  3.4379749  8.3313037  4.8645181  0.9481895  8.1907729  1.2991286
#>  [7]  9.0214365 -0.5624854  9.2656509 10.2523937
#> 
#> [[4]]$beta_mean
#>  [1]  0.443834  2.784132  1.583725 -2.392700  2.778007 -1.957504  4.323497
#>  [8] -4.292634  4.886253  5.356910
#> 
#> [[4]]$tau_median
#> [1] -0.5839221
#> 
#> [[4]]$tau_lower_iqr
#> [1] -4.424293
#> 
#> [[4]]$tau_higher_iqr
#> [1] 2.788672
#> 
#> [[4]]$tau_mean
#> [1] -0.6313193
#> 
#> [[4]]$rhat_pass
#> [1] "PASS"
#> 
#> [[4]]$ess_pass
#> [1] "PASS"
#> 
#> [[4]]$rhat
#>   beta[1]   beta[2]   beta[3]   beta[4]   beta[5]   beta[6]   beta[7]   beta[8] 
#> 0.9990744 0.9997266 0.9997015 0.9995278 0.9993509 1.0013878 1.0002818 1.0011474 
#>   beta[9]  beta[10]  beta[11] 
#> 1.0007296 1.0011203 1.0002239 
#> 
#> [[4]]$ess
#>  beta[1]  beta[2]  beta[3]  beta[4]  beta[5]  beta[6]  beta[7]  beta[8] 
#> 5914.109 4585.647 3525.830 3499.991 3483.391 2623.307 4362.408 2895.566 
#>  beta[9] beta[10] beta[11] 
#> 2417.388 2685.695 2308.284 
#> 
#> [[4]]$model
#> Inference for Stan model: anon_model.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>            mean se_mean   sd   2.5%    25%    50%    75% 97.5% n_eff Rhat
#> beta[1]   -0.63    0.02 1.81  -4.42  -1.74  -0.58   0.54  2.79  5914    1
#> beta[2]    0.44    0.02 1.46  -2.33  -0.55   0.40   1.42  3.44  4586    1
#> beta[3]    2.78    0.04 2.62  -1.96   1.00   2.68   4.36  8.33  3526    1
#> beta[4]    1.58    0.03 1.55  -1.31   0.57   1.52   2.56  4.86  3500    1
#> beta[5]   -2.39    0.03 1.87  -6.57  -3.52  -2.26  -1.12  0.95  3483    1
#> beta[6]    2.78    0.05 2.34  -1.11   1.16   2.54   4.13  8.19  2623    1
#> beta[7]   -1.96    0.03 1.87  -5.95  -3.11  -1.83  -0.67  1.30  4362    1
#> beta[8]    4.32    0.04 2.14   0.83   2.81   4.11   5.59  9.02  2896    1
#> beta[9]   -4.29    0.04 2.19  -9.25  -5.64  -4.08  -2.70 -0.56  2417    1
#> beta[10]   4.89    0.04 1.97   1.71   3.46   4.65   6.05  9.27  2686    1
#> beta[11]   5.36    0.04 2.12   1.95   3.83   5.09   6.60 10.25  2308    1
#> lp__     -14.59    0.09 3.03 -21.20 -16.42 -14.21 -12.34 -9.86  1070    1
#> 
#> Samples were drawn using HMC(diag_e) at Tue Feb  6 11:08:14 2024.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
#> 
#> 
#> [[5]]
#> [[5]]$method
#> [1] "bayesian: joint = TRUE ; tau =   ; chains =  3 ; iter =  500 ; warmup =  100"
#> 
#> [[5]]$beta_median
#>  [1]   0.006524268 -21.244721606  15.412801333  77.130105192  30.838857630
#>  [6] -35.982668060  52.577153789 -47.098517879 131.661832411 -62.605142402
#> [11] 112.158469210 135.556242992
#> 
#> [[5]]$beta_lower_iqr
#>  [1]  4.225685e-04 -8.032600e+01 -3.255176e+01 -1.733734e+01 -1.684637e+01
#>  [6] -7.671598e+01 -1.281618e+01 -9.531378e+01  5.543255e+01 -1.452628e+02
#> [11]  4.979661e+01  4.514558e+01
#> 
#> [[5]]$beta_higher_iqr
#>  [1]   0.04168677  21.56830411  60.01653206 151.23142188  95.74278053
#>  [6]  29.58786157 106.63506719  -7.00188762 208.52815032 -18.17658208
#> [11] 229.12107777 203.54817846
#> 
#> [[5]]$beta_mean
#>  [1]   0.01047889 -22.57973224  14.07313594  74.29929655  31.73075381
#>  [6] -30.78463306  49.80320882 -48.88333408 135.96471895 -67.32547667
#> [11] 122.08430916 131.37422168
#> 
#> [[5]]$rhat_pass
#> [1] "FAIL"
#> 
#> [[5]]$ess_pass
#> [1] "FAIL"
#> 
#> [[5]]$rhat
#>      tau  beta[1]  beta[2]  beta[3]  beta[4]  beta[5]  beta[6]  beta[7] 
#> 1.042249 1.335781 1.263224 2.129662 1.306323 1.226805 1.851544 1.253382 
#>  beta[8]  beta[9] beta[10] beta[11] 
#> 1.938067 1.295102 1.925832 1.649967 
#> 
#> [[5]]$ess
#>        tau    beta[1]    beta[2]    beta[3]    beta[4]    beta[5]    beta[6] 
#> 132.062318   6.872445  16.187986   3.499823  11.863304  10.665574   6.067474 
#>    beta[7]    beta[8]    beta[9]   beta[10]   beta[11] 
#>  15.637190   5.819602  10.238903   4.964817   8.622200 
#> 
#> [[5]]$model
#> Inference for Stan model: anon_model.
#> 3 chains, each with iter=500; warmup=100; thin=1; 
#> post-warmup draws per chain=400, total post-warmup draws=1200.
#> 
#>            mean se_mean    sd    2.5%    25%    50%    75%  97.5% n_eff Rhat
#> tau        0.01    0.00  0.01    0.00   0.00   0.01   0.01   0.04   132 1.04
#> beta[1]  -22.58   10.13 26.56  -80.33 -40.68 -21.24  -1.26  21.57     7 1.34
#> beta[2]   14.07    5.76 23.19  -32.55  -1.68  15.41  29.03  60.02    16 1.26
#> beta[3]   74.30   25.16 47.08  -17.34  38.89  77.13 106.07 151.23     3 2.13
#> beta[4]   31.73    7.56 26.05  -16.85  14.16  30.84  44.95  95.74    12 1.31
#> beta[5]  -30.78    9.04 29.52  -76.72 -53.02 -35.98 -10.47  29.59    11 1.23
#> beta[6]   49.80   12.68 31.23  -12.82  27.89  52.58  73.12 106.64     6 1.85
#> beta[7]  -48.88    6.00 23.71  -95.31 -67.66 -47.10 -31.75  -7.00    16 1.25
#> beta[8]  135.96   18.29 44.13   55.43 106.76 131.66 177.16 208.53     6 1.94
#> beta[9]  -67.33   10.01 32.02 -145.26 -87.21 -62.61 -42.72 -18.18    10 1.30
#> beta[10] 122.08   19.93 44.41   49.80  94.74 112.16 143.63 229.12     5 1.93
#> beta[11] 131.37   14.83 43.55   45.15 104.60 135.56 163.31 203.55     9 1.65
#> lp__      -6.32    0.07  0.86   -8.57  -6.63  -6.12  -5.79  -5.21   176 1.01
#> 
#> Samples were drawn using HMC(diag_e) at Tue Feb  6 11:08:35 2024.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```
