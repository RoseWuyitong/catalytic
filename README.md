
<!-- README.md is generated from README.Rmd. Please edit that file -->

# catalytic

<!-- badges: start -->
<!-- badges: end -->

The goal of catalytic is to …

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
set.seed(1)

## Generate dummy data sets for a logistic regression model
N <- 25 # Number of observations
p <- 10  # Number predictors for the data
beta_true <- c(1, 2, -3, 4, -5, 6, -7, 8, -9, 10, 11)  # true coefficients 

obs_x <- matrix(rnorm(N * p), N, p)
obs_y <- rbinom(N, 1, stats::plogis(obs_x %*% beta_true[-1] + beta_true[1]))

## Initialization
data_ <- initialization(obs_x, obs_y, synthetic_response_generation_family="binomial", regression_type = "logistic")
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
data_ <- frequentist_prior(data_, taus=10)
#> Will only get the point estimate for single tau value for faster process.

## Multiple taus 
data_ <- frequentist_prior(data_, 
                           taus=1:10, 
                           risk_estimate_method = "parametric_bootstrap",
                           discrepancy_function = "logistic_deviance")

data_ <- frequentist_prior(data_,
                           taus = 1:10,
                           risk_estimate_method = "steinian_estimate",
                           discrepancy_function = "logistic_deviance")

data_ <- frequentist_prior(data_,
                           taus = 1:10,
                           risk_estimate_method = "cross_validation",
                           discrepancy_function = "logistic_deviance")

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
#> [1] "frequentist: risk estimate method = parametric_bootstrap ; discrepancy function = logistic_deviance ; bootstrap iteration times = 100 ; tau 0 =  1 ; cv fold num =  5 ; optimal tau =  3"
#> 
#> [[2]]$beta
#> (Intercept)          V1          V2          V3          V4          V5 
#>  -0.1079972   0.1874494   0.9316962   0.6039026  -1.0245888   0.8300029 
#>          V6          V7          V8          V9         V10 
#>  -0.5466316   1.5734161  -1.6575571   1.7190390   1.8815226 
#> 
#> [[2]]$taus
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[2]]$risk
#>  [1] 5.689290 5.205159 5.078246 5.115166 5.943092 6.458870 6.725714 6.980777
#>  [9] 7.249199 7.512305
#> 
#> 
#> [[3]]
#> [[3]]$method
#> [1] "frequentist: risk estimate method = steinian_estimate ; discrepancy function = logistic_deviance ; bootstrap iteration times = 100 ; tau 0 =  1 ; cv fold num =  5 ; optimal tau =  2"
#> 
#> [[3]]$beta
#> (Intercept)          V1          V2          V3          V4          V5 
#>  -0.1179750   0.3057027   1.1470119   0.6804794  -1.2325746   0.9597947 
#>          V6          V7          V8          V9         V10 
#>  -0.6676486   1.9824179  -1.9746820   2.1169091   2.3036149 
#> 
#> [[3]]$taus
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[3]]$risk
#>  [1] 5.210595 4.988431 5.194170 5.528119 5.898986 6.270088 6.625710 6.959717
#>  [9] 7.270623 7.559158
#> 
#> 
#> [[4]]
#> [[4]]$method
#> [1] "frequentist: risk estimate method = cross_validation ; discrepancy function = logistic_deviance ; bootstrap iteration times = 100 ; tau 0 =  1 ; cv fold num =  5 ; optimal tau =  3"
#> 
#> [[4]]$beta
#> (Intercept)          V1          V2          V3          V4          V5 
#>  -0.1079972   0.1874494   0.9316962   0.6039026  -1.0245888   0.8300029 
#>          V6          V7          V8          V9         V10 
#>  -0.5466316   1.5734161  -1.6575571   1.7190390   1.8815226 
#> 
#> [[4]]$taus
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[4]]$risk
#>  [1] 1.921737 2.706657 1.741185 2.317454 2.552712 2.205333 2.661486 3.173324
#>  [9] 2.072023 2.285338
```
