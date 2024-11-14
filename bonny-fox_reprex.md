``` r
devtools::check()
#> ══ Documenting ═════════════════════════════════════════════════════════════════
#> ℹ Updating catalytic documentation
#> ℹ Loading catalytic
#> Writing 'get_cox_gradient.Rd'
#> 
#> ══ Building ════════════════════════════════════════════════════════════════════
#> Setting env vars:
#> • CFLAGS    : -Wall -pedantic
#> • CXXFLAGS  : -Wall -pedantic
#> • CXX11FLAGS: -Wall -pedantic
#> • CXX14FLAGS: -Wall -pedantic
#> • CXX17FLAGS: -Wall -pedantic
#> • CXX20FLAGS: -Wall -pedantic
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/Users/rosewu/catalytic/DESCRIPTION’ ... OK
#> * preparing ‘catalytic’:
#> * checking DESCRIPTION meta-information ... OK
#> * installing the package to build vignettes
#> * creating vignettes ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘catalytic_0.0.0.9000.tar.gz’
#> 
#> ══ Checking ════════════════════════════════════════════════════════════════════
#> Setting env vars:
#> • _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
#> • _R_CHECK_CRAN_INCOMING_                      : FALSE
#> • _R_CHECK_FORCE_SUGGESTS_                     : FALSE
#> • _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
#> • NOT_CRAN                                     : true
#> ── R CMD check ─────────────────────────────────────────────────────────────────
#> * using log directory ‘/private/var/folders/j5/bktz8wt94hv7m_wtcc0k6ymh0000gn/T/RtmphvY1Kv/filed7ed67ef9acf/catalytic.Rcheck’
#> * using R version 4.4.0 (2024-04-24)
#> * using platform: aarch64-apple-darwin20
#> * R was compiled by
#>     Apple clang version 14.0.0 (clang-1400.0.29.202)
#>     GNU Fortran (GCC) 12.2.0
#> * running under: macOS 15.1
#> * using session charset: UTF-8
#> * using options ‘--no-manual --as-cran’
#> * checking for file ‘catalytic/DESCRIPTION’ ... OK
#> * this is package ‘catalytic’ version ‘0.0.0.9000’
#> * package encoding: UTF-8
#> * checking package namespace information ... OK
#> * checking package dependencies ... OK
#> * checking if this is a source package ... OK
#> * checking if there is a namespace ... OK
#> * checking for executable files ... OK
#> * checking for hidden files and directories ... OK
#> * checking for portable file names ... OK
#> * checking for sufficient/correct file permissions ... OK
#> * checking whether package ‘catalytic’ can be installed ... OK
#> * checking installed package size ... OK
#> * checking package directory ... OK
#> * checking for future file timestamps ... NOTE
#> unable to verify current time
#> * checking ‘build’ directory ... OK
#> * checking DESCRIPTION meta-information ... OK
#> * checking top-level files ... NOTE
#> Non-standard files/directories found at top level:
#>   ‘bonny-fox_reprex.R’ ‘bonny-fox_reprex.spin.R’
#>   ‘bonny-fox_reprex.spin.Rmd’
#> * checking for left-over files ... OK
#> * checking index information ... OK
#> * checking package subdirectories ... OK
#> * checking code files for non-ASCII characters ... OK
#> * checking R files for syntax errors ... OK
#> * checking whether the package can be loaded ... OK
#> * checking whether the package can be loaded with stated dependencies ... OK
#> * checking whether the package can be unloaded cleanly ... OK
#> * checking whether the namespace can be loaded with stated dependencies ... OK
#> * checking whether the namespace can be unloaded cleanly ... OK
#> * checking loading without being on the library search path ... OK
#> * checking dependencies in R code ... NOTE
#> Namespaces in Imports field not imported from:
#>   ‘memoise’ ‘pROC’
#>   All declared Imports should be used.
#> * checking S3 generic/method consistency ... OK
#> * checking replacement functions ... OK
#> * checking foreign function calls ... OK
#> * checking R code for possible problems ... OK
#> * checking Rd files ... OK
#> * checking Rd metadata ... OK
#> * checking Rd line widths ... OK
#> * checking Rd cross-references ... OK
#> * checking for missing documentation entries ... OK
#> * checking for code/documentation mismatches ... OK
#> * checking Rd \usage sections ... OK
#> * checking Rd contents ... OK
#> * checking for unstated dependencies in examples ... OK
#> * checking contents of ‘data’ directory ... OK
#> * checking data for non-ASCII characters ... OK
#> * checking LazyData ... OK
#> * checking data for ASCII and uncompressed saves ... OK
#> * checking installed files from ‘inst/doc’ ... OK
#> * checking files in ‘vignettes’ ... OK
#> * checking examples ... [173s/206s] ERROR
#> Running examples in ‘catalytic-Ex.R’ failed
#> The error most likely occurred in:
#> > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
#> > ### Name: get_cox_hessian
#> > ### Title: Compute the Hessian Matrix for Cox Proportional Hazards Model
#> > ### Aliases: get_cox_hessian
#> > 
#> > ### ** Examples
#> > 
#> > hessian_matrix <- get_cox_hessian(
#> +   X = matrix(stats::rnorm(15), ncol = 3),
#> +   time = stats::rexp(5, rate = 0.1),
#> +   status = stats::rbinom(5, 1, 0.5),
#> +   coefs = stats::rnorm(3),
#> +   entry_points = rep(0, 5)
#> + )
#> Error in get_cox_hessian(X = matrix(stats::rnorm(15), ncol = 3), time = stats::rexp(5,  : 
#>   could not find function "get_cox_hessian"
#> Execution halted
#> * checking for unstated dependencies in ‘tests’ ... OK
#> * checking tests ...
#>   Running ‘testthat.R’ [19m/22m]
#>  [19m/22m] OK
#> * checking for unstated dependencies in vignettes ... OK
#> * checking package vignettes ... OK
#> * checking re-building of vignette outputs ... [14m/15m] OK
#> * checking for non-standard things in the check directory ... OK
#> * checking for detritus in the temp directory ... OK
#> * DONE
#> 
#> Status: 1 ERROR, 3 NOTEs
#> See
#>   ‘/private/var/folders/j5/bktz8wt94hv7m_wtcc0k6ymh0000gn/T/RtmphvY1Kv/filed7ed67ef9acf/catalytic.Rcheck/00check.log’
#> for details.
#> 
#> ── R CMD check results ─────────────────────────────── catalytic 0.0.0.9000 ────
#> Duration: 40m 1.8s
#> 
#> ❯ checking examples ... [173s/206s] ERROR
#>   Running examples in ‘catalytic-Ex.R’ failed
#>   The error most likely occurred in:
#>   
#>   > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
#>   > ### Name: get_cox_hessian
#>   > ### Title: Compute the Hessian Matrix for Cox Proportional Hazards Model
#>   > ### Aliases: get_cox_hessian
#>   > 
#>   > ### ** Examples
#>   > 
#>   > hessian_matrix <- get_cox_hessian(
#>   +   X = matrix(stats::rnorm(15), ncol = 3),
#>   +   time = stats::rexp(5, rate = 0.1),
#>   +   status = stats::rbinom(5, 1, 0.5),
#>   +   coefs = stats::rnorm(3),
#>   +   entry_points = rep(0, 5)
#>   + )
#>   Error in get_cox_hessian(X = matrix(stats::rnorm(15), ncol = 3), time = stats::rexp(5,  : 
#>     could not find function "get_cox_hessian"
#>   Execution halted
#> 
#> ❯ checking for future file timestamps ... NOTE
#>   unable to verify current time
#> 
#> ❯ checking top-level files ... NOTE
#>   Non-standard files/directories found at top level:
#>     ‘bonny-fox_reprex.R’ ‘bonny-fox_reprex.spin.R’
#>     ‘bonny-fox_reprex.spin.Rmd’
#> 
#> ❯ checking dependencies in R code ... NOTE
#>   Namespaces in Imports field not imported from:
#>     ‘memoise’ ‘pROC’
#>     All declared Imports should be used.
#> 
#> 1 error ✖ | 0 warnings ✔ | 3 notes ✖
#> Error: R CMD check found ERRORs
```

<sup>Created on 2024-11-14 with [reprex v2.1.1](https://reprex.tidyverse.org)</sup>
