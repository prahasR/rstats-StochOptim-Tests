Tests
================

## Test 1

Prepare R code for three of the unconstrained test functions in [GO Test
Problems](http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm)
that allow dimension greater than 4

-   **Rosenbrock**

``` r
Rosenbrock<-function(x){
    len<-length(x)
    sum((100*(x[-len]^2-x[-1])^2) + (x[-len]-rep(1, len-1))^2)
}

#gradient function
Rosenbrock.g<-function(x){
  len<-length(x)
  g<-rep(NA, len)
  g[1]<- 2*(x[1]-1) + 400*x[1]*(x[1]^2-x[2])
  g[2:(len-1)] <- 2*(x[2:(len-1)]-1) + 400*x[2:(len-1)]*(x[2:(len-1)]^2-x[2:(len-1) + 1]) + 200*(x[2:(len-1)]-x[2:(len-1) - 1]^2)
  g[len] <- 200 * (x[len] - x[len - 1]^2)
  g
}
Rosenbrock(c(1,2,3,4))
```

    ## [1] 2705

-   **Sphere Func.**

``` r
sphere<-function(x){
    sum(x^2)
}

#gradient function
sphere.g<-function(x){
  len<-length(x)
  g<-rep(NA, len)
  g[1:len]<-2*x
  g
}

sphere(c(1,2,3,4))
```

    ## [1] 30

-   **Sum square**

``` r
sum_sq<-function(x){
  len<-length(x)
    sum(seq(1, len, length=len)*(x^2))
}

#gradient function
sum_sq.g<-function(x){
  len<-length(x)
  g<-rep(NA, len)
  g[1:len]<-seq(1, len, length=len)*2*x
  g
}

sum_sq(c(1,2,3,4))
```

    ## [1] 100

-   **Dixon&Price**

``` r
d_and_p<-function(x){
    len<-length(x)
    t1<-(x[1]-1)^2
    t2<-sum(seq(2,len, length=len-1)*(2*x[-1]^2-x[-len])^2)
    t1+t2
}

#gradient function
d_and_p.g<-function(x){
  len<-length(x)
  g<-rep(NA, len)
  g[1]<-2*(x[1]-1) - 4*(2*x[2]^2-x[1])
  g[2:(len-1)]<- 8*seq(2,len-1, length=len-2)*x[2:(len-1)]*(2*x[2:(len-1)]^2-x[1:(len-2)]) - 2*seq(3,len, length=len-2)*(2*x[3:len]^2-x[2:(len-1)])
  g[len]<-8*len*x[len]*(2*x[len]^2-x[len-1])
  g
}
d_and_p(c(1,2,3,4))
```

    ## [1] 4230

## Test 2

Try to minimize these with the base R **optim()** function. Be sure to
document what you do.

``` r
#Rosenbrock
optm_1<-optim(c(1,2,3,4), Rosenbrock, Rosenbrock.g, method="BFGS")
optm_1
```

    ## $par
    ## [1] 1 1 1 1
    ## 
    ## $value
    ## [1] 5.676878e-17
    ## 
    ## $counts
    ## function gradient 
    ##       99       37 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

``` r
#sphere
optm_2<-optim(c(1,2,3,4), sphere, sphere.g, method="CG")
optm_2
```

    ## $par
    ## [1] 2.204415e-07 4.408831e-07 6.613246e-07 8.817661e-07
    ## 
    ## $value
    ## [1] 1.457834e-12
    ## 
    ## $counts
    ## function gradient 
    ##       15        8 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

``` r
#Sum square
optm_3<-optim(c(1,2,3,4), sum_sq, method="SANN")
optm_3
```

    ## $par
    ## [1] -0.05651465  0.08942125  0.04428693 -0.02979915
    ## 
    ## $value
    ## [1] 0.02862218
    ## 
    ## $counts
    ## function gradient 
    ##    10000       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

``` r
#dixon&price
optm_4<-optim(c(1,2,3,4), d_and_p, method="Nelder-Mead")
optm_4
```

    ## $par
    ## [1]  0.9999534  0.7081443  0.5953718 -0.5449882
    ## 
    ## $value
    ## [1] 2.694426e-05
    ## 
    ## $counts
    ## function gradient 
    ##      295       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

## Test 3

Run several (at least 4) solvers at once with the **opm()** function

``` r
require(optimr)
```

    ## Loading required package: optimr

``` r
#Rosenbrock
result1<-opm(c(1,2,3,4), Rosenbrock, Rosenbrock.g, method="ALL", control=list(kkt=FALSE, trace=0))
result1
```

    ##                    p1        p2        p3        p4        value fevals gevals
    ## BFGS        1.0000000 1.0000000 1.0000000 1.0000000 5.676878e-17     99     37
    ## CG          1.1117180 1.2374887 1.5330875 2.3535077 3.534559e-01   3298   1001
    ## Nelder-Mead 0.9999500 0.9995963 0.9987679 0.9973660 3.189882e-05    513     NA
    ## L-BFGS-B    0.9999997 0.9999995 0.9999990 0.9999984 1.207693e-11     43     43
    ## nlm         1.0000000 1.0000000 0.9999999 0.9999999 7.364247e-14     NA     79
    ## nlminb      1.0000000 1.0000000 1.0000000 1.0000000 1.887858e-16     54     40
    ## Rcgmin      1.0000001 1.0000002 1.0000003 1.0000006 1.316622e-13    137     63
    ## Rvmmin      1.0000000 1.0000000 1.0000000 1.0000000 3.722437e-30     65     44
    ## hjn         1.0000006 1.0000012 1.0000024 1.0000047 8.483011e-12   3500     NA
    ##             convergence kkt1 kkt2 xtime
    ## BFGS                  0   NA   NA 0.001
    ## CG                    1   NA   NA 0.016
    ## Nelder-Mead           0   NA   NA 0.001
    ## L-BFGS-B              0   NA   NA 0.000
    ## nlm                   0   NA   NA 0.002
    ## nlminb                0   NA   NA 0.001
    ## Rcgmin                0   NA   NA 0.003
    ## Rvmmin                0   NA   NA 0.004
    ## hjn                   0   NA   NA 0.019

``` r
#sphere
result2<-opm(c(1,2,3,4), sphere, sphere.g, method="ALL", control=list(kkt=FALSE, trace=0))
```

    ## Warning in nlminb(start = spar, objective = efn, gradient = egr, lower =
    ## slower, : NA/NaN function evaluation

``` r
result2
```

    ##                         p1             p2             p3             p4
    ## BFGS         -2.183439e-16  -4.366877e-16   1.221245e-16   3.700743e-16
    ## CG            2.204415e-07   4.408831e-07   6.613246e-07   8.817661e-07
    ## Nelder-Mead  -1.261488e-04   4.478612e-04   9.807649e-05   1.402697e-04
    ## L-BFGS-B     -9.860761e-32  -1.972152e-31  -3.944305e-31  -3.944305e-31
    ## nlm           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
    ## nlminb       3.133255e-167 -1.348723e-167 -1.184606e-167  2.739049e-167
    ## Rcgmin       -6.661338e-16  -1.332268e-15  -1.776357e-15  -2.664535e-15
    ## Rvmmin        1.110223e-16   2.220446e-16  -2.220446e-16  -4.440892e-16
    ## hjn           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
    ##                    value fevals gevals convergence kkt1 kkt2 xtime
    ## BFGS        2.179639e-30      5      3           0   NA   NA 0.000
    ## CG          1.457834e-12     15      8           0   NA   NA 0.000
    ## Nelder-Mead 2.457878e-07    207     NA           0   NA   NA 0.000
    ## L-BFGS-B    3.597681e-61      4      4           0   NA   NA 0.000
    ## nlm         0.000000e+00     NA      1           0   NA   NA 0.000
    ## nlminb      0.000000e+00     30     27           0   NA   NA 0.000
    ## Rcgmin      1.247386e-29      4      2           0   NA   NA 0.001
    ## Rvmmin      3.081488e-31      4      3           2   NA   NA 0.000
    ## hjn         0.000000e+00    128     NA           0   NA   NA 0.001

``` r
#sum_square
result3<-opm(c(1,2,3,4), sum_sq, sum_sq.g, method="ALL", control=list(kkt=FALSE, trace=0))
```

    ## Warning in nlminb(start = spar, objective = efn, gradient = egr, lower =
    ## slower, : NA/NaN function evaluation

``` r
result3
```

    ##                         p1             p2             p3             p4
    ## BFGS         -1.941783e-14  -3.049467e-13  -1.766722e-12   7.956968e-13
    ## CG            5.703672e-07  -7.784541e-08   4.351182e-08   9.908794e-08
    ## Nelder-Mead  -5.818810e-04   3.955017e-05   3.673549e-05  -2.442033e-04
    ## L-BFGS-B      1.892986e-07   3.767959e-06   2.396254e-06   5.000826e-06
    ## nlm           1.623039e-07   4.749345e-08   1.428741e-08  -1.450018e-08
    ## nlminb       3.201704e-168 -5.247493e-169  2.560660e-168 -1.200510e-168
    ## Rcgmin       -2.775558e-17  -1.387779e-17  -1.734723e-17   1.218643e-16
    ## Rvmmin       -6.706204e-16  -6.761823e-16  -2.509268e-16  -2.508378e-15
    ## hjn           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
    ##                    value fevals gevals convergence kkt1 kkt2 xtime
    ## BFGS        1.208497e-23     37     17           0   NA   NA 0.004
    ## CG          3.823921e-13     51     21           0   NA   NA 0.001
    ## Nelder-Mead 5.843035e-07    165     NA           0   NA   NA 0.002
    ## L-BFGS-B    1.456900e-10     10     10           0   NA   NA 0.001
    ## nlm         3.230721e-14     NA     11           0   NA   NA 0.003
    ## nlminb      0.000000e+00     55     49           0   NA   NA 0.002
    ## Rcgmin      6.146199e-32     10      5           0   NA   NA 0.000
    ## Rvmmin      2.672091e-29     20     16           0   NA   NA 0.001
    ## hjn         0.000000e+00    128     NA           0   NA   NA 0.002

``` r
#dixon&price
result4<-opm(c(1,2,3,4), d_and_p, d_and_p.g, method="ALL", control=list(kkt=FALSE, trace=0))
result4
```

    ##                    p1        p2        p3         p4        value fevals gevals
    ## BFGS        1.0000000 0.7071068 0.5946036 -0.5452539 3.899315e-20     57     26
    ## CG          0.9999991 0.7071064 0.5946033 -0.5452537 9.031343e-13    150     55
    ## Nelder-Mead 0.9999534 0.7081443 0.5953718 -0.5449882 2.694426e-05    295     NA
    ## L-BFGS-B    0.9999982 0.7071062 0.5946031  0.5452535 4.912154e-12     23     23
    ## nlm         0.9999999 0.7071068 0.5946035 -0.5452538 9.342096e-14     NA     31
    ## nlminb      1.0000000 0.7071068 0.5946036  0.5452539 6.132918e-22     20     18
    ## Rcgmin      1.0000000 0.7071068 0.5946036  0.5452539 1.338168e-15     48     25
    ## Rvmmin      1.0000000 0.7071068 0.5946036 -0.5452539 3.821045e-29     37     27
    ## hjn         1.0000001 0.7071068 0.5946036  0.5452539 5.463278e-14    660     NA
    ##             convergence kkt1 kkt2 xtime
    ## BFGS                  0   NA   NA 0.017
    ## CG                    0   NA   NA 0.003
    ## Nelder-Mead           0   NA   NA 0.003
    ## L-BFGS-B              0   NA   NA 0.001
    ## nlm                   0   NA   NA 0.001
    ## nlminb                0   NA   NA 0.001
    ## Rcgmin                0   NA   NA 0.001
    ## Rvmmin                0   NA   NA 0.003
    ## hjn                   0   NA   NA 0.009

#### Observations from above Result

-   For problems which has quite significant scaling factor(Eg:
    Rosenbrock function), non gradient methods such as **Nelder-Mead**
    and **hjn** does large function evaluations to arrive at results.

-   For simpler problems such as Sphere function, **Nelder-Mead** and
    **hjn** requires relatively less number of evaluations as compared
    to Rosenbrock problem but still it is larger then as required by
    other gradient optimization methods for the same problem.

-   **CG** method reaches its iteration limit and does not able to find
    the minimum for the Rosenbrock function. Also it can be seen that
    this method takes relatively more function and gradient computations
    compared to other gradient methods.

-   **Rvmmin** method achieves better results compared to other methods
    and for the Rosenbrock problem **Rvmmin** comes up with the best
    accuracy.

-   **BFGS**, **L-BFGS-B** and **Rcgmin** are obtaining good results
    with less number of computaions envolved.

#### Conlusion

-   For purpose of achieving better results we can use **L-BFGS-B**,
    **Rvmmin**, **BFGS** and **Rcgmin** by providing analytic gradient
    to these functions. These methods performs good even for
    large/difficult problems.
