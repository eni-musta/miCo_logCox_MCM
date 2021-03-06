# miCo_logCox_MCM

The logistic/Cox mixture cure model is commonly used to handle survival data in the presence of a cure fraction. It assumes that the population consists of two subpopulations: 
the cured and the noncured ones. Logistic regression is used to model the incidence and a Cox proportional hazards model is used for the latency. Maximum likelihood estimators 
are computed via the Expectation-Maximization algorithm and are implemented in the package smcure (Cai et al. (2012)). Recently, a new two-step estimation method for the 
logistic/Cox  mixture cure model by means of presmoothing was proposed by Musta, Patilea and Van Keilegom (2020). 
However, the presence of mismeasured covariates often results in biased estimators and incorrect inference. Here we propose to correct the estimators for the bias introduced 
by the measurement error by a simex algorithm, which is a very general simulation based method. It essentially estimates this bias by introducing additional error to the data
and then recovers bias corrected estimators through an extrapolation approach. Any estimation method in the absence of measurement error can be used within the simex algorithm.
Here we illustrate the method for the MLE and the presmoothing estimator in the logistic/Cox model but it can be applied to other parametric/semiparametric mixture cure models 
as well requiring only an estimation method in the absence of measurement error. 

We demonstrate the method for prostate cancer data from the Surveillance, Epidemiology and End Results (SEER) database, which is a collection of cancer incidence data from 
population based cancer registries in the US. This dataset ("Data_prostate.txt") consists of 726 observations and the event time is death because of prostate cancer. The 
covariates of interest are the PSA level (prostate-specific antigen) and stage at diagnosis (localized, regional or distant). PSA is often used as a tool to monitor the 
development of the disease but PSA measurements are not error-free because of the inaccuracy of the measuring technique and own fluctuations of the PSA levels. The purpose of 
this study was to analyse the effect of PSA on cure probability and survival while accounting for measurement error and to check how sensitive the results are to the 
measurement error. 

A detailed description of the method, its practical and theoretical behavior and the data application can be found in Musta and Van Keilegom (2021).

One can run the example via 'miCo_logCox-MCM.R'. This produces estimators for the parameters of the incidence, the regression coefficients and the baseline survival of the Cox 
component (latency). The bootstrap procedure is used to estimate the variance of the parameter estimates and the resulting p-values. The R file uses the packages 'smcure' 
for the data and 'np' for bandwidth selection. Additional functions are defined in the file 'functions_miCO_logCox-MCM.R'.




References

Cai, C., Y. Zou, Y. Peng, and J. Zhang (2012). smcure: An r-package for estimating semiparametric mixture cure models. Comput. Meth. Prog. Bio. 108 (3), 1255-1260.

Musta, E., V. Patilea and I. Van Keilegom (2020). A presmoothing approach for estimation in mixture cure models. https://arxiv.org/abs/2008.05338

Musta, E. and I. Van Keilegom (2021). A simulation-extrapolation approach for the mixture cure model with mismeasured covariates. Electronic Journal of Statistics 
Vol. 15 (2021), pp. 3708???3742.
