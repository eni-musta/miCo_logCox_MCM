##########################################################################################
# A simulation-extrapolation approach for the mixture cure model with measurement error  #
##########################################################################################

# logistic/Cox model
# Incidence:  P('cured'|X)=1-phi(gamma_0,X),   phi(gamma,X)=1/(1+exp(-gamma*X)),
#             where gamma_0 is the vector of unknown parameters and X is the vector of covariates
# Latency:  for the uncured subjects, S_u(t|Z)=S_0(t)^exp(beta*Z),
#           where beta_0 is a vectorof unknown parameters and S_0 is the baseline survival function


# Libraries

library(smcure)
library(np)

# Predefined functions

source("functions_miCo_logCox-MCM.R")


# read the prostate-cancer data

Data_prost=read.table("Data_prostate.txt",header = FALSE)
n=dim(Data_prost)[1]    # sample size
Y=Data_prost[,1]        # follow-up time
status=Data_prost[,2]   # censoring status:  censored (status=0), noncensored (status=1)
X=cbind(rep(1,n),Data_prost[,c(3:5)])  # covariates for the incidence: PSA (centered to the mean), Indicator of 'distant' stage, Indicator of 'Regional' stage
Z=Data_prost[,c(3:5)]   # covariates for the incidence: PSA (centered to the mean), Indicator of 'distant' stage, Indicator of 'Regional' stage


# Plot the Kaplan Meier estimator of the survival function

plot(survfit(Surv(Y,status)~1),xlab = 'Time (years)',ylab = 'Survival probability',main='Kaplan-Meier estimator')


################################
# Compute the naive estimators #
################################


# maximum likelihood estimator ('smcure')

Data=data.frame(Y,status,X[,1],X[,2],X[,3],X[,4],Z[,1],Z[,2],Z[,3])
colnames(Data)=c("Y","status","X_1","X_2","X_3","X_4","Z_1","Z_2","Z_3")

MCM_naive=smcure(Surv(Y,status)~Z_1+Z_2+Z_3,cureform=~X_2+X_3+X_4,data=Data, model="ph",Var = FALSE)

gamma_MLE_naive=MCM_naive$b
beta_MLE_naive=MCM_naive$beta


# estimator based on presmoothing

# data is a dataframe containing the follow-up times, censoring status and the covariates used to model the incidence
# In this case all covariates are used for both incidence and latency 
# The continuous covariates are standardized and come before the discrete covariates (as columns of data_log)
# The observations (rows) are ordered according to the follow-up time.

data = data.frame(cbind(Y=Y,Delta=status,X1=X[,2]/sd(X[,2]),X2=X[,3],X3=X[,4]))
ord=order(data[,1],1-data[,2])
data = data[ord,]
g_init = glm(data[,2]~data[,3]+data[,4]+data[,5],family=quasibinomial)$coefficients  # initail value of beta computed by logistic regression

Y_r=max(Y[which(status==1)]) # the largest observed event time
hopt=npcdistbw(formula=Y~X1,data=data,gydat=data[which(data[,1]<=Y_r),1],cxkertype="epanechnikov")$xbw  # compute bandwidth for estimation of  H(t|X)
cmax = 2
hopt=min(hopt,2)  # truncate bandwidth from above


w = weight(data[,3:5],hopt)   # compute the Nadaray-Watson weights for the Beran estimator
pi_hat = beran(data,w)        # nonparametric estimator of cure probabilities
g = optim(par=g_init,logLik_logit,V=cbind(pi_hat,data[,3:5]),method='BFGS',control=list(fnscale=-1,maxit=50))$par  #estimator using logistic-likelihood 
gamma_presmooth_naive = g/c(1,sd(X[,2]),1,1)  # rescale the estimator so that it corresponds to the initial covariates (not standardized ones)


eps = 1e-07     # tolerance criteria for the convergence of the EM algorithm
emmax = 50      # maximum number of iterations for the EM algorithm
beta_init = coxph(Surv(Y, status) ~ Z_1+Z_2+Z_3,data=Data, method = "breslow")$coef # initial estimator for beta as in the classical Cox model without cure

EM_fit_naive = EM(Y, status, X, as.matrix(Z), gamma_presmooth_naive, beta_init, emmax, eps)   # perform the EM algorithm
beta_presmooth_naive = EM_fit_naive$beta



##################
#SIMEX algorithm #
##################


v = 8  # standard deviation of the measurement error


B = 50                     # number of samples for each level of added noise
K = 5                      # number of levels of added naise
Nu = c(0,0.5,1,1.5,2)      # vector of levels of added noise

est_MLE = matrix(0,K,7)           # MLE estimators of beta and gamma for each value of added noise
est_presmooth = matrix(0,K,7)     # estimators of gamma for each value of lambda

# the naive estimators correspond to the '0' level of added noise
est_MLE[1,] = c(gamma_MLE_naive,beta_MLE_naive)
est_presmooth[1,] = c(gamma_presmooth_naive,beta_presmooth_naive)


for (k in 2:K) # for each level of added noise
{
    est_MLE_b = matrix(0,B,7)
    est_presmooth_b = matrix(0,B,7)

    for (b in 1:B) # for each simulation (in one level of added noise)
    {
      #add aditional noise
      U = rnorm(n,0,1)
      X_e = X
      Z_e = Z
      Z_e[,1] = Z_e[,1]+sqrt(Nu[k]*v)*U
      X_e[,2] = X_e[,2]+sqrt(Nu[k]*v)*U

      #estimate ignoring the error

      #smcure

      Data_e = data.frame(Y,status,X_e[,1],X_e[,2],X_e[,3],X_e[,4],Z_e[,1],Z_e[,2],Z_e[,3])
      colnames(Data_e) = c("Y","status","X_1","X_2","X_3","X_4","Z_1","Z_2","Z_3")

      MCM_e = smcure(Surv(Y,status)~Z_1+Z_2+Z_3,cureform=~X_2+X_3+X_4,data=Data_e, model="ph",Var = FALSE)

      est_MLE_b[b,] = c(as.numeric(MCM_e$b),as.numeric(MCM_e$beta))
      

      # estimator based on presmoothing

      data_e = data.frame(cbind(Y=Y,Delta=status,X1=(X_e[,2]-mean(X_e[,2]))/sd(X_e[,2]),X2=X_e[,3],X3=X_e[,4]))
      ord = order(data_e[,1],1-data_e[,2])
      data_e = data_e[ord,]
      g_init = glm(data_e[,2]~data_e[,3]+data_e[,4]+data_e[,5],family=quasibinomial)$coefficients  
     
      w = weight(data_e[,3:5],hopt)  # B.nk(z)
      pi_hat = beran(data_e,w)  # P(T=infinit|x.j) j=1,...,n
      g_e= optim(par=g_init,logLik_logit,V=cbind(pi_hat,data_e[,3:5]),method='BFGS',control=list(fnscale=-1,maxit=15))$par  #estimator using logistic-likelihood, initial value for MLE in the general approach
      gamma_presmooth_e=(g_e-c(g_e[2]*mean(X_e[,2])/sd(X_e[,2]),0,0,0))/c(1,sd(X_e[,2]),1,1)

      beta_init = coxph(Surv(Y, status) ~ Z_1+Z_2+Z_3,data=Data, method = "breslow")$coef # initial estimator for beta as in the classical Cox model without cure
      
      EM_fit_e = EM(Y, status, X_e, as.matrix(Z_e), gamma_presmooth_e, beta_init, emmax, eps)   # perform the EM algorithm
      beta_presmooth_e = EM_fit_e$beta
      
      est_presmooth_b[b,] = c(gamma_presmooth_e,beta_presmooth_e)

    }

    est_MLE[k,] = apply(est_MLE_b, 2, mean)
    est_presmooth[k,] = apply(est_presmooth_b, 2, mean)
  }


# simex maximum likelihood estimator
quad_extrap1 = SimexExtrapolation(estimNui=est_MLE,variableNames=c("gamma_1","gamma_2","gamma_3","gamma_4","beta_1","beta_2","beta_3"),nbParam=7,Nu,orderExtrap=2)
gamma_MLE_simex = quad_extrap1$estim[1:4]
beta_MLE_simex = quad_extrap1$estim[5:7]

# simex estimators based on presmoothing
quad_extrap2 = SimexExtrapolation(estimNui=est_presmooth,variableNames=c("gamma_1","gamma_2","gamma_3","gamma_4","beta_1","beta_2","beta_3"),nbParam=7,Nu,orderExtrap=2)
gamma_presmooth_simex = quad_extrap2$estim[1:4]
beta_presmooth_simex = quad_extrap2$estim[5:7]



#######################################
# Estimate the variance via bootstrap # 
#######################################


nboot = 1000        # number of bootstrap samples
nbeta = dim(Z)[2]  # number of covariates for the latency
ngamma = dim(X)[2] # number of covariates for the incidence

gamma_MLE_naive_boot = matrix(rep(0, nboot * ngamma), nrow = nboot)
beta_MLE_naive_boot = matrix(rep(0, nboot * nbeta), nrow = nboot)
gamma_presmooth_naive_boot = matrix(rep(0, nboot * ngamma), nrow = nboot)
beta_presmooth_naive_boot = matrix(rep(0, nboot * nbeta), nrow = nboot)
gamma_MLE_simex_boot = matrix(rep(0, nboot * ngamma), nrow = nboot)
beta_MLE_simex_boot = matrix(rep(0, nboot * nbeta), nrow = nboot)
gamma_presmooth_simex_boot = matrix(rep(0, nboot * ngamma), nrow = nboot)
beta_presmooth_simex_boot = matrix(rep(0, nboot * nbeta), nrow = nboot)

tempdata = cbind(Y, status, X, Z)
data1 = subset(tempdata, status == 1)
data0 = subset(tempdata, status == 0)
n1 = nrow(data1)
n0 = nrow(data0)

eps=1e-7 # convergence criteria for the EM algorithm  

B=50
K=5
Nu=c(0,0.5,1,1.5,2)      
est_MLE_boot=matrix(0,K,7)     
est_presmooth_boot=matrix(0,K,7)     



  i = 1
  while (i <= nboot) {
    
    # generate a bootstrap sample
        
      id1 = sample(1:n1, n1, replace = TRUE)
      id0 = sample(1:n0, n0, replace = TRUE)
      bootdata = rbind(data1[id1, ], data0[id0, ])
      bootY= bootdata[,1]
      boot_status=bootdata[,2]
      bootX = bootdata[,3:(2+ngamma)]
      bootZ = bootdata[,(3+ngamma):(2+ngamma+nbeta)]
      
    # naive  estimators  
      
    # maximum likelihood
    gamma_init = gamma_MLE_naive
    beta_init = beta_MLE_naive
    em_fit = em(bootY, boot_status, as.matrix(bootX), as.matrix(bootZ), gamma_init, beta_init, emmax, eps)  
    est_MLE_boot[1,]=c(as.numeric(em_fit$gamma),as.numeric(em_fit$beta))
    tau=em_fit$tau
    
    # estimator based on presmoothing
  
    data_boot = data.frame(cbind(Y=bootY,Delta=boot_status,X1=(bootX[,2]-mean(bootX[,2]))/sd(bootX[,2]),X2=bootX[,3],X3=bootX[,4]))
    ord=order(data_boot[,1],1-data_boot[,2])
    data_boot = data_boot[ord,]
    g_init = glm(data_boot[,2]~data_boot[,3]+data_boot[,4]+data_boot[,5],family=quasibinomial)$coefficients  
    
    Y_r=max(bootY[which(boot_status==1)]) 
    hopt=npcdistbw(formula=Y~X1,data=data_boot,gydat=data_boot[which(data_boot[,1]<=Y_r),1],cxkertype="epanechnikov")$xbw  
    cmax = 2
    hopt=min(hopt,2)  
    
    w = weight(data_boot[,3:5],hopt)  
    pi_hat = beran(data_boot,w)  
    g = optim(par=g_init,logLik_logit,V=cbind(pi_hat,data_boot[,3:5]),method='BFGS',control=list(fnscale=-1,maxit=50))$par 
  
    gamma_presmooth_naive_boot[i,]=(g-c(g[2]*mean(bootX[,2])/sd(bootX[,2]),0,0,0))/c(1,sd(bootX[,2]),1,1)
  
    EM_fit_boot = EM(bootY, boot_status, as.matrix(bootX), as.matrix(bootZ), gamma_presmooth_naive_boot[i,], beta_init, emmax, eps)   # perform the EM algorithm
    beta_presmooth_naive_boot[i,] = EM_fit_boot$beta
    tau=max(tau,EM_fit_boot$tau)
    
    est_presmooth_boot[1,]=c(gamma_presmooth_naive_boot[i,],beta_presmooth_naive_boot[i,])
    
    
  if(tau<eps){
      
  #SIMEX algorithm 
  
  for (k in 2:K) 
  {
      est_MLE_b_boot=matrix(0,B,7)
      est_presmooth_b_boot=matrix(0,B,7)
      
      for (b in 1:B) # for each simulation (in one level of added noise)
      { 
        #add aditional noise
        U=rnorm(n,0,1)
        bootX_e=bootX
        bootZ_e=bootZ
        bootZ_e[,1]=bootZ_e[,1]+sqrt(Nu[k]*v)*U
        bootX_e[,2]=bootX_e[,2]+sqrt(Nu[k]*v)*U
        
        # estimate ignoring the error
        
        # maximum likelihood estimator
        Data_e_boot=data.frame(bootY,boot_status,bootX_e[,1],bootX_e[,2],bootX_e[,3],bootX_e[,4],bootZ_e[,1],bootZ_e[,2],bootZ_e[,3])
        colnames(Data_e_boot)=c("Y","status","X_1","X_2","X_3","X_4","Z_1","Z_2","Z_3")
        MCM_e_boot=smcure(Surv(Y,status)~Z_1+Z_2+Z_3,cureform=~X_2+X_3+X_4,data=Data_e_boot, model="ph",Var = FALSE)
        est_MLE_b_boot[b,]=c(as.numeric(MCM_e_boot$b),as.numeric(MCM_e_boot$beta))
        
        # estimator based on presmoothing
        data_e_boot = data.frame(cbind(Y=bootY,Delta=boot_status,X1=(bootX_e[,2]-mean(bootX_e[,2]))/sd(bootX_e[,2]),X2=bootX_e[,3],X3=bootX_e[,4]))
        ord=order(data_e_boot[,1],1-data_e_boot[,2])
        data_e_boot = data_e_boot[ord,]
        g_init = glm(data_e_boot[,2]~data_e_boot[,3]+data_e_boot[,4]+data_e_boot[,5],family=quasibinomial)$coefficients 
        
        w = weight(data_e_boot[,3:5],hopt) 
        pi_hat = beran(data_e_boot,w) 
        g_e= optim(par=g_init,logLik_logit,V=cbind(pi_hat,data_e_boot[,3:5]),method='BFGS',control=list(fnscale=-1,maxit=15))$par  
        gamma_presmooth_boot_e=(g_e-c(g_e[2]*mean(bootX_e[,2])/sd(bootX_e[,2]),0,0,0))/c(1,sd(bootX_e[,2]),1,1)
        
        beta_init = coxph(Surv(Y, status) ~ Z_1+Z_2+Z_3,data=Data_e_boot, method = "breslow")$coef 
        EM_fit_boot_e = EM(bootY, boot_status, as.matrix(bootX_e), as.matrix(bootZ_e), gamma_presmooth_boot_e, beta_init, emmax, eps)   
        beta_presmooth_boot_e = EM_fit_boot_e$beta
        
        est_presmooth_b_boot[b,] = c(gamma_presmooth_boot_e,beta_presmooth_boot_e)
        
      }
      
      est_MLE_boot[k,] = apply(est_MLE_b_boot, 2, mean)
      est_presmooth_boot[k,] = apply(est_presmooth_b_boot, 2, mean)
   }
  
    
    # simex maximum likelihood estimator
    quad_extrap_boot1 = SimexExtrapolation(estimNui=est_MLE_boot,variableNames=c("gamma_1","gamma_2","gamma_3","gamma_4","beta_1","beta_2","beta_3"),nbParam=7,Nu,orderExtrap=2)
    gamma_MLE_simex_boot[i,] = quad_extrap_boot1$estim[1:4]
    beta_MLE_simex_boot[i,] = quad_extrap_boot1$estim[5:7]
    
    # simex estimators based on presmoothing
    quad_extrap_boot2 = SimexExtrapolation(estimNui=est_presmooth_boot,variableNames=c("gamma_1","gamma_2","gamma_3","gamma_4","beta_1","beta_2","beta_3"),nbParam=7,Nu,orderExtrap=2)
    gamma_presmooth_simex_boot[i,] = quad_extrap_boot2$estim[1:4]
    beta_presmooth_simex_boot[i,] = quad_extrap_boot2$estim[5:7]
    
    i=i+1
  }
} 
 
  
  gamma_var_MLE_naive = apply(gamma_MLE_naive_boot, 2, var)
  beta_var_MLE_naive = apply(beta_MLE_naive_boot, 2, var)
  gamma_sd_MLE_naive = sqrt(gamma_var_MLE_naive)
  beta_sd_MLE_naive = sqrt(beta_var_MLE_naive)
  gamma_var_presmooth_naive = apply(gamma_presmooth_naive_boot, 2, var)
  beta_var_presmooth_naive = apply(beta_presmooth_naive_boot, 2, var)
  gamma_sd_presmooth_naive = sqrt(gamma_var_presmooth_naive)
  beta_sd_presmooth_naive = sqrt(beta_var_presmooth_naive)
  
  
  gamma_var_MLE_simex = apply(gamma_MLE_simex_boot, 2, var)
  beta_var_MLE_simex = apply(beta_MLE_simex_boot, 2, var)
  gamma_sd_MLE_simex = sqrt(gamma_var_MLE_simex)
  beta_sd_MLE_simex = sqrt(beta_var_MLE_simex)
  gamma_var_presmooth_simex = apply(gamma_presmooth_simex_boot, 2, var)
  beta_var_presmooth_simex = apply(beta_presmooth_simex_boot, 2, var)
  gamma_sd_presmooth_simex = sqrt(gamma_var_presmooth_simex)
  beta_sd_presmooth_simex = sqrt(beta_var_presmooth_simex)
  
  
  
  log_var=c('Intercept','PSA','Is_Distant','Is_Regional')
  cox_var=c('PSA','Is_Distant','Is_Regional')
  
  
  
  fit_MLE_naive = list()
  fit_MLE_naive$beta = beta_MLE_naive
  fit_MLE_naive$gamma = gamma_MLE_naive
  fit_MLE_naive$gamma_sd = gamma_sd_MLE_naive
  fit_MLE_naive$gamma_zvalue = gamma_MLE_naive/gamma_sd_MLE_naive
  fit_MLE_naive$gamma_pvalue = (1 - pnorm(abs(fit_MLE_naive$gamma_zvalue)))*2
  fit_MLE_naive$beta_sd = beta_sd_MLE_naive
  fit_MLE_naive$beta_zvalue = beta_MLE_naive/beta_sd_MLE_naive
  fit_MLE_naive$beta_pvalue = (1 - pnorm(abs(fit_MLE_naive$beta_zvalue)))*2
  fit_MLE_naive$gammanm =  log_var
  fit_MLE_naive$betanm = cox_var
  print_results(fit_MLE_naive)  
  
  
  fit_MLE_simex = list()
  fit_MLE_simex$beta = beta_MLE_simex
  fit_MLE_simex$gamma = gamma_MLE_simex
  fit_MLE_simex$gamma_sd = gamma_sd_MLE_simex
  fit_MLE_simex$gamma_zvalue = gamma_MLE_simex/gamma_sd_MLE_simex
  fit_MLE_simex$gamma_pvalue = (1 - pnorm(abs(fit_MLE_simex$gamma_zvalue)))*2
  fit_MLE_simex$beta_sd = beta_sd_MLE_simex
  fit_MLE_simex$beta_zvalue = beta_MLE_simex/beta_sd_MLE_simex
  fit_MLE_simex$beta_pvalue = (1 - pnorm(abs(fit_MLE_simex$beta_zvalue)))*2
  fit_MLE_simex$gammanm =  log_var
  fit_MLE_simex$betanm = cox_var
  print_results(fit_MLE_simex)  
  
  
  fit_presmooth_naive = list()
  fit_presmooth_naive$beta = beta_presmooth_naive
  fit_presmooth_naive$gamma = gamma_presmooth_naive
  fit_presmooth_naive$gamma_sd = gamma_sd_presmooth_naive
  fit_presmooth_naive$gamma_zvalue = gamma_presmooth_naive/gamma_sd_presmooth_naive
  fit_presmooth_naive$gamma_pvalue = (1 - pnorm(abs(fit_presmooth_naive$gamma_zvalue)))*2
  fit_presmooth_naive$beta_sd = beta_sd_presmooth_naive
  fit_presmooth_naive$beta_zvalue = beta_presmooth_naive/beta_sd_presmooth_naive
  fit_presmooth_naive$beta_pvalue = (1 - pnorm(abs(fit_presmooth_naive$beta_zvalue)))*2
  fit_presmooth_naive$gammanm =  log_var
  fit_presmooth_naive$betanm = cox_var
  print_results(fit_presmooth_naive)  
  
  
  fit_presmooth_simex = list()
  fit_presmooth_simex$beta = beta_presmooth_simex
  fit_presmooth_simex$gamma = gamma_presmooth_simex
  fit_presmooth_simex$gamma_sd = gamma_sd_presmooth_simex
  fit_presmooth_simex$gamma_zvalue = gamma_presmooth_simex/gamma_sd_presmooth_simex
  fit_presmooth_simex$gamma_pvalue = (1 - pnorm(abs(fit_presmooth_simex$gamma_zvalue)))*2
  fit_presmooth_simex$beta_sd = beta_sd_presmooth_simex
  fit_presmooth_simex$beta_zvalue = beta_presmooth_simex/beta_sd_presmooth_simex
  fit_presmooth_simex$beta_pvalue = (1 - pnorm(abs(fit_presmooth_simex$beta_zvalue)))*2
  fit_presmooth_simex$gammanm =  log_var
  fit_presmooth_simex$betanm = cox_var
  print_results(fit_presmooth_simex)  
  
  
  
  
 