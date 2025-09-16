### FourParModel.R
### R version: R 4.1.0
### R code for estimation of a 4-parameter monotonic time-varying effect model for analysis of recurrent events after
### a single dose or multiple doses of a treatment/exposure 
### Xu et al (2017) Semiparametric estimation of time-varying intervention effects using recurrent event data.
### Statistics in Medicine 2017;36(17):2682-2696. 
### Read the comments below and Example4ParV3.csv for the data format expected

rm(list=ls())

### specify the Workspace
setwd("C:/Users/gmsmxm/Documents/nBox/Stat/Nonlinear/Codes")

### load functions used in the analysis
source("Functions4ParG1.R")

### Data format: one person has multiple rows with the same "id".
### Each record ends with "event" = 1 if the outcome event occurs; otherwise "event" = 0 (censored).
### "treat" is the binary tretment/exposure variable whose effect is to be represented by the 4-parameter function.
### "doe" is the start date of follow-up, in yyyy-mm-dd format.
### "startT" and "endT" are time in months since doe; they define the start and end of each spell of follow-up/at-risk time.
### The data file must include the variables "id", "treat", "startT", "endT" and "event". 
### Covariates can be time-varying (e.g. current age)
### d1, d2, ... are variable names for timing of doses; time is months since doe
### doseX=NULL means a single dose at time "0"; equivalently, specify doseX=c("d1") where d1=0
### Missing doses are represented by a blank cell in the csv file

### Read example data from a csv file
df_E = read.csv("Example4Par.csv",header=TRUE)
df_E[1:10,]

### 4-parameter function g(t)=a*exp(-b*t^c)+d with b,c>0 for representing the treatment effect at time t since time 0
### if there is multiple doses, then G(t)=g(t-d1)+g(t-d2)+g(t-d3)+...
### if someone hadn't received the scheduled dose2, i.e. d2=NA, then g(t-d2)=0 for t>d2
gfun <- function(t,a,b,c,d)  
	 ifelse(is.na(t),0,(a*exp(-b*t^c)+d)*(t>0))

### specify the date of doses of intervention
### E.g., in the example file, "d1","d2","d3","d4"
# single dose
#doseX=NULL
# four doses
doseX=c("d1","d2","d3","d4")

### specify the covariates included in the regression
### time-contant covariate: "gender", "village"
### time-varying covariate: "agegr"
### data has been splitted at the time point where the time-varying covariate changed
X=c("gender", "village","agegr")

#Convert "village" to a factor, as "village" is a categorical covariate 
df_E$village=factor(df_E$village)
levels(df_E$village) 
#here 4 levels,"village1" is the reference group 

#Finally, there are 5 estimators of coefficients for covariates
Covariates=c("gender","village2","village3","village4","agegr")

### Those having incomplete covariates information are excluded in the analysis
df_E = df_E[complete.cases(df_E[,X]),]

### After the estimation, we need to specify the choice of SE
### Robust SE or Naive SE

#----------------------------- The Model  ----------------------------------#

### specify the recurrent event data 
dic.pd <- dicp.dose(data=df_E,doseX=doseX)

### optionally, specify initial values in optimization
a0 <- 0
ln_b0 <- 0
ln_c0 <- 0
d0 <- 0
q <- length(Covariates)
par0 <- c(a0, ln_b0, ln_c0, d0, rep(0,length(Covariates)))

### use the "optim" command to minimize the negative log-likelihood function
fit.pd <- optim(par0, likpd.dose
                ,fX=X
                ,data=df_E
                ,diclst=dic.pd
                ,fun=gfun
                ,method="BFGS"
                ,hessian=TRUE
		    ,control=list(maxit=1000, trace=TRUE))
fit.pd$par;
### compute loglikelihood value 
lik = -likpd.dose(fit.pd$par,fX=X,data=df_E,diclst=dic.pd,fun=gfun)
lik

### show AIC;BIC 
-2*lik+2*length(fit.pd$par)
-2*lik+log(nrow(df_E))*length(fit.pd$par)

### save estimation results
### constraints: b and c are non-negative values
a <- fit.pd$par[1]
b <- exp(fit.pd$par[2])
c <- exp(fit.pd$par[3])
d <- fit.pd$par[4]
H <- fit.pd$hessian

### summary results
coef.pd <- c(a,b,c,d,fit.pd$par[-c(1:4)])

#----------------- Calculate standard errors (SE) ---------------------#

### choose 1 of 2 types of standard errors
### Robust SE with allowing for clustering of observations/events within person or Naive SE

clusterID=unique(df_E$id)
SEs <- Robust.error(fit.pd$par,Hess=H,clusterID=clusterID,data=df_E,fX=X,doseX=doseX,fun=gfun)

## Robust SE allowing for multiple events within cluster/person
error <- sqrt(diag(SEs$ClusterH))

## Naive SE
# error <- sqrt(diag(solve(H)))

### delta method are used for computing SEs of (b,c) 
se.pd <- numeric(length(fit.pd$par))
se.pd[2:3] <- c(b,c)* error[2:3]
se.pd[-c(2:3)] <- error[-c(2:3)]

z.pd <- coef.pd/se.pd
tab.pd <- cbind(coef.pd, se.pd, z.pd)

### labels the outputs
rownames(tab.pd)= c("a","b","c","d",Covariates)

### show results
round(tab.pd,4)

