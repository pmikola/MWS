
rm(list=ls())
library(pracma)
library(ramify)
# Rzucona pinezka upada ostrzem do dolu lub do góry.
# a) Zaproponuj rozklad a priori prawdopodobienstwa p tego, ze pinezka upadnie ostrzem
# do góry.
set.seed(13)
pRand = runif(1,0,1)
Name_of_plot_1 = "A priori of pushpin falling\n with the blade up\n - Normal Distribution"
dist_norm = function(x,mean,sd)
{
  ifelse(x>=0 && x<= 1,(dnorm(x,mean = mean,sd = sd)),0)
}

vdist = Vectorize(dist_norm)

dx = 0.001
x = seq(0,1,dx)
mean = 0.22
sd = 0.08
apriori = vdist(x,mean,sd)
imax = which.max(apriori)
# apriori plot of subjective probability assumption of pushpin falling blade up
plot(x, apriori,type="l",col="red",main= Name_of_plot_1,xlab='p',ylab='f(p)')
GaussMax = x[imax]
points(GaussMax, max(apriori), pch = 19, col = "red")
abline(v = GaussMax, col = "red")
grid(lwd = 1)

# b) Rzuæ pinezki 20 razy (zanotuj wyniki kolejnych rzutów) i na tej podstawie wyznacz
# rozk³ad a posteriori parametru p oraz bayesowski estymator ^p.

ThrowNumb = 20
# 20 throws of pushpin - 0 value means pushpin falling with the blade up
simpleThrowResult = sample(c(0,1), replace=TRUE, size=ThrowNumb)
binomialThrowResult = rbinom(ThrowNumb,1,pRand)
ThResTable = table(binomialThrowResult)
NumOfSucces = as.integer(ThResTable[names(ThResTable) == 0])


likelihood = function(success,trials,prob)
{
  # probability mass function for the binomial distribution
  eq1 = factorial(trials)/(factorial(success)*factorial((trials-success)))
  eq2 = (prob^success)*(1-prob)^(trials-success)
  LL = eq1*eq2
  #LL = dbinom(success,trials,prob)
  return(LL)
}

#Likelihood function and plot
Lik = likelihood(NumOfSucces,ThrowNumb,x)
imax = which.max(Lik)
plot(x, Lik, type = "l",main="Likelihood Estimation", xlab = "p", ylab = "Likelihood")
LL = x[imax]
points(LL, max(Lik), pch = 19, col = "blue")
abline(v = LL, col = "blue")
grid(lwd = 1)

# A Posteriori
aposteriori = Lik*apriori/trapz(Lik*apriori)
plot(x,aposteriori*length(x),type = "h",main="Apriori / Aposteriori", xlab = "p", ylab = "f(p)",ylim=c(0, max(aposteriori*length(x))*1.2),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
grid(lwd=1)
# Alpha and beta parameters estiamtion
apo_mean = mean(aposteriori)
apo_var = var(aposteriori)

alphaBeta = function(mean,var)
{
  alpha = mean**2 / var
  beta = 1/(mean / var)
  return(params = list(alpha = alpha, beta = beta))
}

alpha1 = alphaBeta(apo_mean,apo_var)$alpha
beta1 = alphaBeta(apo_mean,apo_mean)$beta

#Bayes Estimator
bayesEstim1 = alpha1/(alpha1+beta1)
print("bayesestim 1")
print(bayesEstim1)

#c) Rzuc pinezk! jeszcze 20 razy i zanotuj wyniki. Wyznacz rozklad a posteriori oparty
#na wszystkich 40 rzutach i porównaj go z rozkladem uzyskanym po pierwszych 20
#rzutach.

ThrowNumb = 100
# 20 throws of pushpin - 0 value means pushpin falling with the blade up
simpleThrows = simpleThrowResult
simpleThrowResult = c(simpleThrows,sample(c(0,1), replace=TRUE, size=ThrowNumb))
binomialThrow = binomialThrowResult
binomialThrowResult = c(binomialThrow,rbinom(ThrowNumb,1,pRand))
ThResTable = table(binomialThrowResult)
NumOfSucces = as.integer(ThResTable[names(ThResTable) == 0])

#Likelihood function and plot
Lik = likelihood(NumOfSucces,ThrowNumb,x)

# A Posteriori
aposteriori = Lik*apriori/trapz(Lik*apriori)
lines(x,aposteriori*length(x),type = "h", xlab = "p", ylab = "f(p)",col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
lines(x, apriori,type="l",col="blue",xlab='p',ylab='f(p)')
legend( x= "topright", y=0.3, 
        legend=c("aposteriori 20","aposteriori 40", "apriori"), 
        col=c("black", "red", "blue"),
        lty=c(1,1,1))
grid(lwd=1)

# Likelihood plot
imax = which.max(Lik)
plot(x, Lik, type = "l",main="Likelihood Estimation", xlab = "p", ylab = "Likelihood")
LL = x[imax]
points(LL, max(Lik), pch = 19, col = "blue")
abline(v = LL, col = "blue")
grid(lwd = 1)

# Alpha and beta parameters estiamtion
apo_mean = mean(aposteriori)
apo_var = var(aposteriori)

alpha2 = alphaBeta(apo_mean,apo_var)$alpha
beta2 = alphaBeta(apo_mean,apo_mean)$beta

#Bayes Estimator
bayesEstim2 = alpha2/(alpha2+beta2)
print("bayesestim 2")
print(bayesEstim2)
#Z2
#a) Narysuj fgp a priori i wyznacz jej parametry takie jak warto±c oczekiwana, moda,
#odchylenie standardowe; zinterpretuj.

Alpha = 2
Beta = 70
device_trial = 100
broken = 0
x = seq(0, 1, by = dx) 
beta_dist = dbeta(x, 2, 70)
plot(x,beta_dist,type = "l",main="A priori of production line\nalpha = 2\nbeta=70", xlab = "p", ylab = "f(p)")
grid(lwd = 1)

getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Evaluation of beta parameters
beta_mode = getmode(beta_dist)
print("mode 2a")
print(beta_mode)
stdVar = sd(beta_dist)
print("sd 2a")
print(stdVar)
exp_value = mean(beta_dist)
print("mean 2a")
print(exp_value)

#b) Wyznacz analitycznie estymator najwiekszej wiarygodnosci ^thetaML (Maximum Likeli-
# hood) i oblicz jego wartosc dla podanych n; x.
#c) Wyznacz fgp a posteriori i estymator bayesowski ^thetaMAP (MAP = Maximum A Poste-
# riori). W interpretacji odnieœ siê do ^thetaML.

# MaxLikelihood , Likelihood , aposteriori and BayesEstimator
likelihood = likelihood(broken,device_trial,x)
apriori = dbeta(x, 2, 70)
wlikelihood = likelihood*apriori
normconst = trapz(wlikelihood) # equivalent => sum(wlikelihood)
posterior = wlikelihood/normconst
# Plot Data
par(mai=c(1,1,1,1))
plot(x,apriori,ylim=c(0,65),type="l",lwd=1,lty=2,col="blue",main = "FGP / Likelihood",ylab="f(p)",xlab="p")
lines(x,posterior*length(x),type="l",col="blue",lwd=2,lty=1)
lines(x,likelihood*5,type="l",col="red",lwd=1,lty=2,pch=2)
axis(4,at=seq(0,64,by=0.8),labels = seq(0,64,by=0.8))
mtext("Likelihood (Scaled)", side=4, col="red",line=3)
legend( x= "topright", y=0.3, 
        legend=c("aposteriori","apriori", "likelihood"), 
        col=c("blue", "blue", "red"),
        lty=c(1,2,2))
grid(lwd = 1)

# WITHOUT REDEFINITION OF LIKELIHOOD FUN RStudio DOSNT SEE THIS FuNCTION - BUG??
###################################################################################
likelihood = function(success,trials,prob)
{
  # probability mass function for the binomial distribution
  eq1 = factorial(trials)/(factorial(success)*factorial((trials-success)))
  eq2 = (prob^success)*(1-prob)^(trials-success)
  LL = eq1*eq2
  #LL = dbinom(success,trials,prob)
  return(LL)
}
###################################################################################
# Maxlikelihood estimation function
Maxlikelihood = function(prob, trials, success)
{
  Maxlikelihood = sum(likelihood(success,trials,prob))
  #Maxlikelihood = sum(dbinom(x = success, prob = prob, size = trials, log = F))
  return(Maxlikelihood)
}

MaxLik = sapply(x, FUN = Maxlikelihood, trials = device_trial, success = broken)
# set query graphical parameters
par(las = 1, cex.lab = 1.2)
#Likelihood plot
plot(x, MaxLik, type = "l",main="MaxLikelihood Estimation", xlab = "p", ylab = "MaxLikelihood")
# maximum value index location
imax = which.max(MaxLik)
MLE = x[imax]
# Maximum point and vertical check mark of MLE
points(MLE, max(MaxLik), pch = 19, col = "red")
abline(v = MLE, col = "red")
grid(lwd = 1)

# Alpha and beta parameters estiamtion
apo_mean = mean(posterior)
apo_var = var(posterior)

alpha3 = alphaBeta(apo_mean,apo_var)$alpha
beta3 = alphaBeta(apo_mean,apo_mean)$beta

#Bayes Estimator
bayesEstim3 = alpha3/(alpha3+beta3)
print("bayesEstim 2c")
print(bayesEstim3)

###########################################################################################
#d) Zmien parametry rozkladu a priori na alfa = 1 beta = 1 i powtórz obliczenia przywiazujac
#szczególna uwage do interpretacji.
#############################################################################################
# WITHOUT REDEFINITION OF LIKELIHOOD FUN RStudio DOSNT SEE THIS FuNCTION - BUG??
###################################################################################
likelihood = function(success,trials,prob)
{
  # probability mass function for the binomial distribution
  eq1 = factorial(trials)/(factorial(success)*factorial((trials-success)))
  eq2 = (prob^success)*(1-prob)^(trials-success)
  LL = eq1*eq2
  #LL = dbinom(success,trials,prob)
  return(LL)
}
###################################################################################
Alpha = 1
Beta = 1
device_trial = 100
broken = 0
beta_distf = dbeta(x, 1, 1)
plot(x,beta_distf,type = "l",main="A priori of production line \nalpha = 1\nbeta=1", xlab = "p", ylab = "f(p)")
grid(lwd = 1)

# Evaluation of beta parameters
beta_modef = getmode(beta_distf)
print("mode 2d")
print(beta_modef)
stdVar = sd(beta_distf)
print("sd 2d")
print(stdVar)
exp_valuef = mean(beta_distf)
print("mean 2d")
print(exp_valuef)

#bd) Wyznacz analitycznie estymator najwi¦kszej wiarygodno±ci ^thetaML (Maximum Likeli-
# hood) i oblicz jego warto±c dla podanych n; x.
#cd) Wyznacz fgp a posteriori i estymator bayesowski ^thetaMAP (MAP = Maximum A Poste-
# riori). W interpretacji odnie± siê do ^thetaML.

# MaxLikelihood , Likelihood , aposteriori and BayesEstimator
likelihoodf = likelihood(broken,device_trial,x)
apriorif = dbeta(x, 1, 1)
wlikelihoodf = likelihoodf*apriorif
normconstf = trapz(wlikelihoodf) # equivalent => sum(wlikelihood)
posteriorf = wlikelihoodf/normconstf
# Plot Data
par(mai=c(1,1,1,1))
plot(x,apriorif,ylim=c(0,65),type="l",lwd=1,lty=2,col="blue",main = "FGP / Likelihood",ylab="f(p)",xlab="p")
lines(x,posteriorf*length(x),type="l",col="blue",lwd=2,lty=1)
lines(x,likelihoodf*5,type="l",col="red",lwd=1,lty=2,pch=2)
axis(4,at=seq(0,64,by=1.6),labels = seq(0,64,by=1.6))
mtext("Likelihood (Scaled)", side=4, col="red",line=3)
legend( x= "topright", y=0.3, 
        legend=c("aposteriori","apriori", "likelihood"), 
        col=c("blue", "blue", "red"),
        lty=c(1,2,2))
grid(lwd = 1)

MaxLikf = sapply(x, FUN = Maxlikelihood, trials = device_trial, success = broken)
# set query graphical parameters
par(las = 1, cex.lab = 1.2)
#Likelihood plot
plot(x, MaxLikf, type = "l",main="MaxLikelihood Estimation", xlab = "p", ylab = "MaxLikelihood")
# maximum value index location
imaxf = which.max(MaxLikf)
MLEf = x[imaxf]
# Maximum point and vertical check mark of MLE
points(MLEf, max(MaxLikf), pch = 19, col = "red")
abline(v = MLEf, col = "red")
grid(lwd = 1)

# Alpha and beta parameters estiamtion
apo_meanf = mean(posteriorf)
apo_varf = var(posteriorf)

alphaf = alphaBeta(apo_meanf,apo_varf)$alpha
betaf = alphaBeta(apo_meanf,apo_meanf)$beta

#Bayes Estimator
bayesEstimf = alphaf/(alphaf+betaf)
print("bayesEstim 2d")
print(bayesEstimf)


##########################################################################################
# Z3
# Zalóz, ze czas oczekiwania na obsluge w pewnej kolejce jest modelowany rozkladem wy-
# kladniczym z nieznanym parametrem lambda. Rozwaz nastepujacy rozklad a priori parametru lambda:
# rozklad gamma ze srednia 0.4 i wariancja 0.2. Wyznacz numerycznie (przy pomocy reguly
# Bayesa, bez wykorzystywania rozkladów sprzezonych) i narysuj funkcje gestosci rozkladu a
# posteriori uzyskanego po zaobserwowaniu, ze sredni czas oczekiwania w rozwazanej kolejce,
# wyliczony dla losowo wybranych 20 osób, wynosi 5.1 minuty.

##########################################################################################
#Const
gamma_mean = 0.4
gamma_var = 0.2
m20Time = 5.1 # [h] -> beta
peopleNum = 20 # number of events that we are waiting in avg time 5.1min -> alpha
avg4oneCustomer = m20Time/peopleNum
lambda=seq(0,1,dx)
#Gamma apriori distribution
gm_alpha = alphaBeta(gamma_mean,gamma_var)$alpha
gm_beta = alphaBeta(gamma_mean,gamma_var)$beta
print("ALPHA for Gamma Distribution")
print(gm_alpha)
print("BETA for Gamma Distribution")
print(gm_beta)
apriori = dgamma(lambda, shape=gm_alpha, scale=gm_beta)
plot(lambda, apriori,  xlab="lambda", ylab="f(lambda)",type="l",xlim=c(0,1),ylim=c(0,10), main="Apriori / Aposteriori",col="blue")

expDist = function(lambda)
{
  # valid forms
  # expfun = (lambda * exp(-lambda*m20Time))^peopleNum
  # expfun = (lambda^peopleNum)*exp(-lambda*peopleNum*m20Time)
  expfun = dexp(m20Time,lambda)^peopleNum
  return(expfun)
}
estimVal=function(x)
{
  output = expDist(x)*dgamma(x,shape=gm_alpha, scale=gm_beta)
  return(output)
}

# gamma/exp posteriori distribution
estimValInteg = integrate(estimVal,0,Inf)
scale = as.numeric(estimValInteg[1])
aposteriori = expDist(lambda)*apriori/scale
lines(lambda, aposteriori, type="l", col="red",lty=2)
legend("topright", legend = c("apriori","aposteriori"),col=c("blue", "red"),lty=1:2)
grid(lwd = 1)

postEstim=function(x)
{
  estimVal(x)/estimValInteg$value
}

est=function(x)
{
  x*postEstim(x)
}

lambdaEstim=integrate(est, lower=0, upper=Inf)
print("Lambda Estimator")
print(lambdaEstim)
meanexp=1/lambdaEstim$value
print("Mean Value")
print(meanexp)








