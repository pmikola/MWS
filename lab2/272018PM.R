#library(lubridate)
#library(chron)
library(boot)
library(MASS)
library(simpleboot)
rm(list=ls())
graphics.off()
########################################################
# Z1 
# a) Dopasuj do tych danych rozkaad Poissona (wyestymuj parametr lambda). Sprawd1 zgodno±c
# otrzymanego rozkaadu z zaobserwowanym danymi porównuj!c grañcznie rzeczywiste
# (zaobserwowane) i wynikaj!ce z modelu liczby "skrêtów w prawo".
########################################################

# Read Skrety.txt
RealTurns = read.delim("skrety.txt",header=TRUE,sep="")

# Get length and values from skrety.txt
RT_Values = as.integer(RealTurns$X1)
RT_len = length(RT_Values)
#RT_len = lengths(RealTurns)

# All Right Turns in data
sumRT = sum(RT_Values)
RTsum = cat(sprintf('Number of right turns is %s turns in given data\n', sumRT))

# Deffinition of  PoissonDensity function from data for dpois func
PoissDens = function(data)
{
  dens = dpois(0:max(data), mean(data))
  return(dens)
}

dens = PoissDens(RT_Values)

###################################
#format( seq.POSIXt(as.POSIXct(Sys.Date()), as.POSIXct(Sys.Date()+1), by = "3 min"),
#        "%H%M", tz="Poland")

#z = seq.POSIXt(as.POSIXct(Sys.Date()), as.POSIXct(Sys.Date()+1), by = "3 min")

# Time_Interval = seq(ymd_hm('2015-01-01 00:00'),ymd_hm('2016-12-31 23:45'), by = '3 mins')
###################################

# Plot histogram and Poisson distribution from RTvalues
histogram = hist(RT_Values,10,main= 'Right Turns',xlab='Time [min]',ylab='Number of Turns',xaxt='n',freq = F)
axis(side=1, at=seq(0,13,1), labels=seq(0,897,69))
lines(0:max(RT_Values),dens, col = 'red')
legend("topright", c("Data", "Poisson"), col=c("black", "red"), lwd=2)
grid(lwd = 1)


# b) Metod¹ bootstrapu nieparametrycznego oszacuj odchylenie standard. estymatora lambda z daszkiem

#set.seed(1) # fixed sample

# random pois dist probe
x = rpois(1000,mean(RT_Values))

# f for boot function
meanFunc = function(x,i)
{
  mean(x[i])
}
# bootstrap
RT_Values_boot = boot(data = RT_Values, statistic = meanFunc, R = RT_len)
print(RT_Values_boot)

########################################################
# Z2
# a) Narysuj histogram ww. odstêpów. Czy rozk³ad gamma móg³by byc odpowiednim
# modelem dla zarejestrowanych danych?

# Read fotony.txt
PhotonInterSpace = read.delim("fotony.txt",header=TRUE,sep="")

InterSpace = as.integer(PhotonInterSpace$X33.312)
InterSpace_len = length(InterSpace)
#set.seed(1) # fixed sample
# FREQ -> normalization

histogram = hist(InterSpace,50,main= 'Gamma Interspace',xlab='Interspace',ylab='Number of photons',freq = F)
grid(lwd = 1)
# yes

# b) Metod¹ momentów oraz metod¹ najwiêkszej wiarygodnoœci wyznacz estymatory parametr
#  alfa i beta rozk³adu gamma odpowiadaj¹ce zarejestrowanym danym.

#Generate MiniBatch from given data (Sample Variance)
set.seed(1) # fixed sample
GenerateMiniDataset = function(data,n)
{
  data_len = length(data)
  new_data_set = c(1:n)
  for (i in 1:n) {
    rand_index = as.integer(runif(1,1,data_len))
    new_data_set[i] = data[rand_index] 
  }
  return(new_data_set)
}

# Methods of moments for Sample Variance, not population
MethodOfMoments = function(data,n)
{
  new_data_set = GenerateMiniDataset(data,n)
  nds_len = length(new_data_set)
  mean = 1/nds_len * sum(new_data_set)
  variance = (1/(nds_len-1)) * sum((new_data_set-mean)**2)
  standard_deviation = sqrt(variance)
  return_list = list("mean" = mean
                     ,"variance" = variance
                     ,"standard_deviation" = standard_deviation)
  return(return_list)
}

#print(MethodOfMoments(InterSpace,10000))

EstimatorNW = function(sample,alphaCoef,betaCoef)
{
  x = rgamma(sample,shape = alphaCoef,rate = betaCoef)
  fit = fitdistr(x, "gamma", start=list(shape=alphaCoef, scale=1/betaCoef),lower = 0.0001)$estimate
  return(fit)
}

# Gamma Coeficients for fit and comparision from population
GammaCoefficients = function(data)
{
  alpha = (mean(data)**2)/var(data)
  beta = sqrt(alpha/var(data))
  return_list = list("alpha" = alpha,"beta" = beta)
  return(return_list)
}
# Gamma Coeficients for fit and comparision from sample of population
GammaCoefficientsMM = function(data,n)
{
  new_data_set = GenerateMiniDataset(data,n)
  moments = MethodOfMoments(new_data_set,n)
  alpha = (moments$mean**2)/moments$variance
  beta = sqrt(alpha/moments$variance)
  return_list = list("alpha" = alpha,"beta" = beta)
  return(return_list)
}

alphaBeta = function(mean,var)
{
  alpha = mean**2 / var
  beta = 1/(mean / var)
  return(params = list(alpha = alpha, beta = beta))
}

alphaCoef = GammaCoefficients(InterSpace)$alpha
betaCoef = GammaCoefficients(InterSpace)$beta

no_samples = 2000

alphaCoefMM = GammaCoefficientsMM(InterSpace,no_samples)$alpha
betaCoefMM = GammaCoefficientsMM(InterSpace,no_samples)$beta
# c) Narysuj na jednym wykresie histogram z punktu a. oraz funkcje
# gêstoœci rozk³adu gamma o parametrach uzyskanych w punkcie b.
# co mozna powiedzieæ na podstawie tych wykresów?

# set.seed(1) # fixed sample

fit = EstimatorNW(no_samples,alphaCoef,betaCoef)
curve(dgamma( x , shape = fit[1] , scale = fit[2]),col="red",lwd=2,add=T)

#fitMM = EstimatorNW(no_samples,alphaCoefMM,betaCoefMM) #check
curve(dgamma( x , shape = alphaCoefMM , scale = 1/betaCoefMM),col="blue",lwd=1,add=T)
legend("topright", c("Data", "GammaNW","GammaMoM"), col=c("black", "red","blue"), lwd=2)

grid(lwd = 1)

# d) Metod¹ bootstrapu parametrycznego wyznacz odchylenia standardowe estymatorów
# z punktu b. oraz przedzia³y ufnoœci na poziomie ufnoœci 95% dla parametrów alfa oraz beta.
# Porównaj wyniki uzyskane dla estymatorów metody momentów i estymatorów
# najwiêkszej wiarygodnoœci

# bootstrap
InterSpace_boot = boot(data = InterSpace, statistic = meanFunc, R = InterSpace_len)
print(InterSpace_boot)

# Confidence interval for beta for estimated alpha and beta
#x = rbeta(no_samples, alphaCoef, betaCoef)
#x.boot = one.boot(x, meanFunc, R=no_samples)
#print(boot.ci(x.boot, type="bca"))


# Comparison of estimator methods
# Maximum Likelihood estimation
NW_alpha = fit[1]
NW_beta = fit[2]
NW_mean = NW_alpha*(NW_beta)
NW_var = NW_alpha*(NW_beta)**2
NW_std = sqrt(NW_var)


#Method of Moments for sampleset
MM_mean = alphaCoefMM*(1/betaCoefMM)
MM_var = alphaCoefMM*(1/betaCoefMM)**2
MM_std = sqrt(MM_var)

# Build-in R function
mean = alphaCoef*(1/betaCoef)
var = alphaCoef*(1/betaCoef)**2
std = sqrt(MM_var)
cat(sprintf('\n NW mean = %s\n NW var. = %s\n NW std. = %s\n NW_alpha = %s\n NW_beta = %s\n
 MM mean = %s\n MM var. = %s\n MM std. = %s\n MM_alpha = %s\n MM_beta = %s\n
 R  mean = %s\n R  var. = %s\n R  std. = %s\n R alpha = %s\n R beta = %s\n\n'
            , NW_mean,NW_var,NW_std,NW_alpha,NW_beta
            , MM_mean,MM_var,MM_std,alphaCoefMM,1/betaCoefMM
            , mean,var,std,alphaCoef,1/betaCoef))


# Quandile function for distribution data from beta and alpha coef
gammaR = rgamma(no_samples,alphaCoef,1/betaCoef)
gammaMM = rgamma(no_samples,alphaCoefMM,1/betaCoefMM)
gammaNW = rgamma(no_samples,NW_alpha,1/NW_beta)
cat(sprintf('Range: Min %s <-> Max %s\n\n   ', min(gammaR),max(gammaR)))
Tr = (1-0.9)/2
qR = quantile(gammaR,c(Tr, 1-Tr),na.rm = TRUE)
print(qR)

cat(sprintf('Range: Min %s <-> Max %s\n\n   ', min(gammaMM),max(gammaMM)))
qMM = quantile(gammaMM,c(Tr, 1-Tr),na.rm = TRUE)
print(qMM)

cat(sprintf('Range: Min %s <-> Max %s   ', min(gammaNW),max(gammaNW)))
qNW = quantile(gammaNW,c(Tr, 1-Tr),na.rm = TRUE)
print(qNW)



#Computation of the standard error of the MM_mean
sem_MM = MM_std/sqrt(InterSpace_len)
#5% and 95% confidence intervals of the MM_mean
cat(sprintf('\n\n\n\n\n\nConfidence intervals of the MM_mean:\n 5 percent:  %s\n and\n 95 percent: %s\n\n'
            , MM_mean-2*sem_MM,MM_mean+2*sem_MM))

#Computation of the error of the MM_var
err_MM_var = (InterSpace_len-1)*MM_var/qchisq(c(.95,.05), InterSpace_len-1)
cat(sprintf('Confidence intervals of the MM_var:\n 5 percent:  %s\n and\n 95 percent: %s\n\n'
            , err_MM_var[1],err_MM_var[2]))


# Computation of the alphaCoefMM err
err_MM_alphaBeta5 = alphaBeta(MM_mean-2*sem_MM,err_MM_var[1])
err_MM_alphaBeta95 = alphaBeta(MM_mean+2*sem_MM,err_MM_var[2])

cat(sprintf('Confidence intervals of the alphaMM
            :\n 5 percent:  %s and 95 percent: %s\n\n'
            , err_MM_alphaBeta5$alpha,err_MM_alphaBeta95$alpha))
cat(sprintf('Confidence intervals of the betaMM
            :\n 5 percent:  %s and 95 percent: %s\n\n'
            , err_MM_alphaBeta5$beta,err_MM_alphaBeta95$beta))
###############################################################

#Computation of the standard error of the NW_mean
sem_NW = NW_std/sqrt(InterSpace_len)
#5% and 95% confidence intervals of the mean
cat(sprintf('Confidence intervals of the NW_mean:\n 5 percent:  %s\n and\n 95 percent: %s\n\n'
            , NW_mean-2*sem_NW,NW_mean+2*sem_NW))

#Computation of the NW_var err
err_NW_var = (InterSpace_len-1)*NW_var/qchisq(c(.95,.05), InterSpace_len-1)
cat(sprintf('Confidence intervals of the NW_var:\n 5 percent:  %s\n and\n 95 percent: %s\n\n'
            , err_NW_var[1],err_NW_var[2]))

# Computation of the alphaCoefMM err
err_NW_alphaBeta5 = alphaBeta(NW_mean-2*sem_NW,err_NW_var[1])
err_NW_alphaBeta95 = alphaBeta(NW_mean+2*sem_NW,err_NW_var[2])

cat(sprintf('Confidence intervals of the alphaNW
            :\n 5 percent:  %s and 95 percent: %s\n\n'
            , err_NW_alphaBeta5$alpha,err_NW_alphaBeta95$alpha))
cat(sprintf('Confidence intervals of the betaNW
            :\n 5 percent:  %s and 95 percent: %s\n\n'
            , err_NW_alphaBeta5$beta,err_NW_alphaBeta95$beta))
###############################################################


#Computation of the standard error of the mean
sem = std/sqrt(InterSpace_len)
#5% and 95% confidence intervals of the NW_mean
cat(sprintf('Confidence intervals of the mean:\n 5 percent:  %s\n and\n 95 percent: %s\n\n'
            , mean-2*sem,mean+2*sem))

#Computation of the var error
err_var = (InterSpace_len-1)*var/qchisq(c(.95,.05), InterSpace_len-1)
cat(sprintf('Confidence intervals of the var:\n 5 percent:  %s\n and\n 95 percent: %s\n\n'
            , err_var[1],err_var[2]))

# Computation of the alphaCoefMM
err_alphaBeta5 = alphaBeta(mean-2*sem,err_var[1])
err_alphaBeta95 = alphaBeta(mean+2*sem,err_var[2])

cat(sprintf('Confidence intervals of the alpha
            :\n 5 percent:  %s and 95 percent: %s\n\n'
            , err_alphaBeta5$alpha,err_alphaBeta95$alpha))
cat(sprintf('Confidence intervals of the beta
            :\n 5 percent:  %s and 95 percent: %s\n\n'
            , err_alphaBeta5$beta,err_alphaBeta95$beta))
###############################################################
###############################################################
# Z3
# a) Wyestymuj œredni¹ i wariancje "generuj¹cego rozk³adu

Gnorm = c(-27.708919904,-1.667698942,12.989606236
           ,7.443708265,8.790545464
           ,-5.079198171,5.016025220
           ,-19.941320219
           ,-7.641259981
           ,-8.463103348
           ,10.276204434,2.801759173,-13.177411156
           ,2.866651175,-6.657347527,-6.333788150,15.098091896,4.587872107,8.577079846
           ,15.375578463,5.523607716,10.711963309,12.561499866,19.388521240
           ,9.662048222,22.771257610
           ,3.316518780,2.329828257,26.818685160
           ,-4.135963979,2.960030558,15.878139117
)

# Mean and Variance
GN_Mean = mean(Gnorm,na.rm=TRUE)
GN_Variance = var(Gnorm,na.rm=TRUE)
cat(sprintf('\n\n\n Mean = %s\n Var. = %s\n\n\n',GN_Mean,GN_Variance))


# b) Podaj 90%, 95% oraz 99% przedzia³y ufnoœci dla mean

# confidence interval for diffrent confidence value
print(t.test(Gnorm, conf.level = 0.99))
print(t.test(Gnorm, conf.level = 0.95))
print(t.test(Gnorm, conf.level = 0.90))