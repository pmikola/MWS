rm(list=ls())
library("pracma")

# Z1
# Dane przedstawiają liczbę zanotowanych samobójstw w Stanach Zjednoczonych
# w 1970 roku z podziałem na poszczególne miesiące. Zbadaj, czy zamieszczone poniżej
# dane wskazują na sezonową zmienność liczby samobójstw, czy raczej świadczą o stałej
# intensywności badanego zjawiska.

months = seq(1,12,1)
daysInMonths = c(31,28,31,30,31,30,31,31,30,31,30,31)
suicideCount = c(1867,1789,1944,2094,2097,1981,1887,2024,1928,2032,1978,1859)

# Simple Linear Reggresion
plot(months,suicideCount,pch=20,main= 'Suicides per month \n + Regressions',xlab='Months',ylab='Suicide Number',col='red')
abline(lm(suicideCount ~ months))

# Polyfit
fit = lm(suicideCount~poly(months,3,raw=TRUE))
xx = seq(min(months),max(months), length=length(months))
lines(xx, predict(fit, data.frame(X = xx)), col="blue")

# Loess
y.loess = loess(suicideCount ~ months, span=0.7, data.frame(x=months, y=suicideCount))
y.predict = predict(y.loess, data.frame(x=months))
lines(months,y.predict,col="green")

legend("bottomleft", 
       legend = c("Suicides","Linear", "Polyfit","Loess"), 
       col = c(rgb(1,0,0,0.7),
               rgb(0,0,0,0.7),
               rgb(0,0,1,0.7),
               rgb(0,1,0,0.7)),
       pch = c(20,1,1,1), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.3, 0.01))

grid(lwd = 1)


# H0 - stała intensywność / jednostajność (nie zależne)
# czy?
# H1 - sezonowość (zależne)

# H0: T <= c
# H1: T > c

# Test of distribution compliance with Chi-sq
alfa = 0.1
p = 1/sum(daysInMonths)
pMonth = daysInMonths * p
print(cat("Probability Per Month: ", pMonth)) # months

# Deegres of freedom = r-1 where r is number of months
c = qchisq(1 - alfa,length(months)-1)
T = sum((suicideCount-sum(pMonth)*1/sum(daysInMonths))^2/(sum(pMonth)*sum(daysInMonths)))
pValue=1-pchisq(T,length(months) - 1)
print(pValue)
print(cat("T: ", T))
chiTest = chisq.test(pMonth,suicideCount)
# Contribiution in % of a given cell to the total Chi-sq
contrib = 100*chiTest$residuals^2/chiTest$statistic
print(round(contrib, 5))
print(chiTest)
print(chiTest$p.value)
print(T)

if (T <= c) {
        print("H0 WIN - brak zależności")
} else {
        print("H1 WIN - sezonowość")
}

# Z2
# a) Wy estymuj wartość średnią i odchylenie standardowe temperatury ciała, osobno dla
# mężczyzn i kobiet, a następnie wykreśl tzw. wykresy kwantyl-kwantyl dla rozkładu
# normalnego (na płaszczyźnie X-Y narysuj punkty o współrzędnych (x; y), gdzie x jest
# kwantylem rzędu alfa (0; 1) z próby, a y jest kwantylem rzędu alfa rozkładu normal-
# nego o wy estymowanych parametrach). Co możesz powiedzieć na podstawie otrzyma-
# nych wykresów o zgodności rozkładu temperatury ciała mężczyzn/kobiet z rozkładem
# normalnym? (W celu sprawdzenia, jak taki wykres mogłyby wyglądać gdyby dane
# pochodziły z rozkładu normalnego można przeprowadzić odpowiednie symulacje.)

DATA = read.delim("tempciala.txt",header=TRUE,sep=",")

# Indexes
Midx = which(DATA$p == 1)
Kidx = which(DATA$p == 2)

# Temp and Pulse Data
tempM = DATA$temp[Midx]
tempK = DATA$temp[Kidx]
pulseM = DATA$tetno[Midx]
pulseK = DATA$tetno[Kidx]

tMmean = mean(tempM)
tKmean = mean(tempK)
tMsd = sd(tempM)
tKsd = sd(tempK)
pMmean = mean(pulseM)
pKmean = mean(pulseK)
pMsd = sd(pulseM)
pKsd = sd(pulseK)

#qqnorm(tempM,main="Temperature of Men",xlab='Quantiles',ylab='Temperature')
#qqline(tempM, lwd=2, col="red")
#grid(lwd = 1)

#qqnorm(tempK,main="Temperature of Women",xlab='Quantiles',ylab='Temperature')
#qqline(tempK, lwd=2, col="red")
#grid(lwd = 1)

#qqnorm(pulseM,main="Pulse of Men",xlab='Quantiles',ylab='Pulse')
#qqline(pulseM, lwd=2, col="red")
#grid(lwd = 1)

#qqnorm(pulseK,main="Pulse of Women",xlab='Quantiles',ylab='Pulse')
#qqline(pulseK, lwd=2, col="red")
#grid(lwd = 1)

x = seq(0.01, 0.99, 0.005)
tempQM = quantile(tempM, x)
tempQK = quantile(tempK, x)
tMnorm = qnorm(x, mean = tMmean, sd = tMsd)
tKnorm = qnorm(x, mean = tKmean, sd = tKsd)

qqplot(tempQM,tMnorm, xlab = 'Probe Quantile', ylab = 'Quantile Normal Distribution', main = 'Men',col="red")
abline(0,1)
grid(lwd = 1)

qqplot(tempQK,tKnorm, xlab = 'Probe Quantile', ylab = 'Quantile Normal Distribution', main = 'Women',col="red")
abline(0,1)
grid(lwd = 1)

# Simulation Men rnorm
tMSim = rnorm(200, mean = tMmean, sd = tMsd)
tMSimMean = mean(tMSim)
tMSimSd = sd(tMSim)
tMSimQuantile = quantile(tMSim, x);
tMNormSim = qnorm(x, mean = tMSimMean, sd = tMSimSd)
qqplot(tMSimQuantile, tMNormSim, xlab = 'Probe Quantile', ylab = 'Quantile Normal Distribution', main = 'Simulation for Men',col="red")
abline(0,1)
grid(lwd = 1)

# Simulation WOMen rnorm
tKSim = rnorm(200, mean = tKmean, sd = tKsd)
tKSimMean = mean(tKSim)
tKSimSd = sd(tKSim)
tKSimQuantile = quantile(tKSim, x);
tKNormSim = qnorm(x, mean = tKSimMean, sd = tKSimSd)
qqplot(tKSimQuantile, tKNormSim, xlab = 'Probe Quantile', ylab = 'Quantile Normal Distribution', main = 'Simulation for Women',col="red")
abline(0,1)
grid(lwd = 1)



# b) Przeprowadź testy, osobno dla mężczyzn i kobiet, 
# tego że średnia temperatura ciała jest równa 36.6 stopni 
# wobec hipotezy alternatywnej, że ta średnia temperatura jest jednak inna. 
# H0 : equal to 36.6
# H1 : Not equal to 36.6

alfa = 0.01 			
n = length(tempM)
T = ((tMmean - 36.6) / tMsd) * sqrt(n)
c = qt(1 - (alfa/2), n)
if (T < c) {
        print("H0 WIN 4 MEN - EQ")
} else {
        print("H1 WIN 4 MEN - NEQ")
}

m = length(tempK)
T = ((tKmean - 36.6) / tKsd) * sqrt(m)
c = qt(1 - alfa/2, m)
if (T < c) {
        print("H0 WIN 4 WOMEN - EQ")
} else {
        print("H1 WIN 4 WOMEN - NEQ")
}

