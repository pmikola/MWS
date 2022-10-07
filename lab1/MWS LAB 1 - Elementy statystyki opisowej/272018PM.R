
rm(list=ls())
#1)
# Funkcja rozkładu prawdopodobieństwa
par(mfrow=c(1,1))
given_dens = function(x)
{
  ifelse(x>=-2 && x<= -1 || x>=1 && x<=2,0.5,0)
}
dx = 0.001
x = seq(-3,3,dx)
vdens = Vectorize(given_dens)

plot(x,vdens(x),col="red",main= 'Funkcja rozkładu gęstości prawdopodobieństwa',xlab='x',ylab='f(x)')
grid(lwd = 1)
# a ) Wyznacz i narysuj dystrybuante tego rozkłądu
#skalowanie
scaled_dens = function(x)
{
  div = integrate(vdens,-3,3)$value #final estimate of the integral
  new_scaled_vdens = vdens(x) / div
  return(new_scaled_vdens)
}

#dystrybuanta
CDF = function(x,dx)
  {
    CDF = cumsum(scaled_dens(x))*dx
  return(CDF)
}

y = CDF(x,dx)
plot(x,y,main= 'Dystrybuanta zadanego rozkładu',xlab='x',ylab='F(x)')
grid(lwd = 1)

#dystrybuanta odwrotna
ICDF = approxfun(y,x)
plot(ICDF,x,xlim=c(0,1),main= 'Dystrybuanta odwrócona')
grid(lwd = 1)
# b) Wygeneruj 1000-elementową próbe losową z tego rozkładu

#wygenerowanie próbek w zakresie zadanego rozkładu
randSamp = runif(1000,0,1)

#wygenerowanie próby losowej poprzez przemnożenie odwróconej dystrybuanty przez wygenerowane próbki + histogram
ICDFRS = ICDF(randSamp)
histogram = hist(ICDFRS,50,main= 'Zadany rozkład oraz histogram próby losowej',xlab='x',ylab='f(x)',freq=F)

# c) porównanie histogramu z funkcją gęstości prawdopodobieństwa
lines(x,vdens(x),col="red")
grid(lwd = 1)

# d) Dla wygenerowanej w poprzednim punkcie próby losowej narysuj wykres dystrybuanty
# empirycznej i porównaj go z dystrybuant¡ z punktu a.
y1 = ecdf(ICDFRS)
plot(y1,x,main= 'Dystrybuanta empiryczna oraz zadana',xlab='x',ylab='F(x)',col='blue')
lines(x,y,col="red")
grid(lwd = 1)

# 2)
# a) Narysuj jak zmieniaªa si¦ liczba wypadków i ich ofiar śmiertelnych w kolejnych latach
# (dwa osobne wykresy).
# b) Wykresy otrzymane w poprzednim punkcie są dość postrzępione. Spróbuj na każdym
# powyższych dwóch wykresów dorysować gładkę krzywę ilustrujące te same zależno-
# ści, lecz nie rozpraszające uwagi odbiorcy nieistotnymi krótkotrwałymi wahaniami.

# wczytanie csv do zmiennej DATA

DATA = read.csv("katastrofy.csv",header=TRUE,sep=",")

#Wyłuskanie daty i innych danych
allDates = as.Date(DATA$Date, "%m/%d/%Y")
accYears = as.integer(format(allDates, "%Y"))


DATA$Fatalities[is.na(DATA$Fatalities)] = 0
DATA$Aboard[is.na(DATA$Aboard)] = 0

Fatals = as.integer(DATA$Fatalities)
Aboard = as.integer(DATA$Aboard)

#agregacja danych do data frameów
FatalitiesDF = aggregate(x = Fatals,by = list(accYears), FUN = sum)
AboardDF = aggregate(x = Aboard, by = list(accYears), length)

#fit z uzyciem fit
names(AboardDF) <-  c("Year","Accidents")
par(mfrow=c(1,2))

X = AboardDF$Year[1:98]
Y = AboardDF$Accidents[1:98]

fit  = lm(Y~poly(X,3,raw=TRUE))
xx = seq(min(X),max(X), length=length(X))
plot(AboardDF,pch=3,main = "Air Crashes")
lines(xx, predict(fit, data.frame(X = xx)), col="red")
grid(lwd = 1)


# fit z uzyciem loess
names(AboardDF) <-  c("Year","Accidents")
plot(AboardDF,pch=20,main = "Air Crashes")

X = AboardDF$Year[1:98]
Y = AboardDF$Accidents[1:98]

y.loess = loess(Y ~ X, span=0.3, data.frame(x=X, y=Y))
y.predict = predict(y.loess, data.frame(x=X))
lines(X,y.predict,col="red")
grid(lwd = 1)


#dev.new()

#fit z uzyciem fit
names(FatalitiesDF) =  c("Year","Deaths")

X = FatalitiesDF$Year[1:98]
Y = FatalitiesDF$Deaths[1:98]

fit  <- lm(Y~poly(X,3,raw=TRUE))
xx <- seq(min(X),max(X), length=length(X))
plot(FatalitiesDF,pch=3,main = "Air Crashes")
lines(xx, predict(fit, data.frame(X = xx)), col="red")
grid(lwd = 1)


# fit z uzyciem loess
names(FatalitiesDF) =  c("Year","Deaths")
plot(FatalitiesDF,pch=20,main = "Air Crashes")

X = FatalitiesDF$Year[1:98]
Y = FatalitiesDF$Deaths[1:98]

y.loess = loess(Y ~ X, span=0.3, data.frame(x=X, y=Y))
y.predict = predict(y.loess, data.frame(x=X))
lines(X,y.predict,col="red")
grid(lwd = 1)






