rm(list=ls())
# a=c(79.98,80.04,80.02,80.04,80.03,80.03,80.04,79.97,80.05,80.03,80.02,80.00,80.02);
# b=c(80.02,79.94,79.98,79.97,79.97,80.03,79.95,79.97);
# print(mean(a))
# print(mean(b))
# Pewien generator liczb losowych wygenerował 4 liczby z rozkładu normalnego N(u1; sigma2),
# a inny generator liczb losowych wygenerował 5 liczb z rozkładu normalnego N(u2; sigma2) (ta
# sama wariancja w obydwu przypadkach). Pliki z tymi liczbami to los1.txt i los2.txt.

los1 = as.numeric(c(t(read.delim2("los1.txt",header=FALSE,sep=""))))
los2 = as.numeric(c(t(read.delim2("los2.txt",header=FALSE,sep=""))))
#print(los1)
#print(los2)
# Z1
# a. Oszacuj wartości oczekiwane u1 i u2 rozkładów, z których generowane były liczby oraz
# ich różnicę u1 - u2.

u_los1 = mean(los1)
u_los2 = mean(los2)
los_diff = u_los1 - u_los2
print(los_diff)
#print(t.test(los1,los2))

# b. Oszacuj wariancję sigma2.

n1=length(los1)
n2=length(los2)
si1=sum((los1-u_los1)^2)/(n1-1)
si2=sum((los2-u_los2)^2)/(n2-1)
var=((n1-1)*si1+(n2-1)*si2)/(n1+n2-2)
print(var)

# c. Oszacuj odchylenie standardowe błędu wyestymowanej w punkcie a. różnicy u1 - u2.
s=var*((1/n1)+(1/n2))
deviation=sqrt(s)
print(deviation)
# d. Czy w opisanej sytuacji właściwym testem do sprawdzenia równości między średnimi
# będzie test jednostronny, czy dwustronny?

# Dla testu jednostronnego u1 =/= u2 - hipoteza H1
# Dla testu dwustronnego u1 = u2 - Hipoteza H0

# e. Czy hipoteza o równości średnich zostałaby odrzucona przez dwustronny test na po-
# ziomie istotności alfa = 0:1?

ttest1=t.test(los1,los2,var.equal=TRUE, alternative = c("two.sided"),conf.level = 0.9)
print(ttest1) 

T = (u_los1 - u_los2)/(sqrt(var)*sqrt(((1/n1) + 1/n2)));
alfa = 0.1;
c = qt((1 - alfa/2), n1 + n2 - 2);

if(abs(T)<=c){
  message("H0: means are equal\n");
}else{    #H1
  message("H1: means are diffrent\n");
}


# Z2
# Wpliku lozyska.txt podane są czasy (w milionach cykli) pracy (do momentu uszkodzenia)
# łożysk wykonanych z dwóch różnych materiałów.


# a. Przeprowadź test tego, że nie ma różnicy między łożyskami wykonanymi z tych ma-
# teriałów, zakładając, że czas pracy do momentu uszkodzenia opisuje się rozkładem
# normalnym.

data= read.csv(file= "lozyska.txt", sep = ",")
data1=data[,1]
data2=data[,2]

n=length(data1)
m=length(data2)

ttest2=t.test(data1,data2)
print(ttest2)
# b. Przeprowadź analogiczny test, bez zakładania normalności rozkładów,

Wt=wilcox.test(data1,data2)
print(Wt)

# c. Który z powyższych testów jest w rozważanym przypadku odpowiedniejszy?
# test drugi ponieważ założenie rozkładu normalności nie jest pewne

# Z3
# W poniższej tabeli przedstawiona jest długość drogi hamowania y (pewnego pojazdu, na
# pewnego rodzaju drodze) w funkcji prędkości v. Dopasuj liniowe zależności opisujące dłu-
# gość drogi y oraz pierwiastek tej długości
# sqrt(y) w funkcji prędkości v. Która z liniowych
# zależności jest bardziej zgodna z danymi? Dlaczego?

y_kmh = c(0.0,33.0,33.0,49.1,65.2,78.5,93.0)
y_m = c(0.0,4.7,4.1,10.3,22.3,34.4,44.4)

plot(y_kmh,y_m,pch=20,col="red",main= 'Breaking distance vs. Speed',xlab='Speed [km/h]',ylab='Breaking Distance [m]')

# Linear Regression
abline(lm(y_m ~ y_kmh),col="red")

# Sqrt(Lin) Regression
abline(lm(sqrt(y_m) ~ y_kmh),col="red")
# Polyfit
fit = lm(y_m~poly(y_kmh,5,raw=TRUE))
xx = seq(min(y_kmh),max(y_kmh), length=length(y_kmh))
lines(xx, predict(fit, data.frame(X = xx)), col="blue")

# Loess
y.loess = loess(y_m ~ y_kmh, span=1, data.frame(x=y_kmh, y=y_m))
y.predict = predict(y.loess, data.frame(x=y_kmh))
lines(y_kmh,y.predict,col="green")


legend("top", 
       legend = c("Stop Point","y_lin","sqrt(y_lin)", "Polyfit","Loess"), 
       col = c(rgb(1,0,0,0.7),
               rgb(1,0,0,0.7),
               rgb(1,0,0,0.7),
               rgb(0,0,1,0.7),
               rgb(0,1,0,0.7)),
       pch = c(20,1,1,1,1), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.3, 0.01))
grid(lwd = 1)