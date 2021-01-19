
#importuje biblioteki
library (dplyr)
library (moments)
library (ggplot2)
library(tidyr)
library(reshape2)
library(stats4)

#ladowanie danych df - dataset na ktorym bede w wiekszosci pracowal
#          original_df - dataset niezmieniony 

PATHNAME = "C:\\Users\\Adam\\Desktop\\projekt_stat\\diabetes.csv"
originaldf <- read.csv(PATHNAME)
df = originaldf
colnames(df)[7] <- "DPF"


#sprawdzenie ilosci brakujacych danych oraz ilosci zer w kolumnach
summary(df)
colSums(df == 0)


#sprawdzenie skosnosci danych w kolumnach
#uzupelnienie odpowiednich kolumn mediana lub srednia arytmetyczna w zalenosci od skosnosci

skewness(df)
df$Insulin[df$Insulin == 0]<- median(df$Insulin[df$Insulin >0])
df$Glucose[df$Glucose == 0]<- median(df$Glucose[df$Glucose >0])
df$BloodPressure [df$BloodPressure  == 0]<- mean(df$BloodPressure [df$BloodPressure  >0])
df$SkinThickness[df$SkinThickness == 0]<- median(df$SkinThickness[df$SkinThickness >0])
df$BMI[df$BMI == 0]<- median(df$BMI[df$BMI >0])


#stworzenie kolejnego datasetu 
#w ktorym odrzucamy wiersze z wartoscia NA w kolumnie Insulin
cleandf = originaldf
cleandf$Insulin[cleandf$Insulin == 0]<- NA
cleandf = cleandf[!is.na(cleandf$Insulin),]



#rysuje histogramy wszystkich cech
ggplot(gather(df), aes(x = value)) +  
  geom_histogram(bins = 10) + facet_wrap(~key, scales = 'free_x') + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

#sprawdzam ilosc wartosci 0 i 1 w kolumnie Outcome
table(df$Outcome)

#definiuje 3 funckje negative log likelihood

NLLGlucose <- function(theta0,theta1) {
  -sum ( -0.5* log(theta1*2*pi) - 0.5*( df$Glucose- theta0)^2/theta1 )
}

NLLBMI <- function(theta0,theta1) {
  -sum ( -0.5* log(theta1) - 0.5*( df$BMI- theta0)^2/theta1 )
}

NLLBP <- function(theta0,theta1) {
  -sum ( -0.5* log(theta1) - 0.5*( df$BloodPressure- theta0)^2/theta1 )
}

#minimalizuje funkcje NLL za pomoca funkcji mle z pakietu stats4
#wypisuje uzyskane paramtery
Glucoseest <- stats4::mle(minuslog=NLLGlucose, start=list(theta0=100,theta1=900))
coef(Glucoseest)

BMIest <- stats4::mle(minuslog=NLLBMI, start=list(theta0=30,theta1=60))
coef(BMIest)

BPest <- stats4::mle(minuslog=NLLBP, start=list(theta0=70,theta1=375))
coef(BPest)

#rysuje histogramy tych trzech cech razem z wykresem rozkladow o wyznaczonych parametrach

h <- hist(df$Glucose, breaks = 10, density = 10,
          col = "lightgray", xlab = "", ylab = "", main = "Glucose histogram") 
xfit <- seq(min(df$Glucose), max(df$Glucose), length = 40) 
yfit <- dnorm(xfit, mean = coef(Glucoseest)[1], sd = sqrt(coef(Glucoseest)[2])) 
yfit <- yfit * diff(h$mids[1:2]) * length(df$Glucose) 
lines(xfit, yfit, col = "black", lwd = 2)

h <- hist(df$BMI, breaks = 10, density = 10,
          col = "lightgray", xlab = "", ylab = "",main = "BMI histogram") 
xfit <- seq(min(df$BMI), max(df$BMI), length = 40) 
yfit <- dnorm(xfit, mean = coef(BMIest)[1], sd = sqrt(coef(BMIest)[2])) 
yfit <- yfit * diff(h$mids[1:2]) * length(df$BMI) 

lines(xfit, yfit, col = "black", lwd = 2)

h <- hist(df$BloodPressure, breaks = 10, density = 10,
          col = "lightgray", xlab = "",ylab= "",main = "BloodPressure") 
xfit <- seq(min(df$BloodPressure), max(df$BloodPressure), length = 40) 
yfit <- dnorm(xfit, mean = coef(BPest)[1], sd = sqrt(coef(BPest)[2])) 
yfit <- yfit * diff(h$mids[1:2]) * length(df$BloodPressure) 

lines(xfit, yfit, col = "black", lwd = 2)


#sprawdzam rzeczywiste wartosci parametrow i porownuje z wyznaczonymi

mean(df$Glucose)
var(df$Glucose)
coef(Glucoseest)

mean(df$BMI)
var(df$BMI)
coef(BMIest)

mean(df$BloodPressure)
var(df$BloodPressure)
coef(BPest)

#rysuje wykresy pudelkowe

ggplot(stack(data.frame(df$Glucose,df$BloodPressure,df$SkinThickness,df$BMI,df$Age,df$Pregnancies)), aes(x = ind, y = values)) +
  geom_boxplot() +
  coord_flip() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())


qplot(cleandf$Insulin, geom = "boxplot",xlab = "",main = "Insulin boxplot") + 
  theme(plot.title = element_text(hjust = 0.5))

qplot(df$DPF, geom = "boxplot",xlab = "", main = "DPF boxplot") + 
  theme(plot.title = element_text(hjust = 0.5))

#sprawdzam liczbe wartosci odstajacych

length(boxplot(cleandf$Insulin)$out)
length(boxplot(df$DPF)$out)
length(boxplot(df$SkinThickness)$out)


#rysuje heatmap z pomoca wspolczynnika korelacji

ggplot(melt(cor(df)), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())  +
  geom_text(aes(label = round(value, 2)))


#rysuje zaleznosci pomiedzy cechami 
#razem z doposowana prosta za pomoca regresji liniowej

ggplot(df, aes(x = Glucose, y = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") + 
  ggtitle("Zachorowania na cukrzyce w zaleznosci od poziomu glukozy we krwi") + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(cleandf, aes(x = Insulin, y = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") + 
  ggtitle("Zachorowania na cukrzyce w zaleznosci od poziomu insuliny") + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(df, aes(x = BMI, y = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") + 
  ggtitle("Zachorowania na cukrzyce w zaleznosci od wskaznika BMI") + 
  theme(plot.title = element_text(hjust = 0.5))















