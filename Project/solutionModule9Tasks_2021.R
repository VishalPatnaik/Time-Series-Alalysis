
library(TSA)
library(fUnitRoots)
library(forecast)
library(CombMSC)
library(lmtest)
library(fGarch)
library(rugarch)
library(tseries)

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}
residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH", "garch", "fGARCH")[1]){
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "garch"){
    res.model = model$residuals[start:model$n.used]  
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else if (class == "fGARCH"){
    res.model = model@residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}


# --- TASK 1 ---
data("google") # This is already returns series
par(mfrow=c(1,1))
plot(google,type='o',main="Time series plot of daily returns of Google stock")
# There is sign of neither a trend nor seasonality. Observations are bouncing around the mean level. But changing variance is obvious.
mean(google) # Mean is very close to zero

McLeod.Li.test(y=google,main="McLeod-Li Test Statistics for Daily Google Returns")
# McLeod-Li test is significnat at 5% level of significance for all lags. This gives a strong idea about existence of volatiliy clustering.

par(mfrow=c(1,2))
acf(google, main="The ACF plot for return series")
pacf(google, main="The PACF plot for return series")

eacf(google)
# ACF, PACF and EACF all shows pattern of white noise for the correlation structure. However, there is an ARCH effect present in the series.

#So we'll use absolute value and square transformations to figure out this ARCH effect.
abs.google = abs(google)
sq.google = google^2

par(mfrow=c(1,2))
acf(abs.google, ci.type="ma",main="The ACF plot for absolute return series")
pacf(abs.google, main="The PACF plot for absolute return series")

eacf(abs.google)
# After the absolute value transformation, we observe many signficicant lags in 
# both ACF and PACF. Also, EACF do not suggest an ARMA(0,0) model.
# From the EACF, we can identify ARMA(1,1), ARMA(1,2), and ARMA(2,2) models for absolute 
# value series. 
# These models correspond to parameter settings:
# max(p,q) = 1 and q = 1 ==> max(p,1) = 1 and q = 1 ==> p can be either 0 or 1.
# max(p,q) = 1 and q = 2 ==> max(p,2) = 1 and q = 2 ==> There is no model.
# max(p,q) = 2 and q = 2 ==> max(p,2) = 2 and q = 2 ==> p can be 0, 1 and 2.
# So the corresponding tentative GARCH models are GARCH(0,1), GARCH(1,1), GARCH(0,2), GARCH(1,2), and GARCH(2,2).

par(mfrow=c(1,2))
acf(sq.google, ci.type="ma",main="The ACF plot for squared return series")
pacf(sq.google, main="The PACF plot for squared return series")
eacf(sq.google)
# After the square transformation, we observe many significant lags in both ACF and PACF. Also, EACF do not suggest an ARMA(0,0) model.
# From the EACF, we can identify ARMA(1,1), ARMA(1,2), ARMA(2,1), and ARMA(2,2) models for squared series. 
# These models correspond to parameter settings:
# max(p,q) = 1 and q = 1 ==> max(p,1) = 1 and q = 1 ==> p can be either 0 or 1.
# max(p,q) = 1 and q = 2 ==> max(p,2) = 1 and q = 2 ==> There is no model.
# max(p,q) = 2 and q = 1 ==> max(p,1) = 2 and q = 1 ==> p can only be 2.
# max(p,q) = 2 and q = 2 ==> max(p,2) = 2 and q = 2 ==> p can be 0, 1 and 2.
# The corresponding tentative GARCH models are GARCH(0,1), GARCH(1,1), GARCH(2,1), GARCH(0,2), GARCH(1,2), and GARCH(2,2).

# Overall, we have {GARCH(0,1), GARCH(1,1), GARCH(2,1), GARCH(0,2), GARCH(1,2), GARCH(2,2)}.

m.01 = garch(google,order=c(0,1),trace = FALSE)
summary(m.01) 
residual.analysis(model = m.01, class= "garch") # use class= "garch" for garch() function from tseries package

m.11 = garch(google,order=c(1,1),trace = FALSE)
summary(m.11) 
residual.analysis(model = m.11, class= "garch")

m.21 = garch(google,order=c(2,1),trace = FALSE)
summary(m.21)
residual.analysis(model = m.21, class = "garch")
residual.analysis(model = m.21, start = 3, class= "garch")

m.02 = garch(google,order=c(0,2),trace = FALSE)
summary(m.02)
residual.analysis(model = m.02, start = 3, class = "garch")

m.12 = garch(google,order=c(1,2),trace = FALSE)
summary(m.12)
residual.analysis(model = m.12, start = 3, class = "garch")

m.22 = garch(google,order=c(2,2),trace = FALSE)
summary(m.22)
residual.analysis(model = m.22, start = 3, class = "garch")

sc.AIC=AIC(m.01, m.11, m.21, m.02, m.12, m.22)
sc.BIC=AIC(m.01, m.11, m.21, m.02, m.12, m.22, k = log(length(google)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

# GARCH(1,1) model is the best one in terms of residual analysis and both AIC and BIC.

par(mfrow=c(1,1))
plot((fitted(m.11)[,1])^2,type='l',ylab='Conditional Variance',xlab='t',main="Estimated Conditional Variances of the Daily Returns")
# Changes in conditional variance at the beginning of the series and between observations 300 and 400, then the conditional variance settles down. 


m.11F = garchFit(formula = ~garch(1,1), data =google )
fGarch::predict(m.11F, n.ahead=10, trace=FALSE, plot=TRUE) 
# Forecasts for the confidance limits are based on the forecasts of conditional variance.


library(rugarch)
model<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                  distribution.model = "norm")
m.11R<-ugarchfit(spec=model,data=google)
par(mfrow=c(1,1))
plot(m.11R)
plot(m.11R, which = 1)

par(mfrow=c(1,1))
forc = ugarchforecast(m.11R, data = google, n.ahead = 10)
plot(forc)


# --- TASK 2 ---

library(tswge)
data("mm.eq")
wave = mm.eq
wave = ts(wave)
par(mfrow=c(1,1))
plot(wave,main="Time series plot of seismic Lg wave series")
# There is sign of neither a trend nor seasonality. Observations are bouncing around the mean level. But changing variance is obvious especially 
# at the begining of the observationi period.
mean(wave)

McLeod.Li.test(y=wave,main="McLeod-Li Test Statistics for seismic Lg wave series")
# McLeod-Li test is significnat at 5% level of significance for all lags. This gives a strong idea about existence of volatiliy clustering.

par(mfrow=c(1,2))
acf(wave, main="The sample ACF plot for seismic Lg wave series")
pacf(wave, main="The sample PACF plot for seismic Lg wave series")

eacf(wave)
# In the ACF, PACF, and EACF plots we observe significant correlations and there is no sign of 
# a white noise process. However, volatility clustering is obvious in the time series plot. So, 
# we will consider fitting an ARMA+GARCH model.

# par(mfrow=c(1,1))
# wave.positive = wave + min(abs(wave))+0.1
# BC = BoxCox.ar(wave.positive)
# BC$ci
# lambda = BC$lambda[which(max(BC$loglike) == BC$loglike)]
# lambda
# BC.wave = ((wave.positive^lambda)-1)/lambda
# 
# adf.test(BC.wave)
# 
# diff.BC.wave = diff(BC.wave)

wave.positive = wave + min(abs(wave))+0.1
waveReturn = diff(log(wave.positive))

plot(waveReturn,main="Time series plot of returns series for seismic Lg wave series")

adf.test(waveReturn)


par(mfrow=c(1,2))
acf(waveReturn, main="The ACF plot of returns series for seismic Lg wave series")
pacf(waveReturn, main="The PACF plot of returns series for seismic Lg wave series")
par(mfrow=c(1,1))
# Due to the changing variance, we do not have a clear picture in ACF and PACF.
# Referring very high autocorrelations {ARMA(0,3), ARMA(0,4), ARMA(0,5), ARMA(0,6)} can be identified

eacf(waveReturn)
#{ARMA(5,4), ARMA(5,5), ARMA(6,5)}

res = armasubsets(y=waveReturn,nar=14,nma=14,y.name='test',ar.method='ols')
plot(res)
#Additional to the previously identified models, {ARMA(2,4), ARMA(2,6)}

#Overall, 
#{ARMA(0,3), ARMA(0,4), ARMA(0,5), ARMA(0,6), ARMA(5,4), ARMA(5,5), ARMA(6,5), ARMA(2,4), ARMA(2,6)}

model_03_css = arima(waveReturn,order=c(0,0,3),method='CSS')
coeftest(model_03_css)
residual.analysis(model = model_03_css)

model_03_ml = arima(waveReturn,order=c(0,0,3),method='ML')
coeftest(model_03_ml)
residual.analysis(model = model_03_ml)


model_04_css = arima(waveReturn,order=c(0,0,4),method='CSS')
coeftest(model_04_css)
residual.analysis(model = model_04_css)

model_04_ml = arima(waveReturn,order=c(0,0,4),method='ML')
coeftest(model_04_ml)
residual.analysis(model = model_04_ml)


model_05_css = arima(waveReturn,order=c(0,0,5),method='CSS')
coeftest(model_05_css)
residual.analysis(model = model_05_css)

model_05_ml = arima(waveReturn,order=c(0,0,5),method='ML')
coeftest(model_05_ml)
residual.analysis(model = model_05_ml)


model_06_css = arima(waveReturn,order=c(0,0,6),method='CSS')
coeftest(model_06_css)
residual.analysis(model = model_06_css)

model_06_ml = arima(waveReturn,order=c(0,0,6),method='ML')
coeftest(model_06_ml)
residual.analysis(model = model_06_ml)


model_54_css = arima(waveReturn,order=c(5,0,4),method='CSS')
coeftest(model_54_css)
residual.analysis(model = model_54_css)

model_54_ml = arima(waveReturn,order=c(5,0,4),method='ML')
coeftest(model_54_ml)
residual.analysis(model = model_54_ml)


model_55_css = arima(waveReturn,order=c(5,0,5),method='CSS')
coeftest(model_55_css)
residual.analysis(model = model_55_css)

model_55_ml = arima(waveReturn,order=c(5,0,5),method='ML')
coeftest(model_55_ml)
residual.analysis(model = model_55_ml)


model_65_css = arima(waveReturn,order=c(6,0,5),method='CSS')
coeftest(model_65_css)
residual.analysis(model = model_65_css)

model_65_ml = arima(waveReturn,order=c(6,0,5),method='ML')
coeftest(model_65_ml)
residual.analysis(model = model_65_ml)


model_24_css = arima(waveReturn,order=c(2,0,4),method='CSS')
coeftest(model_24_css)
residual.analysis(model = model_24_css)

model_24_ml = arima(waveReturn,order=c(2,0,4),method='ML')
coeftest(model_24_ml)
residual.analysis(model = model_24_ml)


model_26_css = arima(waveReturn,order=c(2,0,6),method='CSS')
coeftest(model_26_css)
residual.analysis(model = model_26_css)

model_26_ml = arima(waveReturn,order=c(2,0,6),method='ML')
coeftest(model_26_ml)
residual.analysis(model = model_26_ml)


sc.AIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, model_65_ml, model_24_ml, model_26_ml)
sc.BIC=AIC(model_03_ml, model_04_ml, model_05_ml, model_06_ml, model_54_ml, model_55_ml, model_65_ml, model_24_ml, model_26_ml, k = log(length(waveReturn)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

# ARMA(0,6) model is the best one in this set. Overfitting models are ARMA(1,6) and ARMA(0,7)

model_16_css = arima(waveReturn,order=c(1,0,6),method='CSS')
coeftest(model_16_css)
residual.analysis(model = model_16_css)

model_16_ml = arima(waveReturn,order=c(1,0,6),method='ML')
coeftest(model_16_ml)
residual.analysis(model = model_16_ml)


model_07_css = arima(waveReturn,order=c(0,0,7),method='CSS')
coeftest(model_07_css)
residual.analysis(model = model_07_css)

model_07_ml = arima(waveReturn,order=c(0,0,7),method='ML')
coeftest(model_07_ml)
residual.analysis(model = model_07_ml)

# We will use ARMA(0,6) for the ARMA part of the model.

m06Residuals = model_06_ml$residuals
abs.res = abs(m06Residuals)
sq.res = m06Residuals^2

par(mfrow=c(1,2))
acf(abs.res, main="The ACF plot for absolute residual series")
pacf(abs.res, main="The PACF plot for absolute residual series")
par(mfrow=c(1,1))
#From ACF/PACF we get {ARMA(2,0)}.

eacf(abs.res)
# From the EACF, we can identify ARMA(3,1), ARMA(3,2) models for absolute residual series. 

# {ARMA(2,0), ARMA(3,1), ARMA(3,2)} models correspond to parameter settings:
# max(p,q) = 2 and q = 0 ==> max(p,0) = 2 and q = 0 ==> p can only be 2.
# max(p,q) = 3 and q = 1 ==> max(p,1) = 3 and q = 1 ==> p can only be 3.
# max(p,q) = 3 and q = 2 ==> max(p,2) = 3 and q = 2 ==> p can only be 3.
# The corresponding tentative GARCH models are GARCH(2,0), GARCH(3,1), GARCH(3,2).


par(mfrow=c(1,2))
acf(sq.res, main="The ACF plot for square residual series")
pacf(sq.res, main="The PACF plot for square residual series")
par(mfrow=c(1,1))
#From ACF/PACF we get {ARMA(3,2)}.

eacf(sq.res)
# From the EACF, we can identify ARMA(0,2), ARMA(0,3), and ARMA(1,3) models for squared residual series. 

# {ARMA(3,2), ARMA(0,2), ARMA(0,3), ARMA(1,3)} models correspond to parameter settings:
# max(p,q) = 3 and q = 2 ==> max(p,2) = 3 and q = 2 ==> p can only be 3.
# max(p,q) = 0 and q = 2 ==> max(p,2) = 0 and q = 1 ==> No models can be identified.
# max(p,q) = 0 and q = 3 ==> max(p,3) = 0 and q = 3 ==> No models can be identified.
# max(p,q) = 1 and q = 3 ==> max(p,3) = 1 and q = 3 ==> No models can be identified.
# The corresponding tentative GARCH model is GARCH(3,2) from the squared residuals.

# Overall, we have {GARCH(2,0), GARCH(3,1), GARCH(3,2)} models.

# When combined with the best model for the ARMA part, we have
# {ARMA(0,6)+GARCH(2,0), ARMA(0,6)+GARCH(3,1), ARMA(0,6)+GARCH(3,2)}

model1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0, 2)), 
                   mean.model = list(armaOrder = c(0, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.06_20<-ugarchfit(spec = model1, data = waveReturn, out.sample = 100)
m.06_20
residual.analysis(model=m.06_20, class = "ARMA-GARCH")


model2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 3)), 
                   mean.model = list(armaOrder = c(0, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.06_31<-ugarchfit(spec = model2, data = waveReturn, out.sample = 100)
m.06_31
residual.analysis(model=m.06_31, class = "ARMA-GARCH")


model3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(0, 6), include.mean = FALSE), 
                   distribution.model = "norm")
m.06_32<-ugarchfit(spec = model3, data = waveReturn, out.sample = 100)
m.06_32
residual.analysis(model=m.06_32, class = "ARMA-GARCH")


res.m.06_32 = m.06_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.06_32, main="The ACF plot for residual series")
pacf(res.m.06_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model4<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(2, 8), include.mean = FALSE), 
                   distribution.model = "norm")
m.28_32<-ugarchfit(spec = model4, data = waveReturn, out.sample = 100)
m.28_32
residual.analysis(model=m.28_32, class = "ARMA-GARCH")

res.m.28_32 = m.28_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.28_32, main="The ACF plot for residual series")
pacf(res.m.28_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model5<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(2, 9), include.mean = FALSE), 
                   distribution.model = "norm")
m.29_32<-ugarchfit(spec = model5, data = waveReturn, out.sample = 100)
m.29_32
residual.analysis(model=m.29_32, class = "ARMA-GARCH")

res.m.29_32 = m.29_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.29_32, main="The ACF plot for residual series")
pacf(res.m.29_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model6<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(3, 9), include.mean = FALSE), 
                   distribution.model = "norm")
m.39_32<-ugarchfit(spec = model6, data = waveReturn, out.sample = 100)
m.39_32
residual.analysis(model=m.39_32, class = "ARMA-GARCH")

res.m.39_32 = m.39_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.39_32, main="The ACF plot for residual series")
pacf(res.m.39_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model7<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(4, 9), include.mean = FALSE), 
                   distribution.model = "norm")
m.49_32<-ugarchfit(spec = model7, data = waveReturn, out.sample = 100)
m.49_32
residual.analysis(model=m.49_32, class = "ARMA-GARCH")

res.m.49_32 = m.49_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.49_32, main="The ACF plot for residual series")
pacf(res.m.49_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))


model8<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(5, 8), include.mean = FALSE), 
                   distribution.model = "norm")
m.58_32<-ugarchfit(spec = model8, data = waveReturn, out.sample = 100)
m.58_32
residual.analysis(model=m.58_32, class = "ARMA-GARCH")

res.m.58_32 = m.58_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.58_32, main="The ACF plot for residual series")
pacf(res.m.58_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))


model9<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(6, 8), include.mean = FALSE), 
                   distribution.model = "norm")
m.68_32<-ugarchfit(spec = model9, data = waveReturn, out.sample = 100)
m.68_32
residual.analysis(model=m.68_32, class = "ARMA-GARCH")

res.m.68_32 = m.68_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.68_32, main="The ACF plot for residual series")
pacf(res.m.68_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model10<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 3)), 
                   mean.model = list(armaOrder = c(7, 8), include.mean = FALSE), 
                   distribution.model = "norm")
m.78_32<-ugarchfit(spec = model10, data = waveReturn, out.sample = 100)
m.78_32
residual.analysis(model=m.78_32, class = "ARMA-GARCH") # Ljung-Box test gave worse results!

res.m.78_32 = m.78_32@fit$residuals
par(mfrow=c(1,2))
acf(res.m.78_32, main="The ACF plot for residual series")
pacf(res.m.78_32, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model11<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
                    mean.model = list(armaOrder = c(7, 8), include.mean = FALSE), 
                    distribution.model = "norm")
m.78_21<-ugarchfit(spec = model11, data = waveReturn, out.sample = 100)
m.78_21
residual.analysis(model=m.78_21, class = "ARMA-GARCH")

res.m.78_21 = m.78_21@fit$residuals
par(mfrow=c(1,2))
acf(res.m.78_21, main="The ACF plot for residual series")
pacf(res.m.78_21, main="The PACF plot for residual series")
par(mfrow=c(1,1))

# Final overfitting models
# ARMA(8,8) + GARCH(2,1) and ARMA(7,9) + GARCH(2,1)
# ARMA(7,8) + GARCH(2,2) and ARMA(7,8) + GARCH(3,1)

model12<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
                    mean.model = list(armaOrder = c(8, 8), include.mean = FALSE), 
                    distribution.model = "norm")
m.88_21<-ugarchfit(spec = model12, data = waveReturn, out.sample = 100)
m.88_21
residual.analysis(model=m.88_21, class = "ARMA-GARCH")

model13<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
                    mean.model = list(armaOrder = c(7, 9), include.mean = FALSE), 
                    distribution.model = "norm")
m.79_21<-ugarchfit(spec = model13, data = waveReturn, out.sample = 100)
m.79_21
residual.analysis(model=m.79_21, class = "ARMA-GARCH")

model14<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), 
                    mean.model = list(armaOrder = c(7, 8), include.mean = FALSE), 
                    distribution.model = "norm")
m.78_22<-ugarchfit(spec = model14, data = waveReturn, out.sample = 100)
m.78_22
residual.analysis(model=m.78_22, class = "ARMA-GARCH")

model15<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 3)), 
                    mean.model = list(armaOrder = c(7, 8), include.mean = FALSE), 
                    distribution.model = "norm")
m.78_31<-ugarchfit(spec = model15, data = waveReturn, out.sample = 100)
m.78_31
residual.analysis(model=m.78_31, class = "ARMA-GARCH")

# Final model that gives best diagnostics is ARMA(7,8) + GARCH(2,1)

AICs = c(infocriteria(m.06_31)[1], infocriteria(m.06_32)[1], infocriteria(m.28_32)[1], infocriteria(m.29_32)[1], 
         infocriteria(m.39_32)[1], infocriteria(m.49_32)[1], infocriteria(m.58_32)[1], infocriteria(m.68_32)[1], infocriteria(m.78_32)[1], 
         infocriteria(m.78_21)[1])
BICs = c(infocriteria(m.06_31)[2], infocriteria(m.06_32)[2], infocriteria(m.28_32)[2], infocriteria(m.29_32)[2], 
         infocriteria(m.39_32)[2], infocriteria(m.49_32)[2], infocriteria(m.58_32)[2], infocriteria(m.68_32)[2], infocriteria(m.78_32)[2], 
         infocriteria(m.78_21)[2])
df = data.frame(AIC = AICs, BIC = BICs)   
rownames(df) = c("m.06_31", "m.06_32", "m.28_32", "m.29_32", "m.39_32", "m.49_32", "m.58_32", "m.68_32", "m.78_32", 
                 "m.78_21")

df[order(df$AIC),] # See sorted by AIC
df[order(df$BIC),] # See sorted by BIC

par(mfrow=c(1,1))
plot(m.78_21)
plot(m.78_21, which = 1)

par(mfrow=c(1,1))
forc = ugarchforecast(m.78_21, data = waveReturn, n.ahead = 10)
plot(forc)


