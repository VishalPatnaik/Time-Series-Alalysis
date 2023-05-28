
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

UNRATE <- read.csv("UNRAT.csv", header = TRUE)

UNRATE$DATE <- as.Date(UNRATE$DATE)

rownames(UNRATE) <- UNRATE$DATE
UNRATE <- within(UNRATE, rm(DATE))

# --- TASK 2 ---
wave = UNRATE
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

UNRATE_TS_1 <-  BoxCox.ar(wave + 1)
title(main = "Log-likelihood vs Lambda values.")

waveReturn = diff(log(wave))

plot(waveReturn,main="Time series plot of returns series for seismic Lg wave series")

adf.test(waveReturn)


par(mfrow=c(1,2))
acf(waveReturn, main="The ACF plot of returns series for seismic Lg wave series")
pacf(waveReturn, main="The PACF plot of returns series for seismic Lg wave series")
par(mfrow=c(1,1))
# Due to the changing variance, we do not have a clear picture in ACF and PACF.
# Referring very high autocorrelations {ARMA(0,1), ARMA(1,1)} can be identified

eacf(waveReturn)
#{ARMA(O,1), ARMA(0,2), ARMA(1,1), ARMA(1,2), ARMA(2,2), ARMA(3,2)}

res = armasubsets(y=waveReturn,nar=14,nma=14,y.name='test',ar.method='ols')
plot(res)
#Additional to the previously identified models, {ARMA(14,14)}

#Overall, 
#{ARMA(O,1), ARMA(0,2), ARMA(1,1), ARMA(1,2), ARMA(2,2), ARMA(3,2), ARMA(14,14)}












model_01_css = arima(UNRATE_diff,order=c(0, 0, 0),method='CSS')
coeftest(model_01_css)
residual.analysis(model = model_01_css)

model_01_ml = arima(UNRATE_diff,order=c(0,0,0),method='ML')
coeftest(model_01_ml)
residual.analysis(model = model_01_ml)


model_02_css = arima(UNRATE_diff,order=c(0,0,1),method='CSS')
coeftest(model_02_css)
residual.analysis(model = model_02_css)

model_02_ml = arima(UNRATE_diff,order=c(0,0,1),method='ML')
coeftest(model_02_ml)
residual.analysis(model = model_02_ml)

model_11_css = arima(UNRATE_diff,order=c(1,0,1),method='CSS')
coeftest(model_11_css)
residual.analysis(model = model_11_css)

model_11_ml = arima(UNRATE_diff,order=c(1, 0, 1),method='ML')
coeftest(model_11_ml)
residual.analysis(model = model_11_ml)

model_12_css = arima(UNRATE_diff,order=c(1,0,0),method='CSS')
coeftest(model_12_css)
residual.analysis(model = model_12_css)

model_12_ml = arima(UNRATE_diff,order=c(1,0,0),method='ML')
coeftest(model_12_ml)
residual.analysis(model = model_12_ml)


model_22_css = arima(UNRATE_diff,order=c(3,0,0),method='CSS')
coeftest(model_22_css)
residual.analysis(model = model_22_css)

model_22_ml = arima(UNRATE_diff,order=c(3,0,0),method='ML')
coeftest(model_22_ml)
residual.analysis(model = model_22_ml)


model_32_css = arima(UNRATE_diff,order=c(3,0,6),method='CSS')
coeftest(model_32_css)
residual.analysis(model = model_32_css)

model_32_ml = arima(UNRATE_diff,order=c(3,0,6),method='ML')
coeftest(model_32_ml)
residual.analysis(model = model_32_ml)

model_1414_css = arima(UNRATE_diff,order=c(3,0,9),method='CSS')
coeftest(model_1414_css)
residual.analysis(model = model_1414_css)

model_1414_ml = arima(UNRATE_diff,order=c(3,0,9),method='ML')
coeftest(model_1414_ml)
residual.analysis(model = model_1414_ml)


sc.AIC=AIC(model_01_ml, model_02_ml, model_11_ml, model_12_ml, model_22_ml, model_32_ml)
sc.BIC=AIC(model_01_css, model_02_css, model_11_css, model_22_ml, model_32_css,  k = (length(UNRATE_diff)))

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic")

# ARMA(3,3) model is the best one in this set. Overfitting models are ARMA(4,0) and ARMA(3,1)

model_11_css = arima(UNRATE_diff,order=c(4,0,0),method='CSS')
coeftest(model_11_css)
residual.analysis(model = model_11_css)

model_11_ml = arima(UNRATE_diff,order=c(4,0,0),method='ML')
coeftest(model_11_ml)
residual.analysis(model = model_11_ml)

model_01_css = arima(UNRATE_diff,order=c(3,0,1),method='CSS')
coeftest(model_01_css)
residual.analysis(model = model_01_css)

model_01_ml = arima(UNRATE_diff,order=c(3,0,1),method='ML')
coeftest(model_01_ml)
residual.analysis(model = model_01_ml)

model_02_css = arima(UNRATE_diff,order=c(3,0,0),method='CSS')
coeftest(model_02_css)
residual.analysis(model = model_02_css)

model_02_ml = arima(UNRATE_diff,order=c(3,0,0),method='ML')
coeftest(model_02_ml)
residual.analysis(model = model_02_ml)

# We will use ARMA(3,0) for the ARMA part of the model.

m11Residuals = model_02_ml$residuals
abs.res = m11Residuals
sq.res = m11Residuals^2

par(mfrow=c(1,2))
acf(abs.res, main="The ACF plot for absolute residual series")
pacf(abs.res, main="The PACF plot for absolute residual series")
par(mfrow=c(1,1))
#From ACF/PACF we get {ARMA(3,5)}.

eacf(abs.res)
# From the EACF, we can identify ARMA(3,3), ARMA(3,4), ARMA(4,3), ARMA(4,4), ARMA(3,5), ARMA(4,5), ARMA(5,5), ARMA(6,5) models for absolute residual series. 

# models correspond to parameter settings:
# max(p,q) = 3 and q = 3 ==> max(p,3) = 3 and q = 3 ==> p can only be 3.
# max(p,q) = 3 and q = 4 ==> max(p,1) = 3 and q = 4 ==> No models can be identified.
# max(p,q) = 4 and q = 3 ==> max(p,3) = 4 and q = 3 ==> p can only be 4.
# max(p,q) = 4 and q = 4 ==> max(p,4) = 4 and q = 4 ==> p can only be 4.
# max(p,q) = 3 and q = 5 ==> max(p,5) = 3 and q = 5 ==> No models can be identified.
# max(p,q) = 4 and q = 5 ==> max(p,5) = 4 and q = 5 ==> No models can be identified.
# max(p,q) = 5 and q = 5 ==> max(p,5) = 5 and q = 5 ==> p can only be 5.
# max(p,q) = 6 and q = 5 ==> max(p,5) = 6 and q = 5 ==> p can only be 6.
# The corresponding tentative GARCH models are GARCH(3,3), GARCH(4,3), GARCH(4,4), GARCH(5,5), GARCH(6,5) .


par(mfrow=c(1,2))
acf(sq.res, main="The ACF plot for square residual series")


pacf(sq.res, main="The PACF plot for square residual series")
par(mfrow=c(1,1))
#From ACF/PACF we get {ARMA(1,1)}.

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

model1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(3, 3)), 
                   mean.model = list(armaOrder = c(1, 1), include.mean = FALSE), 
                   distribution.model = "norm")
m.11_33<-ugarchfit(spec = model1, data = waveReturn, out.sample = 100)
m.11_33
residual.analysis(model=m.11_33, class = "ARMA-GARCH")


model2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(3, 4)), 
                   mean.model = list(armaOrder = c(1, 1), include.mean = FALSE), 
                   distribution.model = "norm")
m.11_43<-ugarchfit(spec = model2, data = waveReturn, out.sample = 100)
m.11_43
residual.analysis(model=m.11_43, class = "ARMA-GARCH")


model3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 4)), 
                   mean.model = list(armaOrder = c(1, 1), include.mean = FALSE), 
                   distribution.model = "norm")
m.11_44<-ugarchfit(spec = model3, data = waveReturn, out.sample = 100)
m.11_44
residual.analysis(model=m.11_44, class = "ARMA-GARCH")


model4<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(5, 5)), 
                   mean.model = list(armaOrder = c(1, 1), include.mean = FALSE), 
                   distribution.model = "norm")
m.11_55<-ugarchfit(spec = model4, data = waveReturn, out.sample = 100)
m.11_55
residual.analysis(model=m.11_55, class = "ARMA-GARCH")


model5<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(6, 5)), 
                   mean.model = list(armaOrder = c(1, 1), include.mean = FALSE), 
                   distribution.model = "norm")
m.11_65<-ugarchfit(spec = model5, data = waveReturn, out.sample = 100)
m.11_65
residual.analysis(model=m.11_65, class = "ARMA-GARCH")



res.m.11_44 = m.11_44@fit$residuals
par(mfrow=c(1,2))
acf(res.m.11_44, main="The ACF plot for residual series")
pacf(res.m.11_44, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model4<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 4)), 
                   mean.model = list(armaOrder = c(1, 2), include.mean = FALSE), 
                   distribution.model = "norm")
m.12_44<-ugarchfit(spec = model4, data = waveReturn, out.sample = 100)
m.12_44
residual.analysis(model=m.12_44, class = "ARMA-GARCH")

res.m.12_44 = m.12_44@fit$residuals
par(mfrow=c(1,2))
acf(res.m.12_44, main="The ACF plot for residual series")
pacf(res.m.12_44, main="The PACF plot for residual series")
par(mfrow=c(1,1))

model5<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 4)), 
                   mean.model = list(armaOrder = c(2, 2), include.mean = FALSE), 
                   distribution.model = "norm")
m.22_44<-ugarchfit(spec = model5, data = waveReturn, out.sample = 100)
m.22_44
residual.analysis(model=m.22_44, class = "ARMA-GARCH")

res.m.22_44 = m.22_44@fit$residuals
par(mfrow=c(1,2))
acf(res.m.22_44, main="The ACF plot for residual series")
pacf(res.m.22_44, main="The PACF plot for residual series")
par(mfrow=c(1,1))


# Final model that gives best diagnostics is ARMA(7,8) + GARCH(2,1)

AICs = c( infocriteria(m.12_44)[1], infocriteria(m.22_44)[1])
BICs = c(infocriteria(m.12_44)[2], infocriteria(m.22_44)[2])
df = data.frame(AIC = AICs, BIC = BICs)
rownames(df) = c("m.12_44", "m.22_44")

df[order(df$AIC),] # See sorted by AIC
df[order(df$BIC),] # See sorted by BIC

par(mfrow=c(1,1))
plot(m.22_44)
plot(m.22_44, which = 1)

par(mfrow=c(1,1))
forc = ugarchforecast(m.22_44, data = waveReturn, n.ahead = 10)
plot(forc)

