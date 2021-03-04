# Load libraries
library(xts)
library(quantmod)
library(forecast)
library(tidyr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(fpp2)
library(tsbox)
library(TSstudio)
library(tseries)
library(forecastHybrid)
library(moments)
library(opera)

# Get S&P500 data
start_date <- as.Date("2009-01-01")
spy <- Ad(getSymbols("SPY", auto.assign = FALSE, from=start_date, warning=FALSE))
names(spy) <- "adjusted"
spy_ts <- ts_ts(spy)
spy_ts <- na.remove(spy_ts)


autoplot(spy$adjusted) +
  xlab("Date") +
  ylab("Adjusted closing price") +
  ggtitle("SPY")



cbind("Price" = spy_ts,
      "Log" = log(spy_ts),
      "Diff(log)" = diff(log(spy_ts))) %>%
  autoplot(facets=TRUE) +
  xlab("Year") + ylab("") +
  ggtitle("SPY")

# Box-Cox transformation to stabilize variance
lambda <- BoxCox.lambda(spy)
autoplot(BoxCox(spy, lambda))
spy_transformed <- BoxCox(spy, lambda)
spy_transformed %>% autoplot()
# Fit ARIMA
spy_transformed %>% diff() %>% ggtsdisplay(main="")
fit <- auto.arima(spy_transformed, seasonal=FALSE)
View(fit$fitted)
# Forecast
fit %>% forecast(h=500) %>% autoplot(include=3000)
checkresiduals(fit)

# Forecast combination


train <- window(spy_ts, end=c(2016,1))

h <- length(spy_ts) - length(train)


arima_comb <- forecast(auto.arima(train, biasadj=TRUE),h=h)
#arima_comb %>% autoplot()
stl_comb <- stlf(train, h=h, biasadj=TRUE, method = "arima")
#stl_comb %>% autoplot()

#combination <- (stl_comb[["mean"]] + arima_comb[["mean"]])/2


autoplot(spy_ts) +
  autolayer(stl_comb, series="STL", PI=TRUE) +
  #autolayer(arima_comb, series="ARIMA", PI=FALSE) +
  #autolayer(combination, series="Combination") +
  xlab("Year") + ylab("$") +
  ggtitle("SPY")


# Start from scratch

arima_comb <- forecast(auto.arima(train),h=h)

arima_acc <- arima_comb %>% accuracy(spy_ts)

stl_comb <- stlf(train, h=h)

stl_acc <- stl_comb %>% accuracy(spy_ts)

combination <- (stl_comb[["mean"]] + arima_comb[["mean"]])/2

comb_acc <- combination %>% accuracy(spy_ts)
comb_acc
stl_acc
arima_acc
print(paste("ARIMA accuracy is", arima_acc, "STL accuracy is", stl_acc, 
            "Combination accuracy is", comb_acc))

checkresiduals(arima_comb)


spy_ts %>%
  stl(t.window=13, s.window="periodic", robust=TRUE) %>%
  autoplot()

# Fit hybrid model
hybrid_fit <- hybridModel(spy_ts, models = "afns")
hmf <- forecast(hybrid_fit, h=300)
hmf %>%
  autoplot() +
  ggtitle("Equal weights - 300-day forecast")



hybrid_acc <- hmf %>% accuracy()
hybrid_acc

df <- cbind(spy_ts, error)
df <- as.data.frame(df)
df$predicted <- df$spy_ts - df$error
View(df)
class(df)

# Cross validation

f_hybrid <- function(x, h) {
  return(forecast(hybridModel(x, models="afns"), h=h))
}

error <- tsCV(spy_ts, f_hybrid, h=1, window = 400)
mse <- mean(error^2, na.rm = TRUE)
print(paste("mean squared error is", mse))
hist(error)
plot(error)


n_error <- tsCV(spy_ts, naive, h=1, window = 400)
View(n_error)
plot(n_error)

n_mse <- mean(n_error^2, na.rm=TRUE)
mse
n_ms












f_hw <- function(x, h) {
  return(forecast(HoltWinters(x, gamma = FALSE), h=h))
}

error <- tsCV(spy_ts, f_hw, h=1, window=400)
mse <- mean(error^2, na.rm = TRUE)
print(paste("mean squared error is", mse))


# check prediction intervals
verify_intervals <- function(ts=spy_ts, n_samples=1500, window_length=260*3, forecast_length=22) {
  total_80 <- 0
  total_95 <- 0
  for (i in 1:n_samples) {
    print(paste("now working on sample", i))
    upper <- length(ts) - window_length - forecast_length 
    start <- sample(1:upper, 1)
    end <- start + window_length - 1
    window <- subset(ts, start = start, end = end)
    fit <- HoltWinters(window, gamma = FALSE)
    fc <- forecast(fit, h=forecast_length)
    actual <- ts[(end+1):(end+forecast_length)]
    lower_80 <- fc$lower[, 1]
    upper_80 <- fc$upper[, 1]
    outside_80 <- sum(!(actual <= upper_80 & actual >= lower_80))
    total_80 <- total_80 + outside_80
    lower_95 <- fc$lower[, 2]
    upper_95 <- fc$upper[, 2]
    outside_95 <- sum(!(actual <= upper_95 & actual >= lower_95))
    total_95 <- total_95 + outside_95 
  }
  
  return(c(total_80/(n_samples*forecast_length), total_95/(n_samples*forecast_length)))
}

prop_outside <- verify_intervals()

# sample plot
ts=spy_ts
window_length=260*3
forecast_length=100
upper <- length(ts) - window_length - forecast_length 


start <- sample(1:upper, 1)
end <- start + window_length - 1
window <- subset(ts, start = start, end = end)
fit <- HoltWinters(window, gamma = FALSE)

fc <- forecast(fit, h=forecast_length)

autoplot(spy_ts) +
  autolayer(fc, series="HW", PI=TRUE) +
  xlab("Year") + ylab("$") +
  ggtitle("SPY")


get_error_distributions <- function(ts=spy_ts, n_samples=50, window_length=260*3, forecast_length=50) {
  error_mat <- matrix(0, nrow = forecast_length, ncol = 0)
  for (i in 1:n_samples) {
    print(paste("now working on sample", i))
    upper <- length(ts) - window_length - forecast_length 
    start <- sample(1:upper, 1)
    end <- start + window_length - 1
    window <- subset(ts, start = start, end = end)
    fit <- HoltWinters(window, gamma = FALSE)
    fc <- forecast(fit, h=forecast_length)
    actual <- ts[(end+1):(end+forecast_length)]
    predicted <- fc$mean
    errors <- as.vector(actual) - predicted
    errors <- as.vector(errors)
    error_mat <- cbind(error_mat, errors)
    
  }
  return(error_mat)
}

error_mat <- get_error_distributions(n_samples = 500)


means <- apply(error_mat, 1, mean)
vars <- apply(error_mat, 1, var)
skew <- apply(error_mat, 1, skewness)
kurt <- apply(error_mat, 1, kurtosis)
positive_prop <- apply(error_mat, 1, function(x) sum(x > 0) / length(x))

stats <- data.frame("mean" = means, 
                    "var" = vars, 
                    "skew" = skew, 
                    "kurt" = kurt,
                    "positive_proportion" = positive_prop)


# stat plots
ggplot(stats) +
  geom_point(aes(x=index(stats),  y=mean)) +
  xlab("")

ggplot(stats) +
  geom_point(aes(x=index(stats),  y=var)) +
  xlab("")

ggplot(stats) +
  geom_point(aes(x=index(stats),  y=skew)) +
  xlab("")

ggplot(stats) +
  geom_point(aes(x=index(stats),  y=kurt)) +
  xlab("")

ggplot(stats) +
  geom_point(aes(x=index(stats),  y=positive_prop)) +
  ylim(0, 1) +
  xlab("")

# histograms
hist(error_mat[1,], breaks=30)
hist(error_mat[5,], breaks=30)
hist(error_mat[10,], breaks=30)
hist(error_mat[20,], breaks=30)
hist(error_mat[50,], breaks=30)

print(paste("outside 80:", round((prop_outside[1]*100), 2), "outside 95:", 
            round((prop_outside[2]*100), 2)))


# Forecast combination
# arima, STL, holtwiners

train <- window(spy_ts, end=c(2016,1))
test <- window(spy_ts, start=c(2016,2))

h <- length(test)

HW <- forecast(HoltWinters(train, gamma=FALSE), h=h)
STL <- stlf(train, biasadj = TRUE, h=h, method = "arima")
ARIMA <- forecast(auto.arima(train), h=h)
X <- cbind(HW=HW$mean, STL=STL$mean, ARIMA=ARIMA$mean)
df <- cbind(spy_ts, X)
colnames(df) <- c("SPY" ,"HW", "STL", "ARIMA")
autoplot(df) +
  xlab("Year")

MLpol0 <- mixture(model="MLpol", loss.type="square")
weights <- predict(MLpol0, X, test, type="weights")
head(weights)
View(weights)


oracle.convex <- oracle(Y = test, experts = X, 
                        loss.type = "square", model = "convex")
plot(oracle.convex)

MLpol0 <- mixture(model="MLpol",loss.type = "square")

MLpol <- MLpol0

for (i in 1:length(test)) {
  MLpol <- predict(MLpol, newexperts=X[i, ], newY=test[i])
}

View(MLpol$weights)
summary(MLpol)
plot(MLpol, pause=TRUE)
class(weights)
2*weights + 4
autoplot(spy_ts) +
  autolayer(STL, series="STL", PI=FALSE) +
  autolayer(ARIMA, series="ARIMA", PI=FALSE) +
  autolayer(HW, series="HW", PI=FALSE) +
  xlab("Year") + ylab("$") +
  ggtitle("SPY")

