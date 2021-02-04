# Load libraries
library(xts)
library(quantmod)
library(forecast)
library(tidyr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(fpp2)
library(TSstudio)

# Get S&P500 data
start_date <- as.Date("2005-01-01")
spy <- Ad(getSymbols("SPY", auto.assign = FALSE, from=start_date, warning=FALSE))
names(spy) <- "adjusted"
autoplot(spy$adjusted) +
  xlab("Date") +
  ylab("Adjusted closing price") +
  ggtitle("SPY")

spy_ts <- xts_to_ts(spy)


cbind("Price" = spy_ts,
      "Log" = log(spy_ts),
      "Diff(log)" = diff(log(spy_ts))) %>%
  autoplot(facets=TRUE) +
  xlab("Year") + ylab("") +
  ggtitle("SPY")

spy %>% autoplot()

# Box-Cox transformation to stabilize variance
lambda <- BoxCox.lambda(spy)
autoplot(BoxCox(spy, lambda))
spy_transformed <- BoxCox(spy, lambda)
spy_transformed %>% autoplot()
# Fit ARIMA
spy_transformed %>% diff() %>% ggtsdisplay(main="")
fit <- auto.arima(spy_transformed, seasonal=FALSE)

# Forecast
fit %>% forecast(h=150) %>% autoplot(include=3000)
checkresiduals(fit)


# Forecast combination

train <- window(spy_ts, end=c(2014, 12))
h <- length(spy_ts) - length(train)

ETS <- forecast(ets(train), h=h)
ets(train)


stl_comb <- stlf(train, lambda=0, h=h, biasadj=TRUE)

arima_comb <- forecast(auto.arima(train, lambda=0, biasadj=TRUE),h=h)

combination <- (ets_comb[["mean"]] + arima_comb[["mean"]])/2

autoplot(spy_ts) +
  autolayer(stl_comb, series="STL", PI=FALSE) +
  autolayer(arima_comb, series="ARIMA", PI=FALSE) +
  autolayer(combination, series="Combination") +
  xlab("Year") + ylab("$") +
  ggtitle("SPY")

View(ARIMA)
class(auscafe)
class(spy_ts)
class(arima_comb)

ETS
train <- window(auscafe, end=c(2015,9))
h <- length(auscafe) - length(train)
ETS <- forecast(ets(train), h=h)
ARIMA <- forecast(auto.arima(train, lambda=0, biasadj=TRUE),
                  h=h)
STL <- stlf(train, lambda=0, h=h, biasadj=TRUE)
NNAR <- forecast(nnetar(train), h=h)
TBATS <- forecast(tbats(train, biasadj=TRUE), h=h)
Combination <- (ETS[["mean"]] + ARIMA[["mean"]] +
                  STL[["mean"]] + NNAR[["mean"]] + TBATS[["mean"]])/5
autoplot(auscafe) +
  autolayer(ETS, series="ETS", PI=FALSE) +
  autolayer(ARIMA, series="ARIMA", PI=FALSE) +
  autolayer(STL, series="STL", PI=FALSE) +
  autolayer(NNAR, series="NNAR", PI=FALSE) +
  autolayer(TBATS, series="TBATS", PI=FALSE) +
  autolayer(Combination, series="Combination") +
  xlab("Year") + ylab("$ billion") +
  ggtitle("Australian monthly expenditure on eating out")