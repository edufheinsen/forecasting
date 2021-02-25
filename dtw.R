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
library(dtw)
library(liqueueR)
library(collections)
devAskNewPage(ask = FALSE)

# Get S&P500 data
start_date <- as.Date("2001-01-01")
spy <- Ad(getSymbols("SPY", auto.assign = FALSE, from=start_date, warning=FALSE))
names(spy) <- "adjusted"
spy_ts <- ts_ts(spy)
spy_ts <- na.remove(spy_ts)


# Steps
# 1 - take diff logs for stationarity
# 2 - walk forward find the top (n_similar) similar periods with appropriate 
# distancing to ensure actually different periods
# - return n_similar most similar periods, with their distance measures
# 3 - prediction options: (need to convert back to untransformed scale)
# - take some weighting using the ranks of the periods (recession paper) and
# get likelihood of passing some threshold
# - approach as regression problem - generate time series forecasts from the
# original periods and combine them

get_similar <- function(dl_ts=diff(log(spy_ts)), interval_length=30, skip_val=10, n_top_similar=10) {
  s_date_ref <- length(dl_ts) - interval_length + 1
  ref <- subset(dl_ts, start=s_date_ref)
  right_start_range <- s_date_ref - 2*interval_length
  q <- priority_queue()
  q_dist <- priority_queue()
  for (i in seq(1, right_start_range, skip_val)) {
    end <- i + interval_length - 1
    query <- subset(dl_ts, start = i, end = end)
    alignment <- dtw(query, ref, keep=TRUE)
    dist <- alignment$normalizedDistance
    q$push(query, priority = -dist)
    q_dist$push(dist, priority = -dist)
  }
  
  top_similar_intervals <- q$as_list()[1:n_top_similar]
  min_distances <- unlist(q_dist$as_list()[1:n_top_similar], use.names = FALSE)
 
  return(list(top_similar_intervals, min_distances, ref))
}

int <- ints[[1]]
end_int <- end(int)
window_length <- 22
dl_ts <- diff(log(spy_ts))
freq <- frequency(int)
start <- end_int+(1/freq)
end <- end_int + (window_length/freq)
window_after <- window(dl_ts, start=start, end = end)
window_after %>% autoplot()
threshold<- 0.01
outside_thresh <- sum(!(window_after <= threshold & window_after >= -threshold))
outside_thresh

# modify get_similar() to make periods non overlapping
# define a diff log threshold
# define a time window
# label the intervals with 0 and 1 depending on whether or not the threshold was crossed within the next time_window days
# use rank ordered centroids for the weights
# find probabilities with dot product

roc <- function(n) {
  x <- 1:n
  y <- 1/x
  
  z <- numeric(n)
  for (i in 1:n) {
    z[i] <- sum(y[i:length(y)])/n
  }
  return(z)
}



prob_exceeding <- function(similar_intervals, dl_ts=diff(log(spy_ts)), threshold=0.05, window_length=22) {
  exceed <- numeric(length(similar_intervals))
  for (i in 1:length(similar_intervals)) {
    interval <- similar_intervals[[i]]
    end_int <- end(interval)
    freq <- frequency(interval)
    start <- end_int+(1/freq)
    end <- end_int + (window_length/freq)
    window_after <- window(dl_ts, start=start, end = end)
    outside_thresh <- sum(!(window_after <= threshold & window_after >= -threshold))
    if (outside_thresh > 0) exceed[i] <- 1
  }
  roc <- roc(length(similar_intervals))
  prob_exceeds <- sum(roc * exceed)
  return(prob_exceeds)
}


ls <- get_similar()
ints <- ls[[1]]
prob_exceeding(ints, threshold=0.03)


numeric(5)

top_similar_intervals
min_distances


spy_ts
spy_ts %>% log() %>% diff() %>% autoplot()
