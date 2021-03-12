#Load libraries
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
library(calculus)
devAskNewPage(ask = FALSE)

# Get S&P500 data
start_date <- as.Date("2009-01-01")
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

get_similar <- function(dl_ts=diff(log(spy_ts)), interval_length=45, skip_val=45, n_top_similar=20) {
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
ref <- ls[[3]]




p <- autoplot(spy_ts)

for (i in 1:length(ints)) {
  interval <- ints[[i]]
  start <- start(interval)
  end <- end(interval)
  window <- window(spy_ts, start=start, end=end)
  p <- p + autolayer(window, colour = TRUE, series = paste(i))
}

p <- p + autolayer(window(spy_ts, start=start(ref), end=end(ref)), colour = TRUE, series="Current")

print(p)

top_int <- ints[[1]]
top_int

s <- start(top_int)
e <- end(top_int)

# dollar progress (during similar ints)
for (interval in ints) {
  curr <- numeric(length(interval) + 1)
  curr[1] <- 1
  for (i in 1:length(interval)) {
    curr[i + 1] <- curr[i] * (1 + interval[i])
  }
  print(curr[length(curr)])
}

# dollar progress (out of sample)
dollar_progress <- function(similar_intervals, dl_ts=diff(log(spy_ts)), threshold=0.05, window_length=22) {
  dollar_performance <- matrix(nrow=(window_length+1), ncol=0)
  for (i in 1:length(similar_intervals)) {
    interval <- similar_intervals[[i]]
    end_int <- end(interval)
    freq <- frequency(interval)
    start <- end_int+(1/freq)
    end <- end_int + (window_length/freq)
    window_after <- window(dl_ts, start=start, end = end)
    
    curr <- numeric(length(window_after) + 1)
    curr[1] <- 1
    for (i in 1:length(window_after)) {
      curr[i + 1] <- curr[i] * (1 + window_after[i])
    }
    
    dollar_performance <- cbind(dollar_performance, curr)
  }
  colnames(dollar_performance) <- 1:20
  return(dollar_performance)
}


returns_mat <- dollar_progress(ints)
matplot(returns_mat, type="l")

ints


mean(returns_mat[nrow(returns_mat), ])
sd(returns_mat[nrow(returns_mat), ])
# forecast with similar
# for each of the top similar intervals - run opera for the duration of the period,
# come up with weights - weight the weights with roc?
spy_ts
p_m_forecast <- function(similar_intervals, curr_interval, test_len=22, num_models=3, ts=spy_ts) {
  roc <- roc(length(similar_intervals))
  weight_mat <- matrix(0, nrow=(test_len), ncol=num_models)
  
  for (i in 1:length(similar_intervals)) {
    interval <- similar_intervals[[i]]
    end_int = end(interval)
    freq <- frequency(interval)
    start <- end_int+(1/freq)
    end <- end_int + (test_len/freq)
    test <- window(ts, start=start, end = end)
    h <- length(test)
    
    HW <- forecast(HoltWinters(train, gamma=FALSE), h=h)
    ETS <- forecast(ets(train), h=h)
    ARIMA <- forecast(auto.arima(train), h=h)
    X <- cbind(HW=HW$mean, ETS=ETS$mean, ARIMA=ARIMA$mean)
    df <- cbind(spy_ts, X)
    colnames(df) <- c("SPY" ,"HW", "ETS", "ARIMA")
    MLpol0 <- mixture(model="MLpol", loss.type="square")
    weights <- predict(MLpol0, X, test, type="weights")
    weight_mat <- weight_mat + weights * roc[i] 
  }
  return(weight_mat)
}
train_len <- 45
h <- 22
train <- subset(spy_ts, start = length(spy_ts) - train_len + 1)



arima_fc <- forecast(auto.arima(train, biasadj=TRUE),h=h)$mean
hw_fc <- forecast(HoltWinters(train, gamma=FALSE), h=h)$mean
ets_fc <- forecast(ets(train), h=h)$mean

fc_mat <- cbind(hw_fc, ets_fc, arima_fc)

w <- p_m_forecast(ints) 

View(fc_mat)
View(w)
ncol(w)
nrow(fc_mat)
ncol(fc_mat)

weighted_fc <- numeric(nrow(fc_mat))
fc_mat
w
for (i in 1:nrow(fc_mat)) {
  weighted_fc[i] <- sum(w[i,]*fc_mat[i,])
}

plot(weighted_fc)
weighted_fc


# do this for rolling x-day windows treated as "recent" - make sure all have 
# past out sample, only use past observations for each 
# get N top similar past intervals
# get NEXT_WINDOW-length periods, for each of these periods
# for ALL intervals (recent, top similar, next-windows following top similar), get
# dollar performance
# variance of returns
# skewness and kurtosis of return distribution
#
#
#
#

end_int <- end(interval)
freq <- frequency(interval)
start <- end_int+(1/freq)
end <- end_int + (window_length/freq)
window_after <- window(dl_ts, start=start, end = end)

curr <- numeric(length(window_after) + 1)


get_rolling_windows <- function(time_series, range_start, range_end, window_length, dist_between_windows) {
  rolling_windows <- list()
  freq <- frequency(time_series)
  start_int <- range_start[1] + range_start[2]/freq
  end_int <- start_int + window_length/freq
  rolling_windows[[1]] <- window(time_series, start=start_int, end=end_int)
  range_end_index <- range_end[1] + range_end[2]/freq
  i <- 1
  while (end_int < range_end_index) {
    i <- i + 1
    start_int <- end_int + dist_between_windows/freq
    end_int <- start_int + window_length/freq 
    new_int <- window(time_series, start=start_int, end=end_int)
    rolling_windows[[i]] <- window(time_series, start=start_int, end=end_int)
  }
  return(rolling_windows)
}
k <- get_rolling_windows(spy_ts, c(2016, 1), c(2020, 5), 22, 22)

x <- 1
while (x < 5) {
  x <- x + 1
  print("ed")
}

l <- list(numeric)
l[[1]] <- 5
l <- append(l, 2)

l
s <- floor(start(spy_ts))

f <- frequency(spy_ts)

s <- c(2009, 2)

w <- window(spy_ts, start=s, end=s+2)

w

k <- window(spy_ts, start=(2009+1/freq))

length(k)
length(w)

start(spy_ts)

class(w)

# call as
# return
# list of window indices - correspond to (2016, 1, 1)-(2016, 2, 1), etc


get_rolling_range <- function(time_series, query_window) {
  start <- start(time_series)
  end <- start(query_window) - 2 * length(query_window) / frequency(time_series)
  return(c(start, end))
}

ex_rolling <- k[[1]]

length(ex_rolling)

vec <- get_rolling_range(spy_ts, ex_rolling)


get_top_similar_intervals <- function(time_series_original, query_window, start_range, end_range, num_top_intervals) {
  time_series <- diff(log(time_series_original))
  interval_length <- length(query_window)
  
  q <- priority_queue()
  q_dist <- priority_queue()
  
  freq <- frequency(time_series)
  start_int <- start_range
  end_int <- start_int + interval_length/freq
  curr_window <- window(time_series, start=start_int, end=end_int)
  alignment <- dtw(curr_window, query_window, keep=TRUE)
  dist <- alignment$normalizedDistance
  q$push(curr_window, priority = -dist)
  q_dist$push(dist, priority = -dist)

  while (end_int < end_range) {
    start_int <- end_int + interval_length/freq
    end_int <- start_int + interval_length/freq 
    curr_window <- window(time_series, start=start_int, end=end_int)
    alignment <- dtw(curr_window, query_window, keep=TRUE)
    dist <- alignment$normalizedDistance
    q$push(curr_window, priority = -dist)
    q_dist$push(dist, priority = -dist)
  
  }
  
  top_similar_intervals <- q$as_list()[1:num_top_intervals]
  min_distances <- unlist(q_dist$as_list()[1:num_top_intervals], use.names = FALSE)
  
  return(list(top_similar_intervals, min_distances, query_window))
}

sim <- get_top_similar_intervals(spy_ts, ex_rolling, vec[1], vec[2], 20)
sim[[1]]


get_following_windows <- function(time_series, intervals) {
  following_windows <- list()
  for (interval in intervals) {
    start <- end(interval) + 1/freq
    end <- start + length(interval)/freq
    following <- window(time_series, start=start, end=end)
    append(following_windows, following)
  }
  return(following_windows)
}


get_returns <- function(window) {
  
}

compute_summary_statistics <- function(returns) {
  
}


# pseudocode

# what is interesting
# can direction be predicted? (if 15/20 ended up at positive, will this one end up as positive too)
# can dispersion be predicted (if past similar were volatile, will returns be volatile right after this one too)

rolling_windows <- get_rolling_windows()

for (qry_window in rolling_windows) {
  following_window_actual <- get_following(qry_window)
  c(start_range, end_range) <- get_rolling_range(qry_window)
  similar_ints <- get_top_similar_intervals(qry_window, start_range, end_range)
  following_windows_similar <- get_following(windows(similar_ints))
  # 1 direction
  # dollar progress, weighted by similarity, of "following windows similar"
  # dollar progress of following window actual
  # save both in list to graph/visualize later
  # save single values as data points to plot, calculate correlation, etc
  # 2 volatility
  # same as above with volatility instead of dollar progress
}















