
##################
#### Overview ####
##################

## This code reproduces, in R, the averaged-ranks case in Fig. S2 of the 
## following paper:

## K.A. McKinnon, I.R. Simpson, A.P. Williams, 
## The pace of change of summertime temperature extremes, 
## Proc. Natl. Acad. Sci. U.S.A.
## 121 (42) e2406143121,
## https://doi.org/10.1073/pnas.2406143121 (2024).

## The code produced and shared by the authors, "synthetic_data_tests.py", 
## was used as a reference when reproducing the code.


## The last section of this code then tries a different test-case for the 
## averaged-ranks case: instead of prescribing the variance to increase 
## linearly over the 65 years, it instead enforces the variance
## to increase exponentially over the 65 years. 
## The key steps are as follows:
# var_seq <- seq(0, 0.5, 0.05) ## Changes in the variance tried
# exp_seq <- c(0, exp(1:(years-1)) / exp((years-1))) ## Sequence of exponentially
                                                     ## spaced-apart values
                                                     ## from 0 to 1.
# expo_sd <- sqrt(1 + exp_seq * var_seq[i]) ## Standard deviation used in the
                                            ## normal distribution to produce
                                            ## the data.

## Using this exponentially-spaced sequence of values for the variance, 
## the power of the test weakens. 



## This code is done by Teng Wei Yeo, 11 Nov 2024.




##############################
#### Reproducing Fig. S2. ####
##############################

samples = 100
N = 50 # no. of spatially independent locations
days = 91
years = 65


set.seed(1234)

mem_avg = rep(0, samples)

for (i in 1:samples){
  long <- rnorm(years*days*N, 0, 1)
  arr <- array(long, dim = c(years, days, N))
  medians <- apply(arr, MARGIN = c(1,3), median) # years x N matrix
  medians_mat <- medians[, rep(1:N, each = days)] # duplicate median for each day 
  medians_arr <- array(medians_mat, dim = c(years,days,N)) # turn back into array
  anomalies <- arr - medians_arr
  max_anoms <- apply(anomalies, MARGIN = c(1,3), max) # years x N matrix
  ranks <- apply(max_anoms, MARGIN = 2, rank) # years x N matrix
  means <- apply(ranks, MARGIN = 1, mean) # 1 x years vector
  year_seq <- seq(1, years, by = 1)
  mem_avg[i] <- coef(lm(means~year_seq))["year_seq"]
}

mem_avg

rejection <- quantile(mem_avg, probs = c(0.025, 0.975))



mean_seq <- seq(0,2,0.25)

data_samples <- 100
is_sig_mat <- matrix(0, nrow = data_samples, ncol = length(mean_seq))

for (i in 1:length(mean_seq)){
  print(paste0("Doing ", mean_seq[i], " increase to mean"))
  for (j in 1:data_samples){
    linearised <- seq(0, mean_seq[i], length.out = years)
    linearised_cntr <- scale(linearised, center = TRUE, scale = FALSE)
    data <- rnorm(years*days*N, linearised_cntr, sd = 1)
    
    arr <- array(data, dim = c(years, days, N))
    medians <- apply(arr, MARGIN = c(1,3), median) # years x N matrix
    medians_mat <- medians[, rep(1:N, each = days)] # same median for each day 
    medians_arr <- array(medians_mat, dim = c(years,days,N))
    anomalies <- arr - medians_arr
    max_anoms <- apply(anomalies, MARGIN = c(1,3), max) # years x N matrix
    ranks <- apply(max_anoms, MARGIN = 2, rank) # years x N matrix
    means <- apply(ranks, MARGIN = 1, mean) # average ranks
    year_seq <- seq(1, years, by = 1)
    beta <- coef(lm(means~year_seq))["year_seq"]
    
    is_sig_mat[j,i] = (beta < rejection[1] | beta > rejection[2]) # TRUE/FALSE
  }
}

proportions <- apply(is_sig_mat, MARGIN=2, sum) / data_samples
proportions

var_seq <- seq(0, 0.5, 0.05)
is_sig_mat_var <- matrix(0, nrow = data_samples, ncol = length(var_seq))


for (i in 1:length(var_seq)){
  print(paste0("Doing +", var_seq[i]*100, "% change to variance"))
  for (j in 1:data_samples){
    linearised_sd <- sqrt(seq(1, 1 + var_seq[i], length.out = years))
    data <- rnorm(years*days*N, mean = 0, sd = linearised_sd) # recycling
    
    arr <- array(data, dim = c(years, days, N))
    medians <- apply(arr, MARGIN = c(1,3), median) # years x N matrix
    medians_mat <- medians[, rep(1:N, each = days)] # same median for each day 
    medians_arr <- array(medians_mat, dim = c(years,days,N))
    anomalies <- arr - medians_arr
    max_anoms <- apply(anomalies, MARGIN = c(1,3), max) # years x N matrix
    ranks <- apply(max_anoms, MARGIN = 2, rank) # years x N matrix
    means <- apply(ranks, MARGIN = 1, mean)
    year_seq <- seq(1, years, by = 1)
    beta <- coef(lm(means~year_seq))["year_seq"]
    
    is_sig_mat_var[j,i] = (beta < rejection[1] | beta > rejection[2]) #TRUE/FALSE
  } 
}

proportions_var <- apply(is_sig_mat_var, MARGIN=2, sum) / data_samples
proportions_var

## Visualisation

library(ggplot2)
df_means <- data.frame(rej_probab_mean = proportions,
                       mean_change = mean_seq)

ggplot(data = df_means, aes(x = mean_change, y = rej_probab_mean))+
  geom_point()+
  geom_line() +
  labs(y = "Rejection probability",
       x = "Additive increase in mean after 65 years",
       title = "Linearly increasing mean over 65 years") +
  ylim(0, 1)


df_var <- data.frame(rej_probab_var = proportions_var,
                     var_change = var_seq,
                     type = "Linearly increasing")

ggplot(data = df_var, aes(x = var_change, y = rej_probab_var))+
  geom_point()+
  geom_line() +
  labs(y = "Rejection probability",
       x = "Percent increase in variance after 65 years",
       title = "Linearly increasing variance over 65 years") +
  scale_x_continuous(labels = scales::percent) + 
  ylim(0,1)




###########################################
#### Quadratically increasing variance ####
###########################################

## This checks whether the test can detect when variance is increasing 
## quadratically over time (instead of linearly).

quad_seq <- seq(0, 1, length.out = years)^2

for (i in 1:length(var_seq)){
  print(paste0("Doing +", var_seq[i]*100, "% change to variance"))
  for (j in 1:data_samples){
    quad_sd <- sqrt(1 + quad_seq * var_seq[i])
    data <- rnorm(years*days*N, mean = 0, sd = quad_sd) # recycling
    
    arr <- array(data, dim = c(years, days, N))
    medians <- apply(arr, MARGIN = c(1,3), median) # years x N matrix
    medians_mat <- medians[, rep(1:N, each = days)] # same median for each day 
    medians_arr <- array(medians_mat, dim = c(years,days,N))
    anomalies <- arr - medians_arr
    max_anoms <- apply(anomalies, MARGIN = c(1,3), max) # years x N matrix
    ranks <- apply(max_anoms, MARGIN = 2, rank) # years x N matrix
    means <- apply(ranks, MARGIN = 1, mean)
    year_seq <- seq(1, years, by = 1)
    beta <- coef(lm(means~year_seq))["year_seq"]
    
    is_sig_mat_var[j,i] = (beta < rejection[1] | beta > rejection[2])#TRUE/FALSE
  } 
}

proportions_var2 <- apply(is_sig_mat_var, MARGIN=2, sum) / data_samples
temp <- data.frame(rej_probab_var = proportions_var2,
                   var_change = var_seq,
                   type = "Quadratically increasing")
df_var2 <- rbind(df_var, temp)
                

ggplot(data = df_var2, aes(x = var_change, 
                           y = rej_probab_var, 
                           colour = type))+
  geom_point()+
  geom_line() +
  labs(y = "Rejection probability",
       x = "Percent increase in variance after 65 years",
       title = paste0("Probability of rejection for different increases in\n",
                      "variance, at different rates over 65 years")) +
  scale_x_continuous(labels = scales::percent) + 
  ylim(0,1)



###########################################
#### Exponentially increasing variance ####
###########################################

## This checks whether the test can detect when variance is increasing 
## exponentially over time (instead of linearly).
## This is relevant in cases where the variance stays somewhat constant for
## many years, then increases a lot in later years.

## Arbitrary exponentially increasing sequence:
exp_seq <- c(0, exp(1:(years-1)) / exp((years-1)))


for (i in 1:length(var_seq)){
  print(paste0("Doing +", var_seq[i]*100, "% change to variance"))
  for (j in 1:data_samples){
    expo_sd <- sqrt(1 + exp_seq * var_seq[i])
    data <- rnorm(years*days*N, mean = 0, sd = expo_sd) # recycling
    
    arr <- array(data, dim = c(years, days, N))
    medians <- apply(arr, MARGIN = c(1,3), median) # years x N matrix
    medians_mat <- medians[, rep(1:N, each = days)] # same median for each day 
    medians_arr <- array(medians_mat, dim = c(years,days,N))
    anomalies <- arr - medians_arr
    max_anoms <- apply(anomalies, MARGIN = c(1,3), max) # years x N matrix
    ranks <- apply(max_anoms, MARGIN = 2, rank) # years x N matrix
    means <- apply(ranks, MARGIN = 1, mean)
    year_seq <- seq(1, years, by = 1)
    beta <- coef(lm(means~year_seq))["year_seq"]
    
    is_sig_mat_var[j,i] = (beta < rejection[1] | beta > rejection[2])#TRUE/FALSE
  } 
}

proportions_var3 <- apply(is_sig_mat_var, MARGIN=2, sum) / data_samples

temp2 <- data.frame(rej_probab_var = proportions_var3,
                   var_change = var_seq,
                   type = "Exponentially increasing")
df_var3 <- rbind(df_var2, temp2)

df_var3$type <- 
  factor(df_var3$type, levels = c("Linearly increasing",
                                  "Quadratically increasing",
                                  "Exponentially increasing"))

ggplot(data = df_var3, aes(x = var_change, 
                           y = rej_probab_var, 
                           colour = type))+
  geom_point()+
  geom_line() +
  labs(y = "Rejection probability",
       x = "Percent increase in variance after 65 years",
       title = paste0("Probability of rejection for different increases in\n",
                      "variance, at different rates over 65 years")) +
  scale_x_continuous(labels = scales::percent) + 
  ylim(0,1)


## The test is weaker for the exponential case.


## This is in part due to the nature of the test in being successful at 
## accounting for the presence of a change in variance, but not for 
## accounting the magnitude of the change in variance. 

## Since the `ramp' of exp() is quite steep (i.e., it takes until the 59th 
## year for exp_seq > 0.001), there is only a small change in variance
## for most of the years. The variance increases by a lot only in the last
## few years. 

exp_seq


