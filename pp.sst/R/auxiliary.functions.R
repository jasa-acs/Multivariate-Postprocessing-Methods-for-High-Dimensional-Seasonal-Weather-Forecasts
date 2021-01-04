
#' crps score for normal distribution, resistant to missing values and 0 standard deviation  
#'
#' @param y vector of observations.
#' @param mean,sd mean and sd of the forecasting distribution
#'                   
#' @return vector of the same length as y containing crps scores.
#'
#' @author Claudio Heinrich
#' 
#' @importFrom scoringRules crps
#' 
#' @export 


crps_na_rm = function(y, mean, sd)
  {
    na_loc = which( is.na(y) | is.na(mean) | is.na(sd) | sd == 0)
    
    if(identical(na_loc,integer(0)))
    {
      x = scoringRules::crps(y, family = "normal", mean = mean, sd = sd)
    }else
    {
      x = rep(0,length(y))
      x[-na_loc] = scoringRules::crps(y[-na_loc], family = "normal", mean = mean[-na_loc], sd = sd[-na_loc])
      x[na_loc] = NA
    }
    return(x)
  }


#'@export

GneitingWeightFct = function(x,L)
  {
    
    ret_value = rep(0,length(x))
    
    ret_value[x == 0] = 1
    
    t = x[x != 0 & x <= L]/L
    ret_value[x != 0 & x <= L] = (1-t)*sin(2*pi*t)/(2*pi*t)+(1-cos(2*pi*t))/(2*pi^2*t)   
    
    ret_value[ret_value < 0] = 0 # for rounding errors
    return(ret_value)
  }





#' truncation function for freezing level of sea water
#' 
#' @author Claudio Heinrich
#' 
#' 
#' @export

trc = function (x,truncation.value = -1.79){ 
  x = truncation.value * (x < truncation.value) + x * (x >= truncation.value)
  return(x)}



#' Run a permutation test of the pairwise difference between two vectors of numbers
#' @param a Vector, the scores from one method
#' @param b Vector, the scores from some other method
#' @param N Integer, the size of the permutation distribution
#' @return A list with the mean of the difference and the permutation distribution of that difference and the p-value
#' @examples
#' N = 1e2
#' trend  = 1:N
#' a = trend + .01 + rnorm(N, .001)
#' b = trend - .01 + rnorm(N, .001)
#' l = permutation_test_difference(a,b)
#' q = sum(l$D <= l$d_bar) / length(l$D)
#' @author Alex,Claudio
#' 
#' @export
#' 
permutation_test_difference = function(a,
                                       b,
                                       N = 5e3){
    n = length(a)
    d = a - b
    d_bar = mean(d)
    D = NULL
    for(i in 1:N){
        swap = rbinom(n,1,0.5)
        w_swap = which(swap == 1)
        d_i = d
        d_i[w_swap] = -d_i[w_swap]
        D[i] = mean(d_i)
    }
    
    p_val = sum(d_bar > sort(D))/N + sum(d_bar == sort(D))/(2*N)
    
    return(list(d_bar = d_bar, D = D,p_val = p_val))
}


#' Run a bootstrap test of the pairwise difference between two vectors of numbers
#' @param sc1 Vector, the scores from one method
#' @param sc2 Vector, the scores from some other method
#' @param N Integer, the size of the bootstrap distribution
#' @param q_prob a vector containing probabilities for which the quantiles should be returned from the bootstrap distribution
#' @return A list as returned by the function boot::boot, with additional values p_val and q
#' @examples
#' N = 1e2
#' trend  = 1:N
#' sc1 = trend + .01 + rnorm(N, .001)
#' sc2 = trend - .01 + rnorm(N, .001)
#' l = bootstrap_difference(sc1,sc2,q_prob = c(0.05,0.95))
#' @author Claudio
#' 
#' @importFrom boot boot
#' 
#' @export
#' 
bootstrap_difference = function(sc1,sc2,q_prob, N = 1e3){
  
  data = sc1-sc2
  
  diff = function(x,inds){return(mean(x[inds]))}
  
  bt = boot::boot(data,diff,R = N)
  
  t0 = bt$t0
  t = bt$t
  
  #p-value:
  bt$p_val = sum(t0 > t)/N + sum(t0 == t)/(2*N)
  
  #quantiles
  bt$q = stats::quantile(t,q_prob)
  
  
  return(bt)
}
