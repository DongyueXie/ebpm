nlm_control_defaults <- function() {
  return(list(ndigit = 8, stepmax = 30, check.analyticals = FALSE))
}

mixsqp_control_defaults <- function() {
  return(list(verbose = F))
}



geom_seq <- function(low, up, m){
  N =  ceiling((log(up) - log(low))/log(m)) + 1
  #if(is.infinite(N)){browser()}
  # while(N > 500){
  #   m = m*(2^0.75)
  #   N =  ceiling((log(up) - log(low))/log(m)) + 1
  # }
  out  = low*m^(seq(1,N, by = 1)-1)
  return(out)
}

# lin_seq <- function(low, up, m){
#   out = seq(low, up, length.out = m)
#   return(out)
# }



#' @export get_uniform_mixture
get_uniform_mixture <- function(x, s, grid_res = NULL, m = 2, low = NULL){
  if(is.null(grid_res)){
    grid_res = select_grid_exponential(x = x, s = s, m = m, low = low)
  }
  shape = grid_res$shape
  scale = grid_res$scale
  n = length(shape)
  pi = replicate(n, 1/n)
  g = gammamix(pi = pi, shape = shape, scale = scale)
  return(g)
}

## compute L matrix from data and selected grid
## L_ik = NB(x_i; a_k, b_k/b_k + s_i)
## but for computation in mixsqr, we can simplyfy it for numerical stability
#' @export compute_L
compute_L <- function(x, s, a, b){
  prob = 1 - s/outer(s,b, "+")
  l = dnbinom_cts_log(x,a,prob = prob) ## 
  l_rowmax  = apply(l,1,max)
  L = exp(l -  l_rowmax)
  return(list(L = L, l_rowmax = l_rowmax))
}

# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log <- function(x, a, prob){
  tmp = x*log(1-prob) 
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  out = t(t(log(prob)) * a) + tmp + lgamma(outer(x, a, "+")) - lgamma(x+1)
  out = t(t(out) - lgamma(a))
  return(out)
}

# dnbinom_cts_log_1d <- function(x, a, prob){
#   tmp = x*log(1-prob)
#   tmp[x == 0] = 0 ## R says 0*-Inf = NaN
#   return(a*log(prob) + tmp + lgamma(x+a) - lgamma(x+1) - lgamma(a))
# }

# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log_vec <- function(x, a, prob){
  #browser()
  if(length(x) > 1 && length(a) == 1){a = replicate(length(x), a)}
  tmp = x*log(1-prob)
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  ## compute lgamma(a + x) - lgamma(a)
  if(any(a > 1e+4) && any(x < 1e+4) && any(a/x > 1e+4)){ ## this condition needs refinement!!
    lgamma_diff = lgamma_diff_taylor(a, x)
  }
  else{lgamma_diff = lgamma(a + x)  - lgamma(a)}
  
  return(a*log(prob) + tmp + lgamma_diff - lgamma(x+1))
}

# it is equivalent to dpois in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dpois_cts_log_vec <- function(x, lam){
  x * log(lam) - lam - lgamma(x+1)
}



## compute lgamma(x + dx) - lgamma(dx), y > x
## when dx << x, and x is very large and dx is small
## use taylor expansion: 
## lgamma(y) - lgamma(x) \approx  
##                        digamma(x) * (y - x) + 
##                        1/2 * psigamma(x, deriv = 1) * (y - x)^2, 
##                        1/6 * psigamma(x, deriv = 2) * (y - x)^3+ 
##                        1/24 * psigamma(c, deriv = 3) * (y - x)^4 where c \in (x, y), and I choose c = x
lgamma_diff_taylor <- function(x, dx){
  c = x
  out = digamma(x) * dx + 
    1/2 * psigamma(c, deriv = 1) * dx^2 
    # 1/6 * psigamma(x, deriv = 2) * dx^3 +
    # 1/24 * psigamma(x, deriv = 3) * dx^4 +
    # 1/120 * psigamma(c, deriv = 4) * dx^5
  return(out)
}


select_shape_gamma <- function(x, s, scale, m = 2, d = NULL, low = NULL){
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min((1/scale) * xprime/s) 
  mu_grid_max = 2*max((1/scale) * x/s)
  if(is.null(m)){
    if(is.null(d)){m = 2}
    else{m = ceiling((mu_grid_max/mu_grid_min)^(1/(d-1)))}
  }
  
  ## some specification of mu_grid_min
  if(!is.null(low)){mu_grid_min = min(low, mu_grid_min)}
  mu_grid_min = max(1e-50, mu_grid_min)
  
  shape = geom_seq(mu_grid_min, mu_grid_max, m)
  return(shape)
}


select_scale_exponential <- function(x, s, m = 2, d = NULL, low = NULL){
  ## mu_grid: mu =  1/b is the exponential mean
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min(xprime/s)
  mu_grid_max = 2*max(x/s)
  if(is.null(m)){
    if(is.null(d)){m = 2}
    else{m = ceiling((mu_grid_max/mu_grid_min)^(1/(d-1)))}
  }
  if(!is.null(low)){mu_grid_min = min(low, mu_grid_min)}
  mu_grid_min = max(1e-50, mu_grid_min)
  scale = geom_seq(mu_grid_min, mu_grid_max, m)
  return(scale)
}



##=====================functions used in  older versions=============================
## select grid for b_k
select_grid_exponential <- function(x, s, m = 2, d = NULL, low = NULL){
  ## mu_grid: mu =  1/b is the exponential mean
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min(xprime/s)
  mu_grid_max = 2*max(x/s)
  if(is.null(m)){
    if(is.null(d)){m = 2}
    else{m = ceiling((mu_grid_max/mu_grid_min)^(1/(d-1)))}
  }
  if(!is.null(low)){mu_grid_min = min(low, mu_grid_min)}
  
  mu_grid_min = max(1e-50, mu_grid_min)
  
  mu_grid = geom_seq(mu_grid_min, mu_grid_max, m)
  a = rep(1, length(mu_grid))
  return(list(shape = a, scale = mu_grid))
}

## select grid for a_k
select_grid_gamma <- function(x, s, m = 2, d = NULL, low = NULL, theta = "one"){
  #if(length(x) != length(s)){browser()}
  ## mu_grid: mu =  1/b is the exponential mean
  
  #browser()
  if(identical(theta, "one")){theta = 1}
  if(identical(theta, "max")){theta = max(x)/s[1]}
  
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min((1/theta) * xprime/s) ## same as x/s_0
  mu_grid_max = 2*max((1/theta) * x/s)
  if(is.null(m)){
    if(is.null(d)){m = 2}
    else{m = ceiling((mu_grid_max/mu_grid_min)^(1/(d-1)))}
  }
  
  ## some specification of mu_grid_min
  if(!is.null(low)){mu_grid_min = min(low, mu_grid_min)}
  mu_grid_min = max(1e-50, mu_grid_min)
  
  shape = geom_seq(mu_grid_min, mu_grid_max, m)
  scale = rep(theta, length(shape))
  return(list(shape = shape, scale = scale))
}


