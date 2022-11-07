#' @title Empirical Bayes Poisson Mean with Point Gamma  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Point Gamma: g()  = pi_0 delta() + (1-pi_0) gamma(shape, scale)
#' 
#' @import stats

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over pi_0, shape, scale)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x vector of Poisson observations.
#' @param s vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param g_init The prior distribution \eqn{g}, of the class \code{point_gamma}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. 
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#' @param pi0 Either  \code{"estimate"} which optimizes over pi0 along with  shape and scale, or a number in \code{[0,1]} that fixes pi0
#' @param control A list of control parameters  to be passed to the optimization function. `nlm` is  used. 

#' 
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, and posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} of class \code{point_gamma}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    out = ebpm_point_gamma(x,s)
#' @export


#trace("pg_nlm_fn", quote(if(any(is.nan(log(par[1])))) { browser() }), at=4, print=F)



# ebpm_point_gamma <- function(x, s, g_init = NULL, seed = 123){
#   if(is.null(g_init)){
#     mus = 10^seq(0,5,1)
#     for(mu in mus){
#       g_init = c(0.5, mu*mean(x[x!=0]/s[x != 0]),mu);
#       try(
#         return(ebpm_point_gamma_helper(x, s, g_init, seed = seed)), silent = T
#       )
#     }
#     stop('ebpm_point_gamma fails  for multiple initialization')
#   }
#   else{
#     try(return(ebpm_point_gamma_helper(x, s, g_init, seed = seed)))
#     stop('ebpm_point_gamma fails  for multiple initialization')
#   }
# }
  



## TODO:
## consider the case where X are all  0s.   If  s are  not all 0, then pi* = 1, which is not reachable after transformation
ebpm_point_gamma <- function(x, s = 1, g_init = NULL, fix_g = F, pi0 = "estimate",control = NULL){
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = nlm_control_defaults()}
  if(is.null(g_init)){g_init = point_gamma(0.5,1,1); fix_g =  F}
  #browser()
  if(!fix_g){
    ## MLE
    if(identical(pi0, "estimate")){
      if(!all(x  > 0)){
        fn_params = list(x = x, s = s)
        opt = do.call(nlm, c(list(pg_nlm_fn, transform_param(g_init)), fn_params, control))
        log_likelihood =  -pg_nlm_fn(opt$estimate, x, s)
        opt_g = transform_param_back(opt$estimate)
      }else{  ## in  this case, optimal pi0 is 0, but is not reachable after a transformation in nlm; so fix pi0 at 0
        pi0 = 0
        fn_params = list(x = x, s = s, pi0 = pi0)
        opt = do.call(nlm, c(list(pg_nlm_fn_fix_pi0, transform_param_fix_pi0(g_init)), fn_params, control))
        log_likelihood =  -pg_nlm_fn_fix_pi0(opt$estimate, x, s, pi0)
        opt_g = c(pi0, transform_param_back_fix_pi0(opt$estimate))
      }
    }else{ ## fix pi0
      pi0 = as.numeric(pi0)
      fn_params = list(x = x, s = s, pi0 = pi0)
      opt = do.call(nlm, c(list(pg_nlm_fn_fix_pi0, transform_param_fix_pi0(g_init)), fn_params, control))
      log_likelihood =  -pg_nlm_fn_fix_pi0(opt$estimate, x, s, pi0)
      opt_g = c(pi0, transform_param_back_fix_pi0(opt$estimate))
    }
  }else{ ## fix_g = T
      opt_g = as.numeric(g_init)
      log_likelihood =  ifelse( g_init$pi0 == 0, -pg_nlm_fn_pi0(transform_param_pi0(opt_g), x, s), -pg_nlm_fn(transform_param(opt_g), x, s))
  }
  fitted_g = point_gamma(pi0 = opt_g[1], shape = opt_g[2], scale  = opt_g[3])
  ## posterior mean
  pi = fitted_g$pi0
  a =  fitted_g$shape
  b =  1/fitted_g$scale
  nb = exp(dnbinom_cts_log_vec(x, a, prob = b/(b+s)))
  
  if(pi == 0){pi_hat = replicate(length(x),0)}
  else{pi_hat = pi*as.integer(x == 0)/(pi*as.integer(x ==  0) + (1-pi)*nb)}
  #pi_hat = pi*as.integer(x ==  0)/(pi*as.integer(x ==  0) + (1-pi)*nb)
  
  lam_pm = (1-pi_hat)*(a+x)/(b+s)
  lam_log_pm =  digamma(a + x) - log(b + s)
  #lam_log_pm[x == 0] = -Inf
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(list(fitted_g = fitted_g, posterior = posterior, log_likelihood = log_likelihood))
}

# pg_nlm_fn_fix_pi0 <- function(par, x, s, pi0){  ## for the case where x  > 0 for all a, and optimal pi0 is 0
#   n = length(x)
#   a  =  exp(par[1])
#   b  =  exp(par[2])
#   d <- exp(dnbinom_cts_log_vec(x, a, b/(b+s)))
#   #if(is.nan(sum(log(pi * c + d)))){browser()}
#   return(-sum(log(d)) - n * log(1-pi0))
# }

# pg_nlm_fn_fix_pi02 <- function(par, x, s, pi0){
#   pi = pi0
#   a  = exp(par[1])
#   b  = exp(par[2])
#   d <- exp(dnbinom_cts_log_vec(x, a, b/(b+s)))
#   c = as.integer(x ==  0) - d
#   ## return log(pi0 + (1-pi0) d) if x = 0;
#   ## else return log(1-pi0) + log(d)
#   return(-sum(log(pi*c + d)))
# }

pg_nlm_fn_fix_pi0 <- function(par, x, s, pi0){
  #browser()
  ## d = NB(x, a, b/(b+s))
  ## return - log(pi0 + (1-pi0) d) if x = 0; 
  ## else return - (log(1-pi0) + log(d))
  a  = exp(par[1])
  b  = exp(par[2])
  d_log <- dnbinom_cts_log_vec(x, a, b/(b+s))
  if(pi0 == 0){
    out = sum(d_log)
  }else{
    out = sum((log(1-pi0) + d_log[x!=0])) + sum(log(pi0) + log1p(((1 - pi0)/pi0) * exp(d_log[x == 0])))
  }
  return(-out)
}



pg_nlm_fn <- function(par, x, s){
  ## d = NB(x, a, b/(b+s))
  ## return - log(pi0 + (1-pi0) d) if x = 0; 
  ## else return - (log(1-pi0) + log(d))
  #browser()
  pi0 = 1/(1+ exp(-par[1]))
  a = exp(par[2])
  b  =  exp(par[3])
  d_log <- dnbinom_cts_log_vec(x, a, b/(b+s))
  #d = exp(dnbinom_cts_log_vec(x, a, b/(b+s)))
  #c = as.integer(x ==  0) - d
  out = sum((log(1-pi0) + d_log[x!=0])) + sum(log(pi0) + log1p(((1 - pi0)/pi0) * exp(d_log[x == 0])))
  #if(is.nan(sum(log(pi0 * c + d)))){browser()}
  return(-out)
}

## par0: pi0,shape, scale
## want: logit(pi0), log(shape), log(1/scale)
transform_param <- function(par0){
  par0 = as.numeric(par0)
  par = rep(0,length(par0))
  par[1] = log(par0[1]/(1-par0[1]))
  par[2] = log(par0[2])
  par[3] = - log(par0[3])
  return(par)
}

transform_param_back <- function(par){
  par0 = rep(0,length(par))
  par0[1] = 1/(1+ exp(-par[1]))
  par0[2] = exp(par[2])
  par0[3] = exp(- par[3])
  return(par0)
}

transform_param_fix_pi0 <- function(par0){
  ## only need to optimize over shape, scale
  par0 = c(par0$shape, par0$scale)
  par = rep(0,length(par0))
  par[1] = log(par0[1])
  par[2] = - log(par0[2])
  return(par)
}

transform_param_back_fix_pi0 <- function(par){
  par0 = rep(0,length(par))
  par0[1] = exp(par[1])
  par0[2] = exp(- par[2])
  return(par0)
}

