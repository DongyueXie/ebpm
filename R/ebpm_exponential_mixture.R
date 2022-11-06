#' @title Empirical Bayes Poisson Mean with Mixture  of Exponential  as Prior
#' @description Uses Empirical Bayes to fit the model \deqn{x_j | \lambda_j ~ Poi(s_j \lambda_j)} with \deqn{lambda_j ~ g()}
#' with Mixture of Exponential: \deqn{g()  = sum_k pi_k gamma(shape = 1, rate = b_k)} 
#' b_k is selected to cover the lambda_i  of interest for all data  points x_i
#' @import mixsqp

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param x A vector of Poisson observations.
#' @param s A vector of scaling factors for Poisson observations: the model is \eqn{y[j]~Pois(s[j]*lambda[j])}.
#' @param scale Either a vector, or "estimate", which finds grid for scale (exponential mean) that makes the range the prior mean cover all \code{x/s}
#' @param g_init The prior distribution \eqn{g}, of the class \code{gammamix}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. 
#'
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#'   
#' @param m multiple coefficient when selectig grid, so the  b_k is of the form {low*m^{k-1}}; must be greater than 1; default is 2
#' @param control A list of control parameters  to be passed to the optimization function. `mixsqp` is  used.
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, and to add posterior log mean).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}}, of class \code{gammamix}} 
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'      }
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    x = rpois(100,beta) # simulate Poisson observations
#'    s = replicate(100,1)
#'    m = 2
#'    out = ebpm::ebpm_exponential_mixture(x,s)
#'    
#' @export

## compute ebpm_exponential_mixture problem
ebpm_exponential_mixture <- function(x,s = 1,  scale = "estimate", g_init = NULL, fix_g = FALSE,m = 2, control =  NULL, low = NULL){
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = mixsqp_control_defaults()}
  if(is.null(g_init)){
    fix_g = FALSE ## then automatically unfix g if specified so
    if(identical(scale, "estimate")){scale <- select_scale_exponential(x,s,m, low)}
    g_init = scale2gammamix_init(shape = replicate(length(scale), 1),scale = scale)
  }
  
  if(!fix_g){ ## need to estimate g_hat
    b = 1/g_init$scale ##  from here use gamma(shape = a, rate = b)  where E = a/b
    a = g_init$shape
    tmp <-  compute_L(x,s,a, b)
    L =  tmp$L
    l_rowmax = tmp$l_rowmax
    fit <- mixsqp(L, x0 = g_init$pi, control = control)
    pi = fit$x
    pi = pi/sum(pi) ## seems that some times pi does not sum to one
  }
  else{
    pi = g_init$pi
    a = g_init$shape
    b = 1/g_init$scale
    ## compute loglikelihood
    tmp <-  compute_L(x,s,a, b)
    L =  tmp$L
    l_rowmax = tmp$l_rowmax
  }
  fitted_g = gammamix(pi = pi, shape = a,  scale  = 1/b)
  log_likelihood = sum(log(exp(l_rowmax) * L %*%  pi))
  
  cpm = outer(x,a,  "+")/outer(s, b, "+")
  Pi_tilde = t(t(L) * pi)
  Pi_tilde = Pi_tilde/rowSums(Pi_tilde)
  lam_pm = rowSums(Pi_tilde * cpm)
  c_log_pm = digamma(outer(x,a,  "+")) - log(outer(s, b, "+"))
  lam_log_pm = rowSums(Pi_tilde * c_log_pm)
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(list(fitted_g = fitted_g, posterior = posterior,log_likelihood = log_likelihood))
}


