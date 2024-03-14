fn_ll1 <- function(x,stick0,...){
  p = length(stick0)
  -loglik1(b_t=x[1:p],theta_t=x[p+1],r_invalid_t=x[(p+2):length(x)],...)
}

pl <- function(x,b_exp,b_out,Sig_exp_inv,Sig_out_inv){
  l = -1/2* t(b_out-x*b_exp) %*% solve(solve(Sig_out_inv) + x^2*solve(Sig_exp_inv)) %*% (b_out-x*b_exp)
  return(-l)
}


loglik1 <- function(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                    b_t,theta_t,r_invalid_t,fix_invalid){
  r_vec_t = rep(0,length(b_out))
  r_vec_t[fix_invalid] = r_invalid_t
  beta_Y_star = b_out - theta_t * b_t - r_vec_t
  beta_X_star = b_exp - b_t
  l = -1/2 * (t(beta_Y_star) %*% Sig_out_inv %*% beta_Y_star +
                t(beta_X_star) %*% Sig_exp_inv %*% beta_X_star)
  return(l)
}

loglik <- function(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                   b_t,theta_t,r_vec_t){
  beta_Y_star = b_out - theta_t * b_t - r_vec_t
  beta_X_star = b_exp - b_t
  l = -1/2 * (t(beta_Y_star) %*% Sig_out_inv %*% beta_Y_star +
                t(beta_X_star) %*% Sig_exp_inv %*% beta_X_star)
  return(l)
}

cML_estimate_cor0 <- function(b_exp,b_out,
                              Sig_exp_inv,Sig_out_inv,
                              K,initial_theta = 0,
                              initial_mu = rep(0,length(b_exp)),
                              maxit = 100)
{
  p = length(b_exp)
  ### initialize
  theta = initial_theta
  mu_vec = initial_mu
  if(K==0){
    ite_ind = 0
    r_vec_tilde = rep(0,p)
    theta_tilde = theta; theta_tilde_old = theta_tilde + 1
    while(abs(theta_tilde-theta_tilde_old)>1e-4 & ite_ind<maxit){
      theta_tilde_old = theta_tilde
      mu_vec_tilde =
        solve(theta_tilde^2 * Sig_out_inv + Sig_exp_inv) %*%
        (theta_tilde * Sig_out_inv %*% (b_out - r_vec_tilde) + Sig_exp_inv %*% b_exp)
      theta_tilde = as.numeric(
        (t(b_out-r_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde) /  (t(mu_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde)
      )
      ite_ind = ite_ind + 1
    }
    conv = (ite_ind<maxit)
    if(conv==1){
      se =   sqrt(1/numDeriv::hessian(func = pl, x = theta_tilde, b_exp = b_exp,b_out = b_out,
                                      Sig_exp_inv=Sig_exp_inv,
                                      Sig_out_inv=Sig_out_inv))
    }else{se=NA}
    return(list(A = which(r_vec_tilde!=0), r_vec=rep(0,p),mu_vec=mu_vec_tilde,theta=theta_tilde, se = se,
                converge=conv, ite_vec = ite_ind,
                l =loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                          b_t = mu_vec_tilde, theta_t = theta_tilde, r_vec_t = r_vec_tilde)))
  }
  v_vec = (b_out - theta * mu_vec)
  r_vec = rep(0,p)
  v_importance1 = v_vec^2
  nonzero_bg_ind1 = sort((order(v_importance1,decreasing = T))[1:K])
  r_vec[nonzero_bg_ind1] = v_vec[nonzero_bg_ind1]
  L = L0 = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                  b_t = mu_vec, theta_t = theta, r_vec_t = r_vec)

  se_out = diag(solve(Sig_out_inv))
  ite_vec = NULL
  ###
  A_m_new = which(r_vec!=0)
  A_m_old = 0
  while(any(!is.element(A_m_old,A_m_new))){
    A_m_old = A_m_new
    ite_ind = 0
    theta_tilde = theta; theta_tilde_old = theta_tilde + 1
    mu_vec_tilde = mu_vec
    #solve cMLE given set A
    while( (abs(theta_tilde_old - theta_tilde) > 1e-5) & (ite_ind<maxit))
    {
      theta_tilde_old = theta_tilde
      ite_ind = ite_ind + 1

      ### first, update v_bg
      r_vec_tilde = rep(0,p)
      r_vec_tilde[A_m_new] = t(b_out - theta_tilde*mu_vec_tilde) %*% Sig_out_inv[,A_m_new] %*% solve(Sig_out_inv[A_m_new,A_m_new])
      mu_vec_tilde =
        solve(theta_tilde^2 * Sig_out_inv + Sig_exp_inv) %*%
        (theta_tilde * Sig_out_inv %*% (b_out - r_vec_tilde) + Sig_exp_inv %*% b_exp)
      theta_tilde = as.numeric(
        (t(b_out-r_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde) /  (t(mu_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde)
      )
      ite_ind = ite_ind + 1
    }
    ite_vec = c(ite_vec,ite_ind)
    # if not converge, next
    if(ite_ind>=maxit){next;}
    # if improve likelihood, update set A
    if(loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
              b_t = mu_vec_tilde, theta_t = theta_tilde, r_vec_t = r_vec_tilde) > L){
      L = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                 b_t = mu_vec_tilde, theta_t = theta_tilde, r_vec_t = r_vec_tilde)
      r_vec = r_vec_tilde
      mu_vec = mu_vec_tilde
      theta = theta_tilde
      #update set A
      v_vec = (b_out - theta * mu_vec)
      v_importance1 = v_vec^2/se_out
      A_m_new = sort((order(v_importance1,decreasing = T))[1:K])
    }

  }

  if(theta==initial_theta){
    conv = 0
    return(list(theta = theta,
                b_vec = mu_vec,
                r_vec = r_vec,
                se = NA,
                ite_vec = ite_vec,
                converge=conv,
                A = A_m_new,
                l = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                           b_t = mu_vec, theta_t = theta, r_vec_t = r_vec)))
  }else{
    conv = 1
    ## ONE MORE STEP -
    valid_ind = setdiff(1:p,A_m_new)

    opt_cor = optim(c(mu_vec, theta, r_vec[-valid_ind]),
                    fn_ll1,
                    stick0=b_exp,
                    b_out=b_out,
                    b_exp=b_exp,
                    Sig_exp_inv=Sig_exp_inv,
                    Sig_out_inv=Sig_out_inv,
                    fix_invalid=sort(A_m_new),
                    method='BFGS',
                    hessian=TRUE,control = list(maxit=50000))
    theta = opt_cor$par[p+1]
    se = try(sqrt(diag(solve(opt_cor$hessian)))[p+1])
    if(class(se)=='try-error'){se=NA}
    r_vec = rep(0,p)
    r_vec[-valid_ind] = opt_cor$par[(p+2):length(opt_cor$par)]
    b_vec = opt_cor$par[1:p]


    return(list(theta = theta,
                b_vec = mu_vec,
                r_vec = r_vec,
                se = se,
                ite_vec = ite_vec,
                converge=conv,
                A = A_m_new,
                l = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                           b_t = mu_vec, theta_t = theta, r_vec_t = r_vec)))
  }

}

cML_estimate0 <- function(b_exp,b_out,
                          Sig_exp_inv,Sig_out_inv,
                          K,initial_theta = 0,
                          initial_mu = rep(0,length(b_exp)),
                          maxit = 100)
{
  p = length(b_exp)
  ### initialize
  theta = initial_theta
  mu_vec = initial_mu
  if(K==0){
    ite_ind = 0
    r_vec_tilde = rep(0,p)
    theta_tilde = theta; theta_tilde_old = theta_tilde + 1
    while(abs(theta_tilde-theta_tilde_old)>1e-4 & ite_ind<maxit){
      theta_tilde_old = theta_tilde
      mu_vec_tilde =
        solve(theta_tilde^2 * Sig_out_inv + Sig_exp_inv) %*%
        (theta_tilde * Sig_out_inv %*% (b_out - r_vec_tilde) + Sig_exp_inv %*% b_exp)
      theta_tilde = as.numeric(
        (t(b_out-r_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde) /  (t(mu_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde)
      )
      ite_ind = ite_ind + 1
    }
    conv = (ite_ind<maxit)
    if(conv==1){
      se =   sqrt(1/numDeriv::hessian(func = pl, x = theta_tilde, b_exp = b_exp,b_out = b_out,
                                      Sig_exp_inv=Sig_exp_inv,
                                      Sig_out_inv=Sig_out_inv))
    }else{se=NA}
    return(list(A = which(r_vec_tilde!=0), r_vec=rep(0,p),mu_vec=mu_vec_tilde,theta=theta_tilde, se = se,
                converge=conv, ite_vec = ite_ind,
                l =loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                          b_t = mu_vec_tilde, theta_t = theta_tilde, r_vec_t = r_vec_tilde)))
  }
  v_vec = (b_out - theta * mu_vec)
  r_vec = rep(0,p)
  v_importance1 = v_vec^2
  nonzero_bg_ind1 = sort((order(v_importance1,decreasing = T))[1:K])
  r_vec[nonzero_bg_ind1] = v_vec[nonzero_bg_ind1]
  L = L0 = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                  b_t = mu_vec, theta_t = theta, r_vec_t = r_vec)

  se_out = diag(solve(Sig_out_inv))
  ite_vec = NULL
  ###
  A_m_new = which(r_vec!=0)
  ite_ind = 0
  theta_tilde = theta; theta_tilde_old = theta_tilde + 1
  mu_vec_tilde = mu_vec

  while( (abs(theta_tilde_old - theta_tilde) > 1e-5) & (ite_ind<maxit))
  {
    theta_tilde_old = theta_tilde
    ite_ind = ite_ind + 1

    ### first, update v_bg
    v_importance1 = (b_out - theta_tilde * mu_vec_tilde)^2/se_out
    A_m_new = sort((order(v_importance1,decreasing = T))[1:K])

    r_vec_tilde = rep(0,p)
    r_vec_tilde[A_m_new] = t(b_out - theta_tilde*mu_vec_tilde) %*% Sig_out_inv[,A_m_new] %*% solve(Sig_out_inv[A_m_new,A_m_new])

    mu_vec_tilde =
      solve(theta_tilde^2 * Sig_out_inv + Sig_exp_inv) %*%
      (theta_tilde * Sig_out_inv %*% (b_out - r_vec_tilde) + Sig_exp_inv %*% b_exp)
    theta_tilde = as.numeric(
      (t(b_out-r_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde) /  (t(mu_vec_tilde) %*% Sig_out_inv %*% mu_vec_tilde)
    )
    ite_ind = ite_ind + 1
  }

  # if not converge, next
  if(ite_ind>=maxit){
    conv = 0
    return(list(theta = theta,
                b_vec = mu_vec,
                r_vec = r_vec,
                se = NA,
                ite_vec = ite_vec,
                converge=conv,
                A = A_m_new,
                l = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                           b_t = mu_vec, theta_t = theta, r_vec_t = r_vec)))
  }else{
    conv = 1
    ## ONE MORE STEP -
    valid_ind = setdiff(1:p,A_m_new)

    opt_cor = optim(c(mu_vec, theta, r_vec[-valid_ind]),
                    fn_ll1,
                    stick0=b_exp,
                    b_out=b_out,
                    b_exp=b_exp,
                    Sig_exp_inv=Sig_exp_inv,
                    Sig_out_inv=Sig_out_inv,
                    fix_invalid=sort(A_m_new),
                    method='BFGS',
                    hessian=TRUE,control = list(maxit=50000))
    theta = opt_cor$par[p+1]
    se = try(sqrt(diag(solve(opt_cor$hessian)))[p+1])
    if(class(se)=='try-error'){se=NA}
    r_vec = rep(0,p)
    r_vec[-valid_ind] = opt_cor$par[(p+2):length(opt_cor$par)]
    b_vec = opt_cor$par[1:p]


    return(list(theta = theta,
                b_vec = mu_vec,
                r_vec = r_vec,
                se = se,
                ite_vec = ite_vec,
                converge=conv,
                A = A_m_new,
                l = loglik(b_exp,b_out,Sig_exp_inv,Sig_out_inv,
                           b_t = mu_vec, theta_t = theta, r_vec_t = r_vec)))
  }

}







cML_estimate_random <- function(b_exp, b_out,
                                Sig_exp_inv,Sig_out_inv,
                                K,random_start = 0,
                                maxit = 100,
                                min_theta_range = -0.5,
                                max_theta_range = 0.5)
{
  p = length(b_exp)

  theta_v_RandomCandidate = NULL
  sd_v_RandomCandidate = NULL
  l_v_RandomCandidate = NULL
  invalid_RandomCandidate = NULL
  conv_v_RandomCandidate = NULL
  Sig_exp = solve(Sig_exp_inv)
  for(random_ind in 1:(1+random_start))
  {
    #      ptm <- proc.time()
    if(random_ind == 1)
    {
      initial_theta = 0
      initial_mu = rep(0,length(b_exp))
      #initial_mu = b_exp
    } else {
      initial_theta = runif(1,min = min_theta_range,max = max_theta_range)
      initial_mu = MASS::mvrnorm(1, mu = b_exp, Sigma = Sig_exp)
    }
    #  print(proc.time() - ptm)
    #    ptm <- proc.time()


    MLE_result =
      cML_estimate_cor0(b_exp,b_out,
                        Sig_exp_inv,Sig_out_inv,
                        K = K,initial_theta = initial_theta,
                        initial_mu = initial_mu,
                        maxit = maxit)


    Neg_l = -MLE_result$l

    theta_v_RandomCandidate = c(theta_v_RandomCandidate,MLE_result$theta)
    sd_v_RandomCandidate = c(sd_v_RandomCandidate,MLE_result$se)
    l_v_RandomCandidate = c(l_v_RandomCandidate,Neg_l)
    invalid_RandomCandidate = rbind(invalid_RandomCandidate,
                                    as.numeric(MLE_result$r_vec))
    conv_v_RandomCandidate = c(conv_v_RandomCandidate,MLE_result$converge)


    #  print(proc.time() - ptm)
  }
  keep_conv = which(conv_v_RandomCandidate==1)
  if(length(keep_conv)==0){
    return(list(theta = NA,
                se = NA,
                l = NA,
                r_est = rep(NA,p),
                converge = 0))
  }
  l_v_RandomCandidate = l_v_RandomCandidate[keep_conv]
  theta_v_RandomCandidate = theta_v_RandomCandidate[keep_conv]
  sd_v_RandomCandidate = sd_v_RandomCandidate[keep_conv]
  invalid_RandomCandidate = invalid_RandomCandidate[keep_conv,,drop=FALSE]
  conv_v_RandomCandidate = conv_v_RandomCandidate[keep_conv]

  min_neg_l = which.min(l_v_RandomCandidate)
  theta_est = theta_v_RandomCandidate[min_neg_l]
  sd_est = sd_v_RandomCandidate[min_neg_l]
  l_est = l_v_RandomCandidate[min_neg_l]
  r_est = invalid_RandomCandidate[min_neg_l,]
  conv_est = conv_v_RandomCandidate[min_neg_l]


  return(list(theta = theta_est,
              se = sd_est,
              l = -l_est,
              r_est = r_est,
              converge = conv_est
  )
  )
}

cismr_cML <- function(b_exp,b_out,
                      Sig_exp_inv,Sig_out_inv,
                      K_vec = 0:(length(b_exp) - 2),
                      random_start = 0,
                      maxit = 100,
                      random_seed = 0,
                      n,
                      min_theta_range = -0.5,
                      max_theta_range = 0.5)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }
  rand_theta = NULL
  rand_sd = NULL
  rand_l = NULL
  rand_conv = NULL
  invalid_mat = NULL
  for(K_value in K_vec)
  {
    rand_res = cML_estimate_random(b_exp = b_exp,
                                   b_out = b_out,
                                   Sig_exp_inv = Sig_exp_inv,
                                   Sig_out_inv = Sig_out_inv,
                                   K = K_value,
                                   random_start = random_start,
                                   maxit = maxit,
                                   min_theta_range = min_theta_range,
                                   max_theta_range = max_theta_range
    )
    rand_theta = c(rand_theta,rand_res$theta)
    rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    rand_conv= c(rand_conv,rand_res$converge)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
  }

  ### get result
  theta_v = rand_theta
  sd_v = rand_sd
  l_v = rand_l

  # cML-BIC
  BIC_vec = log(n) * K_vec - 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec,na.rm = T)
  min_ind = which.min(BIC_vec)
  BIC_theta = theta_v[min_ind]
  BIC_se = sd_v[min_ind]
  BIC_p = pnorm(-abs(BIC_theta/BIC_se))*2
  BIC_invalid = which(invalid_mat[min_ind,]!=0)


  return(list(
    BIC_theta = BIC_theta,
    BIC_se = BIC_se,
    BIC_p = BIC_p,
    BIC_invalid = BIC_invalid,
    l_vec = l_v,
    BIC_vec = log(n) * K_vec - 2 * l_v,
    converge_vec = rand_conv)
  )
}

#' Main function to run cisMR-cML with data perturbation
#'
#' @param b_exp A vector of conditional genetic effect sizes on the exposure.
#' @param b_out A vector of conditional genetic effect sizes on the outcome.
#' @param Sig_exp_inv The inverse of the covariance matrix of \code{b_exp}.
#' @param Sig_out_inv The inverse of the covariance matrix of \code{b_out}.
#' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
#' @param random_start Number of random start points for cisMR-cML-BIC, default is 0.
#' @param random_start_pert Number of random start points for cisMR-cML with data perturbation, default is 0.
#' @param maxit Maximum number of iterations for each optimization, default is 200.
#' @param num_pert Number of perturbation, default is 100.
#' @param random_seed Random seed, an integer. Default is
#' 0, which does not set random seed; user could specify a positive integer
#' as random seed to get replicable results.
#' @param min_theta_range The lower bound for the randomly generated initial causal effect size, default is -0.5.
#' @param max_theta_range The upper bound for the randomly generated initial causal effect size, default is 0.5.
#' @param n Sample size of exposure or outcome GWAS.
#'
#' @return A list
#' \describe{
#' \item{BIC_theta}{Estimated causal effect from cisMR-cML-BIC}
#' \item{BIC_se}{Estimate standard error for \code{BIC_theta}}
#' \item{BIC_p}{P-value for \code{BIC_theta}}
#' \item{BIC_invalid}{Invalid IVs selected by cisMR-cML-BIC}
#' \item{BIC_DP_theta}{Estimated causal effect from cisMR-cML-DP }
#' \item{BIC_DP_se}{Estimate standard error for \code{BIC_DP_theta}}
#' \item{BIC_DP_p}{P-value for \code{BIC_DP_theta}}
#' }
#'
#' @import MASS numDeriv
#' @importFrom stats optim pnorm runif sd
#' @export
cismr_cML_DP <- function(b_exp,b_out,
                         Sig_exp_inv,Sig_out_inv,
                         K_vec = 0:(length(b_exp) - 2),
                         random_start = 0,
                         random_start_pert = 0,
                         maxit = 200,
                         num_pert = 100,
                         random_seed = 0,
                         min_theta_range = -0.5,
                         max_theta_range = 0.5,
                         n)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }
  theta_v = theta_MA_v = NULL
  cML_res = cismr_cML(b_exp = b_exp, b_out = b_out, Sig_exp_inv = Sig_exp_inv, Sig_out_inv = Sig_out_inv,
                      K_vec = K_vec, random_start = random_start, random_seed = random_seed, n = n, maxit = maxit,
                      min_theta_range = min_theta_range, max_theta_range = max_theta_range)
  BIC_invalid_mat = matrix(0,nrow=length(b_exp),ncol=num_pert)
  p = length(b_exp)
  Sig_exp = solve(Sig_exp_inv)
  Sig_out = solve(Sig_out_inv)
  for(pt_ind in 1:num_pert){
    b_exp_new = MASS::mvrnorm(1, mu = b_exp, Sigma = Sig_exp)
    b_out_new = MASS::mvrnorm(1, mu = b_out, Sigma = Sig_out)
    cML_res_b = cismr_cML(b_exp = b_exp_new, b_out = b_out_new, Sig_exp_inv = Sig_exp_inv, Sig_out_inv = Sig_out_inv,
                          K_vec = K_vec, random_start = random_start_pert, n = n, maxit = maxit,
                          min_theta_range = min_theta_range, max_theta_range = max_theta_range)
    #    print(cML_res_b$BIC_theta)
    theta_v = c(theta_v,cML_res_b$BIC_theta)
    BIC_invalid_mat[cML_res_b$BIC_invalid,pt_ind] = 1
  }



  # cML-BIC-DP
  BIC_DP_theta = mean(theta_v,na.rm=T)
  BIC_DP_se = sd(theta_v,na.rm=T)
  BIC_DP_p = pnorm(-abs(BIC_DP_theta/BIC_DP_se))*2


  return(list(
    BIC_theta = cML_res$BIC_theta,
    BIC_se = cML_res$BIC_se,
    BIC_p = cML_res$BIC_p,
    BIC_invalid = cML_res$BIC_invalid,
    BIC_DP_theta = BIC_DP_theta,
    BIC_DP_se = BIC_DP_se,
    BIC_DP_p = BIC_DP_p,
    BIC_DP_invalid = BIC_invalid_mat
  ))
}
