
# handle flip

label_flip = function(niter, res_prop, eta_true) {
ids = (niter/2 + 1):niter
eta_post_mode = apply(res_prop$eta_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] }) # estimated eta, posterior mode 
qq = mean(res_prop$q_tk[ids]) 
if (qq > 0.5) {
eta_post = 1 - eta_post_mode
b_mean = mean(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids])
b_sd = sd(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids])
ab_mean = mean(res_prop$beta_tk[ids])
ab_sd = sd(res_prop$beta_tk[ids])
bci = quantile(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids], probs=c(0.025,0.975), na.rm=T)
}
if (qq <=0.5) {
eta_post = eta_post_mode
b_mean = mean(res_prop$beta_tk[ids])
b_sd = sd(res_prop$beta_tk[ids])
ab_mean = mean(res_prop$alpha_tk[ids] + res_prop$beta_tk[ids])
ab_sd = sd(res_prop$alpha_tk[ids] + res_prop$beta_tk[ids])
bci = quantile(res_prop$beta_tk[ids], probs=c(0.025,0.975), na.rm=T)
}
eta_match = mean(eta_post == eta_true)
eta_est = eta_post

return(list(qq = qq, eta_match = eta_match, eta_est = eta_est, b_mean = b_mean, b_sd = b_sd, bci = bci))

}

# handle flip

label_flip_thin = function(niter, res_prop, eta_true, thin = 10) {
ids = seq((niter/2), niter, thin)
eta_post_mode = apply(res_prop$eta_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] }) # estimated eta, posterior mode 
qq = mean(res_prop$q_tk[ids]) 
if (qq > 0.5) {
eta_post = 1 - eta_post_mode
b_mean = mean(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids])
b_sd = sd(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids])
ab_mean = mean(res_prop$beta_tk[ids])
ab_sd = sd(res_prop$beta_tk[ids])
}
if (qq <=0.5) {
eta_post = eta_post_mode
b_mean = mean(res_prop$beta_tk[ids])
b_sd = sd(res_prop$beta_tk[ids])
ab_mean = mean(res_prop$alpha_tk[ids] + res_prop$beta_tk[ids])
ab_sd = sd(res_prop$alpha_tk[ids] + res_prop$beta_tk[ids])
}
eta_match = mean(eta_post == eta_true)
eta_est = eta_post

return(list(qq = qq, eta_match = eta_match, eta_est = eta_est, b_mean = b_mean, b_sd = b_sd))

}

# handle flip - joint

label_flip_joint = function(niter, res, eta_true_1, eta_true_2) {
  ids = (niter/2 + 1):niter
  p00_post = res$pst_tk[ids, 1]
  p01_post = res$pst_tk[ids, 2]
  p10_post = res$pst_tk[ids, 3]
  p11_post = res$pst_tk[ids, 4]
  eta_1_post_mode = apply(res$eta_1_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] })
  eta_2_post_mode = apply(res$eta_2_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] })
  qq1 = mean(p10_post + p11_post) 
  qq2 = mean(p01_post + p11_post)
  if (qq1 > 0.5) {
  eta_1_post = 1 - eta_1_post_mode
  b1_mean = mean(res$beta_1_tk[ids] + res$alpha_1_tk[ids])
  b1_sd = sd(res$beta_1_tk[ids] + res$alpha_1_tk[ids])
  bci1 = quantile(res$beta_1_tk[ids] + res$alpha_1_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  if (qq1 <=0.5) {
  eta_1_post = eta_1_post_mode
  b1_mean = mean(res$beta_1_tk[ids])
  b1_sd = sd(res$beta_1_tk[ids])
  bci1 = quantile(res$beta_1_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  if (qq2 > 0.5) {
  eta_2_post = 1 - eta_2_post_mode
  b2_mean = mean(res$beta_2_tk[ids] + res$alpha_2_tk[ids])
  b2_sd = sd(res$beta_2_tk[ids] + res$alpha_2_tk[ids])
  bci2 = quantile(res$beta_2_tk[ids] + res$alpha_2_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  if (qq2 <=0.5) {
  eta_2_post = eta_2_post_mode
  b2_mean = mean(res$beta_2_tk[ids])
  b2_sd = sd(res$beta_2_tk[ids])
  bci2 = quantile(res$beta_2_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  
  eta_1_match = mean(eta_1_post == eta_true_1)
  eta_2_match = mean(eta_2_post == eta_true_2)
  eta_1_est = eta_1_post
  eta_1_est = eta_1_post
  
  return(list(qq1 = qq1, qq2 = qq2, eta_1_match = eta_1_match, eta_2_match = eta_2_match, 
  eta_1_est = eta_1_post, eta_2_est = eta_2_post, b1_mean = b1_mean, b1_sd = b1_sd, b2_mean = b2_mean, b2_sd = b2_sd, bci1 = bci1, bci2 = bci2))
  
}

# handle flip - thin - joint

label_flip_joint_thin = function(niter, res, eta_true_1, eta_true_2, thin = 10) {
  ids = seq((niter/2), niter, thin)
  p00_post = res$pst_tk[ids, 1]
  p01_post = res$pst_tk[ids, 2]
  p10_post = res$pst_tk[ids, 3]
  p11_post = res$pst_tk[ids, 4]
  eta_1_post_mode = apply(res$eta_1_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] })
  eta_2_post_mode = apply(res$eta_2_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] })
  qq1 = mean(p10_post + p11_post) 
  qq2 = mean(p01_post + p11_post)
  if (qq1 > 0.5) {
    eta_1_post = 1 - eta_1_post_mode
    b1_mean = mean(res$beta_1_tk[ids] + res$alpha_1_tk[ids])
    b1_sd = sd(res$beta_1_tk[ids] + res$alpha_1_tk[ids])
  }
  if (qq1 <=0.5) {
    eta_1_post = eta_1_post_mode
    b1_mean = mean(res$beta_1_tk[ids])
    b1_sd = sd(res$beta_1_tk[ids])
  }
  if (qq2 > 0.5) {
    eta_2_post = 1 - eta_2_post_mode
    b2_mean = mean(res$beta_2_tk[ids] + res$alpha_2_tk[ids])
    b2_sd = sd(res$beta_2_tk[ids] + res$alpha_2_tk[ids])
  }
  if (qq2 <=0.5) {
    eta_2_post = eta_2_post_mode
    b2_mean = mean(res$beta_2_tk[ids])
    b2_sd = sd(res$beta_2_tk[ids])
  }
  
  eta_1_match = mean(eta_1_post == eta_true_1)
  eta_2_match = mean(eta_2_post == eta_true_2)
  eta_1_est = eta_1_post
  eta_1_est = eta_1_post
  
  return(list(qq1 = qq1, qq2 = qq2, eta_1_match = eta_1_match, eta_2_match = eta_2_match, 
              eta_1_est = eta_1_post, eta_2_est = eta_2_post, b1_mean = b1_mean, b1_sd = b1_sd, b2_mean = b2_mean, b2_sd = b2_sd))
  
}



