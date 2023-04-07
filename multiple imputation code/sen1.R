load("sen_res_dir/sen1_res1.Rdata")
sen_1m = sen1_a
load("sen_res_dir/sen1_res2.Rdata")
sen_2m = sen1_a
load("sen_res_dir/sen1_res3.Rdata")
sen_3m = sen1_a
load("sen_res_dir/sen1_res4.Rdata")
sen_4m = sen1_a
load("sen_res_dir/sen1_res5.Rdata")
sen_5m = sen1_a

mean(sen_1m[[5]],sen_2m[[5]],sen_3m[[5]],sen_4m[[5]],sen_5m[[5]])

###########################
#model 5.1
rubin <- function(beta_vec,sigma_vec,lambda_vec){
  theta_hat = mean(beta_vec)
  #coef[1] is the point estiamte of beta (coef)
  Vw = mean(sigma_vec^2)
  #coef[3] is se(coef)
  Vb = 1/(5-1)*((beta1-theta_hat)^2+(beta2-theta_hat)^2+(beta3-theta_hat)^2+
                  (beta4-theta_hat)^2+ (beta5-theta_hat)^2)
  Vt = Vw + Vb + Vb/5
  SE_pool = sqrt(Vt)
  
  adj_lambda = (Vb + Vb/5)/Vt
  df_old = (5-1)/(adj_lambda^2)
  n=444375
  k=23
  df_obs = (n-k+1)/(n=k+3)*(n-k)*(1-adj_lambda)
  df_adj = df_old*df_obs/(df_old+df_obs)
  Wald_pool = sqrt(theta_hat^2/Vt)
  if(Wald_pool > 0){
    p_value = 2*(1-pt(Wald_pool,df_adj))
  }else{
    p_value = 2*pt(Wald_pool,df_adj)
  }
  
  t_value = qt(0.975,df_adj)
  upper_b = theta_hat + t_value *SE_pool
  lower_b = theta_hat - t_value *SE_pool
  
  HR_hat = exp(theta_hat)
  HR_up = exp(upper_b)
  HR_low=exp(lower_b)
  
  lambda = mean(lambda_vec)
  HR_hat_dil = exp(theta_hat/lambda)
  HR_up_dil = exp(upper_b/lambda)
  HR_low_dil=exp(lower_b/lambda)
  
  res = c(HR_hat,HR_low,HR_up,HR_hat_dil,HR_low_dil,HR_up_dil,p_value)
  res_pr = c(paste(sprintf("%0.2f", HR_hat)," (",
                   sprintf("%0.2f", HR_low),",",
                   sprintf("%0.2f", HR_up),")",sep=""),
             
             paste(sprintf("%0.2f", HR_hat_dil)," (",
                   sprintf("%0.2f", HR_low_dil),",",
                   sprintf("%0.2f", HR_up_dil),")",sep=""),
             
             sprintf("%0.3f", p_value))
  names(res) = c("HR","lower","upper","HR_adj","lower_adj","upper_adj","p")
  a = list(res,res_pr)
  return(a)
}

ori = rep(NA,4)
adj = rep(NA,4)
p = rep(NA,4)
for(i in c(3,4)){
  k=24
  s1=summary(sen_1m[[i]])$coefficients
  beta1 = s1[k,1]
  sigma1 = s1[k,3]
  s1=summary(sen_2m[[i]])$coefficients
  beta2 = s1[k,1]
  sigma2 = s1[k,3]
  s1=summary(sen_3m[[i]])$coefficients
  beta3 = s1[k,1]
  sigma3 = s1[k,3]
  s1=summary(sen_4m[[i]])$coefficients
  beta4 = s1[k,1]
  sigma4 = s1[k,3]
  s1=summary(sen_5m[[i]])$coefficients
  beta5 = s1[k,1]
  sigma5 = s1[k,3]
  
  beta_vec = c(beta1,beta2,beta3,beta4,beta5)
  sigma_vec = c(sigma1,sigma2,sigma3,sigma4,sigma5)
  lambda_vec = c(sen_1m[[5]],sen_2m[[5]],sen_3m[[5]],sen_4m[[5]],sen_5m[[5]])
  temp = rubin(beta_vec,sigma_vec,lambda_vec)
  ori[i] = temp[[2]][1]
  adj[i] = temp[[2]][2]
  p[i] = temp[[2]][3]
}

ori_ord = rep(NA,2)
adj_ord = rep(NA,2)
p_ord = rep(NA,2)
for(j in c(1,2)){
  k=c(24,25)
  i=1
  s1=summary(sen_1m[[i]])$coefficients
  beta1 = s1[k[j],1]
  sigma1 = s1[k[j],3]
  s1=summary(sen_2m[[i]])$coefficients
  beta2 = s1[k[j],1]
  sigma2 = s1[k[j],3]
  s1=summary(sen_3m[[i]])$coefficients
  beta3 = s1[k[j],1]
  sigma3 = s1[k[j],3]
  s1=summary(sen_4m[[i]])$coefficients
  beta4 = s1[k[j],1]
  sigma4 = s1[k[j],3]
  s1=summary(sen_5m[[i]])$coefficients
  beta5 = s1[k[j],1]
  sigma5 = s1[k[j],3]
  
  beta_vec = c(beta1,beta2,beta3,beta4,beta5)
  sigma_vec = c(sigma1,sigma2,sigma3,sigma4,sigma5)
  lambda_vec = c(sen_1m[[5]],sen_2m[[5]],sen_3m[[5]],sen_4m[[5]],sen_5m[[5]])
  temp = rubin(beta_vec,sigma_vec,lambda_vec)
  ori_ord[j] = temp[[2]][1]
  adj_ord[j] = temp[[2]][2]
  p_ord[j] = temp[[2]][3]
}

ori_cat = rep(NA,3)
adj_cat = rep(NA,3)
p_cat = rep(NA,3)
for(j in c(1,2,3)){
  k=c(24,25,26)
  i=2
  s1=summary(sen_1m[[i]])$coefficients
  beta1 = s1[k[j],1]
  sigma1 = s1[k[j],3]
  s1=summary(sen_2m[[i]])$coefficients
  beta2 = s1[k[j],1]
  sigma2 = s1[k[j],3]
  s1=summary(sen_3m[[i]])$coefficients
  beta3 = s1[k[j],1]
  sigma3 = s1[k[j],3]
  s1=summary(sen_4m[[i]])$coefficients
  beta4 = s1[k[j],1]
  sigma4 = s1[k[j],3]
  s1=summary(sen_5m[[i]])$coefficients
  beta5 = s1[k[j],1]
  sigma5 = s1[k[j],3]
  
  beta_vec = c(beta1,beta2,beta3,beta4,beta5)
  sigma_vec = c(sigma1,sigma2,sigma3,sigma4,sigma5)
  lambda_vec = c(sen_1m[[5]],sen_2m[[5]],sen_3m[[5]],sen_4m[[5]],sen_5m[[5]])
  temp = rubin(beta_vec,sigma_vec,lambda_vec)
  ori_cat[j] = temp[[2]][1]
  adj_cat[j] = temp[[2]][2]
  p_cat[j] = temp[[2]][3]
}

ori
adj
p

ori_ord
adj_ord
p_ord

ori_cat
adj_cat
p_cat


