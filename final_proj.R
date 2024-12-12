heart_dat <- read.table("~/Desktop/biostat885/project/raw_dat.txt", quote="\"", comment.char="")
# Convert transplant to factor
heart_dat = heart_dat[-which(heart_dat$age < 18), ]
heart_dat$fustat = 1 - heart_dat$fustat
heart_dat$transplant <- factor(heart_dat$transplant)
setwd("/Users/ziyuliu/Desktop/biostat885/project/")
glm_mod = glm(factor(fustat) ~  age * transplant, family = binomial, data = heart_dat)
glm_pred_val_1 = predict(glm_mod, data.frame(age = heart_dat$age, transplant = factor(rep(1, nrow(heart_dat)))), type = "link", se.fit = TRUE)
glm_fit_1 =  plogis(glm_pred_val_1$fit)
glm_pred_val_0 = predict(glm_mod, data.frame(age = heart_dat$age, transplant = factor(rep(0, nrow(heart_dat)))), type = "link", se.fit = TRUE)
glm_fit_0 =  plogis(glm_pred_val_0$fit)

gaussian_kernel <- function(x, x_center, h) {
  exp(-0.5 * ((x - x_center) / h)^2)
}



loocv = function(kernel, h_val, heart_dat) {
  log_lik_val = c()
  for (i in 1:nrow(heart_dat)) {
    trn_dat = heart_dat[-i,]
    tst_dat = heart_dat[i, , drop = FALSE]
    age_val = tst_dat$age
    # print(age_val)
    if (kernel == "uniform") {
      within_age = trn_dat$age >= age_val - h_val & trn_dat$age <= age_val + h_val
      tmp_dat = trn_dat[within_age, ]
      fitted_mod = glm(factor(fustat) ~ age * transplant, family = binomial, data = tmp_dat)
    } else if (kernel == "gaussian") {
      weights = gaussian_kernel(trn_dat$age, age_val, h_val)
      weights = weights / sum(weights)
      fitted_mod = glm(factor(fustat) ~ age * transplant, family = binomial, data = trn_dat, weights = weights)
    }
    pred_prob = predict(fitted_mod, newdata = tst_dat, type = "response")
    response_val = tst_dat$fustat
    log_lik_val = c(log_lik_val, -response_val * log(pred_prob) - (1 - response_val) * log(1 - pred_prob))
  }
  loocv_ll = mean(log_lik_val, na.rm = TRUE)
  return(loocv_ll)
}
library(pROC)
loocv = function(kernel, h_val, heart_dat) {
  actual_vals = c()
  pred_probs = c()
  for (i in 1:nrow(heart_dat)) {
    trn_dat = heart_dat[-i,]
    tst_dat = heart_dat[i, , drop = FALSE]
    age_val = tst_dat$age
    # print(age_val)
    if (kernel == "uniform") {
      within_age = trn_dat$age >= age_val - h_val & trn_dat$age <= age_val + h_val
      tmp_dat = trn_dat[within_age, ]
      fitted_mod = glm(factor(fustat) ~ age * transplant, family = binomial, data = tmp_dat)
    } else if (kernel == "gaussian") {
      weights = gaussian_kernel(trn_dat$age, age_val, h_val)
      weights = weights / sum(weights)
      fitted_mod = glm(factor(fustat) ~ age * transplant, family = binomial, data = trn_dat, weights = weights)
    }
    pred_prob = predict(fitted_mod, newdata = tst_dat, type = "response")
    pred_probs = c(pred_probs, pred_prob)
    actual_vals = c(actual_vals, tst_dat$fustat)
  }
  auc_val = auc(actual_vals, pred_probs)
  return(auc_val)
}

bandwidth_seq_unif = 10:100
auc_vals_unif = sapply(bandwidth_seq_unif, function(i) loocv("uniform", i, heart_dat))
uniform_h_max = bandwidth_seq_unif[which.max(auc_vals_unif)]
uniform_h = bandwidth_seq_unif[which.max(auc_vals_unif[auc_vals_unif != max(auc_vals_unif)])]
bandwidth_seq_gau = 1:100
gaussian_h = bandwidth_seq_gau[which.max(sapply(bandwidth_seq_gau, function(i) loocv("gaussian", i, heart_dat)))]
kern_fit = function(heart_dat, kernel, h_val) {
  # Initialize storage vectors
  pred_val_1 <- c()
  lowerbound_1 <- c()
  upperbound_1 <- c()
  pred_val_0 <- c()
  lowerbound_0 <- c()
  upperbound_0 <- c()
  if (kernel == "uniform") {
    for (i in 1:nrow(heart_dat)) {
      age_val <- heart_dat$age[i]
      transplant_val <- heart_dat$transplant[i]
      
      # Subset within bandwidth
      within_age <- heart_dat$age >= age_val - h_val & heart_dat$age <= age_val + h_val
      tmp_dat <- heart_dat[within_age, ]
      
      
      # Fit GLM
      fitted_mod <- glm(factor(fustat) ~ age * transplant, family = binomial, data = tmp_dat)
      
      # Predict for current values
      new_data_1 <- data.frame(age = age_val, transplant = factor(1))
      pred_surv_1 <- predict(fitted_mod, new_data_1, type = "link", se.fit = TRUE)
      
      new_data_0 = data.frame(age = age_val, transplant = factor(0))
      pred_surv_0 = predict(fitted_mod, new_data_0, type = "link", se.fit = TRUE)
      
      # Store predictions
      pred_val_1 <- c(pred_val_1, plogis(pred_surv_1$fit))
      lowerbound_1 <- c(lowerbound_1, plogis(pred_surv_1$fit - 1.96 * pred_surv_1$se.fit))
      upperbound_1 <- c(upperbound_1, plogis(pred_surv_1$fit + 1.96 * pred_surv_1$se.fit))
      
      pred_val_0 <- c(pred_val_0, plogis(pred_surv_0$fit))
      lowerbound_0 <- c(lowerbound_0, plogis(pred_surv_0$fit - 1.96 * pred_surv_0$se.fit))
      upperbound_0 <- c(upperbound_0, plogis(pred_surv_0$fit + 1.96 * pred_surv_0$se.fit))
    }
  } else if (kernel == "gaussian") {
    for (i in 1:nrow(heart_dat)) {
      age_val <- heart_dat$age[i]
      # print(age_val)
      transplant_val <- heart_dat$transplant[i]
      weights <- gaussian_kernel(heart_dat$age, age_val, h_val)
      fitted_mod <- glm(factor(fustat) ~ age * transplant, family = binomial, data = heart_dat, weights = weights)
      new_data_1 <- data.frame(age = age_val, transplant = factor(1))
      pred_surv_1 <- predict(fitted_mod, new_data_1, type = "link", se.fit = TRUE)
      new_data_0 = data.frame(age = age_val, transplant = factor(0))
      pred_surv_0 = predict(fitted_mod, new_data_0, type = "link", se.fit = TRUE)
      pred_val_1 <- c(pred_val_1, plogis(pred_surv_1$fit))
      lowerbound_1 <- c(lowerbound_1, plogis(pred_surv_1$fit - 1.96 * pred_surv_1$se.fit))
      upperbound_1 <- c(upperbound_1, plogis(pred_surv_1$fit + 1.96 * pred_surv_1$se.fit))
      pred_val_0 <- c(pred_val_0, plogis(pred_surv_0$fit))
      lowerbound_0 <- c(lowerbound_0, plogis(pred_surv_0$fit - 1.96 * pred_surv_0$se.fit))
      upperbound_0 <- c(upperbound_0, plogis(pred_surv_0$fit + 1.96 * pred_surv_0$se.fit))
    }
  }
  return(list(pred_val_1 = pred_val_1, lowerbound_1 = lowerbound_1, upperbound_1 = upperbound_1, pred_val_0 = pred_val_0, lowerbound_0 = lowerbound_0, upperbound_0 = upperbound_0))
}

kern_fit(heart_dat, "gaussian", gaussian_h)
pred_check = kern_fit(heart_dat, "gaussian", gaussian_h)
library(ggplot2)
res_gen = function(kern_res, heart_dat, kern) {
  # Sort predictions and plot
  od_list <- order(heart_dat$age)
  heart_dat <- heart_dat[od_list, ]
  s_pred_val_1 <- kern_res$pred_val_1[od_list]
  s_lowerbound_1 <- kern_res$lowerbound_1[od_list]
  s_upperbound_1 <- kern_res$upperbound_1[od_list]
  # glm_fit_1 = glm_fit_1[od_list]
  s_pred_val_0 <- kern_res$pred_val_0[od_list]
  s_lowerbound_0 <- kern_res$lowerbound_0[od_list]
  s_upperbound_0 <- kern_res$upperbound_0[od_list]
  df = data.frame(x = rep(heart_dat$age, 2), y = c(s_pred_val_1, s_pred_val_0), group = rep(c("Transplant", "No transplant"), each = nrow(heart_dat)), ymin = c(s_lowerbound_1, s_lowerbound_0), ymax = c(s_upperbound_1, s_upperbound_0))
  my_plot = ggplot(df, aes(x = x, y = y, color = group, group = group)) + geom_line(size =1) + geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = group), alpha= 0.2) + theme_minimal() + labs(x = "Age", y = "Survival probability")
  ggsave(paste0("kernel_", kern, "_CI.png"), plot = my_plot)
}
pred_check_unif_max = kern_fit(heart_dat, "uniform",uniform_h_max)
pred_check_unif = kern_fit(heart_dat, "uniform", uniform_h)
# Loop through rows
uniform_h
gaussian_h
res_gen(pred_check, heart_dat, "gaussian")
res_gen(pred_check_unif, heart_dat, "uniform")
res_gen(pred_check_unif_max, heart_dat, "uniform")
 library(mgcv)

od_list <- order(heart_dat$age)
heart_dat <- heart_dat[od_list,]

bic_list = c()
k_list = c()
# Fit GAM model
for (k1 in 1:10) {
  for (k2 in 1:10) {
    gam_mod <- gam(factor(fustat) ~ s(age, k = k1) + s(age, by = transplant, k = k2), 
                   family = binomial(link = "logit"), 
                   data = heart_dat)
    bic_list = c(bic_list, BIC(gam_mod))
    k_list = rbind(k_list, c(k1, k2))
  }
}
which.min(bic_list)
k_list[which.min(bic_list),]
gam_mod <- gam(factor(fustat) ~ s(age, k = 4) + s(age, by = transplant, k = 1), 
               family = binomial(link = "logit"), 
               data = heart_dat)


# Predict values for transplant = 1
s_pred_val_1 <- predict(gam_mod, 
                        newdata = data.frame(age = heart_dat$age, transplant = factor(rep(1, nrow(heart_dat)))), 
                        type = "link", se.fit = TRUE)
s_lowerbound_1 <- plogis(s_pred_val_1$fit - 1.96 * s_pred_val_1$se.fit)
s_upperbound_1 <- plogis(s_pred_val_1$fit + 1.96 * s_pred_val_1$se.fit)

# Predict values for transplant = 0
s_pred_val_0 <- predict(gam_mod, 
                        newdata = data.frame(age = heart_dat$age, transplant = factor(rep(0, nrow(heart_dat)))), 
                        type = "link", se.fit = TRUE)
s_lowerbound_0 <- plogis(s_pred_val_0$fit - 1.96 * s_pred_val_0$se.fit)
s_upperbound_0 <- plogis(s_pred_val_0$fit + 1.96 * s_pred_val_0$se.fit)
df = data.frame(x = rep(heart_dat$age, 2), y = c(plogis(s_pred_val_1$fit), plogis(s_pred_val_0$fit)), group = rep(c("Transplant", "No Transplant"), each = nrow(heart_dat)), ymin = c(s_lowerbound_1, s_lowerbound_0), ymax = c(s_upperbound_1, s_upperbound_0))

my_plot = ggplot(df, aes(x = x, y = y, color = group, group = group)) + geom_line(size =1) + geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = group), alpha= 0.2) + theme_minimal() + labs(x = "Age", y = "Survival probability")
ggsave("gam_bic.png", plot = my_plot)

glm_mod = glm(factor(fustat) ~ age * transplant, family = binomial(link = "logit"), data = heart_dat)
s_pred_val_1 = predict(glm_mod, 
                       newdata = data.frame(age = heart_dat$age, transplant = factor(rep(1, nrow(heart_dat)))), 
                       type = "link", se.fit = TRUE)
s_lowerbound_1 <- plogis(s_pred_val_1$fit - 1.96 * s_pred_val_1$se.fit)
s_upperbound_1 <- plogis(s_pred_val_1$fit + 1.96 * s_pred_val_1$se.fit)
s_pred_val_0 <- predict(glm_mod, 
                        newdata = data.frame(age = heart_dat$age, transplant = factor(rep(0, nrow(heart_dat)))), 
                        type = "link", se.fit = TRUE)
s_lowerbound_0 <- plogis(s_pred_val_0$fit - 1.96 * s_pred_val_0$se.fit)
s_upperbound_0 <- plogis(s_pred_val_0$fit + 1.96 * s_pred_val_0$se.fit)
df = data.frame(x = rep(heart_dat$age, 2), y = c(plogis(s_pred_val_1$fit), plogis(s_pred_val_0$fit)), group = rep(c("Transplant", "No Transplant"), each = nrow(heart_dat)), ymin = c(s_lowerbound_1, s_lowerbound_0), ymax = c(s_upperbound_1, s_upperbound_0))
glm_plot = ggplot(df, aes(x = x, y = y, color = group, group = group)) + geom_line(size =1) + geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = group), alpha= 0.2) + theme_minimal() + labs(x = "Age", y = "Survival probability")
ggsave("glm_CI.png", plot = glm_plot)
summary(glm_mod)
BIC(glm_mod)
min(bic_list)
mean(heart_dat$age)
sd(heart_dat$age)


sum(heart_dat$fustat == 1)
sum(heart_dat$fustat == 0)

sum(heart_dat$transplant == 1)
sum(heart_dat$transplant == 0)

summary(glm(factor(fustat) ~ transplant, family = binomial(link = "logit"), data = heart_dat))
summary(glm_mod)
