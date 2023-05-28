remove(list=ls())

wd <- "~path to root directory"
setwd(wd)

## Path to output folder
path = "1- Sample Input Data/Analysis files/"

source("4- Model/Programs_Feb_2020.R")
library(tidyverse)

marker <- "Overweight"
date <- paste0("NS_imp_Mar23")
date_file <- paste0("Mar23")

B <- 20
boot_vals <- 1:B
plot_data <- NULL
SE_par <- NULL

for(j in 1:B){
  
  
  all_data <- readRDS(paste0("1- Sample Input Data/Over_data_w_cov_",date_file,"_multiple_impute.rds")) %>% 
    filter(.imp == j | .imp == 0) %>% 
    mutate(
      Sex = case_when(
        Sex == "Both" ~ 0,
        Sex == "Female" ~ 1,
        Sex == "Male" ~ -1
      ),
      Region = factor(Region), 
      "SMART" = case_when(
        ShortSource=="SMART" ~ 1,
        TRUE ~ 0
      ),
      "Surveillance" = case_when(
        ShortSource=="Surveillance" ~ 1,
        TRUE ~ 0
      )
    ) %>% 
    select(c("country","year", "Point.Estimate.NS", "Point.Estimate.Imp","SE_val", "ShortSource",
             "Region","SEV", "Sex", "SMART", "Surveillance", "MCI_5_yr", "SDI")) %>% 
    rename("Y" = "Point.Estimate.NS", 
           "Y_all" = "Point.Estimate.Imp",
           "SE_var" = "SE_val") %>% 
    arrange(country, year) #%>% filter(!is.na(Y))
  
  cat(j,length(unique(all_data$country)), 
      length(c(all_data$Y)),
      length(c(all_data$Y[!is.na(all_data$Y)])), 
      table(all_data$Region),"\n") 
  
  all_data$P_Region <- all_data$Region
  levels(all_data$P_Region)
  referent_lev <- 2
  all_data$P_Region <- relevel(all_data$P_Region,ref = referent_lev)
  levels(all_data$P_Region)
  
  ## Transforming outcome ##
  data_w_out <- all_data %>% 
    select(c("country", "year", "Y", "SE_var")) %>% 
    mutate(
      SE_pred = 0, 
      SE_var = SE_var*((1/Y) + 1/(1-Y)), 
      Y = log(Y/(1-Y))
    )
  
  
  Pcov_data <- model.matrix(~P_Region,data=all_data,contrasts.arg = list(P_Region=diag(nlevels(all_data$P_Region))))[,-1]
  
  TRANS=FALSE
  ##Number of penalized and random splines.  Equally spaced throughout the B.knots.
  DF_P <- seq(1993,2022,5)
  DF_R <- NULL #quantile(data_w_out$year,probs = 0.5)#
  B.knots <- range(data_w_out$year[!is.na(data_w_out$Y)])
  B.knots[1] <- B.knots[1] - 1
  B.knots[2] <- B.knots[2] + 10
  q.order = 2
  plots=FALSE
  slope = TRUE
  
  #Specifying the covariance matrix of the random effects.  
  #   - cov_mat="CS" (default) -> compound symmetric covariance matrix
  #   - cov_mat="UN"  -> unstructured symmetric covariance matrix
  #   - cov_mat="VC"  -> Diagonal only covariance matrix
  cov_mat <- "VC"
  
  
  cov_data <- as.matrix(data.frame(model.matrix(
    ~ SMART + Sex + P_Region + MCI_5_yr, 
    data = all_data)[,-1]
  ))
  zero_covs <- colnames(cov_data)[1]
  
  
  ##################### Covariate analysis with multiple penalized functions #################################
  t1 <- try(Estimation <- cmnpe(data_w_out, DF_P, DF_R, B.knots,q.order , cov_data=cov_data, Pcov_data = Pcov_data, cov_mat = cov_mat,plots=plots,TRANS=TRANS,zero_covs=zero_covs,slope = slope))
  
  if(is.null(attr(t1,"class"))){
  cat(j,2*Estimation$df - 2*c(summary(Estimation$result$model)$logLik) + summary(Estimation$result$model)$AIC+2*c(summary(Estimation$result$model)$logLik),"\n")
  
  t_plot_data <- Estimation$plot_data
  saveRDS(t_plot_data, file =  paste0(path,"/Plot data for ",marker," imputation ",j,".rds"))
  if(j < B){
    remove(Estimation)
    remove(all_data)
  }
  }else{
    boot_vals <- boot_vals[boot_vals != j]
    cat(j,"Model Error.\n")
  }
}


### Summarizing imputed samples:
Point.Est_mat <- NULL
Point.Est_fixed_mat <- NULL
Point.Est_fixpen_mat <- NULL
Var_mean_mat <- NULL
Var_pred_mat <- NULL
B <- length(boot_vals)
for(j in boot_vals){
  plot_data <- read_rds( paste0(path,"/Plot data for ",marker," imputation ",j,".rds"))
  cat(dim(plot_data),"\n")
  Point.Est_mat <- cbind(Point.Est_mat,plot_data$pred)
  Point.Est_fixed_mat <- cbind(Point.Est_fixed_mat,plot_data$pred_fixed)
  Point.Est_fixpen_mat <- cbind(Point.Est_fixpen_mat,plot_data$pred_fixpen)
  Var_mean_mat  <- cbind(Var_mean_mat ,plot_data$sigma_T_est^2)
  Var_pred_mat  <- cbind(Var_pred_mat ,plot_data$sigma_Y_est^2)
}

plot_data$pred <- apply(Point.Est_mat,1,mean)
plot_data$pred_fixed <- apply(Point.Est_fixed_mat,1,mean)
plot_data$pred_fixpen <- apply(Point.Est_fixpen_mat,1,mean)
B_Var <- apply(Point.Est_mat,1,var)
W_mean <- apply(Var_mean_mat,1,mean)
W_pred <- apply(Var_pred_mat,1,mean)
plot_data$SE_mean_pred <- sqrt(W_mean + (B+1)/B*B_Var)
plot_data$SE_pred <- sqrt(W_pred + (B+1)/B*B_Var)
plot_data$lower_CI <- plot_data$pred - 1.96*plot_data$SE_mean_pred
plot_data$upper_CI <- plot_data$pred + 1.96*plot_data$SE_mean_pred
plot_data$lower_CI2 <- plot_data$pred - 1.96*plot_data$SE_pred
plot_data$upper_CI2 <- plot_data$pred + 1.96*plot_data$SE_pred

remove(list = c("Point.Est_mat","Point.Est_fixed_mat", 
                "Point.Est_fixpen_mat", "Var_mean_mat",
                "Var_pred_mat"))

Estimation$plot_data <- plot_data

data_one <- readRDS(paste0("1- Sample Input Data/Over_data_w_cov_",date_file,"_multiple_impute.rds")) %>% 
  filter(.imp == j | .imp == 0) %>% 
  select(-".imp")


### Outputting the data to "path" folder.
output_function(Estimation, all_data, data_one, marker, date, path)
