# Title: Pediatric Donor Risk Index 
# Author: Dave Merola
# Date Created: 06-19-2024
# Description: Construction and testing of prediction models for graft failure
# and mortality in pediatric liver transplant patients using the UNOS data.

rm(list = ls())

# Constants ----
# Specify directories
DIR_RAW_DATA <- "~/OneDrive/Research/UNOS Data/STAR_STATA/SAS Export to STATA 202403/"
DIR_FORMATS <- "~/OneDrive/Research/UNOS Data/STAR_STATA/CODE DICTIONARY - FORMATS 202403/Liver/"
DIR_DATA <- "~/OneDrive/Research/Transplant DRI/Data/"
DIR_OUTPUT <- "~/OneDrive/Research/Transplant DRI/Output/"

# Load Libraries ----
library(tableone)
library(performance)
library(randomForest)
library(pROC)
library(mice)
library(glmnet)
library(haven)
library(tidyverse)

# Functions ----
perf_metrics <- function(obs, pred, cutoff) {
  
  # Description: This function generates performance metrics based on observed
  #              and predicted values. Inputs must have the same length with 
  #              0 indicating negative and 1 indicating positive event status.
  # 
  # Arguments: obs: A binary numeric vector of observed values (0, 1)
  #            pred: A numeric vector of predicted probabilities [0:1]
  #            cutoff: A scalar indicating the cutoff in which to classify the 
  #                    predicted values as positive or negative
  # 
  # Output: Returns a vector containing the cutoff, sensitivity, specificity,
  #         PPV, NPV, F1 score, and MCC
  
  # A quality control step
  if((length(obs) == length(pred)) == FALSE) {
    print("ERROR: OBSERVED & PREDICTED VECTORS HAVE DIFFERENT LENGTHS")
    stop()
  } 
  
  # Create a standardized dataset of observed and predicted values
  data_perf <- as.data.frame(cbind(as.numeric(obs), as.numeric(pred))) %>% 
    `colnames<-`(c("observed", "predicted"))
  
  # Assess sensitivity
  sensitivity <- data_perf %>% 
    filter(observed == 1) %>% 
    summarize(sensitivity = sum(predicted >= cutoff) / n()) %>% 
    pull()
  
  # Assess specificity
  specificity <- data_perf %>% 
    filter(observed == 0) %>% 
    summarize(specificity = sum(predicted < cutoff) / n()) %>% 
    pull()
  
  # Assess PPV
  ppv <- data_perf %>% 
    filter(predicted >= cutoff) %>% 
    summarize(ppv = sum(observed == 1) / n()) %>% 
    pull()
  
  # Assess NPV
  npv <- data_perf %>% 
    filter(predicted < cutoff) %>% 
    summarize(npv = sum(observed == 0) / n()) %>% 
    pull()
  
  # Assess F1 Score
  f1_score <- sensitivity / (sensitivity + (0.5*((1-specificity)+(1-sensitivity))))
  
  # Matthews Correlation Coefficient
  mcc <- ((sensitivity * specificity) - ((1-specificity) * (1-sensitivity))) / sqrt((sensitivity + (1-specificity)) * (sensitivity + (1-sensitivity)) * (specificity + (1-specificity)) * (specificity + (1-sensitivity)))
  
  # Brier Score
  brier <- sum((data_perf$predicted-data_perf$observed)^2) / nrow(data_perf)
  
  # Mean predicted
  mean_pred <- mean(data_perf$predicted)
  
  performance_summary <- cbind(cutoff, sensitivity, specificity, ppv, npv, f1_score, mcc, brier, mean_pred)
  
  return(performance_summary)
}

miss_vals_dims <- function(mat) {
  
  # Description: This function takes in a matrix or data frame and returns 
  #              the indices (column and row numbers) of missing values.
  # 
  # Arguments: An object of class matrix or data frame.
  # 
  # Output: A matrix containing the row and column numbers of missing values.
  
  missing_indices <- which(is.na(mat), arr.ind = TRUE)
  
  return(missing_indices)
}

marginal_prob <- function(data, model, vars_marg, vars_sum) {
  
  # Description: This function estimates the marginalized probability for each 
  #              observation of a dataset. Variables in "vars_marg" argument
  #              are 'fixed' and variables in "vars_sum" argument are summed over.
  # 
  # Arguments: data: an object of class data.frame containing all variables of interest
  #            model: an object of class glm that has been fit to the data frame 
  #                   specified in the 'data' argument
  #            vars_marg: a character vector containing the variable name(s) in 
  #                       'data' in which marginal effects are desired
  #            vars_sum: a character vector containing the variable name(s) in 
  #                       'data' in which we wish to sum (marginalize) over
  # 
  # Output: A numeric vector of probabilities is returned that has a length  
  #         corresponding to the number of observations in the data provided to 
  #         the function.
  
  # Loop through each observation
  
  for (i in 1:nrow(data)){
    
    # Print progress indicator
    print(paste0("Row No. ", i, " of ", nrow(data)))
    
    # Create a new data frame that contains 1) fixed values for the variables  
    # whose marginal effect we are interested in and 2) values of the variables
    # we will be summing over
    data_fixed_i <- data
    data_fixed_i[vars_marg] <- data_fixed_i[vars_marg][i, ]
    
    # Predict the probabilities using the original logistic regression model and 
    # take mean
    mean_pred_i <- mean(predict(model, newdata = data_fixed_i, type = "response"))
    
    # Store mean of the predicted values in vector
    if (i == 1){
      mean_pred <- rep(NA, nrow(data))
    }
    
    mean_pred[i] <- mean_pred_i
    
  }
  
  return(mean_pred)
}

# Load Data ----

# Convert files to lower-memory .RData format
# data_liver <- read_dta(paste0(DIR_RAW_DATA, "Liver/LIVER_DATA.DTA"))
# save(data_liver, file = paste0(DIR_DATA, "liver_data.RData"))
# data_age_months <- read_sas(paste0(DIR_RAW_DATA, "Liver/data0014198_Merola/liver_non_std.sas7bdat"))
# save(data_age_months, file = paste0(DIR_DATA, "age_months.RData"))

# Load .RData files
load(file = paste0(DIR_DATA, "liver_data.RData"))
load(file = paste0(DIR_DATA, "age_months.RData"))

formats <- read.delim(paste0(DIR_FORMATS, "LIVER_FORMATS_FLATFILE.DAT"), 
                      header = FALSE, 
                      col.names = c("label", "fmtname", "type", "code"))

# Data Management ----

## Variable & Cohort Selection ----
data_liver_cohort <- data_liver %>% 
  
  # Select variables of interest and tag as outcome, eligibility, or covariate
  select(trr_id_code = "TRR_ID_CODE", 
         donor_id = "DONOR_ID", 
         pt_code = "PT_CODE", 
         out_gstatus = "GSTATUS", 
         out_gtime = "GTIME", 
         out_pstatus = "PSTATUS", 
         out_ptime = "PTIME", 
         out_composite_death_date = "COMPOSITE_DEATH_DATE",
         elig_tx_date = "TX_DATE", 
         elig_age = "AGE", 
         elig_age_don = "AGE_DON", 
         elig_val_dt_ddr = "VAL_DT_DDR", 
         elig_val_dt_tcr = "VAL_DT_TCR", 
         elig_val_dt_trr = "VAL_DT_TRR", 
         elig_val_dt_ldr = "VAL_DT_LDR",
         elig_multiorg = "MULTIORG", 
         cov_abo_mat = "ABO_MAT", 
         cov_academic_level_trr = "ACADEMIC_LEVEL_TRR", 
         cov_academic_prg_trr = "ACADEMIC_PRG_TRR", 
         cov_albumin_tx = "ALBUMIN_TX", 
         cov_alcohol_heavy_don = "ALCOHOL_HEAVY_DON", 
         cov_amis = "AMIS", 
         cov_antihype_don = "ANTIHYPE_DON", 
         cov_arginine_don = "ARGININE_DON", 
         cov_artificial_li_trr = "ARTIFICIAL_LI_TRR", 
         cov_ascites_tx = "ASCITES_TX", 
         cov_blood_inf_conf_don = "BLOOD_INF_CONF_DON", 
         cov_blood_inf_don = "BLOOD_INF_DON", 
         cov_bmi_tcr = "BMI_TCR", 
         cov_bmi_don_calc = "BMI_DON_CALC", 
         cov_bmis = "BMIS", 
         cov_bun_don = "BUN_DON", 
         cov_cardarrest_neuro = "CARDARREST_NEURO", 
         cov_cdc_risk_hiv_don = "CDC_RISK_HIV_DON", 
         cov_clin_infect_don = "CLIN_INFECT_DON", 
         cov_cod_cad_don = "COD_CAD_DON", 
         cov_cold_isch = "COLD_ISCH",
         cov_creat_don = "CREAT_DON", 
         cov_creat_tx = "CREAT_TX", 
         cov_ddavp_don = "DDAVP_DON", 
         cov_death_circum_don = "DEATH_CIRCUM_DON", 
         cov_death_mech_don = "DEATH_MECH_DON", 
         cov_diab = "DIAB", 
         cov_diabetes_don = "DIABETES_DON", 
         cov_diag = "DIAG", 
         cov_dial_tx = "DIAL_TX", 
         cov_drmis = "DRMIS", 
         cov_ebv_igg_cad_don = "EBV_IGG_CAD_DON", 
         cov_ebv_igm_cad_don = "EBV_IGM_CAD_DON", 
         cov_ebv_serostatus = "EBV_SEROSTATUS", 
         cov_enceph_tx = "ENCEPH_TX", 
         cov_ethcat = "ETHCAT", 
         cov_ethcat_don = "ETHCAT_DON", 
         cov_func_stat_trr = "FUNC_STAT_TRR", 
         cov_gender = "GENDER", 
         cov_gender_don = "GENDER_DON", 
         cov_hist_cancer_don = "HIST_CANCER_DON", 
         cov_hist_cig_don = "HIST_CIG_DON", 
         cov_hist_cocaine_don = "HIST_COCAINE_DON", 
         cov_hist_diabetes_don = "HIST_DIABETES_DON", 
         cov_hist_hypertens_don = "HIST_HYPERTENS_DON", 
         cov_hist_insulin_dep_don = "HIST_INSULIN_DEP_DON", 
         cov_hist_oth_drug_don = "HIST_OTH_DRUG_DON", 
         cov_history_mi_don = "HISTORY_MI_DON", 
         cov_hlamis = "HLAMIS", 
         cov_inotrop_support_don = "INOTROP_SUPPORT_DON", 
         cov_inr_tx = "INR_TX", 
         cov_insulin_don = "INSULIN_DON", 
         cov_intracranial_cancer_don = "INTRACRANIAL_CANCER_DON", 
         cov_life_sup_trr = "LIFE_SUP_TRR", 
         cov_lv_eject_don = "LV_EJECT_DON", 
         cov_malig_trr = "MALIG_TRR", 
         cov_med_cond_trr = "MED_COND_TRR", 
         cov_meld_peld_lab_score = "MELD_PELD_LAB_SCORE", 
         cov_non_hrt_don = "NON_HRT_DON", 
         cov_on_vent_trr = "ON_VENT_TRR", 
         cov_other_inf_conf_don = "OTHER_INF_CONF_DON", 
         cov_other_inf_don = "OTHER_INF_DON", 
         cov_portal_vein_trr = "PORTAL_VEIN_TRR", 
         cov_prev_ab_surg_trr = "PREV_AB_SURG_TRR", 
         cov_protein_urine = "PROTEIN_URINE", 
         cov_pt_diuretics_don = "PT_DIURETICS_DON", 
         cov_pt_steroids_don = "PT_STEROIDS_DON", 
         cov_pt_t3_don = "PT_T3_DON", 
         cov_pt_t4_don = "PT_T4_DON", 
         cov_pulm_inf_conf_don = "PULM_INF_CONF_DON", 
         cov_pulm_inf_don = "PULM_INF_DON", 
         cov_recov_out_us = "RECOV_OUT_US", 
         cov_recuscit_dur = "RESUSCIT_DUR", 
         cov_sgot_don = "SGOT_DON", 
         cov_sgpt_don = "SGPT_DON", 
         cov_tattoos = "TATTOOS", 
         cov_tbili_don = "TBILI_DON", 
         cov_tbili_tx = "TBILI_TX", 
         cov_tipss_trr = "TIPSS_TRR", 
         cov_tx_procedur_ty = "TX_PROCEDUR_TY", 
         cov_urine_inf_conf_don = "URINE_INF_CONF_DON", 
         cov_urine_inf_don = "URINE_INF_DON", 
         cov_vasodil_don = "VASODIL_DON", 
         cov_share_ty = "SHARE_TY", 
         cov_exc_case = "EXC_CASE", 
         cov_cmv_don = "CMV_DON") %>% 
  
  # Apply eligibility criteria (nrow = 364,895 non-unique patient observations)
  
  # Select all transplants (n = 199,921 unique patients)
  filter(trr_id_code != "") %>% 
  
  # Inclusion 1: Initial transplant between 6/30/2004 - 3/31/2023 (n = 130,328)
  group_by(pt_code) %>% 
  arrange(elig_tx_date) %>% 
  mutate(i = seq(1, n(), 1)) %>% 
  ungroup() %>% 
  filter(i == 1 & 
           elig_tx_date >= as.Date("2004-06-30") & 
           elig_tx_date <= as.Date("2023-03-31")) %>% 
  
  # Inclusion 2: Recipient aged <18 years at time of transplant (n = 9,565)
  filter(elig_age < 18) %>% 
  
  # Inclusion 3: Donor aged <18 years (n = 6,498)
  # Note: 10 NAs for elig_age_don
  filter(elig_age_don < 18) %>% 
  
  # Exclusion 1: Multi-organ transplant recipient (n = 5,613)
  filter(elig_multiorg != "Y") %>% 
  
  # Exclusion 2: Transplant and candidate registrations not validated by UNOS (n = 5,596)
  filter(is.na(elig_val_dt_ddr) == FALSE & 
           is.na(elig_val_dt_tcr) == FALSE & 
           is.na(elig_val_dt_trr) == FALSE) 

## Variable Coding ---- 
data_liver_clean <- data_liver_cohort %>% 
  mutate(outcome = as.logical(case_when(out_gtime > 365 ~ 0,
                                        out_gtime <= 365 & out_gstatus == 1 ~ 1,
                                        out_gtime <= 365 & out_gstatus == 0 & 
                                          is.na(out_composite_death_date) == FALSE ~ 1,
                                        is.na(out_gtime) == TRUE & out_ptime <= 365 & 
                                          is.na(out_composite_death_date) == FALSE ~ 1,
                                        out_gtime <= 365 & out_gstatus == 0 & 
                                          is.na(out_composite_death_date) == TRUE ~ NA)), 
         time = case_when(out_gtime <= 365 ~ out_gtime,
                          out_gtime > 365 ~ 365,
                          is.na(out_gtime) == TRUE ~ out_ptime), # 1 pt missing gtime
         cov_abo_mat = as.factor(case_when(cov_abo_mat %in% c("1", "2") ~ "Compatible",
                                           cov_abo_mat == 3 ~ "Incompatible",
                                           TRUE ~ NA)),
         cov_academic_level_trr = as.factor(case_when(cov_academic_level_trr == "1" ~ "Full Academic Load",
                                                      cov_academic_level_trr %in% c("2", "3", "4") ~ "Reduced Academic Load",
                                                      cov_academic_level_trr == "996" | elig_age <= 5 ~ "Not School Age",
                                                      TRUE ~ NA_character_)),
         cov_academic_prg_trr = as.factor(case_when(cov_academic_level_trr == "Not School Age" ~ "Not School Age",
                                                    cov_academic_prg_trr == "1" ~ "Within One Grade Level of Peers",
                                                    cov_academic_prg_trr %in% c("2", "3") ~ "Delayed or Special Ed",
                                                    TRUE ~ NA_character_)),
         cov_albumin_tx = factor(case_when(cov_albumin_tx < 3.5 ~ "<3.5",
                                           cov_albumin_tx >= 3.5 & cov_albumin_tx <= 5.0 ~ "3.5-5.0",
                                           cov_albumin_tx > 5.0 ~ ">5.0",
                                           TRUE ~ NA_character_),
                                 levels = c("<3.5", "3.5-5.0", ">5.0")),
         cov_alcohol_heavy_don = as.factor(case_when(cov_alcohol_heavy_don == "Y" ~ 1,
                                                      cov_alcohol_heavy_don == "N" ~ 0,
                                                      TRUE ~ NA)),
         cov_amis = as.factor(cov_amis),
         cov_antihype_don = as.factor(case_when(cov_antihype_don == "Y" ~ 1,
                                                cov_antihype_don == "N" ~ 0,
                                                TRUE ~ NA)),
         cov_arginine_don = as.factor(case_when(cov_arginine_don == "Y" ~ 1,
                                                cov_arginine_don == "N" ~ 0,
                                                TRUE ~ NA)),
         cov_artificial_li_trr = as.logical(cov_artificial_li_trr),
         cov_ascites_tx = factor(case_when(cov_ascites_tx == 1 ~ "Absent",
                                              cov_ascites_tx == 2 ~ "Slight",
                                              cov_ascites_tx == 3 ~ "Moderate",
                                              cov_ascites_tx == 4 ~ "N/A",
                                              TRUE ~ NA_character_),
                                 levels = c("Absent", "Slight", "Moderate", "N/A")),
         cov_blood_inf_conf_don = as.logical(case_when(cov_blood_inf_conf_don == "Y" ~ 1,
                                            cov_blood_inf_conf_don == "N" ~ 0,
                                            TRUE ~ NA)),
         cov_blood_inf_don = as.logical(cov_blood_inf_don),
         cov_bmi_tcr = factor(case_when(cov_bmi_tcr < 18.5 ~ "<18.5",
                                        cov_bmi_tcr >= 18.5 & cov_bmi_tcr <= 25 ~ "18.5-25",
                                        cov_bmi_tcr > 25 ~ ">25",
                                         TRUE ~ NA_character_),
                               levels = c("<18.5","18.5-25",">25")), 
         cov_bmi_don_calc = factor(case_when(cov_bmi_don_calc < 18.5 ~ "<18.5",
                                             cov_bmi_don_calc >= 18.5 & cov_bmi_don_calc <= 25 ~ "18.5-25",
                                             cov_bmi_don_calc > 25 ~ ">25",
                                             TRUE ~ NA_character_),
                                   levels = c("<18.5","18.5-25",">25")), 
         cov_bmis = factor(cov_bmis),
         cov_bun_don = factor(case_when(cov_bun_don < 5 ~ "<5",
                                        cov_bun_don >= 5 & cov_bun_don <= 18 ~ "5-18",
                                        cov_bun_don > 18 ~ ">18",
                                        TRUE ~ NA_character_),
                              levels = c("<5", "5-18", ">18")),
         cov_cardarrest_neuro = as.logical(case_when(cov_cardarrest_neuro == "Y" ~ 1,
                                                     cov_cardarrest_neuro == "N" ~ 0,
                                                     TRUE ~ NA)),
         cov_cdc_risk_hiv_don = as.logical(case_when(cov_cdc_risk_hiv_don == "Y" ~ 1,
                                                     cov_cdc_risk_hiv_don == "N" ~ 0,
                                                     TRUE ~ NA)),
         cov_clin_infect_don = as.logical(case_when(cov_clin_infect_don == "Y" ~ 1,
                                                    cov_clin_infect_don == "N" ~ 0,
                                                     TRUE ~ NA)),
         cov_cod_cad_don = factor(case_when(cov_cod_cad_don == "1" ~ "Anoxia",
                                            cov_cod_cad_don == "2" ~ "Stroke",
                                            cov_cod_cad_don == "3" ~ "Head Trauma",
                                            cov_cod_cad_don == "4" ~ "CNS Tumor",
                                            cov_cod_cad_don == "999" ~ "Other",
                                            TRUE ~ NA_character_)),
         cov_cold_isch = factor(case_when(cov_cold_isch < 8 ~ "<8",
                                          cov_cold_isch >= 8 ~ "≥8",
                                          TRUE ~ NA_character_)),
         cov_creat_don = factor(case_when(cov_creat_don < 0.3 ~ "<0.3",
                                          cov_creat_don >= 0.3 & cov_creat_don < 0.7 ~ "0.3-0.7",
                                          cov_creat_don >= 0.7 ~ "≥0.7",
                                          TRUE ~ NA_character_),
                                levels = c("<0.3","0.3-0.7","≥0.7")),
         cov_creat_tx = factor(case_when(cov_creat_tx < 0.3 ~ "<0.3",
                                         cov_creat_tx >= 0.3 & cov_creat_tx < 0.7 ~ "0.3-0.7",
                                         cov_creat_tx >= 0.7 ~ "≥0.7",
                                         TRUE ~ NA_character_),
                               levels = c("<0.3","0.3-0.7","≥0.7")),
         cov_ddavp_don = as.logical(case_when(cov_ddavp_don == "Y" ~ 1,
                                          cov_ddavp_don == "N" ~ 0,
                                          TRUE ~ NA)),
         cov_death_circum_don = factor(case_when(cov_death_circum_don == 1 ~ "Motor Vehicle Accident",
                                                 cov_death_circum_don == 2 ~ "Suicide",
                                                 cov_death_circum_don == 3 ~ "Homicide",
                                                 cov_death_circum_don == 4 ~ "Child Abuse",
                                                 cov_death_circum_don == 5 ~ "Accident (Non-MVA)",
                                                 cov_death_circum_don == 6 ~ "Natural Causes",
                                                 cov_death_circum_don == 997 ~ "Other",
                                                 TRUE ~ NA)),
         cov_death_mech_don = factor(case_when(cov_death_mech_don == 1 ~ "Drowning",
                                               cov_death_mech_don == 2 ~ "Seizure",
                                               cov_death_mech_don == 3 ~ "Drug Intoxication",
                                               cov_death_mech_don == 4 ~ "Asphyxiation",
                                               cov_death_mech_don == 5 ~ "Cardiovascular",
                                               cov_death_mech_don %in% c(995, 7, 8) ~ "Gunshot or Stab Wound",
                                               cov_death_mech_don == 9 ~ "Blunt Injury",
                                               cov_death_mech_don == 10 ~ "SIDS",
                                               cov_death_mech_don == 11 ~ "Intrcranial Hemorrhage",
                                               cov_death_mech_don == 12 ~ "Natural Causes",
                                               cov_death_mech_don %in% c(6, 997) ~ "Other Causes",
                                               TRUE ~ NA)),
         cov_diab = as.logical(case_when(cov_diab == "1" ~ 0,
                                         cov_diab %in% c("2", "3", "4", "5") ~ 1,
                                         TRUE ~ NA)),
         cov_diabetes_don = as.logical(case_when(cov_diabetes_don == "Y" ~ 1,
                                                 cov_diabetes_don == "N" ~ 0,
                                                 TRUE ~ NA)),
         cov_diag = factor(case_when(cov_diag %in% c("4100", "4101", "4102", 
                                                     "4103", "4104", "4105", 
                                                     "4106", "4107", "4108", 
                                                     "4110", "4200") ~ "Acute Liver Failure",
                                     cov_diag %in% c("4300", "4301", "4302", 
                                                     "4303", "4304", "4305", 
                                                     "4306", "4307", "4308", 
                                                     "4315") ~ "Metabolic Liver Disease",
                                     cov_diag %in% c("4230", "4231", "4235", 
                                                     "4240", "4241", "4242", 
                                                     "4245", "4250", "4255", 
                                                     "4260", "4264") ~ "Cholestatic Liver Disease",
                                     cov_diag %in% c("4270", "4271", "4272", 
                                                     "4275") ~ "Biliary Atresia",
                                     cov_diag %in% c("4400", "4401", "4402", 
                                                     "4403", "4404", "4405", 
                                                     "4410", "4420", "4430", 
                                                     "4450", "4451", "4455") ~ "Liver Cancers",
                                     cov_diag %in% c("4201", "4202", "4203", 
                                                     "4204", "4205", "4206", 
                                                     "4207", "4208", "4209", 
                                                     "4210", "4212", "4213", 
                                                     "4214", "4215", "4216", 
                                                     "4217", "4218", "4219", 
                                                     "4220", "4230", "4231", 
                                                     "4235", "4240", "4241", 
                                                     "4242", "4245", "4250", 
                                                     "4255", "4260", "4264", 
                                                     "4265", "4270", "4271", 
                                                     "4272", "4275", "4280", 
                                                     "4285", "4290", "4300", 
                                                     "4301", "4302", "4303", 
                                                     "4304", "4305", "4306", 
                                                     "4307", "4308", "4315", 
                                                     "4400", "4401", "4402", 
                                                     "4403", "4404", "4405", 
                                                     "4410", "4420", "4430", 
                                                     "4450", "4451", "4455", 
                                                     "4500", "4510", "4520", 
                                                     "4592", "4593") ~ "Viral Hepatitis",
                                     cov_diag %in% c("4210", "4212", "4213", 
                                                     "4214", "4215", "4216", 
                                                     "4217", "4218", "4219", 
                                                     "4220", "4230", "4231", 
                                                     "4235", "4240", "4241", 
                                                     "4242", "4245") ~ "Autoimmune Disease",
                                     cov_diag %in% c("4214", "4215", "4216", 
                                                     "4217", "4218", "4219") ~ "NASH/Alcoholic Liver Disease",
                                     cov_diag %in% c("4208", "4213", "4265", 
                                                     "4280", "4285", "4290", 
                                                     "4500", "4510", "4520", 
                                                     "4597", "4598", "999") ~ "Other")),
         cov_dial_tx = as.logical(case_when(cov_dial_tx == "Y" ~ 1,
                                            cov_dial_tx == "N" ~ 0,
                                            TRUE ~ NA)),
         cov_drmis = factor(cov_drmis),
         cov_ebv_serostatus = as.logical(case_when(cov_ebv_igg_cad_don == "P" | 
                                                      cov_ebv_igm_cad_don == "P" | 
                                                      cov_ebv_serostatus == "P" ~ 1,
                                                   cov_ebv_igg_cad_don == "N" | 
                                                     cov_ebv_igm_cad_don == "N" | 
                                                     cov_ebv_serostatus == "N" ~ 0,
                                                    TRUE ~ NA)),
         cov_enceph_tx = factor(case_when(cov_enceph_tx == 1 ~ "None",
                                          cov_enceph_tx == 2 ~ "1-2",
                                          cov_enceph_tx == 3 ~ "3-4",
                                          cov_enceph_tx == 4 ~ "N/A",
                                          TRUE ~ NA_character_),
                                levels = c("N/A", "None", "1-2", "3-4")),
         cov_ethcat = factor(case_when(cov_ethcat == 1 ~ "White, Non-Hispanic",
                                       cov_ethcat == 2 ~ "Black, Non-Hispanic",
                                       cov_ethcat == 4 ~ "Hispanic Latino",
                                       cov_ethcat == 5 ~ "Asian, Non-Hispanic",
                                       cov_ethcat == 6 ~ "Amer Indian/Alaska Native, Non-Hispanic",
                                       cov_ethcat == 7 ~ "Native Hawaiian/Pacific Islander, Non-Hispanic",
                                       cov_ethcat == 9 ~ "Multiracial, Non-Hispanic",
                                       TRUE ~ NA_character_)),
         cov_ethcat_don = factor(case_when(cov_ethcat_don == 1 ~ "White, Non-Hispanic",
                                           cov_ethcat_don == 2 ~ "Black, Non-Hispanic",
                                           cov_ethcat_don == 4 ~ "Hispanic Latino",
                                           cov_ethcat_don == 5 ~ "Asian, Non-Hispanic",
                                           cov_ethcat_don == 6 ~ "Amer Indian/Alaska Native, Non-Hispanic",
                                           cov_ethcat_don == 7 ~ "Native Hawaiian/Pacific Islander, Non-Hispanic",
                                           cov_ethcat_don == 9 ~ "Multiracial, Non-Hispanic",
                                           TRUE ~ NA_character_)),
         cov_func_stat_trr = factor(case_when(cov_func_stat_trr == 996 | elig_age < 1 ~ "N/A (age < 1 year)",
                                              cov_func_stat_trr %in% c(3, 4010, 4020, 4030) ~ "10%-30%",
                                              cov_func_stat_trr %in% c(2, 4040, 4050) ~ "40%-60%",
                                              cov_func_stat_trr %in% c(1, 4060, 4070, 4080, 4090, 4100) ~ "70%-100%",
                                              TRUE ~ NA_character_),
                                    levels = c("10%-30%",
                                               "40%-60%",
                                               "70%-100%",
                                               "N/A (age < 1 year)")),
         cov_gender = factor(cov_gender),
         cov_gender_don = factor(cov_gender_don),
         cov_hist_cancer_don = as.logical(case_when(cov_hist_cancer_don == "Y" ~ 1,
                                                cov_hist_cancer_don == "N" ~ 0,
                                                TRUE ~ NA)),
         cov_hist_cig_don = as.logical(case_when(cov_hist_cig_don == "Y" ~ 1,
                                                 cov_hist_cig_don == "N" ~ 0,
                                                 TRUE ~ NA)),
         cov_hist_cocaine_don = as.logical(case_when(cov_hist_cocaine_don == "Y" ~ 1,
                                                     cov_hist_cocaine_don == "N" ~ 0,
                                                     TRUE ~ NA)),
         cov_hist_diabetes_don = as.logical(case_when(cov_diab == "1" ~ 0,
                                                      cov_diab %in% c("2", "3", "4", "5") ~ 1,
                                                      TRUE ~ NA)),
         cov_hist_hypertens_don = as.logical(case_when(cov_hist_hypertens_don == "Y" ~ 1,
                                                       cov_hist_hypertens_don == "N" ~ 0,
                                                       TRUE ~ NA)),
         cov_hist_insulin_dep_don = as.logical(case_when(cov_hist_insulin_dep_don == "Y" ~ 1,
                                                         cov_hist_insulin_dep_don == "N" ~ 0,
                                                         TRUE ~ NA)),
         cov_hist_oth_drug_don = as.logical(case_when(cov_hist_oth_drug_don == "Y" ~ 1,
                                                      cov_hist_oth_drug_don == "N" ~ 0,
                                                      TRUE ~ NA)),
         cov_history_mi_don = as.logical(case_when(cov_history_mi_don == "Y" ~ 1,
                                                   cov_history_mi_don == "N" ~ 0,
                                                   TRUE ~ NA)),
         cov_hlamis = factor(case_when(cov_hlamis %in% c(0:2) ~ "0-2",
                                       cov_hlamis %in% 3 ~ "3",
                                       cov_hlamis %in% 4 ~ "4",
                                       cov_hlamis %in% 5 ~ "5",
                                       cov_hlamis %in% 6 ~ "6",
                                       TRUE ~ NA_character_)),
         cov_inotrop_support_don = as.logical(case_when(cov_inotrop_support_don == "Y" ~ 1,
                                                        cov_inotrop_support_don == "N" ~ 0,
                                                        TRUE ~ NA)),
         cov_inr_tx = factor(case_when(cov_inr_tx < 2 ~ "<2",
                                       cov_inr_tx >= 2 & cov_inr_tx <= 3 ~ "2-3",
                                       cov_inr_tx > 3 ~ ">3",
                                       TRUE ~ NA_character_),
                             levels = c("<2", "2-3",">3")),
         cov_insulin_don = as.logical(case_when(cov_insulin_don == "Y" ~ 1,
                                                cov_insulin_don == "N" ~ 0,
                                                TRUE ~ NA)),
         cov_intracranial_cancer_don = as.logical(case_when(cov_intracranial_cancer_don == "Y" ~ 1,
                                                            cov_intracranial_cancer_don == "N" ~ 0,
                                                            TRUE ~ NA)),
         cov_life_sup_trr = as.logical(case_when(cov_life_sup_trr == "Y" ~ 1,
                                                 cov_life_sup_trr == "N" ~ 0,
                                                 TRUE ~ NA)),  
         cov_lv_eject_don = factor(case_when(cov_lv_eject_don < 50 ~ "<50%",
                                             cov_lv_eject_don >= 50 & cov_lv_eject_don <= 75 ~ "50%-75%",
                                             cov_lv_eject_don > 75 ~ "75%",
                                             TRUE ~ NA_character_),
                                   levels = c("<50%","50%-75%","75%")),
         cov_malig_trr = as.logical(case_when(cov_malig_trr == "Y" ~ 1,
                                              cov_malig_trr == "N" ~ 0,
                                              TRUE ~ NA)),  
         cov_med_cond_trr = factor(case_when(cov_med_cond_trr == 1 ~ "Hospitalized - ICU",
                                             cov_med_cond_trr == 2 ~ "Hospitalized - Non-ICU",
                                             cov_med_cond_trr == 3 ~ "Not Hospitalized",
                                             TRUE ~ NA_character_)),
         cov_meld_peld_lab_score = factor(case_when(cov_meld_peld_lab_score < 12 ~ "<12",
                                                    cov_meld_peld_lab_score %in% c(12:30) ~ "12-30",
                                                    cov_meld_peld_lab_score > 30 ~ ">30",
                                                    TRUE ~ NA_character_),
                                          levels = c("<12", "12-30", ">30")),
         cov_non_hrt_don = as.logical(case_when(cov_non_hrt_don == "Y" ~ 1,
                                                cov_non_hrt_don == "N" ~ 0,
                                                TRUE ~ NA)),
         cov_on_vent_trr = as.logical(cov_on_vent_trr),
         cov_other_inf_conf_don = as.logical(case_when(cov_other_inf_conf_don == "Y" ~ 1,
                                                       cov_other_inf_conf_don == "N" ~ 0,
                                                       TRUE ~ NA)),
         cov_other_inf_don = as.logical(cov_other_inf_don),
         cov_portal_vein_trr = as.logical(case_when(cov_portal_vein_trr == "Y" ~ 1,
                                                    cov_portal_vein_trr == "N" ~ 0,
                                                    TRUE ~ NA)),
         cov_prev_ab_surg_trr = as.logical(case_when(cov_prev_ab_surg_trr == "Y" ~ 1,
                                                     cov_prev_ab_surg_trr == "N" ~ 0,
                                                     TRUE ~ NA)),
         cov_protein_urine = as.logical(case_when(cov_protein_urine == "Y" ~ 1,
                                                  cov_protein_urine == "N" ~ 0,
                                                  TRUE ~ NA)),
         cov_pt_diuretics_don = as.logical(case_when(cov_pt_diuretics_don == "Y" ~ 1,
                                                     cov_pt_diuretics_don == "N" ~ 0,
                                                     TRUE ~ NA)),
         cov_pt_steroids_don = as.logical(case_when(cov_pt_steroids_don == "Y" ~ 1,
                                                    cov_pt_steroids_don == "N" ~ 0,
                                                    TRUE ~ NA)),
         cov_pt_t3_don = as.logical(case_when(cov_pt_t3_don == "Y" ~ 1,
                                              cov_pt_t3_don == "N" ~ 0,
                                              TRUE ~ NA)),
         cov_pt_t4_don = as.logical(case_when(cov_pt_t4_don == "Y" ~ 1,
                                              cov_pt_t4_don == "N" ~ 0,
                                              TRUE ~ NA)),
         cov_pulm_inf_conf_don = as.logical(case_when(cov_pulm_inf_conf_don == "Y" ~ 1,
                                                      cov_pulm_inf_conf_don == "N" ~ 0,
                                                      TRUE ~ NA)),
         cov_pulm_inf_don = as.logical(cov_pulm_inf_don),
         cov_recov_out_us = as.logical(case_when(cov_recov_out_us == "Y" ~ 1,
                                                 cov_recov_out_us == "N" ~ 0,
                                                 TRUE ~ NA)),
         cov_recuscit_dur = factor(case_when(cov_recuscit_dur %in% c(0, NA) ~ "0",
                                             cov_recuscit_dur %in% c(1:15) ~ "1-15",
                                             cov_recuscit_dur %in% c(16:30) ~ "16-30",
                                             cov_recuscit_dur %in% c(31:45) ~ "31-45",
                                             cov_recuscit_dur %in% c(46:120) ~ ">45",
                                             TRUE ~ NA_character_),
                                   levels = c("0", "1-15", "16-30", "31-45", ">45")),
         cov_sgot_don = factor(case_when(cov_sgot_don <= 40 ~ "≤40",
                                         cov_sgot_don > 40 & cov_sgot_don <= 80 ~ "41-80",
                                         cov_sgot_don > 80 & cov_sgot_don <= 120 ~ "81-120",
                                         cov_sgot_don > 120 & cov_sgot_don <= 160 ~ "121-160",
                                         cov_sgot_don > 161 ~ ">160",
                                         TRUE ~ NA_character_),
                               levels = c("≤40","41-80","81-120","121-160",">160")),
         cov_sgpt_don = factor(case_when(cov_sgpt_don <= 40 ~ "≤40",
                                         cov_sgpt_don > 40 & cov_sgpt_don <= 80 ~ "41-80",
                                         cov_sgpt_don > 80 & cov_sgpt_don <= 120 ~ "81-120",
                                         cov_sgpt_don > 120 & cov_sgpt_don <= 160 ~ "121-160",
                                         cov_sgpt_don > 161 ~ ">160",
                                         TRUE ~ NA_character_),
                               levels = c("≤40","41-80","81-120","121-160",">160")),
         cov_tattoos = as.logical(case_when(cov_tattoos == "Y" ~ 1,
                                            cov_tattoos == "N" ~ 0,
                                            TRUE ~ NA)),
         cov_tbili_don = factor(case_when(cov_tbili_don <= 2 ~ "≤2",
                                          cov_tbili_don > 2 ~ ">2",
                                          TRUE ~ NA_character_),
                                levels = c("≤2",">2")),
         cov_tbili_tx = factor(case_when(cov_tbili_tx <= 2 ~ "≤2",
                                         cov_tbili_tx > 2 ~ ">2",
                                          TRUE ~ NA_character_),
                                levels = c("≤2",">2")),
         cov_tipss_trr = as.logical(case_when(cov_tipss_trr == "Y" ~ 1,
                                              cov_tipss_trr == "N" ~ 0,
                                              TRUE ~ NA)),
         cov_tx_procedur_ty = factor(case_when(cov_tx_procedur_ty %in% c(701, 704) ~ "Whole",
                                               cov_tx_procedur_ty %in% c(702) ~ "Partial",
                                               cov_tx_procedur_ty %in% c(703) ~ "Split")),
         cov_urine_inf_conf_don = as.logical(case_when(cov_urine_inf_conf_don == "Y" ~ 1,
                                                       cov_urine_inf_conf_don == "N" ~ 0,
                                                       TRUE ~ NA)),
         cov_urine_inf_don = as.logical(cov_urine_inf_don),
         cov_vasodil_don = as.logical(case_when(cov_vasodil_don == "Y" ~ 1,
                                                cov_vasodil_don == "N" ~ 0,
                                                TRUE ~ NA)),
         cov_share_ty = factor(case_when(cov_share_ty == 3 ~ "Local",
                                         cov_share_ty == 4 ~ "Regional",
                                         cov_share_ty == 5 ~ "National")),
         cov_exc_case = as.logical(case_when(cov_exc_case == "Yes" ~ 1,
                                             cov_exc_case == "No" ~ 0,
                                             TRUE ~ NA)),
         cov_cmv_don = as.logical(case_when(cov_cmv_don == "P" ~ 1,
                                            cov_cmv_don == "N" ~ 0,
                                            TRUE ~ NA)),
         cov_age = factor(case_when(elig_age %in% c(0:2) ~ "0-2",
                                    elig_age %in% c(3:5) ~ "3-5",
                                    elig_age %in% c(6:8) ~ "6-8",
                                    elig_age %in% c(9:11) ~ "9-11",
                                    elig_age %in% c(12:14) ~ "12-14",
                                    elig_age %in% c(15:17) ~ "15-17",
                                    TRUE ~ NA_character_),
                          levels = c("0-2","3-5","6-8","9-11","12-14","15-17")),
         cov_age_don = factor(case_when(elig_age_don %in% c(0:2) ~ "0-2",
                                        elig_age_don %in% c(3:5) ~ "3-5",
                                        elig_age_don %in% c(6:8) ~ "6-8",
                                        elig_age_don %in% c(9:11) ~ "9-11",
                                        elig_age_don %in% c(12:14) ~ "12-14",
                                        elig_age_don %in% c(15:17) ~ "15-17",
                                        TRUE ~ NA_character_),
                              levels = c("0-2","3-5","6-8","9-11","12-14","15-17"))) %>% 
  
  # Remove variables not needed for analysis
  select(-c(cov_ebv_igg_cad_don, cov_ebv_igm_cad_don, i, out_gstatus, out_gtime,
         out_pstatus, out_ptime, out_composite_death_date, donor_id,
         trr_id_code, cov_academic_prg_trr, cov_academic_level_trr, cov_cmv_don),
         -starts_with("elig_")) %>% 
  
  # Reorder variables
  select(pt_code, outcome, time, everything())

## Missing Data Exploration ----

# Assess percent missing values in each variable
pct_missing <- unlist(lapply(data_liver_clean, function(x) sum(is.na(x))))/nrow(data_liver_clean)
sort(pct_missing[pct_missing > 0], decreasing = TRUE)

# Identify variables with >30% missing values
sort(pct_missing[pct_missing > 0.3], decreasing = TRUE)
VARS_EXCL_MISS <- names(sort(pct_missing[pct_missing > 0.3], decreasing = TRUE))

# Remove variables with >30% missing values
data_liver_clean_miss <- data_liver_clean %>% 
  select(-all_of(VARS_EXCL_MISS)) 

# Identify remaining variables with any missingness
pct_missing <- unlist(lapply(data_liver_clean_miss, function(x) sum(is.na(x))))/nrow(data_liver_clean_miss)
VARS_INCL_MISS <- names(sort(pct_missing[pct_missing > 0], decreasing = TRUE))

## Multiple Imputation ----

### Parameters & Specification ---- 
# Run mice algorithm with no iterations
imp <- mice(data_liver_clean_miss, maxit = 0)

# Extract predictorMatrix and methods of imputation 
predM <- imp$predictorMatrix
meth <- imp$method

# Set patient id variable to zero in predictor matrix because it is non-informative.
predM[, "pt_code"] <- 0

# Specify that the outcome not be included in the imputation model for a set of
# variables to facilitate model convergence.
predM[c(6,16,17,21,25,28,29,33,40,41,45,47,52,53,54,67,72), "outcome"] <- 0

# Specify random forest as the method for imputation of all variables except 
# outcome, which will not be imputed.
meth[!meth == ""] <- "rf"
meth["outcome"] <- ""

### Run Imputation Algorithm -----
# data_imputed <- mice(data_liver_clean_miss,
#                      m = 10,
#                      maxit = 10,
#                      predictorMatrix = predM,
#                      method = meth,
#                      seed = 777,
#                      visitSequence = "monotone")

# Save/Load mids object
# save(data_imputed, file = paste0(DIR_DATA, "data_imputed.Rdata"))
load(file = paste0(DIR_DATA, "data_imputed.Rdata"))

### Algorithm Convergence Check -----
# Review for any problems that occurred during the imputation procedure
data_imputed$loggedEvents

### Inspection of Imputed Values ----
# QC - Evaluate convergence of MICE algorithm by plotting:
# mean + SD of each imputed variable (y-axis) vs. iteration (x-axis)
# for every dataset (separate lines)
# 
# Inspection reveals no major trends in means and that the curves mix well.
# A lack of trend in the means suggests the algorithm converged.
# When inspecting the variability (SD) in the plots, a high degree of variation  
# between the datasets (lines) vs within the datasets indicates that there is a 
# relatively high increase in variance due to missing data. This is not evident.
# Note: Rare binary variables with few missing values (e.g., cig, diabetes, cancer)
# have very low variability, which is to be expected.
plot(data_imputed, layout = c(10,10))

# Extract first imputed dataset and evaluate dimensions of missing values
data_imputed_1 <- complete(data_imputed, 1)
miss_vals_dims(data_imputed_1)

# List variable names for which outcome had to be excluded as a predictor
colnames(predM[,c(6,16,17,21,25,28,29,33,40,41,45,47,52,53,54,67,72)])

# List variable names that contain at least one missing value
colnames(data_liver_clean_miss)[apply(data_liver_clean_miss, 2, anyNA)]

# Analysis ----

### Partition Data for Training & Testing ----
# Only analyze complete cases with no missing outcome
data_imputed_1_cc <- data_imputed_1[complete.cases(data_imputed_1),]

# Create random row index to select from analytic dataset
set.seed(2)
random_index <- sample(x = 1:nrow(data_imputed_1_cc), 
                       size = (0.75 * nrow(data_imputed_1_cc)), 
                       replace = FALSE)

# Select rows from analytic file based on random index to establish each dataset
data_train <- data_imputed_1_cc[random_index, ]
data_test <- data_imputed_1_cc[-random_index, ]

# Create candidate predictor matrix and outcome vector for TRAINING dataset
data_train_x <- model.matrix( ~ ., select(data_train, starts_with("cov_")))[,-1]
data_train_y <- as.numeric(pull(select(data_train, "outcome")))

# Create candidate predictor matrix and outcome vector for TESTING dataset
data_test_x <- model.matrix( ~ ., select(data_test, starts_with("cov_")))[,-1]
data_test_y <- as.numeric(pull(select(data_test, "outcome")))

## LASSO Logistic Model ----

#### Model Fit with Cross-Validation ----
model_lasso_logit <- cv.glmnet(x = data_train_x, 
                               y = data_train_y, 
                               nfolds = 10,
                               type.measure = "deviance", # Use deviance for cross-validation
                               alpha = 1,                 # Specify LASSO penalty
                               family = "binomial", 
                               nlambda = 100)

lambda_min_lasso_logit <- model_lasso_logit$lambda.min

# Plot model deviance as a function of lambda
# plot(model_lasso_logit)

# Use best lambda to predict TEST data
pred_lasso_logit <- predict(model_lasso_logit, 
                            s = "lambda.min",
                            type = 'response',
                            newx = data_test_x) 

#### Extract Coefficients ----
coef_lasso_logit <- coef(model_lasso_logit, s = "lambda.min")

#### Assess Performance -----
# Specify metrics object and values to search
metrics_lasso_logit <- perf_metrics(obs = data_test_y, pred = pred_lasso_logit, cutoff = NA)

cutoff_vals <- seq(min(pred_lasso_logit), max(pred_lasso_logit), 0.001)

# Grid search for best cutoff 
for (i in 1:length(cutoff_vals)) {
  
  if (i == 1){
    metrics_lasso_logit <- perf_metrics(obs = data_test_y, pred = pred_lasso_logit, cutoff = cutoff_vals[i])
  }
  
  metrics_lasso_logit <- rbind(metrics_lasso_logit, perf_metrics(obs = data_test_y, pred = pred_lasso_logit, cutoff = cutoff_vals[i]))
}

# Obtain cutoff that optimizes F1 score
metrics_lasso_logit_optimal <- metrics_lasso_logit %>%
  as.data.frame() %>% 
  arrange(-f1_score) %>% 
  slice(1) 

# AUROC
auc <- roc(response = data_test_y, predictor = as.numeric(pred_lasso_logit))$auc

# Gather all performance metrics for model in single object
model <- "LASSO Logit"
perf_lasso_logit <- cbind(model, metrics_lasso_logit_optimal, as.data.frame(auc))

## Ridge Logistic Model ----

#### Model Fit with Cross-Validation ----
model_ridge_logit <- cv.glmnet(x = data_train_x, 
                               y = data_train_y, 
                               nfolds = 10,
                               type.measure = "deviance", # Use deviance for cross-validation
                               alpha = 0,                 # Specify penalty
                               family = "binomial", 
                               nlambda = 100)

lambda_min_ridge_logit <- model_ridge_logit$lambda.min

# Plot model deviance as a function of lambda
# plot(model_ridge_logit)

# Use best lambda to predict TEST data
pred_ridge_logit <- predict(model_ridge_logit, 
                            s = "lambda.min",
                            type = 'response',
                            newx = data_test_x) 

#### Extract Coefficients ----
coef_ridge_logit <- coef(model_ridge_logit, s = "lambda.min")

#### Assess Performance -----
# Specify metrics object and values to search
metrics_ridge_logit <- perf_metrics(obs = data_test_y, pred = pred_ridge_logit, cutoff = NA)

cutoff_vals <- seq(min(pred_ridge_logit), max(pred_ridge_logit), 0.001)

# Grid search for best cutoff 
for (i in 1:length(cutoff_vals)) {
  
  if (i == 1){
    metrics_ridge_logit <- perf_metrics(obs = data_test_y, pred = pred_ridge_logit, cutoff = cutoff_vals[i])
  }
  
  metrics_ridge_logit <- rbind(metrics_ridge_logit, perf_metrics(obs = data_test_y, pred = pred_ridge_logit, cutoff = cutoff_vals[i]))
}

# Obtain cutoff that optimizes F1 score
metrics_ridge_logit_optimal <- metrics_ridge_logit %>%
  as.data.frame() %>% 
  arrange(-f1_score) %>% 
  slice(1) 

# AUROC - Warning: object 'auc' has same name as prior object(s)
auc <- roc(response = data_test_y, predictor = as.numeric(pred_ridge_logit))$auc

# Gather all performance metrics for model in single object
model <- "Ridge Logit"
perf_ridge_logit <- cbind(model, metrics_ridge_logit_optimal, as.data.frame(auc))

## Elastic Net Logistic Model ----

#### Model Fit with Cross-Validation ----
model_elas_logit <- cv.glmnet(x = data_train_x, 
                              y = data_train_y, 
                              nfolds = 10,
                              type.measure = "deviance",  # Use deviance for cross-validation
                              alpha = 0.25,               # Specify penalty
                              family = "binomial", 
                              nlambda = 100)

lambda_min_elas_logit <- model_elas_logit$lambda.min

# Plot model deviance as a function of lambda
# plot(model_elas_logit)

# Use best lambda to predict TEST data
pred_elas_logit <- predict(model_elas_logit, 
                            s = "lambda.min",
                            type = 'response',
                            newx = data_test_x) 

#### Extract Coefficients ----
coef_elas_logit <- coef(model_elas_logit, s = "lambda.min")

#### Assess Performance -----
# Specify metrics object and values to search
metrics_elas_logit <- perf_metrics(obs = data_test_y, pred = pred_elas_logit, cutoff = NA)

cutoff_vals <- seq(min(pred_elas_logit), max(pred_elas_logit), 0.001)

# Grid search for best cutoff 
for (i in 1:length(cutoff_vals)) {
  
  if (i == 1){
    metrics_elas_logit <- perf_metrics(obs = data_test_y, pred = pred_elas_logit, cutoff = cutoff_vals[i])
  }
  
  metrics_elas_logit <- rbind(metrics_elas_logit, perf_metrics(obs = data_test_y, pred = pred_elas_logit, cutoff = cutoff_vals[i]))
}

# Obtain cutoff that optimizes F1 score
metrics_elas_logit_optimal <- metrics_elas_logit %>%
  as.data.frame() %>% 
  arrange(-f1_score) %>% 
  slice(1) 

# AUROC - Warning: object 'auc' has same name as prior object(s)
auc <- roc(response = data_test_y, predictor = as.numeric(pred_elas_logit))$auc

# Gather all performance metrics for model in single object
model <- "Elastic Net Logit"
perf_elas_logit <- cbind(model, metrics_elas_logit_optimal, as.data.frame(auc))

## Adaptive LASSO Logistic Model ----

#### Model Fit with Cross-Validation ----

# Ascertain Ridge coefficients from above
best_ridge_coef <- as.numeric(coef_ridge_logit)[-1]

# Perform adaptive LASSO with 10-fold cross validation
model_adap_logit <- cv.glmnet(x = data_train_x, 
                              y = data_train_y, 
                              nfolds = 10,
                              type.measure = "deviance", # Use deviance for cross-validation
                              alpha = 1,                 # Specify penalty
                              family = "binomial", 
                              penalty.factor = 1 / abs(best_ridge_coef),
                              keep = TRUE)

lambda_min_adap_logit <- model_adap_logit$lambda.min

# Plot model deviance as a function of lambda
# plot(model_adap_logit)

# Use best lambda to predict TEST data
pred_adap_logit <- predict(model_adap_logit, 
                           s = "lambda.min",
                           type = 'response',
                           newx = data_test_x) 

#### Extract Coefficients ----
coef_adap_logit <- coef(model_adap_logit, s = "lambda.min")

#### Assess Performance -----
# Specify metrics object and values to search
metrics_adap_logit <- perf_metrics(obs = data_test_y, pred = pred_adap_logit, cutoff = NA)

cutoff_vals <- seq(min(pred_adap_logit), max(pred_adap_logit), 0.001)

# Grid search for best cutoff 
for (i in 1:length(cutoff_vals)) {
  
  if (i == 1){
    metrics_adap_logit <- perf_metrics(obs = data_test_y, pred = pred_adap_logit, cutoff = cutoff_vals[i])
  }
  
  metrics_adap_logit <- rbind(metrics_adap_logit, perf_metrics(obs = data_test_y, pred = pred_adap_logit, cutoff = cutoff_vals[i]))
}

# Obtain cutoff that optimizes F1 score
metrics_adap_logit_optimal <- metrics_adap_logit %>%
  as.data.frame() %>% 
  arrange(-f1_score) %>% 
  slice(1) 

# AUROC - Warning: object 'auc' has same name as prior object(s)
auc <- roc(response = data_test_y, predictor = as.numeric(pred_adap_logit))$auc

# Gather all performance metrics for model in single object
model <- "Adaptive LASSO Logit"
perf_adap_logit <- cbind(model, metrics_adap_logit_optimal, as.data.frame(auc))

## MLE Logit Model ----
# Extract LASSO-selected variable names and format to facilitate placement in
# glm() formula
coef_lasso_logit
coef_names <- rownames(coef_lasso_logit)[coef_lasso_logit[,1] != 0]
paste(coef_names, collapse = " + ")

#### Model Fit ----
# Note: Removed 'commented out' variables due to sparsity. If one dummy from 
# a factor was selected by LASSO, then entire factor was included in this model.
model_mle_logit <- glm(outcome ~ cov_arginine_don + #cov_artificial_li_trr + 
                         cov_ascites_tx + cov_blood_inf_don + 
                         cov_bmi_tcr + cov_bmi_don_calc + 
                         cov_bun_don + #cov_cod_cad_don + 
                         cov_cold_isch + cov_death_circum_don + 
                         cov_death_mech_don + 
                         cov_diag + 
                         cov_ebv_serostatus + cov_enceph_tx + 
                         cov_ethcat + 
                         cov_func_stat_trr + 
                         cov_gender + #cov_hist_cocaine_don + 
                         cov_hist_hypertens_don + #cov_history_mi_don + 
                         cov_inr_tx + cov_life_sup_trr + 
                         cov_med_cond_trr + 
                         cov_meld_peld_lab_score + 
                         cov_protein_urine + cov_pt_steroids_don + 
                         cov_pt_t4_don + cov_pulm_inf_don + 
                         cov_recuscit_dur + 
                         cov_sgot_don + 
                         cov_sgpt_don + cov_tattoos + cov_tbili_don + 
                         cov_tx_procedur_ty + cov_urine_inf_don + 
                         cov_vasodil_don + cov_share_ty + cov_age + 
                         cov_age_don,
                       family = binomial(link = 'logit'),
                       data = data_train)

summary(model_mle_logit)

# Predict TEST data
pred_mle_logit <- predict(model_mle_logit, 
                            type = 'response',
                            newdata = data_test) 

#### Extract Coefficients ----
coef_mle_logit <- as.matrix(coef(model_mle_logit))

#### Assess Performance -----
# Specify metrics object and values to search
metrics_mle_logit <- perf_metrics(obs = data_test_y, pred = pred_mle_logit, cutoff = NA)

cutoff_vals <- seq(min(pred_mle_logit), max(pred_mle_logit), 0.001)

# Grid search for best cutoff 
for (i in 1:length(cutoff_vals)) {
  
  if (i == 1){
    metrics_mle_logit <- perf_metrics(obs = data_test_y, pred = pred_mle_logit, cutoff = cutoff_vals[i])
  }
  
  metrics_mle_logit <- rbind(metrics_mle_logit, perf_metrics(obs = data_test_y, pred = pred_mle_logit, cutoff = cutoff_vals[i]))
}

# Obtain cutoff that optimizes F1 score
metrics_mle_logit_optimal <- metrics_mle_logit %>%
  as.data.frame() %>% 
  arrange(-f1_score) %>% 
  slice(1) 

# AUROC
auc <- roc(response = data_test_y, predictor = as.numeric(pred_mle_logit))$auc

# Hosmer-Lemeshow calibration
performance_hosmer(model_mle_logit, n_bins = 10)

# Gather all performance metrics for model in single object
model <- "MLE Logit"
perf_mle_logit <- cbind(model, metrics_mle_logit_optimal, as.data.frame(auc))


## AUROC Plot for Final Model ----
fig_auc_final_logit <-ggplot(data = as.data.frame(metrics_lasso_logit), 
       aes(x = 1 - specificity, y = sensitivity)) +
  geom_smooth(size = 0.5, col = 'gray', se = FALSE) + 
  geom_point(size = .2, alpha = 0.5) + 
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  geom_abline(col = 'red') +
  ggtitle("Area Under the ROC Curve - LASSO Logit Model") + 
  xlab("1 - Specificity") + ylab("Sensitivity") + 
  annotate("text", x = 0.1, y = 0.95, 
           label = paste0("AUC: ", round(perf_lasso_logit$auc, 4)), size = 3, 
           color = "black", hjust = 0.5, vjust = 0.5)

## Precision-Recall Plot for Best Model ----
fig_prerec_final_logit <- ggplot(data = as.data.frame(metrics_lasso_logit), 
       aes(x = sensitivity, y = ppv)) +
  geom_smooth(size = 0.5, col = 'gray', se = FALSE) + 
  geom_point(size = .2, alpha = 0.5) + 
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  ggtitle("Precision vs. Recall Curve - LASSO Logit Model") + 
  xlab("Sensitivity") + ylab("Positive Predictive Value")

## Donor Risk Index ----
#### Marginalized Donor Risk Model ----
# This section sums over recipient factors in the MLE logit model so
# that donor factor coefficients reflect the population as a whole

# Create vector of donor factor variable names
MLE_LOGIT_VARS <- names(model_mle_logit$model)
MLE_LOGIT_DONOR_FACTORS <- c(grep("don", MLE_LOGIT_VARS, value = TRUE),
                             "cov_cold_isch", # addt'l relevant donor factors
                             "cov_share_ty", 
                             "cov_tx_procedur_ty")

MLE_LOGIT_RECIP_FACTORS <- names(model_mle_logit$model)[!names(model_mle_logit$model) %in% 
                                                          c(MLE_LOGIT_DONOR_FACTORS, "outcome")]

# Obtain vector of the means of expected values for each observation in the 
# training dataset
set.seed(576)
marginal_probs <- marginal_prob(data = model_mle_logit$data, 
                                model = model_mle_logit, 
                                vars_marg = MLE_LOGIT_DONOR_FACTORS, 
                                vars_sum = MLE_LOGIT_RECIP_FACTORS)

# Use the above vector of probabilities to create a new binary response variable
data_train_outcome_marg <- data_train %>% 
  mutate(outcome_marg = rbinom(length(marginal_probs), 1, marginal_probs)) %>% 
  select(pt_code, outcome, outcome_marg, everything())

# Specify new model formula to obtain 'marginal' donor coefficients
formula_marg_mle_logit <- as.formula(paste("outcome_marg ~", 
                                           paste(MLE_LOGIT_DONOR_FACTORS, 
                                                 collapse = " + ")))

# Regress the new response variables on the donor factors only in training data
model_mle_logit_marg <- glm(formula_marg_mle_logit, 
                            data = data_train_outcome_marg, 
                            family = binomial(link = 'logit'))

# Summary of the new model
summary(model_mle_logit_marg)

##### Extract Coefficients ----
coef_mle_logit_marg <- as.matrix(coef(model_mle_logit_marg))

##### Estimate Donor Risk Score ----
# Estimate donor risk score of patients in the test dataset
donor_risk_index_test_data <- data_test %>% 
  mutate(pred_probs = predict(model_mle_logit_marg, type = 'response', newdata = data_test),
         pred_link = predict(model_mle_logit_marg, type = 'link', newdata = data_test),
         dri = round(pred_link / 0.5),
         n = 1,
         pred_probs_quantile = ntile(pred_probs, 5)) %>% 
  select(pt_code, n, outcome, pred_probs, pred_probs_quantile, pred_link, dri)

# Estimate observed & predicted risk by quantile of predicted risk
donor_risk_index_quantile <- donor_risk_index_test_data %>% 
  group_by(pred_probs_quantile) %>% 
  summarize(n = sum(n),
            dri_min = min(dri),
            dri_max = max(dri),
            outcome_risk_pred = mean(pred_probs),
            outcome_risk_obs = mean(outcome)) %>% 
  ungroup()

# Estimate observed and predicted risk by donor risk score
donor_risk_index_points <- donor_risk_index_test_data %>% 
  group_by(dri) %>% 
  summarize(n = sum(n),
            outcome_risk_pred = mean(pred_probs),
            outcome_risk_obs = mean(outcome)) %>% 
  ungroup()

# Estimate correlation between DRI-predicted and observed risk by QUANTILE of 
# DRI-predicted risk
corr_spearman <- cor(donor_risk_index_quantile$outcome_risk_pred, 
                     donor_risk_index_quantile$outcome_risk_obs, 
                     method = "spearman")

corr_pearson <- cor(donor_risk_index_quantile$outcome_risk_pred, 
                    donor_risk_index_quantile$outcome_risk_obs, 
                    method = "pearson")

dri_corr_risk <- as.data.frame(cbind(corr_pearson, corr_spearman))
dri_corr_risk

# Histogram of donor risk index 
fig_dri_histogram_test <- ggplot(data = donor_risk_index_test_data) + 
  geom_histogram(aes(x = dri), binwidth = 0.5) + 
  xlab("Donor Risk Index") + ylab("Frequency") +
  scale_x_continuous(limits = c(-11, 12), breaks = seq(-11, 12, 1)) 

summary(donor_risk_index_test_data$dri)

## Plots: Select Donor Factors vs. Predicted Risk -----
# Plot relationship between predicted risk (using DRI model) vs: year of age, 
# cause of death, organ type, organ share type in ENTIRE cohort (training + testing)

# Create data table containing select donor factors and the outcome
data_final_figs <- data_liver_cohort %>% 
  
  # Merge age_in_months variable
  left_join(data_age_months, by = join_by("trr_id_code" == "TRR_ID_CODE")) %>% 
  
  # Extract numeric age variable from (uncleaned) cohort file
  filter(pt_code %in% data_imputed_1_cc$pt_code) %>% 
  select(pt_code, age = elig_age, age_months = AGE_IN_MONTHS) %>% 
  
  # Add variables of interest
  mutate(predicted_outcome = predict(model_mle_logit_marg, 
                                     type = 'response', 
                                     newdata = data_imputed_1_cc),
         cov_death_mech_don = data_imputed_1_cc$cov_death_mech_don,
         cov_tx_procedur_ty = data_imputed_1_cc$cov_tx_procedur_ty,
         cov_share_ty = data_imputed_1_cc$cov_share_ty,
         age_months = factor(case_when(age_months %in% c(0:5) ~ "0-5",
                                age_months %in% c(6:11) ~ "6-11",
                                age_months %in% c(12:17) ~ "12-17",
                                age_months %in% c(18:23) ~ "18-23",
                                age_months %in% c(24:29) ~ "24-29",
                                age_months %in% c(30:35) ~ "30-35",
                                TRUE ~ NA_character_),
                             levels = c("0-5","6-11","12-17","18-23","24-29",
                                        "30-35"))) %>% 
  select(pt_code, predicted_outcome, everything())
     
#### Age (Years) vs. Predicted Risk -----
fig_age_v_pred_risk <- ggplot(data_final_figs, aes(x = factor(age), y = 100 * predicted_outcome)) +
  stat_summary(fun = mean, geom = "point", size = 1) +  # Plot mean
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # Plot standard error
  labs(title = "",
       x = "Age (Years)",
       y = "Graft Failure or Mortality (%)",
       caption = "Note: Estimates shown are means and standard errors calculated in entire cohort.") +
  geom_hline(yintercept = 100 * mean(data_final_figs$predicted_outcome), color = "red", linetype = "dashed", size = 0.5) + 
  theme_minimal() + 
  theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))

#### Age (Months) vs. Predicted Risk -----
x_axis_counts <- data_final_figs %>%
  filter(!is.na(age_months)) %>%
  group_by(age_months) %>%
  summarize(n = n())

# Create a custom function get sample sizes in labels
custom_labels <- x_axis_counts %>%
  mutate(label = paste0(age_months, "\n(n=", n, ")"))

# Merge custom labels into the plot
fig_age_mon_v_pred_risk <- ggplot(filter(data_final_figs, !is.na(age_months)),
                                  aes(x = age_months, y = 100 * predicted_outcome)) +
  stat_summary(fun = mean, geom = "point", size = 1) +  # Plot mean
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # Plot standard error
  labs(title = "",
       x = "Age (Months)",
       y = "Graft Failure or Mortality (%)",
       caption = "Note: Estimates shown are means and standard errors calculated in entire cohort.") +
  geom_hline(yintercept = 100 * mean(data_final_figs$predicted_outcome, na.rm = TRUE),
             color = "red", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic")) +
  scale_x_discrete(breaks = x_axis_counts$age_months,  # Ensure all age_months are on the axis
                   labels = custom_labels$label)  # Add custom labels with sample size

#### Death Mechanism vs. Predicted Risk -----
fig_death_mech_v_pred_risk <- ggplot(data_final_figs, aes(x = reorder(cov_death_mech_don, predicted_outcome, FUN = mean), y = 100 * predicted_outcome)) +
  stat_summary(fun = mean, geom = "point", size = 1) +  # Plot mean
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # Plot standard error
  labs(title = "",
       x = "Cause of Death Mechanism",
       y = "Graft Failure or Mortality (%)",
       caption = "Note: Estimates shown are means and standard errors calculated in entire cohort.") +
  geom_hline(yintercept = 100 * mean(data_final_figs$predicted_outcome), color = "red", linetype = "dashed", size = 0.5) + 
  theme_minimal() + 
  theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

#### Liver Organ Type vs. Predicted Risk ----
fig_tx_proc_ty_v_pred_risk <- ggplot(data_final_figs, aes(x = reorder(cov_tx_procedur_ty, predicted_outcome, FUN = mean), y = 100 * predicted_outcome)) +
  stat_summary(fun = mean, geom = "point", size = 1) +  # Plot mean
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # Plot standard error
  labs(title = "",
       x = "Liver Organ Type",
       y = "Graft Failure or Mortality (%)",
       caption = "Note: Estimates shown are means and standard errors calculated in entire cohort.") +
  geom_hline(yintercept = 100 * mean(data_final_figs$predicted_outcome), color = "red", linetype = "dashed", size = 0.5) + 
  theme_minimal() + 
  theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

#### Organ Share Type vs. Predicted Risk ----
fig_share_ty_v_pred_risk <- ggplot(data_final_figs, aes(x = reorder(cov_share_ty, predicted_outcome, FUN = mean), y = 100 * predicted_outcome)) +
  stat_summary(fun = mean, geom = "point", size = 1) +  # Plot mean
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # Plot standard error
  labs(title = "",
       x = "Organ Share Type",
       y = "Graft Failure or Mortality (%)",
       caption = "Note: Estimates shown are means and standard errors calculated in entire cohort.") +
  geom_hline(yintercept = 100 * mean(data_final_figs$predicted_outcome), color = "red", linetype = "dashed", size = 0.5) + 
  theme_minimal() + 
  theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Output ----

## DRI by Quartile ----
donor_risk_index_quantile_clean <- donor_risk_index_quantile %>% 
  rename(`DRI-Predicted Risk Quartile` = pred_probs_quantile,
         `Mean Predicted Risk` = outcome_risk_pred,
         `Mean Observed Risk` = outcome_risk_obs) 

write.csv(donor_risk_index_quantile_clean, 
          file = paste0(DIR_OUTPUT, "Predicted and Observed Risk by Quantile.csv"),
          row.names = FALSE)

## DRI by Points ----
donor_risk_index_points_clean <- donor_risk_index_points %>% 
  rename(`Mean Predicted Risk` = outcome_risk_pred,
         `Mean Observed Risk` = outcome_risk_obs,
         `Donor Risk Index` = dri)

write.csv(donor_risk_index_points_clean, 
          file = paste0(DIR_OUTPUT, "Predicted and Observed Risk by DRI Points.csv"),
          row.names = FALSE)

## DRI Correlation to Observed Risk ----
write.csv(dri_corr_risk, 
          file = paste0(DIR_OUTPUT, "Correlation of DRI-Predicted and Observed Risk.csv"),
          row.names = FALSE)

## DRI Histogram in Test Dataset -----
ggsave(paste0(DIR_OUTPUT, "Figure - DRI Histogram in Test Dataset.pdf"),
       plot = fig_dri_histogram_test)

## Select Factors vs. Predicted Risk of Graft Failure Figures -----
ggsave(paste0(DIR_OUTPUT, "Figure - Age vs Predicted Risk.pdf"),
       plot = fig_age_v_pred_risk, height = 4.26, width = 6)

ggsave(paste0(DIR_OUTPUT, "Figure - Age (Months) vs Predicted Risk.pdf"),
       plot = fig_age_mon_v_pred_risk, height = 4.26, width = 6)

ggsave(paste0(DIR_OUTPUT, "Figure - COD Mechanism vs Predicted Risk.pdf"),
       plot = fig_death_mech_v_pred_risk, height = 4.26, width = 6)

ggsave(paste0(DIR_OUTPUT, "Figure - Organ Type vs Predicted Risk.pdf"),
       plot = fig_tx_proc_ty_v_pred_risk, height = 4.26, width = 6)

ggsave(paste0(DIR_OUTPUT, "Figure - Organ Share Type vs Predicted Risk.pdf"),
       plot = fig_share_ty_v_pred_risk, height = 4.26, width = 6)


## Coefficients ----

# Clean coefficient tables and convert to data frame
coef_lasso_logit_clean <- as.data.frame(as.matrix(coef_lasso_logit)) %>% 
  rownames_to_column() %>% 
  rename(coef_lasso_logit = s1,
         variable = rowname)

coef_ridge_logit_clean <- as.data.frame(as.matrix(coef_ridge_logit)) %>% 
  rownames_to_column() %>% 
  rename(coef_ridge_logit = s1,
         variable = rowname) 

coef_elas_logit_clean <- as.data.frame(as.matrix(coef_elas_logit)) %>% 
  rownames_to_column() %>%
  rename(coef_elas_logit = s1,
         variable = rowname) 

coef_adap_logit_clean <- as.data.frame(as.matrix(coef_adap_logit)) %>% 
  rownames_to_column() %>%
  rename(coef_adap_logit = s1,
         variable = rowname) 

coef_mle_logit_clean <- as.data.frame(coef_mle_logit) %>% 
  rownames_to_column() %>%
  rename(coef_mle_logit = V1,
         variable = rowname)

# Combine all coefficient tables into a single data frame
coefficients_clean <- coef_lasso_logit_clean %>% 
  left_join(coef_ridge_logit_clean, by = 'variable') %>% 
  left_join(coef_elas_logit_clean, by = 'variable') %>% 
  left_join(coef_adap_logit_clean, by = 'variable') %>% 
  left_join(coef_mle_logit_clean, by = 'variable') %>% 
  mutate(across(2:6, ~ na_if(., 0)))

write.csv(coefficients_clean, 
          file = paste0(DIR_OUTPUT, "Model Coefficients.csv"),
          row.names = FALSE)

## Marginalized Donor Factors in MLE Logit Model -----
coef_mle_logit_marg_clean <- as.data.frame(coef_mle_logit_marg) %>% 
  rownames_to_column() %>% 
  rename(coef_mle_logit_marg = V1,
         variable = rowname)

write.csv(coef_mle_logit_marg_clean, 
          file = paste0(DIR_OUTPUT, "Marginalized Model Coefficients for DRI.csv"),
          row.names = FALSE)

## All Model Performance Metrics ----
performance_clean <- as.data.frame(rbind(perf_lasso_logit, 
                                         perf_ridge_logit, 
                                         perf_elas_logit, 
                                         perf_adap_logit, 
                                         perf_mle_logit))

write.csv(performance_clean, 
          file = paste0(DIR_OUTPUT, "All Model Performance Metrics.csv"),
          row.names = FALSE)

# Final Model Performance by Cutoff
write.csv(metrics_lasso_logit, 
          file = paste0(DIR_OUTPUT, "Final Model Performance Metrics.csv"),
          row.names = FALSE)

# ROC and Precision-Recall Plots 
ggsave(paste0(DIR_OUTPUT, "Figure - ROC Curve Final Model.pdf"),
       plot = fig_auc_final_logit)

ggsave(paste0(DIR_OUTPUT, "Figure - Precision Recall Curve Final Model.pdf"),
       plot = fig_prerec_final_logit)

## Baseline Tables & Variable Summary----
### After Imputation ----
# Specify variables to be summarized
VARS_CAND_PREDICTORS <- data_imputed_1_cc %>% 
  select(starts_with("cov_")) %>% 
  names()

# Generate baseline table in analytic dataset (after imputation)
baseline_table <- CreateTableOne(vars = VARS_CAND_PREDICTORS, data = data_imputed_1_cc)
baseline_table_print <- print(baseline_table, showAllLevels = TRUE)

# Export tables to a CSV file
write.csv(baseline_table_print, 
          file = paste0(DIR_OUTPUT, "Baseline Table - Analytic Cohort.csv"))

### Before Imputation (Missing Values) ----
# Extract dataset prior to imputation and select same patients as in analytic  
# file for comparability
data_imputed_0 <- complete(data_imputed, 0)
data_imputed_0_cc <- filter(data_imputed_0, pt_code %in% data_imputed_1_cc$pt_code)

# Generate baseline table in data prior to imputation 
# Note: Percentages for covariates shown are calculated AFTER removing missing values
baseline_table_miss <- CreateTableOne(vars = VARS_CAND_PREDICTORS, data = data_imputed_0_cc)
baseline_table_miss_print <- print(baseline_table_miss, showAllLevels = TRUE, missing = TRUE)

# Export tables to a CSV file
write.csv(baseline_table_miss_print, 
          file = paste0(DIR_OUTPUT, "Baseline Table with Missing Value Distribution.csv"))


### Complete vs. Missing Outcome -----
# Generate baseline table in analytic dataset (after imputation), stratified by 
# missing indicator of outcome 
data_imputed_1_cc_ind <- data_imputed_1 %>% 
  mutate(`Missing Outcome` = case_when(is.na(outcome) == TRUE ~ 1,
                                       TRUE ~ 0))
baseline_table_1_cc_ind <- CreateTableOne(vars = VARS_CAND_PREDICTORS, 
                                          data = data_imputed_1_cc_ind,
                                          strata = "Missing Outcome")
baseline_table_1_cc_ind_print <- print(baseline_table_1_cc_ind, showAllLevels = TRUE,
                                       test = TRUE,
                                       smd = TRUE)

# Export tables to a CSV file
write.csv(baseline_table_1_cc_ind_print, 
          file = paste0(DIR_OUTPUT, "Baseline Table Stratified by Missing Outcome Indicator.csv"))

### Training vs. Testing Dataset -----
# Generate baseline table stratified by training vs testing dataset
data_imputed_1_cc_tt_ind <- data_imputed_1_cc %>% 
  mutate(Training = case_when(pt_code %in% data_train$pt_code ~ 1,
                              TRUE ~ 0)) 
baseline_table_1_cc_tt_ind <- CreateTableOne(vars = VARS_CAND_PREDICTORS, 
                                          data = data_imputed_1_cc_tt_ind,
                                          strata = "Training")
baseline_table_1_cc_tt_ind_print <- print(baseline_table_1_cc_tt_ind, showAllLevels = TRUE,
                                       test = TRUE,
                                       smd = TRUE)

# Export tables to a CSV file
write.csv(baseline_table_1_cc_tt_ind_print, 
          file = paste0(DIR_OUTPUT, "Baseline Table Stratified by Training Dataset Indicator.csv"))


# Testing, QC, and Notes ----

