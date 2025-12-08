#################################################
################## FUNCTIONS ####################
#################################################



################################################################################################################################
################################################################################################################################
#                                                   0. DRF & RESAMPLING                                                        #
################################################################################################################################
################################################################################################################################

# FUNCTION interpolate_rr :
interpolate_rr <- function(df, disease, metric) {
  
  column <- paste0(disease, "_", metric)
  
  # Select the corresponding RR column
  df_sub <- df %>%
    filter(disease == !!disease) %>%
    select(step, rr = all_of(column))
  
  # Complete steps
  df_complete <- df_sub %>% complete(step = seq(0, 12000, by = 10))
  
  # Existing points
  x <- df_complete$step[!is.na(df_complete$rr)]
  y <- df_complete$rr[!is.na(df_complete$rr)]
  
  # Interpolation
  interp <- case_when(
    # Quadratic model
    disease %in% c("mort", "cvd") ~ spline(x, y, xout = df_complete$step, method = "fmm")$y,
    # Cubic model
    disease == "dem"              ~ spline(x, y, xout = df_complete$step, method = "natural")$y,
    # Linear model
    TRUE                          ~ approx(x, y, xout = df_complete$step, method = "linear", rule = 2)$y
  )
  
  # Add the interpolated column
  df_complete <- df_complete %>%
    mutate(rr_interpolated = interp,
           disease = disease,
           metric = metric)
  
  return(df_complete)
}



# FUNCTION generate_RR : Generate random RR (dw) values in a normal distribution based on existing RR and their IC (Monte-Carlo)
#set.seed()
generate_RR_distrib = function (RR, low, sup, N) {          # N : number of random values
  lRR <- log(RR)                                            # Conversion in log scale
  l_low <- log(low)
  l_sup <- log(sup)
  
  sd1 <- (lRR - l_low) / qnorm(1-0.05/2)
  sd2 <- (l_sup - lRR) / qnorm(1-0.05/2) 
  sd <- mean( c(sd1, sd2))                                  # Estimation of standard deviation assuming symmetrical confidence intervals
  
  distr_RR <- exp(rnorm(N, lRR, sd))                        # Generation of log-normal distribution (random samples)
                               
  distr_RR[distr_RR<0]=0                                    # just need to truncat values
  
  return(distr_RR)                                          # Return simulated RR value
}


# FUNCTION graph_sim_DRF : Graphical representation of the n possible DRF from normal distributions
graph_sim_DRF = function (dis, data) {
  graph_drf_sim_dis <- ggplot(data %>% 
                                 filter(disease == dis,
                                        step %in% 0:12000),
                               aes(x = step,
                                   y = rr_interpolated,
                                   group = simulation_id,
                                   color = disease))+
    scale_color_manual(values = colors_disease[dis])+
    geom_line(na.rm = TRUE,
              alpha = 0.05)+
    labs(title = names_disease[dis],
         x = "Steps per day",
         y = "RR")+
    theme(legend.position = "none")
  
  return(graph_drf_sim_dis)
}


# FUNCTION graph_DRF : Graphical representation of the mean DRF (with IC95)
graph_DRF <- function(dis, data, rr_mean, rr_lci, rr_uci) {
  graph_drf_dis <- data %>%
    filter(disease == dis, step %in% 0:12000) %>%
    ggplot(aes(x = step)) +
    geom_ribbon(aes(
      ymin = !!sym(rr_lci),
      ymax = !!sym(rr_uci)
    ),
    fill = colors_disease[dis],
    alpha = 0.3
    ) +
    geom_line(aes(y = !!sym(rr_mean)),
              color = colors_disease[dis],
              linewidth = 1) +
    labs(
      title = names_disease[dis],
      x = "Steps per day",
      y = "RR"
    ) +
    theme(legend.position = "none")
  
  return(graph_drf_dis)
}




################################################################################################################################
################################################################################################################################
#                                               1. HEALTH IMPACT ASSESSMENT                                                    #
################################################################################################################################
################################################################################################################################

##############################################################
#                    DISEASE PARAMETERS                      #
##############################################################

# FUNCTION dis_setting : Get the corresponding parameters for each disease
dis_setting = function (dis) {
  rr_women_mid <-  get(paste0("rr_", dis, "_women"))
  rr_women_low <-  get(paste0("rr_", dis, "_women_low"))
  rr_women_up <-  get(paste0("rr_", dis, "_women_up"))
  rr_men_mid <- get(paste0("rr_", dis, "_men"))
  rr_men_low <-  get(paste0("rr_", dis, "_men_low"))
  rr_men_up <-  get(paste0("rr_", dis, "_men_up"))
  ref_women <-  get(paste0("ref_", dis, "_w"))
  ref_men <- get(paste0("ref_", dis, "_m"))
  return(data.frame("rr_men_mid" = rr_men_mid, "rr_men_low" = rr_men_low, "rr_men_up" = rr_men_up,
                    "rr_women_mid" = rr_women_mid, "rr_women_low" = rr_women_low, "rr_women_up" = rr_women_up,
                    "ref_women" = ref_women, "ref_men" = ref_men))
}


##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################

## -------------------------------------------------------
## CALCULATE RISK REDUCTION FOR BREAST CANCER)
## -------------------------------------------------------

# FUNCTION log_reduction_risk : Calculate the disease risk reduction percentage for each individual with a log linear regression
  # (% of decrease in disease risk comparing to the baseline : 2000 steps per day)

log_reduction_risk = function(data, dis, rr_women, rr_men, ref_women, ref_men, week_base) {
  
  # --- For women ---
  RR_baseline_w <- exp(log(rr_women) * (week_base / ref_women))
  RR_obs_w    <- exp(log(rr_women) * ((data$week_time + week_base) / ref_women))
  RR_rel_w      <- RR_obs_w / RR_baseline_w
  
  # --- For men ---
  RR_baseline_m <- exp(log(rr_men) * (week_base / ref_men))
  RR_obs_m    <- exp(log(rr_men) * ((data$week_time + week_base)/ ref_men))
  RR_rel_m      <- RR_obs_m / RR_baseline_m
  
  # --- Risk reduction compared to baseline ---
  data[["reduction_risk"]] <- ifelse(
    data$sex == "Female",
    1 - RR_rel_w,
    1 - RR_rel_m
  )
  
  return(data)
}



## -------------------------------------------------------
## ASSOCIATION RISK REDUCTIONS TO INDIVIDUALS
## -------------------------------------------------------
reduction_risk <- function(data_list, rr_table, week_base, dis_vec, bound_vec) {
  
  result_list <- list()
  
  for(bound in bound_vec) {
    bound_list <- list()
    
    for(dis in dis_vec) {
      dis_data <- data_list[[bound]][[dis]]
      
      if(dis == "bc") {
        params <- dis_setting(dis)
        
        dis_data <- log_reduction_risk(
          data      = dis_data,
          dis       = dis,
          rr_women  = params[[paste0("rr_women_", bound)]],
          rr_men    = params[[paste0("rr_men_", bound)]],
          ref_women = params$ref_women,
          ref_men   = params$ref_men,
          week_base = week_base
        )
        
      } else {
        dis_rr <- rr_table %>%
          filter(disease == dis) %>%
          select(step, reduction_risk = !!sym(paste0("reduction_risk_", bound)))
        
        dis_data <- dis_data %>%
          left_join(dis_rr, by = "step")
      }
      
      bound_list[[dis]] <- dis_data
    }
    
    result_list[[bound]] <- bound_list
  }
  
  return(result_list)
}




##############################################################
#                           CASES                            #
##############################################################

# FUNCTION reduc_incidence : Calculate the reduced disease incidence (number of prevented new cases)
reduc_incidence <- function(data) {
 data <- data %>% 
   mutate(reduc_incidence = reduction_risk * rate)
  
  return(data)
}




##############################################################
#                           DALY                             #
##############################################################
# Goal : To know the number of sick or death years prevented for each individual by walking

# FUNCTION daly : Calculate DALY (Disability-Adjusted Life Years) for each disease
daly = function(data) { 
  data <- data %>% 
    mutate(daly = years_remaining * dw * reduc_incidence)

  return(data) 
}




##############################################################
#                      ECONOMIC IMPACT                       #
##############################################################

## MEDICAL COSTS ----

# FUNCTION medic_costs : Calculate the medical costs associated with the reduced disease incidence for each individual
medic_costs = function(data, dis) {
  data [[paste0("medic_costs")]] <- get(paste0(dis, "_cost")) * data[[paste0("reduc_incidence")]] 
  
  return(data)
}



##############################################################
#                        CALCULATE HIA                       #
##############################################################

## -------------------------------------------------------
## CENTRAL VALUE ANALYSIS: Cases prevented
## -------------------------------------------------------
# FUNCTION calc_cases: Calculate the disease reduction percentage, reduced incidence, DALY and medical costs prevented for each individual
calc_cases <- function(data_list, rr_table, week_base, dis_vec, bound_vec) {
  
  # 1. Reduction risk
  risk_list <- reduction_risk(data_list, rr_table, week_base, dis_vec, bound_vec)
  
  # 2. Number of prevented cases
  cases_list <- list()
  
  for(bound in bound_vec) {
    bound_cases <- list()
    
    for(dis in dis_vec) {
      dis_data <- risk_list[[bound]][[dis]]
      
      dis_data <- reduc_incidence(dis_data)
      
      bound_cases[[dis]] <- dis_data
    }
    
  
  # 3. Gather results per bound
    cases_list[[bound]] <- bind_rows(bound_cases, .id = "disease")
  }
    
  return(cases_list)
}




## -------------------------------------------------------
## RESAMPLING ANALYSIS
## -------------------------------------------------------
# FUNCTION calc_HIA




##############################################################
#                        HIA OUTCOMES                        #
##############################################################
## -------------------------------------------------------
## CENTRAL VALUE ANALYSIS: Cases prevented
## -------------------------------------------------------

# FUNCTION cases_prevented : Total of prevented cases
cases_prevented <- function(data, bound_vec, dis_vec, group) {
  
  result_list <- list()
  
  for(bound in bound_vec) {
    health_walkers_bound <- data[[bound]]
    
    # Survey design
    health_walkers_svy <- health_walkers_bound %>%
      as_survey_design(ids = ident_ind,
                       strata = c(age_grp10, sex),
                       weights = pond_indc,
                       nest = TRUE)
    
    
    cases_prev_bound <- data.frame()
    for(dis in dis_vec) {
      
      dis_cases <- health_walkers_svy %>%
        filter(disease == dis)
      
        dis_cases <- dis_cases %>%
          group_by(across(all_of(group))) %>%
          summarise(tot_cases = survey_total(reduc_incidence, na.rm = TRUE)) %>%
          ungroup()
      
      # Add disease
      dis_cases <- dis_cases %>%
        mutate(disease = dis)
      
      # Results for all diseases
      cases_prev_bound <- bind_rows(cases_prev_bound, dis_cases)
    }
    
    # Results for all diseases per bound
    result_list[[bound]] <- cases_prev_bound
  }
  return(result_list)
}




## -------------------------------------------------------
## RESAMPLING ANALYSIS
## -------------------------------------------------------

# FUNCTION burden_prevented : Total of prevented cases, DALY and saved medical costs, for each disease
burden_prevented <- function(survey_data, dis, group){
  
  dis_burden <- survey_data %>%
    filter(disease == dis) %>%
    group_by(across(all_of(group))) %>%
    summarise(
      tot_cases     = survey_total(!!sym("cases_mid"), na.rm = TRUE),
      tot_cases_low = survey_total(!!sym("cases_low"), na.rm = TRUE),
      tot_cases_up  = survey_total(!!sym("cases_up"), na.rm = TRUE),
      
      tot_daly      = survey_total(!!sym("daly_mid"), na.rm = TRUE),
      tot_daly_low  = survey_total(!!sym("daly_low"), na.rm = TRUE),
      tot_daly_up   = survey_total(!!sym("daly_up"), na.rm = TRUE),
      
      tot_medic_costs     = survey_total(!!sym("medic_costs_mid"), na.rm = TRUE),
      tot_medic_costs_low = survey_total(!!sym("medic_costs_low"), na.rm = TRUE),
      tot_medic_costs_up  = survey_total(!!sym("medic_costs_up"), na.rm = TRUE)
    ) %>%
    mutate(disease = dis)
  
  return(dis_burden)
}




################################################################################################################################
################################################################################################################################
#                                                       2. RESAMPLING                                                          #
################################################################################################################################
################################################################################################################################


##############################################################
#                        MONTE-CARLO                         #
##############################################################

# FUNCTION calc_replicate_IC : Calculate interval of confidence by combining replications obtained with generated RR samples to generate a posterior distribution for each outcome
#set.seed()
calc_replicate_IC = function(data, outcome){
  vec = c()
  se_name = paste0(outcome, "_se")
  
  for (i in 1:nrow(data)){
    sam = rnorm(n=200, mean = as.numeric(data[i,outcome]), sd = as.numeric(data[i,se_name]) ) # Generation of samples : uncertainty estimation
    vec = c(vec, sam)
  }
  IC = quantile(vec, probs = c(0.025, 0.5, 0.975))
  return(IC)  
}


##############################################################
#                       RUBIN'S RULE                         #
##############################################################

# FUNCTION calc_IC_Rubin : Calculate interval of confidence with using Rubin's rules
calc_IC_Rubin = function(data, outcome){
  zq <- qnorm(1-0.05/2)
  se_name = paste0(outcome, "_se")
  
  theta =  sum(data[,outcome])/nrow(data)                             # Pooled mean differences
  V_w = sum((data[,se_name] )^2)/nrow(data)                           # Within imputation variance
  V_b = sum((data[,outcome] - theta )^2) / (nrow(data)-1)             # Between imputation variance
  
  V_tot = V_w + V_b +  V_b / (nrow(data))                             # Total variance    
  
  IC = (c(theta-zq*sqrt(V_tot), theta,theta+zq*sqrt(V_tot)))     # Confidence interval
  return(IC)  
}



##############################################################
#                 UNCERTAINTIES ANALYSIS                     #
##############################################################

# FUNCTION HIA_burden_age_IC : Get a table with HIA outcomes and IC 
  # set.seed() if use calc_replicate_IC
HIA_burden_IC = function(data, dis_vec, age_vec, sex_vec, outcome_vec, IC_func){
  
  HIA_burden <- data.frame()
  
  # If age_vec is NULL → no loop for age
  if (is.null(age_vec)) age_vec <- NA
  
  # If sex_vec is NULL → no loop for sex
  if (is.null(sex_vec)) sex_vec <- NA
  
  for (dis in dis_vec){
    for (age in age_vec){
      for (sx in sex_vec){
        
        data_sub <- data %>%
          filter(
            disease == dis,
            if (!is.na(age)) age_grp10 == age else TRUE,
            if (!is.na(sx)) sex == sx else TRUE
          )
        
        if (nrow(data_sub) == 0) next
        
        row <- data.frame(
          disease = dis,
          age_grp10 = if (!is.na(age)) age else NA,
          sex = if (!is.na(sx)) sx else NA
        )
        
        for (out in outcome_vec){
          IC <- IC_func(data_sub, out)
          row[[out]] <- round(IC[2], 3)                         # Median
          row[[paste0(out, "_low")]] <- round(IC[1], 3)         # Low bound
          row[[paste0(out, "_sup")]] <- round(IC[3], 3)         # Up bound
        }
        
        HIA_burden <- bind_rows(HIA_burden, row)
      }
    }
  }
  
  return(HIA_burden)
}


