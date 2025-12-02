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
  rr_women <-  get(paste0("rr_", dis, "_women"))
  rr_women_low <-  get(paste0("rr_", dis, "_women_low"))
  rr_women_up <-  get(paste0("rr_", dis, "_women_up"))
  rr_men <- get(paste0("rr_", dis, "_men"))
  rr_men_low <-  get(paste0("rr_", dis, "_men_low"))
  rr_men_up <-  get(paste0("rr_", dis, "_men_up"))
  ref_women <-  get(paste0("ref_", dis, "_w"))
  ref_men <- get(paste0("ref_", dis, "_m"))
  return(data.frame("rr_men" = rr_men, "rr_men_low" = rr_men_low, "rr_men_up" = rr_men_up,
                    "rr_women" = rr_women, "rr_women_low" = rr_women_low, "rr_women_up" = rr_women_up,
                    "ref_women" = ref_women, "ref_men" = ref_men))
}


##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################

## DISEASE RISK REDUCTION ----


# FUNCTION log_reduction_risk : Calculate the disease risk reduction percentage for each individual with a log linear regression
# (% of decrease in disease risk comparing to the baseline : if people did not walk)
log_reduction_risk = function(data, dis, rr_women, rr_women_low, rr_women_up, rr_men, rr_men_low, rr_men_up, ref_women, ref_men) {
  data[["relative_rr_mid"]] <- ifelse(                                       # Calculate risk reduction percentage
    data$sex == "Female",
    1-(exp(log(rr_women) * data$week_time / ref_women)),                     # for women % of decrease for this disease risk
    1-(exp(log(rr_men) * data$week_time / ref_men)) )                        # for men % of decrease for this disease risk
  
  data[["relative_rr_low"]] <- ifelse(                                       # Low bound    
    data$sex == "Female",
    1-(exp(log(rr_women_low) * data$week_time / ref_women)),                       
    1-(exp(log(rr_men_low) * data$week_time / ref_men)) )                          
  
  data[["relative_rr_up"]] <- ifelse(                                       # Upper bound    
    data$sex == "Female",
    1-(exp(log(rr_women_up) * data$week_time / ref_women)),                       
    1-(exp(log(rr_men_up) * data$week_time / ref_men)) )                         
  
  return(data)
}





## REDUCED DISEASE INCIDENCE ----

# FUNCTION reduc_incidence : Calculate the reduced disease incidence (number of prevented new cases)
reduc_incidence <- function(data) {
  
 data <- data %>% 
   mutate(cases_mid = rate * relative_rr_mid,
          cases_low = rate * relative_rr_low,
          cases_up = rate * relative_rr_up)
  
  return(data)
}




##############################################################
#                           DALY                             #
##############################################################
# Goal : To know the number of sick or death years prevented for each individual by walking

# FUNCTION daly : Calculate DALY (Disability-Adjusted Life Years) for each disease
daly = function(data) { 
  data[["daly"]] <- data[["years_remaining"]] * data[["dw"]] * data[["reduc_incidence"]]
  return(data) }




##############################################################
#                      ECONOMIC IMPACT                       #
##############################################################

## MEDICAL COSTS ----

# FUNCTION medic_costs : Calculate the medical costs associated with the reduced disease incidence for each individual
medic_costs = function(data, dis, bound) {
  data [[paste0("medic_costs_", bound)]] <- get(paste0(dis, "_cost")) * data[[paste0("cases_", bound)]] 
  return(data)
}



##############################################################
#                        CALCULATE HIA                       #
##############################################################

# FUNCTION calc_HIA : Calculate the reduced incidence, DALY and medical costs prevented for each individual



#calc_HIA = function(data, bound, data_rr, dis_vec){
  
 




##############################################################
#                        HIA OUTCOMES                        #
##############################################################
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







