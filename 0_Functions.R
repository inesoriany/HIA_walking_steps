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



# FUNCTION generate_RR : Generate random RR values in a normal distribution based on existing RR and their IC (Monte-Carlo)
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
#                DISEASE REDUCTION RISK                      #
##############################################################

## DISEASE RISK REDUCTION ----
# FUNCTION reduction_risk : Calculate the disease risk reduction percentage for each individual with a linear regression
# (% of decrease in disease risk comparing to the baseline : if people did not walk)

reduction_risk = function(data, dis, bound, data_rr) {
    rr_obs <- sym(paste0(dis, "_rr"))  
    
    rr2000 <- 1
    
      # Risk reduction
    data <- data %>%
      mutate(
        !!paste0(dis, "_reduction_risk") := rr2000 / !!rr_obs
      )
      
    return(data)
}
    



## REDUCED DISEASE INCIDENCE ----
# FUNCTION reduc_incidence : Calculate the reduced disease incidence (number of prevented new cases)
reduc_incidence = function (data, incidence_rate, reduction_risk, dis) {
  for (i in 1:nrow(data)) {
    data[i, paste0(dis,"_reduc_incidence")] <- data[i, incidence_rate] * data[i, reduction_risk]
    
    if (!is.na(data[i, paste0(dis,"_reduc_incidence")]) & data[i, paste0(dis,"_reduc_incidence")] > 0.30) {             
      data[i, paste0(dis,"_reduc_incidence")] <-  0.30                                    # cap reduction to 30%
    }
  }  
  return(data)
}



##############################################################
#                           DALY                             #
##############################################################
# Goal : To know the number of sick or death years prevented for each individual by walking

#FUNCTION daly : Calculate DALY (Disability-Adjusted Life Years) for each disease
daly = function(data, dis, bound) { 
  data[[paste0(dis, "_daly")]] <- data$years_remaining * get((paste0(dis, "_dw_", bound))) * data[[paste0(dis, "_reduc_incidence")]] 
  return(data) }




##############################################################
#                      ECONOMIC IMPACT                       #
##############################################################

## MEDICAL COSTS ----
# FUNCTION medic_costs : Calculate the medical costs associated with the reduced disease incidence for each individual
medic_costs = function(data, dis) {
  data [[paste0(dis, "_medic_costs")]] <- get(paste0(dis, "_cost")) * data[[paste0(dis, "_reduc_incidence")]]
  return(data)
}



##############################################################
#                        CALCULATE HIA                       #
##############################################################

# FUNCTION calc_HIA : Calculate the disease reduction percentage, reduced incidence, DALY and medical costs prevented for each individual
calc_HIA = function(data, bound, data_rr, dis_vec){
  
  for (dis in dis_vec) {
    
    # 1. Reduction risk
    data <- reduction_risk (data, dis, bound, data_rr)
    
    # 2. Reduced incidence
    dis_incidence_rate <- paste0(dis, "_rate")
    dis_reduction_risk <- paste0(dis, "_reduction_risk")
    
    data <- reduc_incidence(data, dis_incidence_rate, dis_reduction_risk, dis) 
    
    # 3. DALY prevented  
    data <- daly(data, dis, bound)
    
    # 4. Medical costs prevented
    data <- medic_costs(data, dis)
  }
  
  return(data)
}




##############################################################
#                        HIA OUTCOMES                        #
##############################################################
# FUNCTION burden_prevented : Total of prevented cases, DALY and saved costs, for each disease
burden_prevented = function(data, dis, group){
  
  dis_burden <- data %>% 
    group_by(across(all_of(group))) %>% 
    summarise(tot_cases = survey_total(!!sym(paste0(dis, "_reduc_incidence")), na.rm = TRUE),       # Total of prevented cases per disease
              tot_daly = survey_total(!!sym(paste0(dis, "_daly")), na.rm = TRUE),                 # Total of prevented DALY per disease
              tot_medic_costs = survey_total(!!sym(paste0(dis, "_medic_costs")), na.rm = TRUE)    # Total of saved medical costs per disease
    ) %>%             
    mutate(disease = dis)
  
  return(dis_burden)
}







