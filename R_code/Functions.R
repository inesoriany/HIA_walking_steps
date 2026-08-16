#################################################
################## FUNCTIONS ####################
#################################################

################################################################################################################################
################################################################################################################################
#                                                         0. DATASET                                                           #
################################################################################################################################
################################################################################################################################

##############################################################
#                        INCIDENCE                           #
##############################################################
# FUNCTION age10_categ : Create 10-year age categories from 5-year age categories
age10_categ <- function(age_grp) {
  age_min <- as.numeric(sub("-.*", "", age_grp))
  paste0(floor(age_min / 10) * 10, "-", floor(age_min / 10) * 10 + 9)
}


##############################################################
#                      WALK DATASET                          #
##############################################################

# FUNCTION age_grp10 : 10 year category
age_grp_10 = function(age){
  age_grp = gsub(",", "-",
                 gsub("\\[|\\]|\\(|\\)", "",
                      cut( age, breaks = seq(0,150, by = 10), include.lowest = T, right = F)))
  post = sub(".*-","",age_grp)
  age_grp = paste0(sub("-.*", "", age_grp),
                   "-", as.numeric(post)-1)
  return(age_grp)
}


# FUNCTION walk_dataset : Creating a dataset combining the emp dataset with diseases incidence and mortality data, and calculating walking exposure.
walk_dataset <- function(data, diseases_10, insee, dis_vec, 
                            walk_dist_var = "nbkm_tot_walking", 
                            walk_dist_jour_var = "nbkm_tot_walking_jour", 
                            step_length, walk_speed) {
  
  # Re-write sexe as female and male and convert as factors
  data<- data %>% 
    mutate(sexe = as.character(sexe)) %>%                                 
    mutate(sexe = fct_recode(sexe, "Male" = "1", "Female" = "2")) %>%     
    rename(sex = sexe)  %>% 
  
  # Daily steps
    mutate(step_commute = .data[[walk_dist_var]] / step_length)  %>% 
    mutate(step_commute_jour = .data[[walk_dist_jour_var]] / step_length)  %>% 
  
  # Day time spent walking (min)
    mutate(day_time = .data[[walk_dist_var]] * 60 / walk_speed)  %>% 
    mutate(day_time_jour = .data[[walk_dist_jour_var]] * 60 / walk_speed)  %>% 
  
  # Create age categories
    mutate(age_grp10 = age_grp_10(age)) %>%
    filter(age >= 20 & age < 90)  %>%                       
  
  # Add population counts per sex and age group
    left_join(diseases_10 %>% 
        select(pop_age_grp10, sex, age_grp10),     
      by = c("sex", "age_grp10"))  %>% 
  
  # Add diseases incidences (include all *_incidence_mid columns present in diseases_10)
    left_join(diseases_10 %>% 
        select(sex, age_grp10, matches("_incidence_mid$")),    
      by = c("sex", "age_grp10")                                                                                               ) 
  
  # Add mortality rates
  insee <- insee %>% 
    rename(sex = sexe)
  
  data <- data %>% 
    left_join(
      insee %>% select(MR, sex, age),       
      by = c("sex", "age")                  
    ) %>% 
    rename(mort_rate = MR) %>%
  
  # Calculate death incidence
    mutate(mort_incidence_mid = mort_rate * pop_age_grp10)
  
  # Calculate disease incidence rates
  for (dis in dis_vec) {
    data <- data %>%
      mutate(
        !!paste0(dis, "_rate") := if_else(
          !is.na(pop_age_grp10),
          .data[[paste0(dis, "_incidence_mid")]] / pop_age_grp10,
          NA_real_
        )
      )
  }
  
  # Add life-expectancy for each sex
  data <- data %>%
    mutate(life_exp = if_else(sex == "Female", 85.99324, 79.59503)) %>%
    # Add the years of life remaining, potentially affected by diseases or premature death
    mutate(years_remaining = pmax(life_exp - age, 0))
  
  return(data)
}



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
generate_RR_distrib = function (RR, low, up, N) {          # N : number of random values
  lRR <- log(RR)                                            # Conversion in log scale
  l_low <- log(low)
  l_up <- log(up)
  
  sd1 <- (lRR - l_low) / qnorm(1-0.05/2)
  sd2 <- (l_up - lRR) / qnorm(1-0.05/2) 
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
#                                                       I. MONTE CARLO                                                         #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                                 1. Health impact assessment                                                  #
################################################################################################################################

##############################################################
#                    DISEASE INCIDENCE                       #
##############################################################

# FUNCTION incidence_replicate : Randomly associate incidence to individuals
  # set.seed()
incidence_replicate = function(data_list, incidence_distrib_table, dis_vec) {
  
  for (dis in dis_vec) {
    dis_data <- data_list[[dis]]
    
    if (dis == "mort") {
      dis_data[[dis]] <- dis_data 

    } else {

      # Randomly associate incidence rate for each individual based on age_grp10 and sex
      rates_age_sex <- incidence_distrib_table %>%
        filter(disease == dis) %>%
        group_by(age_grp10, sex) %>%
        summarise(rates = list(rate), .groups = "drop")

      dis_data <- dis_data %>%
          left_join(rates_age_sex, by = c("age_grp10", "sex")) %>%
          mutate(rate = map_dbl(rates, ~ if (length(.x) > 0) sample(.x, 1) else NA_real_)) %>%
          select(-rates)
      
      data_list[[dis]] <- dis_data
    }
  }
  return(data_list)
}



##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################
# FUNCTION reduction_risk_replicate : Randomly associate risk reductions to individuals
  # set.seed()
reduction_risk_replicate <- function(data_list, reduction_risk_distrib_table, dis_vec) {

  for (dis in dis_vec) {
    dis_data <- data_list[[dis]]

      # Randomly associate incidence rate for each individual based on step
      reduction_risk_step <- reduction_risk_distrib_table %>%
        filter(disease == dis) %>%
        group_by(step) %>%
        summarise( reduction_risk_values = list(reduction_risk), .groups = "drop")

      dis_data <- dis_data %>%
        left_join(reduction_risk_step, by = c("step")) %>%
        mutate(reduction_risk = map_dbl(reduction_risk_values, ~ if(length(.x) > 0)sample(.x, 1)
          else
          NA_real_)) %>%
        select(-reduction_risk_values)
      
  
      if (dis == "mort") {
        dis_data <- dis_data %>%
          mutate(reduction_risk = if_else(reduction_risk > (1-0.45), (1-0.45), reduction_risk))    # Cap mortality reduction to 45%
      }


      data_list[[dis]] <- dis_data
    }

  return(data_list)
}


##############################################################
#                           CASES                            #
##############################################################

# FUNCTION reduc_incidence : Calculate the reduced disease incidence (number of prevented new cases)
reduc_incidence = function(data) {
 data <- data %>% 
   mutate(cases = reduction_risk * rate)
  
  return(data)
}



##############################################################
#                   RANDOM ASSOCIATION                       #
##############################################################
# FUNCTION value_associate : Randomly associate values to individuals
  # set.seed()
value_associate = function(data_list, distrib_table, dis_vec, variable) {
  
  for (dis in dis_vec) {
    
    dis_val <- distrib_table %>%
      filter(disease == dis) %>%
      pull(paste0("simulated_", variable))
    
    data_list[[dis]] <- data_list[[dis]] %>%
      mutate(!!sym(variable) := sample(dis_val, size = n(), replace = TRUE))
  }
  return(data_list)
}



##############################################################
#                           DALY                             #
##############################################################
# Goal : To know the number of sick or death years prevented for each individual by walking

# FUNCTION daly : Calculate DALY (Disability-Adjusted Life Years) for each disease
daly = function(data, dep_duration_table , dis, prop_relapse, duration_recovery) { 
  if (dis == "dep") {

    # Randomly associate depression duration for each individual and compute DALY
    duration_age_sex <- dep_duration_table %>%
      filter(disease == dis) %>%
      group_by(age_grp10, sex) %>%
      summarise(duration_values = list(simulated_duration_dep), .groups = "drop")

    data <- data %>%
      left_join(duration_age_sex, by = c("age_grp10", "sex")) %>%
      mutate(duration_dep = map_dbl(duration_values, ~ if (length(.x) > 0) sample(.x, 1) else NA_real_)) %>%
      select(-duration_values) %>%
      mutate(daly = cases * dw * prop_relapse * pmin(years_remaining, duration_dep) / 
                    (duration_recovery + duration_dep) * years_remaining)
    
      
  } else {
    data <- data %>%
      mutate(daly = years_remaining * dw * cases)
  }
  return(data) 
}




##############################################################
#                      ECONOMIC IMPACT                       #
##############################################################
# FUNCTION medic_costs : Calculate the medical costs associated with the reduced disease incidence for each individual
medic_costs = function(data, dis) {
  data [[paste0("medic_costs")]] <- get(paste0(dis, "_cost")) * data[[paste0("cases")]] 
  
  return(data)
}



##############################################################
#                        CALCULATE HIA                       #
##############################################################
# FUNCTION calc_HIA_replicate : Calculate the disease reduction percentage, cases, DALY and medical costs prevented for 1 run
  # set.seed()
calc_HIA_replicate = function(data_list, incidence_distrib_table, dep_duration_table, reduction_risk_distrib_table, dw_distrib_table, dis_vec, 
                              prop_relapse, duration_recovery, vsl) {
  

  # 1. Disease incidence
  if (!is.null(incidence_distrib_table)) {
    data_list <- incidence_replicate(data_list = data_list,
                                    incidence_distrib_table = incidence_distrib_table,
                                    dis_vec = dis_vec)
  }
  
  # 2. Disease reduction risk 
  data_list <- reduction_risk_replicate(data_list = data_list,
                                        reduction_risk_distrib_table = reduction_risk_distrib_table,
                                        dis_vec = dis_vec)
  
 # 3. Disability weights
  data_list <- value_associate(data_list = data_list,
                               distrib_table = dw_distrib_table,
                               dis_vec = dis_vec,
                               variable = "dw")
  
  for (dis in dis_vec) {
    
    dis_data <- data_list[[dis]]

    # 4. Cases prevented
    dis_data <- reduc_incidence(dis_data)
    
    # 5. DALY
    dis_data <- daly(dis_data, dep_duration_table, dis, prop_relapse, duration_recovery)
    
    # 6. Economic impact
    dis_data <- medic_costs(dis_data, dis)
    
    dis_data <- dis_data %>%
      mutate(soc_costs = daly * vsl)
    
    data_list[[dis]] <- dis_data
  }
  
  return(data_list)
}


##############################################################
#                        HIA OUTCOMES                        #
##############################################################
# FUNCTION burden_prevented_replicate : Prevented cases, DALY and saved medical costs
burden_prevented_replicate = function(data_list, dis_vec, group) {
  
burden_run <- data.frame()

  for (dis in dis_vec) {
    dis_data <- data_list[[dis]]
    dis_burden_replicate <- list() 
    
      # Survey design
      surv_dis_replicate <- dis_data %>% 
        as_survey_design(ids = ident_ind, weights = pond_indc)
      
      # Burden for 1 disease
      dis_burden_replicate <- surv_dis_replicate %>% 
        group_by(across(all_of(group))) %>%
        summarise(
          tot_cases = survey_total(cases, na.rm = TRUE),
          tot_daly = survey_total(daly, na.rm = TRUE),
          tot_medic_costs = survey_total(medic_costs, na.rm = TRUE),
          tot_soc_costs = survey_total(soc_costs, na.rm = TRUE)) %>% 
        mutate(disease = dis)
    
    burden_run <- bind_rows(burden_run, dis_burden_replicate)         # Results for all diseases for 1 run
  }
  
  return(burden_run)
}


# FUNCTION HIA_burden_total : HIA calculations and Total of prevented cases, DALY and saved medical costs for N simulations
  # set.seed()
HIA_burden_total = function(data_list, function_calc_HIA, incidence_distrib_table, dep_duration_table, reduction_risk_distrib_table, dw_distrib_table, 
                            dis_vec, 
                            prop_relapse, duration_recovery, vsl, group, N, show_progress = TRUE) {
  
  burden_total <- data.frame()
  burden_total_list <- list()
  
  # Progress bar setup
  progress_bar <- NULL
  if (show_progress && interactive() && N > 0) {
    progress_bar <- progress::progress_bar$new(
      format = "  Simulation :current/:total [:bar] :percent ETA: :eta",
      total = N,
      clear = FALSE,
      width = 60
    )
  }
  
  for (i in 1:N) {
    data_list_replicate <- function_calc_HIA(data_list, incidence_distrib_table, dep_duration_table, reduction_risk_distrib_table, dw_distrib_table, dis_vec, prop_relapse, duration_recovery, vsl)
    burden_total_list[[i]] <- burden_prevented_replicate(data_list_replicate, dis_vec, group)  %>% 
      mutate(simulation_id = i)
    
    # Update progress bar
    if (!is.null(progress_bar)) progress_bar$tick()
  }
  burden_total <- bind_rows(burden_total_list)

  return(burden_total)
}




################################################################################################################################
#                                                  2. Uncertainty analysis                                                     #
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

# FUNCTION HIA_burden_IC : Get a table with HIA outcomes and IC 
  # set.seed() if use calc_replicate_IC
  HIA_burden_IC = function(data, dis_vec, outcome_vec, IC_func){
  
  HIA_burden <- data.frame()
  
  presence_age <- "age_grp10" %in% names(data)
  presence_sex <- "sex" %in% names(data)
  presence_area <- "area_type" %in% names(data)
  

  age_vec <- if (presence_age) unique(data$age_grp10) else NA
  sex_vec <- if (presence_sex) unique(data$sex) else NA
  area_vec <- if (presence_area) unique(data$area_type) else NA
  
  for (dis in dis_vec){
    for (age in age_vec){
      for (sexe in sex_vec){
        for(area in area_vec) {
        
          data_sub <- data
          
          # Filter disease
          data_sub <- filter(data_sub, disease == dis)
          
          # Filter age (if column exists)
          if (presence_age && !is.na(age)) {
            data_sub <- filter(data_sub, age_grp10 == age)
          }
          
          # Filter sex (if column exists)
          if (presence_sex && !is.na(sexe)) {
            data_sub <- filter(data_sub, sex == sexe)
          }
          
          # Filtre area_type (if column exists)
          if (presence_area && !is.na(area)) {
            data_sub <- filter(data_sub, area_type == area)
          }

          if (nrow(data_sub) == 0) next
          
          # Construire ligne dynamiquement
          row <- data.frame(disease = dis)
          
          if (presence_age) row$age_grp10 <- age
          if (presence_sex) row$sex <- if (!is.na(sexe)) sexe else NA
          if (presence_area) row$area_type <- area
          
          for (out in outcome_vec){
            IC <- IC_func(data_sub, out)
            
            row[[out]] <- round(IC[2], 3)
            row[[paste0(out, "_low")]] <- round(IC[1], 3)
            row[[paste0(out, "_up")]] <- round(IC[3], 3)
          }
          
          HIA_burden <- dplyr::bind_rows(HIA_burden, row)
        }
      }
    }
  }
  
  return(HIA_burden)
}



################################################################################################################################
#                                                3. Economic unit value (€)                                                    #
################################################################################################################################

# FUNCTION unit_value : Calculate the economic value of 1 km walked
unit_value = function(km, km_low, km_up, euro, euro_low, euro_up, N = 1000) {
  km_sd1 <- (km - km_low) / qnorm(1-0.05/2)
  km_sd2 <- (km_up - km) / qnorm(1-0.05/2)
  km_sd <- mean(c(km_sd1, km_sd2))
  distr_km <- rnorm(N, km, km_sd)
  
  euro_sd1 <- (euro - euro_low) / qnorm(1-0.05/2)
  euro_sd2 <- (euro_up - euro) / qnorm(1-0.05/2)
  euro_sd <- mean(c(euro_sd1, euro_sd2))
  distr_euro <- rnorm(N, euro, euro_sd)
  
  distr_unit <- distr_euro / distr_km
  
  return(distr_unit)
}


# FUNCTION euro_km_value : Calculate distance walked to save 1€
euro_km_unit = function(km, km_low, km_up, euro, euro_low, euro_up, N = 1000) {
  
  km_sd1 <- (km - km_low) / qnorm(1-0.05/2)
  km_sd2 <- (km_up - km) / qnorm(1-0.05/2)
  km_sd <- mean(c(km_sd1, km_sd2))
  distr_km <- rnorm(N, km, km_sd)
  
  euro_sd1 <- (euro - euro_low) / qnorm(1-0.05/2)
  euro_sd2 <- (euro_up - euro) / qnorm(1-0.05/2)
  euro_sd <- mean(c(euro_sd1, euro_sd2))
  distr_euro <- rnorm(N, euro, euro_sd)
  
  distr_unit <- distr_km / distr_euro
  
  return(distr_unit)
}


# FUNCTION euro_step_value : Calculate number of steps to save 1€
euro_step_unit = function(step, step_low, step_up, euro, euro_low, euro_up, N = 1000) {
  
  step_sd1 <- (step - step_low) / qnorm(1-0.05/2)
  step_sd2 <- (step_up - step) / qnorm(1-0.05/2)
  step_sd <- mean(c(step_sd1, step_sd2))
  distr_step <- rnorm(N, step, step_sd)
  
  euro_sd1 <- (euro - euro_low) / qnorm(1-0.05/2)
  euro_sd2 <- (euro_up - euro) / qnorm(1-0.05/2)
  euro_sd <- mean(c(euro_sd1, euro_sd2))
  distr_euro <- rnorm(N, euro, euro_sd)
  
  distr_unit <- distr_step / distr_euro
  
  return(distr_unit)
}




################################################################################################################################
################################################################################################################################
#                                                  4. SENSITIVITY ANALYSIS                                                     #
################################################################################################################################
################################################################################################################################
# FUNCTION dis_setting : Get the corresponding parameters for each disease
dis_setting = function (dis) {
  rr_women <-  get(paste0("rr_", dis, "_women"))
  rr_women_lb <-  get(paste0("rr_", dis, "_women_lb"))
  rr_women_ub <-  get(paste0("rr_", dis, "_women_ub"))
  rr_men <- get(paste0("rr_", dis, "_men"))
  rr_men_lb <-  get(paste0("rr_", dis, "_men_lb"))
  rr_men_ub <-  get(paste0("rr_", dis, "_men_ub"))
  ref_women <-  get(paste0("ref_", dis, "_w"))
  ref_men <- get(paste0("ref_", dis, "_m"))
  return(data.frame("rr_men" = rr_men, "rr_men_lb" = rr_men_lb, "rr_men_ub" = rr_men_ub,
                    "rr_women" = rr_women, "rr_women_lb" = rr_women_lb, "rr_women_ub" = rr_women_ub,
                    "ref_women" = ref_women, "ref_men" = ref_men))
}




# FUNCTION calc_HIA_replicate : Calculate the disease reduction percentage, cases, DALY and medical costs prevented for 1 run
  # set.seed()
calc_alt_HIA <- function(data_list, incidence_distrib_table, dep_duration_table, rr_distrib_table, dw_distrib_table, 
                          dis_vec,
                          prop_relapse, duration_recovery, vsl) {
  
  # 1. Disease incidence
  if (!is.null(incidence_distrib_table)) {
    data_list <- incidence_replicate(data_list = data_list,
                                     incidence_distrib_table = incidence_distrib_table,
                                     dis_vec = dis_vec)}
  
  # 2. Relative risks
  data_list <- value_associate(data_list = data_list,
                               distrib_table = rr_distrib_table,
                               dis_vec = dis_vec,
                               variable = "rr")
  
  # 3. Disability weights
  data_list <- value_associate(data_list = data_list,
                               distrib_table = dw_distrib_table,
                               dis_vec = dis_vec,
                               variable = "dw")
  
  
  for (dis in dis_vec) {
    
    dis_data <- data_list[[dis]]
    
    # 4. Disease reduction risk
    dis_data <- dis_data %>%
      mutate(reduction_risk = 1 - exp(log(rr) * week_time / ref))
    
    # Cap mortality reduction to 45%
    if (dis == "mort") {
      dis_data <- dis_data %>%
        mutate(reduction_risk = if_else(reduction_risk > (1 - 0.45), 1 - 0.45, reduction_risk))}
    
    # 5. Cases prevented
    dis_data <- reduc_incidence(dis_data) %>%
      mutate(cases = ifelse(cases > 0.40, 0.40, cases))
    

    # 6. DALY
    dis_data <- daly(dis_data, dep_duration_table, dis, prop_relapse, duration_recovery)

    
    # 7. Economic impact
    dis_data <- medic_costs(dis_data, dis)
    
    dis_data <- dis_data %>%
      mutate(soc_costs = daly * vsl)
    
    data_list[[dis]] <- dis_data
  }
  
  return(data_list)
}