#################################################
############ HEALTH IMPACT ASSESSMENT ###########
#################################################


###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                       HIA - 7 000 STEPS                                                                                 #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  srvyr,        # Survey
  survey,
  ggplot2       # Data visualization
)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Walkers dataset
emp_step <- import(here("data_clean", "EMP_dis_walkers.xlsx"))

# Risk reductions
reduc_rr_table <- import(here("data_clean", "DRF", "reduction_risk_central.xlsx"))

# Disability weights
dw_table <- import(here("data", "dw_table.xlsx"))

# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "bc", "cvd", "cancer", "diab2", "dem", "dep")

# Bound
bound_vec <- c("mid", "low", "up")


################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_step <- emp_step %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = 7000)


# EMP Dataset per disease
RECO_walkers_list <- list()

for(bound in bound_vec) {
  bound_list <- list()
  
  for(dis in dis_vec) {
    dis_walkers <- emp_step %>%
      filter(disease == dis)
    
    bound_list[[dis]] <- dis_walkers
  }
  
  RECO_walkers_list[[bound]] <- bound_list
}



################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

RECO_HIA_list <- calc_HIA(data_list = RECO_walkers_list,
                          rr_table = reduc_rr_table,
                          dw_table = dw_table,
                          week_base = week_base,
                          dis_vec = dis_vec,
                          bound_vec = bound_vec)




################################################################################################################################
#                                   5. HIA OUTCOMES : Total of prevented cases, DALY, costs                                    #
################################################################################################################################

##############################################################
#                        PER DISEASE                         #
##############################################################
RECO_burden <- burden_prevented(data_list = RECO_HIA_list, 
                                dis_vec = dis_vec,
                                bound_vec,
                                group = NULL)


##############################################################
#                          PER SEX                           #
##############################################################
RECO_burden_sex <- burden_prevented(data_list = RECO_HIA_list, 
                                    dis_vec = dis_vec,
                                    bound_vec,
                                    group = "sex")




# Re-organize by decreasing order
RECO_burden_sex_order <- RECO_burden_sex %>% 
  left_join(RECO_burden %>% select(disease, TOTAL_mixed = tot_cases_mid), by = "disease") %>% 
  arrange(desc(tot_cases_mid)) %>%                      
  mutate(disease = factor(disease, levels = unique(disease))) 




## -------------------------------------------------------
## ADDITIONNAL GAINS
## -------------------------------------------------------
# Import 2019 cases prevented
burden_sex_2019 <- import(here("output", "Tables", "2019", "cases_prev_2019_sex.xlsx"))

# Data preparation
burden_2019 <- burden_sex_2019 %>% 
  rename_with(.fn = ~ paste0(.x, "_2019"),
              .cols = -c (disease, sex))

RECO_burden <- RECO_burden_sex %>% 
  rename_with(.fn = ~ paste0(.x, "_RECO"),
              .cols = -c (disease, sex))


# Additional prevented cases for each disease according to sex
add_RECO_burden_sex <- burden_2019 %>%
  left_join(RECO_burden, by = c("sex", "disease"), suffix = c("_2019", "_RECO")) %>%
  mutate(across(
    ends_with("_RECO"),
    ~ . - get(sub("_RECO$", "_2019", cur_column())),
    .names = "{.col}_diff"
  ))




## -------------------------------------------------------
## DECREASING ORDER
## -------------------------------------------------------
order_disease <- burden_sex_2019 %>% 
  group_by(disease) %>% 
  summarise(tot = sum(tot_cases_mid)) %>% 
  arrange(desc(tot)) %>% 
  pull(disease)

IC_2019_cases_sex$disease <- factor(IC_2019_cases_sex$disease,
                                    levels = order_disease)

IC_RECO_cases_sex$disease <- factor(IC_RECO_cases_sex$disease,
                                    levels = order_disease)



##############################################################
#                          GLOBAL                            #
##############################################################
# Total of prevented cases
RECO_cases <- cases_prevented (data = RECO_cases_list,
                                   bound_vec = c("low", "mid", "up"),
                                   dis_vec = dis_vec,
                                   group = NULL)



# Gather results with IC
IC_RECO_cases <- RECO_cases$mid %>% 
  mutate(tot_cases_low = RECO_cases$low[,"tot_cases"], 
         tot_cases_up = RECO_cases$up[,"tot_cases"]) %>% 
  select(disease, tot_cases, tot_cases_se, tot_cases_low, tot_cases_up)



## -------------------------------------------------------
## ADDITIONNAL GAINS
## -------------------------------------------------------
# Import 2019 cases prevented
IC_2019_cases <- import(here("output", "Tables", "2019", "cases_prev_2019.xlsx"))

# Data preparation
IC_2019 <- IC_2019_cases %>% 
  rename_with(.fn = ~ paste0(.x, "_2019"),
              .cols = -c (disease))

IC_RECO <- IC_RECO_cases %>% 
  rename_with(.fn = ~ paste0(.x, "_RECO"),
              .cols = -c (disease))


# Additional prevented cases for each disease according to sex
IC_add_RECO_cases <- IC_2019 %>%
  left_join(IC_RECO, by = "disease", suffix = c("_2019", "_RECO")) %>%
  mutate(across(
    ends_with("_RECO"),
    ~ . - get(sub("_RECO$", "_2019", cur_column())),
    .names = "{.col}_diff"
  ))


## -------------------------------------------------------
## MORBIDITY
## -------------------------------------------------------
IC_RECO_cases %>% 
  filter(disease != "mort") %>% 
  summarise(tot_cases = sum(tot_cases, na.rm = TRUE),
            tot_cases_low = sum(tot_cases_low, na.rm = TRUE),
            tot_cases_up = sum(tot_cases_up, na.rm = TRUE))



## -------------------------------------------------------
## ALL DISEASES
## -------------------------------------------------------
IC_RECO_cases %>% 
  summarise(tot_cases = sum(tot_cases, na.rm = TRUE),
            tot_cases_low = sum(tot_cases_low, na.rm = TRUE),
            tot_cases_up = sum(tot_cases_up, na.rm = TRUE))



################################################################################################################################
#                                                     6. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Cases prevented had the 2019 population walked 7000 steps according to sex, compared to 2019 levels
plot_RECO_cases_prev <- 
  ggplot() +
  geom_bar(data = IC_2019_cases_sex %>%  filter(disease != "bc"),
           mapping = aes(x = disease, y = tot_cases, fill = sex, alpha = "2019 baseline"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  
  geom_errorbar(data = IC_2019_cases_sex %>%  filter(disease != "bc"),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "2019 baseline"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_fill_manual(values = colors_sex) +
  
  
  geom_bar(data = IC_RECO_cases_sex %>%  filter(disease != "bc"), 
           mapping = aes(x = disease, y = tot_cases, fill = sex, alpha = "7000 steps recommendation"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(name   = "Scenario",
                     values = c("2019 baseline" = 1, "7000 steps recommendation" = 0.4)) +
  
  geom_errorbar(data = IC_RECO_cases_sex %>%  filter(disease != "bc"),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "7000 steps recommendation"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_x_discrete(labels = names_disease) + 
  ylab("Cases prevented") +
  xlab("Disease") +
  theme_minimal()

plot_RECO_cases_prev 






################################################################################################################################
#                                                      8. EXPORT DATA                                                          #
################################################################################################################################
# Tables
export(RECO_burden_sex, here("output", "Tables", "7000 steps", "cases_prev_7000steps.xlsx"))
export(add_RECO_burden_sex, here("output", "Tables", "7000 steps", "add_cases_prev_7000steps.xlsx"))

# Plot 
ggsave(here("output", "Plots", "7000 steps", "cases_prev_7000steps.png"), plot = plot_RECO_cases_prev)










