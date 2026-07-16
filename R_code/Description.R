##############################################
################ DESCRIPTION #################
##############################################

# GOALS : Description
  # HEALTH :
    # Health event incidence distribution
  # WALKING : 
    # Mean walking distance by age group and sex, by area type, by revenue
    # Proportion of people by distance walked
    # Rate of deadly accidents
  # DRIVING :
    # Proportion of people reporting any short (<2km) car trip
    # Mean length of short car travel <2km 
    # Distribution of people reporting any short trips (0.5 - 1 - 1.5 - 2km) (mutually exclusive)



################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################

pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data manipulation
  tidyr,        # Data manipulation
  epikit,       # Age categories creation
  survey,       # Survey management
  srvyr,        # Survey management
  ggplot2,      # Data visualization
  rlang
)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# EMP 2019 : distances for walking
emp <- import(here("data", "emp_dataset_walk_individual.xlsx")) 

# EMP 2019 subset for walkers
emp_walkers <- import(here("data_clean", "EMP_walkers.xlsx"))

# EMP 2019 subset for walking trips
emp_walk_trips <- import(here("data_clean", "EMP_walking_trips.xlsx"))


# EMP 2019 subset for car trips
emp_car_trips <- import(here("data_clean", "EMP_car_trips.xlsx"))


################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("R_code", "Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "cancer", "diab2", "dem", "dep")


################################################################################################################################
################################################################################################################################
#                                                      4. DESCRIPTION                                                          #
################################################################################################################################
################################################################################################################################


# Total population
emp_20_89 <-  emp %>% 
  filter(age >= 20 & age <90)

pop_tot <- sum(emp_20_89$pond_indc)
pop_tot


##############################################################
#                        SURVEY DESIGNS                      #
##############################################################

# Survey design ponderated by day
jour_walkers <- emp_walkers %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)

jour_walk_trips <- emp_walk_trips %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10, area_type),
                   nest = TRUE)


# Survey design ponderated by individual
indiv_walkers <- emp_walkers %>% 
  filter (pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE)


indiv_walk_trips <- emp_walk_trips %>% 
  filter (pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10, area_type),
                   nest = TRUE)


################################################################################################################################
#                                                           HEALTH                                                             #
################################################################################################################################

## Mortality rates distribution per sex and age group
mortality <- indiv_walkers %>% 
  group_by(sex, age_grp10) %>% 
  summarise(mort_rate_mean = survey_mean(mort_rate, proportion = TRUE, na.rm = TRUE))

mean_mortality_rates <- ggplot(mortality, aes(x = age_grp10, y = mort_rate_mean, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(0.7), stat = "identity") +
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  labs (y = "Mean mortality rate",
        x ="Age group") +
  theme_minimal()
plot(mean_mortality_rates)

# Export plot
ggsave(here("output", "Plots", "Description", "Diseases", "plot_mean_mortality_rates.png"), plot = mean_mortality_rates)



## Incidence distribution per age and sex
for(dis in dis_vec) {
  dis_distrib <- ggplot() +
    geom_point(data = emp_walkers, 
               mapping = aes(x = .data[["age_grp10"]],
                             y = .data[[paste0(dis, "_incidence")]],
                             color = .data[["sex"]]), 
               size = 1, alpha = 0.5) +
    
    geom_line(data = emp_walkers, 
                mapping = aes(x = .data[["age_grp10"]],
                              y = .data[[paste0(dis, "_incidence")]],
                              group = .data[["sex"]],
                              color = .data[["sex"]]), 
                size = 1) +
    
    scale_color_manual(name = "Sex",
                       values = c("Female" = "darkorange1",
                                  "Male" = "chartreuse4")) +
    
    labs(title = names_disease[[dis]],
         y = "Incidence",
         x = "Age group") +

    theme_minimal() +
    theme(legend.position = "top")
  
  assign(paste0(dis, "_distrib"), dis_distrib)
  print(dis_distrib)
}

# Export
for (dis in dis_vec){
ggsave(here("output", "Plots", "Description", "Diseases", paste0("plot_",dis, "_incidence.png")), plot = get(paste0(dis,"_distrib")))
}





################################################################################################################################
#                                                          WALKING                                                             #
################################################################################################################################


##############################################################
#                  TOTAL WALKED DISTANCE                     #
##############################################################

## Total walked distance in 2019
km_total_2019 <- as.numeric(svytotal(~nbkm_tot_walking, jour_walkers)) *365.25/7                              # Total km per year
km_total_2019_IC <- as.numeric(confint(svytotal(~nbkm_tot_walking, jour_walkers) *365.25/7 ))                 # Confidence interval

km_total_2019 * 1e-9 # billion km
km_total_2019_IC * 1e-9

# In time Per person (min/pers/year)
km_total_2019 / (pop_tot * walk_speed)


## Total walked distance per day in 2019
km_total_day <- svytotal(~nbkm_tot_walking, jour_walkers)                   # Total km per day
km_total_day_IC <- as.numeric(confint(km_total_day))
km_total_day *1e-6
km_total_day_IC * 1e-6


# Total walking distances per day, by age group
svyby(~nbkm_tot_walking, by = ~age_grp10, jour_walkers, svytotal, na.rm = T)  



# Creation of walking distances categories (in case but for now no need)
emp_walkers <- emp_walkers %>% 
  mutate(dist_grp = case_when(
    nbkm_tot_walking < 1                           ~  "0-1 km",
    nbkm_tot_walking >= 1 & nbkm_tot_walking < 2       ~  "1-2 km",
    nbkm_tot_walking >= 2 & nbkm_tot_walking < 5       ~  "2-5 km",
    nbkm_tot_walking >= 5 & nbkm_tot_walking < 10      ~  "5-10 km",
    nbkm_tot_walking >= 10                         ~  "10 km +" 
  ))%>% 
  mutate(dist_grp = as.factor(dist_grp))


## Plot : Proportion of people by distance walked
proportion_km <- emp_walkers %>% 
  filter (pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE) %>% 
  group_by(dist_grp = factor(dist_grp, levels = c("0-1 km", "1-2 km", "2-5 km", "5-10 km", "10 km +"))) %>%   # Put distance range in the right order 
  summarise(total_pondere = survey_total (1, na.rm = TRUE)) %>%         # Sum of ponderation of individual for each distance group
  mutate(proportion = total_pondere / sum(total_pondere))               # Proportion

prop_walkers_km <- ggplot(proportion_km, aes(x = dist_grp, y = proportion)) +
  geom_col() +
  labs(title = "Proportion of people by distance walked",
       x = "Distance range",
       y = "Proportion") +
  theme_minimal()
prop_walkers_km

# Export plot
ggsave(here("output", "Plots", "Description", "Walk", "plot_prop_walkers_km.png"), plot = prop_walkers_km)




# Proportion of distances walked by each sex
prop_sex <-  emp_walkers %>% 
  filter(pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strat = sex,
                   nest = TRUE) %>% 
  group_by (sex = sex) %>% 
  summarise(tot_km = survey_total(nbkm_tot_walking)) %>% 
  mutate(proportion = tot_km / sum(tot_km))






##############################################################
#                    MEAN WALKED DISTANCE                    #
##############################################################
## Mean km walked per day
km_mean <- svymean(~nbkm_tot_walking, jour_walkers, na.rm = TRUE)         # Mean km per day
km_mean

km_mean_IC <- confint(km_mean)
km_mean_IC


# In exposure time (min)
km_mean * 60 / walk_speed



# Mean walking distances per day, by age group
svyby(~nbkm_tot_walking, by = ~age_grp10, jour_walkers, svymean, na.rm = T)



#################################################################
## Mean walking distances per day, per age group and per sex ----
mean_distance_people <- svyby(~nbkm_tot_walking, by = ~sex + age_grp10, jour_walkers, svymean, na.rm = T)


# Plot : Mean walking distance by age group and sex
zq <- qnorm(1-0.05/2)      # Level of confidence at 95%

mean_distance_people <- jour_walkers %>% 
  group_by(sex , age_grp10) %>% 
  summarise(mean_distance = survey_mean(nbkm_tot_walking, na.rm = TRUE)) %>% 
  rename(sex = sex) 



mean_km_walkers = ggplot(mean_distance_people, aes(x = age_grp10, y = mean_distance,
                                              ymin = mean_distance - zq*mean_distance_se, ymax = mean_distance + zq*mean_distance_se,
                                              fill = sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  labs(x = "Age group",
         y = "Mean distance walked (km per day)")
plot(mean_km_walkers)

  # Export plot
ggsave(here("output", "Plots", "Description", "Walk", "plot_mean_km_walkers.png"), plot = mean_km_walkers)



  # Test Anova: Age difference
anova_age <- svyglm(nbkm_tot_walking ~ age_grp10, jour_walkers)
summary(anova_age)

regTermTest(anova_age, ~ age_grp10)
      # p_value = 3.5064e-06                Highly significant (p<0.001)


  # T-test: Sex difference 
svyttest(nbkm_tot_walking ~ sex, jour_walkers)
      # p-value = 0.00208                   Statistically significant (<0.05) 




  # T-test: Sex difference for each age category
# test_t_sex_age FUNCTION: Perform a T test for a given age category 
test_t_sex_age <- function(age_cat, design) {
  sub_design <- subset(design, age_grp10 == age_cat)
  test <- svyttest(nbkm_tot_walking ~ sex, design = sub_design)          # T-test between sex 
  
  data.frame(
    age_grp10 = age_cat,
    statistic = test$statistic,
    p_value = test$p.value
  )
}

age_cat <- unique(jour_walkers$variables$age_grp10)
# T-test on sex for each age category
test_sex_per_age <- do.call(bind_rows, lapply(age_cat, test_t_sex_age, design = jour_walkers))





#####################################################
## Walking levels per area ----
mean_distance_area <- jour_walk_trips %>% 
  group_by(area_type) %>% 
  summarise(mean_distance = survey_mean(nbkm_tot_walking, na.rm = TRUE)) %>% 
  mutate(area_type = factor(area_type, levels = c("rural", "periurban", "urban"))) %>% 
  mutate(
    ic_lower = mean_distance - zq * mean_distance_se,
    ic_upper = mean_distance + zq * mean_distance_se
  )


mean_km_area = ggplot(mean_distance_area, aes(x = area_type, y = mean_distance,
                                                   ymin = mean_distance - zq*mean_distance_se, ymax = mean_distance + zq*mean_distance_se,
                                                   fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(name = "Municipality density degree",
                    labels = c("rural" = "sparsely and very sparsely populated",
                              "periurban" = "intermediate density",
                              "urban"= "densely populated"),
                    values = c(
                      "rural" = "palegreen3",
                      "periurban" = "darkorange",
                      "urban" = "slateblue"
                    )) +
  scale_x_discrete(labels = c("rural" = "Rural",
                              "periurban" = "Periurban", 
                              "urban" = "Urban")) +
  labs(x = "Area type",
       y = "Mean distance walked (km per day)") +
  theme_minimal()
plot(mean_km_area)

  # Export plot 
ggsave(here("output", "Plots", "Description", "Walk", "plot_mean_km_area.png"), plot = mean_km_area)


  # Test Anova
anova_area <- svyglm(nbkm_tot_walking ~ area_type, jour_walk_trips)
summary(anova_area)

regTermTest(anova_area, ~ area_type)
# p_value = 5.0434e-09                 Highly significant (p<0.0001)




##############################################################
#                       INTERMODALITY                        #
##############################################################

# Mean distance walked per day for people doing intermodal walking trips vs people doing only walking trips
mean_walk_trips <- jour_walk_trips %>% 
  group_by(age_grp10, sex) %>%
  summarise(mean_main_mid = survey_mean(nbkm_main_walk, na.rm = TRUE),
            mean_intermodal_mid = survey_mean(nbkm_intermodal_walk, na.rm = TRUE),
            mean_total_mid = survey_mean(nbkm_tot_walking, na.rm = TRUE)) %>% 
  mutate(mean_main_low = mean_main_mid - zq * mean_main_mid_se,
         mean_main_up = mean_main_mid + zq * mean_main_mid_se,
         mean_intermodal_low = mean_intermodal_mid - zq * mean_intermodal_mid_se,
         mean_intermodal_up = mean_intermodal_mid + zq * mean_intermodal_mid_se,
         mean_total_low = mean_total_mid - zq * mean_total_mid_se, 
         mean_total_up = mean_total_mid + zq * mean_total_mid_se)


# Distribution marche seule et marche intermodale par classe d'âge et sexe 
plot_mean_walk <- ggplot() +
  geom_bar(data = mean_walk_trips,
           mapping = aes(x = age_grp10, y = mean_main_mid, fill = sex, alpha = "Main walk"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  
  geom_errorbar(data = mean_walk_trips,
                mapping = aes(x = age_grp10, ymin = mean_main_low, ymax = mean_main_up, group = sex, alpha = "Main walk"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_fill_manual(values = colors_sex) +
  
  
  geom_bar(data = mean_walk_trips, 
           mapping = aes(x = age_grp10, y = mean_total_mid, fill = sex, alpha = "Total walk adjusted"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(name   = "Scenario",
                     values = c("2019 baseline" = 1, "Total walk adjusted" = 0.4)) +
  
  geom_errorbar(data = mean_walk_trips,
                mapping = aes(x = age_grp10, ymin = mean_total_low, ymax = mean_total_up, group = sex, alpha = "Total walk adjusted"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  ylab("Distance walked (km)") +
  xlab("Age group") +
  theme_minimal() 


plot_mean_walk



# Mean intermodal walking distance associated to different modes of transport







##############################################################
#                  RATE OF DEADLY ACCIDENTS                  #
##############################################################

# Rate of deadly accidents per km from 2019 levels
deaths_per_km_walked <- 483 / km_total_2019                           # Number of dead walkers per km in 2019 (ONISR 2020 - Bilan 2019)
deaths_per_km_walked








################################################################################################################################
#                                                       WALKING STEPS                                                         #
################################################################################################################################


##############################################################
#                       TOTAL STEPS                          #
##############################################################

## Total steps in 2019
step_total_2019 <- as.numeric(svytotal(~step_commute, jour_walkers)) *365.25/7                                      # Total steps per year
step_total_2019_IC <- as.numeric(confint(svytotal(~step_commute, jour_walkers) *365.25/7 ))                         # Confidence interval

step_total_2019 * 1e-9 # billion steps
step_total_2019_IC * 1e-9

# In time Per person (min/pers/year)
step_total_2019 / (pop_tot * walk_speed)


## Total steps per day in 2019
step_total_day <- svytotal(~step_commute, jour_walkers)                   # Total steps per day
step_total_day_IC <- as.numeric(confint(step_total_day))
step_total_day *1e-6
step_total_day_IC * 1e-6


# Total steps per day, by age group
svyby(~step_commute, by = ~age_grp10, jour_walkers, svytotal, na.rm = T)  



# Creation of walking distances categories (in case but for now no need)
emp_walkers_step <- emp_walkers %>%
  mutate(
    step_grp = cut(
      step_commute,
      breaks = c(0, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000),
      labels =  c("0-2000", "2000-3000", "3000-4000", "4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000-11000", "11000-12000"),
      right = FALSE
    )
  )


## Plot : Proportion of people by steps
proportion_step <- emp_walkers_step %>% 
  filter(!is.na(pond_indc)) %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE) %>% 
  group_by(step_grp = factor(step_grp, levels = c("0-2000", "2000-3000", "3000-4000", "4000-5000", "5000-6000", "6000-7000", "7000-8000", "8000-9000", "9000-10000", "10000-11000", "11000-12000"))) %>%   # Put distance range in the right order 
  summarise(total_pondere = survey_total (1, na.rm = TRUE)) %>%         # Sum of ponderation of individual for each step group
  mutate(proportion = total_pondere / sum(total_pondere))               # Proportion

prop_walkers_step <- ggplot(proportion_step, aes(x = step_grp, y = proportion)) +
  geom_col() +
  labs(title = "Distribution of people by step walked per day",
       x = "Number of steps per day",
       y = "Proportion") +
  theme_minimal()
prop_walkers_step

# Export plot
ggsave(here("output", "Plots", "Description", "Step", "plot_prop_walkers_step.png"), plot = prop_walkers_step)




# Proportion of number of steps for each sex
prop_sex_step <-  emp_walkers_step %>% 
  filter(pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strat = sex,
                   nest = TRUE) %>% 
  group_by (sex = sex) %>% 
  summarise(tot_step = survey_total(step_commute)) %>% 
  mutate(proportion = tot_step / sum(tot_step))







##############################################################
#                    MEAN NUMBER OF STEPS                    #
##############################################################
## Mean number of steps per day
step_mean <- svymean(~step_commute, jour_walkers, na.rm = TRUE)
step_mean

step_mean_IC <- confint(step_mean)
step_mean_IC


# In exposure time (min)
step_mean * 60 / (walk_speed * 1e3 / 0.686)



# Mean walking distances per day, by age group
svyby(~step_commute, by = ~age_grp10, jour_walkers, svymean, na.rm = T)



#################################################################
## Mean walking distances per day, per age group and per sex ----
mean_step_people <- svyby(~step_commute, by = ~sex + age_grp10, jour_walkers, svymean, na.rm = T)


# Plot : Mean walking distance by age group and sex
zq <- qnorm(1-0.05/2)      # Level of confidence at 95%

mean_step_people <- jour_walkers %>% 
  group_by(sex , age_grp10) %>% 
  summarise(mean_step = survey_mean(step_commute, na.rm = TRUE)) %>% 
  rename(sex = sex) 



mean_step_walkers = ggplot(mean_step_people, aes(x = age_grp10, y = mean_step,
                                                   ymin = mean_step - zq*mean_step_se, ymax = mean_step + zq*mean_step_se,
                                                   fill = sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  labs(x = "Age group",
       y = "Mean number of commute steps (steps per day)")
plot(mean_step_walkers)

# Export plot
ggsave(here("output", "Plots", "Description", "Step",  "plot_mean_step_walkers.png"), plot = mean_step_walkers)



# Test Anova: Age difference
anova_age_step <- svyglm(step_commute ~ age_grp10, jour_walkers)
summary(anova_age_step)

regTermTest(anova_age_step, ~ age_grp10)
# p_value = 5.5593e-10                 Highly significant (p<0.001)


# T-test: Sex difference 
svyttest(step_commute ~ sex, jour_walkers)
# p-value = 0.08037                  Significant (<0.05) 






#####################################################
## Walking levels per area ---- TO MODIFY !!
mean_step_area <- jour_walkers %>% 
  group_by(area_type) %>% 
  summarise(mean_step = survey_mean(step_commute, na.rm = TRUE)) %>% 
  mutate(area_type = factor(area_type, levels = c("rural", "semi_urban", "urban"))) %>% 
  mutate(
    ic_lower = mean_step - zq * mean_step_se,
    ic_upper = mean_step + zq * mean_step_se
  )


mean_step_area = ggplot(mean_step_area, aes(x = area_type, y = mean_step,
                                              ymin = mean_step - zq*mean_step_se, ymax = mean_step + zq*mean_step_se,
                                              fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(name = "Urban unit slices",
                    labels = c("rural" = "<5,000 inhabitants",
                               "semi_urban" = "<50,000 inhabitants",
                               "urban"= "≥ 50,000 inhabitants"),
                    values = c(
                      "rural" = "palegreen3",
                      "semi_urban" = "darkorange",
                      "urban" = "slateblue"
                    )) +
  scale_x_discrete(labels = c("rural" = "Rural",
                              "semi_urban" = "Semi-urban", 
                              "urban" = "Urban")) +
  labs(x = "Area type",
       y = "Mean number of commute steps (steps per day)") +
  theme_minimal()
plot(mean_step_area)

# Export plot 
ggsave(here("output", "Plots", "Description", "plot_mean_step_area.png"), plot = mean_step_area)


# Test Anova
anova_area_step <- svyglm(step_commute ~ area_type, jour_walk_trip)
summary(anova_area_step)

regTermTest(anova_area_step, ~ area_type)
# p_value = 3.6767e-12                Highly significant (p<0.0001)






################################################################################################################################
#                                                          DRIVING                                                             #
################################################################################################################################

##############################################################
#                     SHORT TRIPS (<2km)                     #
##############################################################

# French adult reporting any short (<2km) car trip in the past day according to sex and age
# emp_car_trip is trip-level data, so count each individual once by ident_ind.
drivers_2km <- emp_car_trips %>% 
  filter(!is.na(pond_jour), nbkm_car > 0) %>% 
  group_by(ident_ind, sex, age_grp10) %>% 
  summarise(short_trip = any(nbkm_car <= 2),
            pond_jour = first(pond_jour),
            .groups = "drop") %>% 
  filter(short_trip) %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour) %>% 
  group_by(sex, age_grp10) %>% 
  summarise(total = survey_total(short_trip, na.rm = TRUE)) %>% 
  rename(Sex = sex)

zq <- qnorm(1-0.05/2)

nb_drivers_2km <- ggplot(drivers_2km, aes(x = age_grp10, y = total,
                                            ymin = total - zq*total_se, ymax = total + zq*total_se, fill = Sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  ylab ("Number of drivers driving <2km in the past day") +
  xlab("Age group") +
  theme_minimal()
plot(nb_drivers_2km)

  # Export plot
ggsave(here("output", "Plots", "Description", "Drivers", "plot_drivers_2km.png"), plot = nb_drivers_2km)



## Proportion of the French adult population reporting any short (<2km) car trip in the past day according to sex and age
mean_drivers_2km <- emp_car_trips %>% 
  filter(!is.na(pond_jour), nbkm_car > 0) %>% 
  group_by(ident_ind, sex, age_grp10) %>% 
  summarise(short_trip = any(nbkm_car <= 2),
            pond_jour = first(pond_jour),
            .groups = "drop") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour) %>% 
  group_by(sex, age_grp10) %>% 
  summarise(perc = 100 * survey_mean(short_trip, na.rm = TRUE)) %>% 
  rename(Sex = sex)

zq <- qnorm(1-0.05/2)

perc_drivers_2km <- ggplot(mean_drivers_2km, aes(x = age_grp10, y = perc,
                                            ymin = perc - zq*perc_se, ymax = perc + zq*perc_se, fill = Sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  ylab ("Proportion of drivers driving <2km in the past day (%)") +
  xlab("Age group") +
  theme_minimal()
plot(perc_drivers_2km)

  # Export plot
ggsave(here("output", "Plots", "Description", "Drivers", "plot_prop_drivers_2km.png"), plot = perc_drivers_2km)



# Age distribution of people reporting any short (<2km) car trip in the past day
drivers_pyramid <- drivers_2km %>%
  mutate(total = ifelse(Sex == "Male", -total, total))

pyramid_drivers_2km <- ggplot(drivers_pyramid, aes(x = age_grp10, y = total, fill = Sex)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = function(x) abs(x)) +
  scale_fill_manual(values = c("Female" = "darkorange1", "Male" = "chartreuse4")) +
  labs(title = "Age pyramid of drivers that reported a short car trip (<2km) in the past day",
       y = "Number of drivers", x = "Age group") +
  theme_minimal()

plot(pyramid_drivers_2km)
  # Export plot
ggsave(here("output", "Plots", "Description", "Drivers", "pyramid_drivers_2km.png"), plot = pyramid_drivers_2km, width = 8, height = 6)




# Survey design for statistical test

# Test Anova: Age difference
anova_age_drivers <- svyglm(perc ~ age_grp10, )
summary(anova_age_drivers)

regTermTest(anova_age_drivers, ~ age_grp10)
# p_value =                 Highly significant (p<0.001)


# T-test: Sex difference 
svyttest(perc ~ Sex, )
# p-value =                    Not statistically significant (>0.05) 




##############################################################
#                MEAN DRIVEN DISTANCE (<2km)                 #
##############################################################

# Mean distance driven (km) in the past day among those reporting short car trips <2km 
zq <- qnorm(1-0.05/2)

mean_short_trips <- emp_car_trips %>% 
  filter(pond_jour != "NA", nbkm_car > 0, nbkm_car <= 2) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_jour) %>% 
  summarise(day_mean = survey_mean(nbkm_car, na.rm = TRUE)) %>% 
  mutate(ic_lower = day_mean - zq * day_mean_se,
         ic_upper = day_mean + zq * day_mean_se)


# Mean distance driven (km) in the past day among those reporting short car trips <2km according to sex and age
mean_drivers_2km <- emp_car_trips %>% 
  filter(pond_jour != "NA", nbkm_car > 0, nbkm_car <= 2) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_jour) %>% 
  group_by(sex, age_grp10) %>% 
  summarise(day_mean = survey_mean(nbkm_car, na.rm = TRUE)) %>% 
  rename(Sex = sex)


mean_km_drivers_2km <- ggplot(mean_drivers_2km, aes(x = age_grp10, y = day_mean,
                                            ymin = day_mean - zq*day_mean_se, ymax = day_mean + zq*day_mean_se, fill = Sex)) +
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  ylab ("Mean length of short car travel <2km (km)") +
  xlab("Age group") +
  theme_minimal() 
plot(mean_km_drivers_2km)

  # Export plot
ggsave(here("output", "Plots", "Description", "Drivers", "plot_mean_drivers_2km.png"), plot = mean_km_drivers_2km)



