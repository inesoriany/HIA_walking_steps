##############################################
################ DESCRIPTION #################
##############################################




################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################

pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data manipulation
  tidyr,        # Data manipulation
  purrr,        # Loop
  epikit,       # Age categories creation
  survey,       # Survey management
  srvyr,        # Survey management
  ggplot2,      # Data visualization
  scales,       # Format visualization
  patchwork     # Graphs combination
)



################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# EMP 2019 trips
emp_trip <- import(here("data_clean", "EMP_walking_trips.xlsx")) 

# EMP 2019 subset for walkers
emp_walkers <- import(here("data_clean", "EMP_walkers.xlsx"))

# EMP 2019 subset for car trips
emp_car_trips <- import(here("data_clean", "EMP_car_trips.xlsx"))

# Diseases incidence
diseases_incidence <- import(here("data_clean", "Diseases", "incidence_table.xlsx"))


################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("R_code", "Parameters.R"))

# Diseases considered
dis_vec = c("mort", "bc", "cc", "cvd", "cancer", "diab2", "dem", "dep")
morbi_vec = c("bc", "cc", "cvd", "cancer", "diab2", "dem", "dep")


################################################################################################################################
#                                                    4. SURVEY DESIGNS                                                         #
################################################################################################################################
# Survey design ponderated by day
jour_walkers <- emp_walkers %>% 
  filter(!is.na(pond_jour)) %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)

jour_trips <- emp_trip %>% 
  filter(!is.na(pond_jour)) %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10, area_type),
                   nest = TRUE)

jour_car <- emp_car_trips %>% 
  filter(!is.na(pond_jour)) %>%
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



# Survey design ponderated by individual
indiv_walkers <- emp_walkers %>% 
  filter(!is.na(pond_indc)) %>%
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE)


indiv_trips <- emp_trip %>% 
  filter(!is.na(pond_indc)) %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10, area_type),
                   nest = TRUE)

indiv_car <- emp_car_trips %>% 
  filter(!is.na(pond_indc)) %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



################################################################################################################################
#                                                        5. HEALTH                                                             #
################################################################################################################################
# Death incidence for 100,000
mort_incidence <- indiv_walkers  %>% 
  summarise (rate_mean = survey_mean(mort_rate, vartype = "ci")) %>% 
  mutate(rate_mean*100000, rate_mean_low*100000, rate_mean_upp*100000)


# Incidence distribution per age and sex
list_incidence <- lapply(morbi_vec, function(dis) {
  
  ggplot(
    diseases_incidence %>% 
      filter(disease == dis, measure == "incidence"),
    aes(x = age_grp10, y = mid, color = sex, group = sex)
  ) +
    geom_line(show.legend = FALSE) +
    geom_errorbar(
      aes(ymin = low, ymax = up),
      width = 0.2,
      show.legend = FALSE
    ) +
    geom_point(size = 2) +
    scale_color_manual(values = colors_sex) +
    scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
    labs(
      title = names_disease[[dis]],
      y = "Incidence",
      x = "Age group",
      color = "Sex"
    ) +
    theme_minimal() +
    theme(legend.position = "top") +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
})

combined_plot_incidence <- wrap_plots(list_incidence, ncol = 3)

print(combined_plot_incidence)



################################################################################################################################
#                                                    6. TRANSPORT MODE                                                         #
################################################################################################################################

##############################################################
#                          DISTANCE                          #
##############################################################
# Total walking including intermodal trips
walking <- indiv_trips  %>% 
    summarise(tot_km = survey_total(nbkm_tot_walking_jour, na.rm = TRUE))  %>% 
    mutate(mode = "tot_walk")


# Transport mode
transport_mode <- indiv_trips  %>% 
    group_by(mode)  %>% 
    summarise(tot_km = survey_total(mdisttot_fin*pond_jour/(pond_indc*7), na.rm = TRUE))  %>% 
    filter(mode != "walk")

# Modal share in distance
modal_share <- bind_rows(walking, transport_mode) %>%
  mutate(share_perc = tot_km * 100 / sum(tot_km))



##############################################################
#                           TIME                             #
##############################################################
# Total walking including intermodal trips and walking in private spaces
walking_time <- indiv_trips  %>% 
    summarise(tot_min = survey_total((nbkm_tot_walking_jour + baseline_step*step_length)*60/walk_speed, na.rm = TRUE))  %>% 
    mutate(mode = "tot_walk")


# Transport mode
transport_mode_time <- indiv_trips  %>% 
    filter(mode != "walk")  %>% 
    # Associate transport speeds to corresponding mode and area
    mutate (speed = case_when(
      area_type == "urban" & mode == "bike" ~ bike_urban_speed,
      area_type == "periurban" & mode == "bike" ~ bike_periurban_speed,
      area_type == "rural" & mode == "bike" ~ bike_rural_speed,
      area_type == "urban" & mode == "public_transport" ~ public_transport_urban_speed,
      area_type == "periurban" & mode == "public_transport" ~ public_transport_periuban_speed,
      area_type == "rural" & mode == "public_transport" ~ public_transport_rural_speed,
      area_type == "urban" & mode == "car" ~ car_urban_speed,
      area_type == "periurban" & mode == "car" ~ car_periurban_speed,
      area_type == "rural" & mode == "car" ~ car_rural_speed,
      TRUE ~ NA_real_
    ))  %>% 
    group_by(mode)  %>% 
    summarise(tot_min = survey_total((mdisttot_fin*pond_jour/(pond_indc*7))*60/speed, na.rm = TRUE))  


# Modal share in distance
modal_share_time <- bind_rows(walking_time, transport_mode_time) %>%
  mutate(share_perc = tot_min * 100 / sum(tot_min))




################################################################################################################################
#                                                    7. WALKING - distance                                                     #
################################################################################################################################

# Total population
pop_tot <- sum(emp_walkers$pond_indc)
pop_tot



##############################################################
#                  TOTAL WALKED DISTANCE                     #
##############################################################

## Total walked distance in 2019
km_total_2019 <- as.numeric(svytotal(~nbkm_tot_walking, jour_walkers)) *365.25/7                              # Total km per year
km_total_2019_IC <- as.numeric(confint(svytotal(~nbkm_tot_walking, jour_walkers) *365.25/7 ))                 # Confidence interval

km_total_2019 * 1e-9 # billion km
km_total_2019_IC * 1e-9



## Total walked distance per day in 2019
km_total_day <- svytotal(~nbkm_tot_walking, jour_walkers)/7                   # Total km per day
km_total_day_IC <- as.numeric(confint(km_total_day))/7
km_total_day *1e-6
km_total_day_IC * 1e-6


intermodal_km_total_day <- (svytotal(~nbkm_intermodal_walk, jour_walkers))                          # Total km per day
intermodal_km_total_day_IC <- as.numeric(confint(intermodal_km_total_day))
intermodal_km_total_day
intermodal_km_total_day_IC 
intermodal_km_total_day / step_length * 1e-9
intermodal_km_total_day_IC /  step_length * 1e-9
intermodal_km_total_day / km_total_day                    # Share of intermodal walk


# Total walking distances per day, by age group
svyby(~nbkm_tot_walking, by = ~age_grp10, jour_walkers, svytotal, na.rm = T)  


# Proportion of distances walked by each sex
prop_sex <-  indiv_walkers  %>% 
  group_by (sex) %>% 
  summarise(tot_km = survey_total(nbkm_tot_walking_jour, na.rm = TRUE)) %>% 
  mutate(proportion = tot_km / sum(tot_km))



#############################################################
#                    MEAN WALKED DISTANCE                    #
##############################################################
# Mean walked distance per day for an average French person
mean_distance_jour <- indiv_walkers %>%
    summarise(mean_km = survey_mean(nbkm_tot_walking_jour, na.rm = TRUE, vartype = "ci"))                 # 1.1608 (1.1109-1.2107) km per day


## EMP METHODOLOGY
day <- emp_walkers  %>% 
    mutate(km_pond = nbkm_tot_walking_jour * pond_jour/7, 
          km_main_pond = nbkm_main_walk * pond_jour/7,
          km_inter_pond = nbkm_intermodal_walk * pond_jour/7)

  # Exclusive + intermodal 
  km_mean <- sum(day$km_pond, na.rm = TRUE) / sum(day$pond_indc, na.rm = TRUE)  # Total walk (exclusive + intermodal)
  km_mean                             # 1.132096 km per day
  km_mean / step_length               # 1583.351 steps per day  
  km_mean*60 / walk_speed             # 14.1512 minutes per day


  # Exclusive walking
  main_km_mean <- sum(day$km_main_pond, na.rm = TRUE) / sum(day$pond_indc, na.rm = TRUE)     # Mean km per day
  main_km_mean                             # 0.7246853 km per day
  main_km_mean / step_length               # 1013.546 steps per day
  main_km_mean*60 / walk_speed             # 9.058567 minutes per day


  # Intermodal walk
  intermodal_km_mean <- sum(day$km_inter_pond, na.rm = TRUE) / sum(day$pond_indc, na.rm = TRUE)         # Mean km per day 
  intermodal_km_mean                       # 0.4074109 km per day
  intermodal_km_mean / step_length         # 569.8055 steps per day  



# Plot : Mean walking distance by age group and sex

    # Total walk including intermodal
    mean_distance_people <- indiv_walkers %>% 
    group_by(sex , age_grp10) %>% 
    summarise(mean_km = survey_mean(nbkm_tot_walking_jour, na.rm = TRUE, vartype = "ci"))

    # Main walk
    main_mean_distance <- indiv_walkers  %>% 
        group_by(sex, age_grp10)  %>% 
        summarise(mean_km = survey_mean(nbkm_main_walk_jour, na.rm = TRUE, vartype = "ci"))


plot_mean_km_walkers <- 
  ggplot() +
  geom_bar(data = main_mean_distance,
    mapping = aes(x = age_grp10, y = mean_km, fill = sex, alpha = "Exclusively walking"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(data = main_mean_distance,
    mapping = aes(x = age_grp10, ymin = mean_km_low, ymax = mean_km_upp, group = sex, alpha = "Exclusively walking"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_fill_manual(values = colors_sex) +

  geom_bar(data = mean_distance_people, 
    mapping = aes(x = age_grp10, y = mean_km, fill = sex, alpha = "Total walking including intermodal walk"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(
    data = mean_distance_people,
    mapping = aes(x = age_grp10, ymin = mean_km_low, ymax = mean_km_upp, group = sex, alpha = "Total walking including intermodal walk"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_alpha_manual(
    values = c("Exclusively walking" = 1, "Total walking including intermodal walk" = 0.4),
    guide = "none") +
    
  scale_x_discrete(labels = names_disease) + 
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +  
  ylab("Mean distance walked (km per day)") +
  xlab("Age group") +
  theme_minimal()

plot_mean_km_walkers



# Test Anova: Age difference
anova_age <- svyglm(nbkm_tot_walking_jour ~ age_grp10, indiv_walkers)
summary(anova_age)

regTermTest(anova_age, ~ age_grp10)
      # p_value = 0.037478                Highly significant (p<0.001)


  # T-test: Sex difference 
svyttest(nbkm_tot_walking_jour ~ sex, indiv_walkers)
      # p-value = 0.1822                   Statistically significant (<0.05) 




  # T-test: Sex difference for each age category
# test_t_sex_age FUNCTION: Perform a T test for a given age category 
test_t_sex_age <- function(age_cat, design) {
  sub_design <- subset(design, age_grp10 == age_cat)
  test <- svyttest(nbkm_tot_walking_jour ~ sex, design = sub_design)          # T-test between sex 
  
  data.frame(
    age_grp10 = age_cat,
    statistic = test$statistic,
    p_value = test$p.value
  )
}

age_cat <- unique(indiv_walkers$variables$age_grp10)
# T-test on sex for each age category
test_sex_per_age <- do.call(bind_rows, lapply(age_cat, test_t_sex_age, design = indiv_walkers))




#############################################################
#                        MEAN PER AREA                      #
#############################################################
mean_distance_area <- indiv_walkers %>% 
  group_by(area_type) %>% 
  summarise(mean_km = survey_mean(nbkm_tot_walking_jour, na.rm = TRUE, vartype = "ci")) %>% 
  mutate(area_type = factor(area_type, levels = c("urban", "periurban", "rural"))) 


plot_mean_km_area <- ggplot(mean_distance_area, aes(x = area_type, y = mean_km,
                                                   ymin = mean_km_low, ymax = mean_km_upp,
                                                   fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(name = "Municipality density degree",
                    values = colors_area) +
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  labs(x = "Area type",
       y = "Mean distance walked (km per day)") +
  theme_minimal()
plot(plot_mean_km_area)



  # Test Anova
anova_area <- svyglm(nbkm_tot_walking_jour ~ area_type, indiv_walkers)
summary(anova_area)

regTermTest(anova_area, ~ area_type)
# p_value = < 2.22e-16                 Highly significant (p<0.0001)


##############################################################
#                  RATE OF DEADLY ACCIDENTS                  #
##############################################################

# Rate of deadly accidents per km from 2019 levels
deaths_per_km_walked <- 483 / km_total_2019                           # Number of dead walkers per km in 2019 (ONISR 2020 - Bilan 2019)
deaths_per_km_walked



################################################################################################################################
#                                                     8. WALKING - steps                                                       #
################################################################################################################################

##############################################################
#                       TOTAL STEPS                          #
##############################################################

## Total steps in 2019
step_total_2019 <- as.numeric(svytotal(~step_commute, jour_walkers)) *365.25/7                            # Total steps per year
step_total_2019_IC <- as.numeric(confint(svytotal(~step_commute, jour_walkers) *365.25/7 ))               # Confidence interval

step_total_2019 * 1e-9 # billion steps
step_total_2019_IC * 1e-9


## Total steps per day 
step_total_day <- as.numeric(svytotal(~step_commute, jour_walkers))                                     # Total steps per day
step_total_day_IC <- as.numeric(confint(svytotal(~step_commute, jour_walkers)))                         # Confidence interval

step_total_day * 1e-9 # billion steps
step_total_day_IC * 1e-9


#############################################################
#                          MEAN STEPS                       #
#############################################################

## Mean number of steps per day
mean_step_jour <- indiv_walkers %>%
    summarise(mean_step = survey_mean(step_commute_jour, na.rm = TRUE, vartype = "ci"))              # 1623.4999 (1553.6724 - 1693.3275) steps per day



# Plot : Mean steps by age group and sex (+ baseline steps)

    # Total walk including intermodal
    mean_step_people <- indiv_walkers %>% 
    group_by(sex , age_grp10) %>% 
    summarise(mean_step = survey_mean(step_commute_jour + baseline_step, na.rm = TRUE, vartype = "ci"))

    # Main walk
    main_mean_step <- indiv_walkers  %>% 
        group_by(sex, age_grp10)  %>% 
        summarise(mean_step = survey_mean(nbkm_main_walk_jour/step_length + baseline_step, na.rm = TRUE, vartype = "ci"))


plot_mean_steps_walkers <- 
  ggplot() +
  geom_bar(data = main_mean_step,
    mapping = aes(x = age_grp10, y = mean_step, fill = sex, alpha = "Exclusively walking"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(data = main_mean_step,
    mapping = aes(x = age_grp10, ymin = mean_step_low, ymax = mean_step_upp, group = sex, alpha = "Exclusively walking"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_fill_manual(values = colors_sex) +

  geom_bar(data = mean_step_people, 
    mapping = aes(x = age_grp10, y = mean_step, fill = sex, alpha = "Total walking including intermodal walk"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(
    data = mean_step_people,
    mapping = aes(x = age_grp10, ymin = mean_step_low, ymax = mean_step_upp, group = sex, alpha = "Total walking including intermodal walk"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_alpha_manual(
    values = c("Exclusively walking" = 1, "Total walking including intermodal walk" = 0.4),
    guide = "none") +
  scale_x_discrete(labels = names_disease) + 
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  ylab("Mean steps walked (steps per day)") +
  xlab("Age group") +
  theme_minimal()

plot_mean_steps_walkers



#############################################################
#                        MEAN PER AREA                      #
#############################################################
mean_step_area <- indiv_walkers %>% 
  group_by(area_type) %>% 
  summarise(mean_step = survey_mean(step_commute + baseline_step, na.rm = TRUE, vartype = "ci")) %>% 
  mutate(area_type = factor(area_type, levels = c("urban", "periurban", "rural"))) 


plot_mean_step_area = ggplot(mean_step_area, aes(x = area_type, y = mean_step,
                                              ymin = mean_step_low, ymax = mean_step_upp,
                                              fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(name = "Municipality density degree",
                    labels = c("urban", "periurban", "rural"),
                    values = colors_area) +
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  labs(x = "Area type",
       y = "Mean steps per day per area") +
  theme_minimal()
plot_mean_step_area


# Test Anova
anova_area_step <- svyglm(step_commute_jour ~ area_type, indiv_walkers)
summary(anova_area_step)

regTermTest(anova_area_step, ~ area_type)
# p_value < 2.22e-16                 Highly significant (p<0.0001)



################################################################################################################################
#                                                         9.DRIVING                                                            #
################################################################################################################################
emp_short_drivers <- emp_car_trips %>% 
  filter(!is.na(pond_jour), nbkm_car > 0) %>% 
  group_by(ident_ind, sex, age_grp10) %>%           # emp_car_trip is trip-level data, so count each individual once by ident_ind
  summarise(short_trip = any(nbkm_car <= 2),
            pond_indc = first(pond_indc),
            .groups = "drop") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc)


# Proportion of French adult reporting any short (<2km) car trip in the past day
prop_short_drivers <- emp_short_drivers  %>% 
  summarise(perc = 100 * survey_mean(short_trip, na.rm = TRUE,  vartype = "ci")) 


##############################################################
#                     SHORT TRIPS (<2km)                     #
##############################################################
# French adult reporting any short (<2km) car trip in the past day according to sex and age
drivers_2km <- emp_short_drivers %>% 
  group_by(sex, age_grp10) %>% 
  summarise(nb_trip = survey_total(short_trip, na.rm = TRUE, vartype = "ci"))


plot_nb_drivers_2km <- ggplot(drivers_2km, aes(x = age_grp10, y = nb_trip,
                                            ymin = nb_trip_low, ymax = nb_trip_upp, fill = sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = colors_sex) +
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  ylab ("Number of drivers driving <2km in the past day") +
  xlab("Age group") +
  theme_minimal()
plot_nb_drivers_2km



# Proportion of the French adult population reporting any short (<2km) car trip in the past day according to sex and age
prop_drivers_2km <- emp_short_drivers %>% 
  group_by(sex, age_grp10) %>% 
  summarise(perc = 100 * survey_mean(short_trip, na.rm = TRUE, vartype = "ci"))


plot_perc_drivers_2km <- ggplot(prop_drivers_2km, aes(x = age_grp10, y = perc,
                                            ymin = perc_low, ymax = perc_upp, fill = sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = colors_sex) +
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  ylab ("Proportion of drivers driving <2km in the past day (%)") +
  xlab("Age group") +
  theme_minimal()
plot_perc_drivers_2km




# Age distribution of people reporting any short (<2km) car trip in the past day
drivers_pyramid <- drivers_2km %>%
  mutate(nb_trip = ifelse(sex == "Male", -nb_trip, nb_trip))

plot_pyramid_drivers_2km <- ggplot(drivers_pyramid, aes(x = age_grp10, y = nb_trip, fill = sex)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = function(x) label_comma(big.mark = ",", decimal.mark = ".")(abs(x))) +
  scale_fill_manual(values = colors_sex) +
  labs(title = "Age pyramid of drivers that reported a short car trip (<2km) in the past day",
       y = "Number of drivers", x = "Age group") +
  theme_minimal()

plot_pyramid_drivers_2km





# Test Anova: Age difference
anova_age_drivers <- svyglm(nbkm_car_jour ~ age_grp10, indiv_car)
summary(anova_age_drivers)

regTermTest(anova_age_drivers, ~ age_grp10)
# p_value = < 2.22e-16                 Highly significant (p<0.001)


# T-test: Sex difference 
svyttest(nbkm_car_jour ~ sex, indiv_car )
# p-value = 8.731e-08                   Statistically significant (<0.05) 
# Mean difference observed :  3.42 [2.17 - 4.68] km




##############################################################
#                MEAN DRIVEN DISTANCE (<2km)                 #
##############################################################

# Mean distance driven (km) in the past day among those reporting short car trips <2km 
mean_short_trips <- emp_car_trips %>% 
  filter(!is.na(pond_indc), nbkm_car_jour > 0, nbkm_car_jour <= 2) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_indc) %>% 
  summarise(day_mean = survey_mean(nbkm_car_jour, na.rm = TRUE, vartype = "ci"))


# Mean distance driven (km) in the past day among those reporting short car trips <2km according to sex and age
mean_drivers_2km <- emp_car_trips %>% 
  filter(!is.na(pond_indc), nbkm_car_jour > 0, nbkm_car_jour <= 2) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_indc) %>% 
  group_by(sex, age_grp10) %>% 
  summarise(day_mean = survey_mean(nbkm_car_jour, na.rm = TRUE, vartype = "ci"))


plot_mean_km_drivers_2km <- ggplot(mean_drivers_2km, aes(x = age_grp10, y = day_mean,
                                            ymin = day_mean_low, ymax = day_mean_upp, fill = sex)) +
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  ylab ("Mean distance of short car travel <2km (km)") +
  xlab("Age group") +
  theme_minimal() 
plot_mean_km_drivers_2km





################################################################################################################################
#                                                     10. EXPORT DATA                                                          #
################################################################################################################################

# HEALTH
    # Incidence
    ggsave(here("output", "Plots", "Description", "Diseases", "morbi_incidence.png"), plot = combined_plot_incidence)

# WALKING
    # Sex proportion
    export(prop_sex, here("output", "Tables", "Description", "Walk", "sex_proportion.xlsx"))
    # Mean walk
    ggsave(here("output", "Plots", "Description", "Walk", "plot_mean_km_walkers.png"), plot = plot_mean_km_walkers)
    ggsave(here("output", "Plots", "Description", "Steps", "plot_mean_steps_walkers.png"), plot = plot_mean_steps_walkers)
    ggsave(here("output", "Plots", "Description", "Walk", "plot_mean_km_area.png"), plot = plot_mean_km_area)
    ggsave(here("output", "Plots", "Description", "Steps", "plot_mean_step_area.png"), plot = plot_mean_step_area)



# DRIVING
    # Short drivers proportion
    export(prop_short_drivers, here("output", "Tables", "Description", "Drivers", "short_drivers_proportion.xlsx"))
    # Drivers profile
    ggsave(here("output", "Plots", "Description", "Drivers", "plot_drivers_2km.png"), plot = plot_nb_drivers_2km)
    ggsave(here("output", "Plots", "Description", "Drivers", "plot_prop_drivers_2km.png"), plot = plot_perc_drivers_2km)
    ggsave(here("output", "Plots", "Description", "Drivers", "pyramid_drivers_2km.png"), plot = plot_pyramid_drivers_2km, width = 8, height = 6)
    # Mean
    export(mean_short_trips, here("output", "Tables", "Description", "Drivers", "mean_km_short_drivers.xlsx"))
    ggsave(here("output", "Plots", "Description", "Drivers", "plot_mean_km_drivers_2km.png"), plot = plot_mean_km_drivers_2km)
