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
  epikit,       # Age categories creation
  survey,       # Survey management
  srvyr,        # Survey management
  ggplot2,       # Data visualization
  patchwork           # graphs combination
)



################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# EMP 2019
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
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)

jour_trips <- emp_trip %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10, area_type),
                   nest = TRUE)

jour_car <- emp_car_trips %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



# Survey design ponderated by individual
indiv_walkers <- emp_walkers %>% 
  filter (pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE)


indiv_trips <- emp_trip %>% 
  filter (pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10, area_type),
                   nest = TRUE)

indiv_car <- emp_car_trips %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
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

combined_plot_incidence <- 
  reduce(list_incidence, `+`) + 
  plot_layout(ncol = 3)

print(combined_plot_incidence)



################################################################################################################################
#                                                    6. TRANSPORT MODE                                                         #
################################################################################################################################

# Total walking including intermodal trips
walking <- indiv_trips  %>% 
    summarise(tot_km = survey_total(nbkm_tot_walking, na.rm = TRUE))  %>% 
    mutate(mode = "tot_walk")


# Transport mode
transport_mode <- indiv_trips  %>% 
    group_by(mode)  %>% 
    summarise(tot_km = survey_total(mdisttot_fin, na.rm = TRUE))  %>% 
    filter(mode != "walk")

# Modal share
modal_share <- bind_rows(walking, transport_mode) %>%
  mutate(share_perc = tot_km * 100 / sum(tot_km))



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
km_total_day <- svytotal(~nbkm_tot_walking, jour_walkers)                   # Total km per day
km_total_day_IC <- as.numeric(confint(km_total_day))
km_total_day *1e-6
km_total_day_IC * 1e-6

intermodal_km_total_day <- (svytotal(~nbkm_intermodal_walk, jour_walkers))                          # Total km per year
intermodal_km_total_day_IC <- as.numeric(confint(intermodal_km_total_day))
intermodal_km_total_day
intermodal_km_total_day_IC 
intermodal_km_total_day / step_length * 1e-9
intermodal_km_total_day_IC /  step_length * 1e-9
intermodal_km_total_day / km_total_day                    # Share of intermodal walk


# Total walking distances per day, by age group
svyby(~nbkm_tot_walking, by = ~age_grp10, jour_walkers, svytotal, na.rm = T)  


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



#############################################################
#                    MEAN WALKED DISTANCE                    #
##############################################################
## Mean km walked per day
km_mean <- svymean(~nbkm_tot_walking, jour_walkers, na.rm = TRUE)         # Mean km per day
km_mean                             # 1.3468 (1.290339-1.403312) km per day
km_mean / step_length               # 1883.7 (1804.67-1962.675) steps per day  


km_mean_IC <- confint(km_mean)
km_mean_IC
km_mean_IC/step_length


# Intermodal walk
intermodal_km_mean <- svymean(~nbkm_intermodal_walk, jour_walkers, na.rm = TRUE)         # Mean km per day 
intermodal_km_mean                       # 0.48774 (0.4564521-0.519022) km per day
intermodal_km_mean / step_length         # 682.15 (638.3945-725.9048) steps per day  

intermodal_km_mean_IC <- confint(intermodal_km_mean)
intermodal_km_mean_IC
intermodal_km_mean_IC/step_length


# Plot : Mean walking distance by age group and sex
zq <- qnorm(1-0.05/2)      # Level of confidence at 95%

    # Total walk including intermodal
    mean_distance_people <- jour_walkers %>% 
    group_by(sex , age_grp10) %>% 
    summarise(mean_distance = survey_mean(nbkm_tot_walking, na.rm = TRUE))

    # Main walk
    main_mean_distance <- jour_walkers  %>% 
        group_by(sex, age_grp10)  %>% 
        summarise(mean_distance = survey_mean(nbkm_main_walk, na.rm = TRUE))


plot_mean_km_walkers <- 
  ggplot() +
  geom_bar(data = main_mean_distance,
    mapping = aes(x = age_grp10, y = mean_distance, fill = sex, alpha = "Exclusively walking"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(data = main_mean_distance,
    mapping = aes(x = age_grp10, ymin = mean_distance - zq*mean_distance_se, ymax = mean_distance + zq*mean_distance_se, group = sex, alpha = "Exclusively walking"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_fill_manual(values = colors_sex) +

  geom_bar(data = mean_distance_people, 
    mapping = aes(x = age_grp10, y = mean_distance, fill = sex, alpha = "Total walking including intermodal walk"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(
    data = mean_distance_people,
    mapping = aes(x = age_grp10, ymin = mean_distance - zq*mean_distance_se, ymax = mean_distance + zq*mean_distance_se, group = sex, alpha = "Total walking including intermodal walk"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_alpha_manual(
    name = "Walking type",
    values = c("Exclusively walking" = 1, "Total walking including intermodal walk" = 0.4)) +
  scale_x_discrete(labels = names_disease) + 
  ylab("Mean distance walked (km per day)") +
  xlab("Age group") +
  theme_minimal()

plot_mean_km_walkers



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




#############################################################
#                        MEAN PER AREA                      #
#############################################################
mean_distance_area <- jour_walkers %>% 
  group_by(area_type) %>% 
  summarise(mean_distance = survey_mean(nbkm_tot_walking, na.rm = TRUE)) %>% 
  mutate(area_type = factor(area_type, levels = c("urban", "periurban", "rural"))) %>% 
  mutate(
    ic_lower = mean_distance - zq * mean_distance_se,
    ic_upper = mean_distance + zq * mean_distance_se
  )


plot_mean_km_area = ggplot(mean_distance_area, aes(x = area_type, y = mean_distance,
                                                   ymin = mean_distance - zq*mean_distance_se, ymax = mean_distance + zq*mean_distance_se,
                                                   fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(name = "Municipality density degree",
                    values = colors_area) +
  labs(x = "Area type",
       y = "Mean distance walked (km per day)") +
  theme_minimal()
plot(plot_mean_km_area)



  # Test Anova
anova_area <- svyglm(nbkm_tot_walking ~ area_type, jour_walkers)
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
step_total_day_IC <- as.numeric(confint(svytotal(~step_commute, jour_walkers)))              # Confidence interval

step_total_day * 1e-9 # billion steps
step_total_day_IC * 1e-9


#############################################################
#                          MEAN STEPS                       #
#############################################################
## Mean number of steps per day
step_mean <- svymean(~step_commute, jour_walkers, na.rm = TRUE)
step_mean

step_mean_IC <- confint(step_mean)
step_mean_IC                        # 1883.7 (1804.67-1962.675) steps per day  


# Plot : Mean steps by age group and sex
zq <- qnorm(1-0.05/2)      # Level of confidence at 95%

    # Total walk including intermodal
    mean_step_people <- jour_walkers %>% 
    group_by(sex , age_grp10) %>% 
    summarise(mean_step = survey_mean(step_commute, na.rm = TRUE))

    # Main walk
    main_mean_step <- jour_walkers  %>% 
        group_by(sex, age_grp10)  %>% 
        summarise(mean_step = survey_mean(nbkm_main_walk/step_length, na.rm = TRUE))


plot_mean_steps_walkers <- 
  ggplot() +
  geom_bar(data = main_mean_step,
    mapping = aes(x = age_grp10, y = mean_step, fill = sex, alpha = "Exclusively walking"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(data = main_mean_step,
    mapping = aes(x = age_grp10, ymin = mean_step - zq*mean_step_se, ymax = mean_step + zq*mean_step_se, group = sex, alpha = "Exclusively walking"),
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
    mapping = aes(x = age_grp10, ymin = mean_step - zq*mean_step_se, ymax = mean_step + zq*mean_step_se, group = sex, alpha = "Total walking including intermodal walk"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_alpha_manual(
    name = "Walking type",
    values = c("Exclusively walking" = 1, "Total walking including intermodal walk" = 0.4)) +
  scale_x_discrete(labels = names_disease) + 
  ylab("Mean steps walked (steps per day)") +
  xlab("Age group") +
  theme_minimal()

plot_mean_steps_walkers



#############################################################
#                        MEAN PER AREA                      #
#############################################################
mean_step_area <- jour_walkers %>% 
  group_by(area_type) %>% 
  summarise(mean_step = survey_mean(step_commute, na.rm = TRUE)) %>% 
  mutate(area_type = factor(area_type, levels = c("urban", "periurban", "rural"))) %>% 
  mutate(
    ic_lower = mean_step - zq * mean_step_se,
    ic_upper = mean_step + zq * mean_step_se,
    share = mean_step / sum(mean_step)
  ) 


plot_mean_step_area = ggplot(mean_step_area, aes(x = area_type, y = mean_step,
                                              ymin = mean_step - zq*mean_step_se, ymax = mean_step + zq*mean_step_se,
                                              fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4)) +
  geom_errorbar(position = position_dodge(.7), width = .25) + 
  scale_fill_manual(name = "Municipality density degree",
                    labels = c("urban", "periurban", "rural"),
                    values = colors_area) +
  labs(x = "Area type",
       y = "Mean steps per day per area") +
  theme_minimal()
plot_mean_step_area


# Test Anova
anova_area_step <- svyglm(step_commute ~ area_type, jour_walkers)
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
  summarise(total = survey_total(short_trip, na.rm = TRUE)) %>% 
  rename(Sex = sex)

zq <- qnorm(1-0.05/2)

plot_nb_drivers_2km <- ggplot(drivers_2km, aes(x = age_grp10, y = total,
                                            ymin = total - zq*total_se, ymax = total + zq*total_se, fill = Sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = colors_sex) +
  ylab ("Number of drivers driving <2km in the past day") +
  xlab("Age group") +
  theme_minimal()
plot_nb_drivers_2km



# Proportion of the French adult population reporting any short (<2km) car trip in the past day according to sex and age
prop_drivers_2km <- emp_short_drivers %>% 
  group_by(sex, age_grp10) %>% 
  summarise(perc = 100 * survey_mean(short_trip, na.rm = TRUE)) %>% 
  rename(Sex = sex)

zq <- qnorm(1-0.05/2)

plot_perc_drivers_2km <- ggplot(prop_drivers_2km, aes(x = age_grp10, y = perc,
                                            ymin = perc - zq*perc_se, ymax = perc + zq*perc_se, fill = Sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = colors_sex) +
  ylab ("Proportion of drivers driving <2km in the past day (%)") +
  xlab("Age group") +
  theme_minimal()
plot_perc_drivers_2km




# Age distribution of people reporting any short (<2km) car trip in the past day
drivers_pyramid <- drivers_2km %>%
  mutate(total = ifelse(Sex == "Male", -total, total))

plot_pyramid_drivers_2km <- ggplot(drivers_pyramid, aes(x = age_grp10, y = total, fill = Sex)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = function(x) abs(x)) +
  scale_fill_manual(values = colors_sex) +
  labs(title = "Age pyramid of drivers that reported a short car trip (<2km) in the past day",
       y = "Number of drivers", x = "Age group") +
  theme_minimal()

plot_pyramid_drivers_2km





# Test Anova: Age difference
anova_age_drivers <- svyglm(nbkm_car ~ age_grp10, jour_car)
summary(anova_age_drivers)

regTermTest(anova_age_drivers, ~ age_grp10)
# p_value = 1.4346e-09                Highly significant (p<0.001)


# T-test: Sex difference 
svyttest(nbkm_car ~ sex, jour_car )
# p-value = 5.81e-05                   Statistically significant (<0.05) 
# Mean difference observed :  2,90 [1,49 - 4,32] km




##############################################################
#                MEAN DRIVEN DISTANCE (<2km)                 #
##############################################################

# Mean distance driven (km) in the past day among those reporting short car trips <2km 
zq <- qnorm(1-0.05/2)

mean_short_trips <- emp_car_trips %>% 
  filter(pond_jour != "NA", nbkm_car > 0, nbkm_car <= 2) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_jour) %>% 
  summarise(day_mean = survey_mean(nbkm_car, na.rm = TRUE, vartype = "ci"))


# Mean distance driven (km) in the past day among those reporting short car trips <2km according to sex and age
mean_drivers_2km <- emp_car_trips %>% 
  filter(pond_jour != "NA", nbkm_car > 0, nbkm_car <= 2) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_jour) %>% 
  group_by(sex, age_grp10) %>% 
  summarise(day_mean = survey_mean(nbkm_car, na.rm = TRUE)) %>% 
  rename(Sex = sex)


plot_mean_km_drivers_2km <- ggplot(mean_drivers_2km, aes(x = age_grp10, y = day_mean,
                                            ymin = day_mean - zq*day_mean_se, ymax = day_mean + zq*day_mean_se, fill = Sex)) +
  scale_fill_manual(values = c("Female" = "darkorange1",
                               "Male" = "chartreuse4")) +
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
    ggsave(here("output", "Plots", "Description", "Disease", "morbi_incidence.png"), plot = combined_plot_incidence)

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
