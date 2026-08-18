##############################################
################ DESCRIPTION #################
##############################################



# Reprendre les données de l'EMP : la marche est-elle le premier mode de transport ?

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

# Survey design ponderated by individual
indiv_walkers <- emp_walkers %>% 
  filter (pond_indc != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



################################################################################################################################
#                                                        5. HEALTH                                                             #
################################################################################################################################
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
emp_20_89 <-  emp %>% 
  filter(age >= 20 & age <90)

pop_tot <- sum(emp_20_89$pond_indc)
pop_tot



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
km_mean
km_mean/step_length                # Mean number of steps per day


km_mean_IC <- confint(km_mean)
km_mean_IC
km_mean_IC/step_length


# In exposure time (min)
km_mean * 60 / walk_speed


# Mean walking distances per day, by age group
svyby(~nbkm_tot_walking, by = ~age_grp10, jour_walkers, svymean, na.rm = T)


## Mean walking distances per day, per age group and per sex ----
mean_distance_people <- svyby(~nbkm_tot_walking, by = ~sex + age_grp10, jour_walkers, svymean, na.rm = T)


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



################################################################################################################################
#                                                     8. WALKING - steps                                                       #
################################################################################################################################



################################################################################################################################
#                                                      5. EXPORT DATA                                                          #
################################################################################################################################

# HEALTH
    # Incidence
    ggsave(here("output", "Plots", "Description", "Disease", "morbi_incidence.png"), plot = combined_plot_incidence)

# WALKING
    # Sex proportion
    export(prop_sex, here("output", "Tables", "Description", "sex_proportion.xlsx"))
    # Distance
    ggsave(here("output", "Plots", "Description", "Walk", "plot_mean_km_walkers.png"), plot = mean_km_walkers)
