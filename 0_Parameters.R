#################################################
#################  PARAMETERS ###################
#################################################


## STEP LENGTH (Murtagh et al., 2020)
step_length <- 0.715*1e-3   # km


## WALKING SPEED (Barban et al, 2022) ----
walk_speed <- 4.8  # km/h


## DRIVING SPEED (Kahlmeier, Götschi et al, 2017) ----
paris_car_speed <- 31  # km/h 
urban_car_speed <-  32
rural_car_speed <- 60 


## BASELINE: 2000 steps daily
baseline_step = 2000
  # Number of minimum weekly steps 
week_base <- 7* 2000 / (80/0.715)     




##############################################################
#                          Diseases                          #
##############################################################

# Disease name
names_disease <- c(
  "cancer" = "Cancer",
  "cvd" = "CVD",
  "dem" = "Dementia",
  "diab2" = "T2 diabetes",
  "dep" = "Depressive episodes",
  "mort" = "All-cause mortality"
)

# Disease colour
colors_disease <- c(
  "cancer" = "firebrick2",      
  "cvd" = "gold" ,
  "dem" = "pink" ,
  "diab2" = "palegreen3",
  "dep" = "slateblue",
  "mort" = "steelblue"
)

# Sex colour
colors_sex <- c(
  "Female" = "darkorange1",
  "Male" = "chartreuse4")



## DEPRESSION ----
  # Duration of a depression episode
  duration_dep <- 2      # months



## MEDICAL COSTS ----
bc_cost <- 44087
cancer_cost <- 14807
cvd_cost <- 55702
dem_cost <- 16839
diab2_cost <- 75201
dep_cost <- NA
mort_cost <- NA




##############################################################
#                       Social cost                          #
##############################################################

## VALUE OF A STATISTICAL LIFE YEAR FOR 2019 FRANCE ----
vsl <- 133000



##############################################################
#                      CO2 emissions                         #
##############################################################

# CO2 emissions per km driven
CO2_emit <- 124                    # 124g CO2 per km








