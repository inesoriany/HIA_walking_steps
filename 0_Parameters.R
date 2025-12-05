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
  # Number of minimum weekly steps 
week_base <- 7* 2000 / (80/0.715)     

##############################################################
#                          Diseases                          #
##############################################################

# Disease name
names_disease <- c(
  "bc" = "Breast cancer incidence",
  "cancer" = "Cancer incidence",
  "cvd" = "Cardiovascular disease incidence",
  "dem" = "Dementia",
  "diab2" = "Type 2 diabetes",
  "dep" = "Depressive episodes",
  "mort" = "All-cause mortality"
)

# Disease colour
colors_disease <- c(
  "bc" = "firebrick2",
  "cancer" = "darkorange",      
  "cvd" = "gold" ,
  "dem" = "pink" ,
  "diab2" = "palegreen3",
  "dep" = "slateblue",
  "mort" = "steelblue"
)



## BREAST CANCER ----
# Relative risks (Monninkhof et al., 2007)
ref_bc_w <- 60 
rr_bc_women_low <- .92 
rr_bc_women <- .94
rr_bc_women_up <-.97

ref_bc_m <- NA
rr_bc_men_up <-NA
rr_bc_men <- NA
rr_bc_men_low <- NA



## MEDICAL COSTS ----
cc_cost <- 26716
dem_cost <- 22748
bc_cost <- 46968
cancer_cost <- 14808
cvd_cost <- 20938
diab2_cost <- 36514
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








