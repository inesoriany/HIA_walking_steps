#################################################
#################  PARAMETERS ###################
#################################################


## STEP LENGTH (Murtagh et al., 2020)
step_length <- 0.715*1e-3   # km


## WALKING SPEED (Barban et al, 2022) ----
walk_speed <- 4.8  # km/h


## BASELINE: 2000 steps daily
baseline_step = 2000



##############################################################
#                          Diseases                          #
##############################################################

# Disease name
names_disease <- c(
  "cancer" = "Cancer",
  "cvd" = "CVD",
  "dem" = "Dementia",
  "diab2" = "T2 diabetes",
  "dep" = "MDD",
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
  # Proportion of adults that relapse
  prop_relapse <- 0.49

  # Duration of recovery
  duration_recovery <- 5      # years


## MEDICAL COSTS ----
bc_cost <- 46242
cc_cost <- 25126
cancer_cost <- 14807
cvd_cost <- 55702
dem_cost <- 16839
diab2_cost <- 75201
dep_cost <- NA
mort_cost <- NA


##############################################################
#              Sensitivity Analysis - Diseases               #
##############################################################

# Colon cancer (Harris et al., 2009)
ref_cc_m <- 168*30.1/11.25                                  # Disease reference volume men (in min)
rr_cc_men_lb <-.67                                          # Disease relative risk for men (upper bound)
rr_cc_men <- .80                                            # Disease relative risk for men
rr_cc_men_ub <- .96                                         # Disease relative risk for men (lower bound)
ref_cc_w <- 168*30.9/11.25                                  # Disease reference volume women
rr_cc_women_lb <- .76                                       # Disease relative risk for women (upper bound)
rr_cc_women <-  .86                                         # Disease relative risk for women 
rr_cc_women_ub <-.98                                        # Disease relative risk for women (lower bound)

# Dementia (Hamer & Chida, 2009)
ref_dem = ref_dem_m = ref_dem_w <- 168*33/11.25             # same for men and women
rr_dem_lb = rr_dem_men_lb =rr_dem_women_lb <-.6 
rr_dem = rr_dem_men = rr_dem_women <- .72
rr_dem_ub = rr_dem_men_ub =rr_dem_women_ub<-.86 

# Breast cancer (Monninkhof et al., 2007)
ref_bc_w <- 60 
rr_bc_women_lb <- .92 
rr_bc_women <- .94
rr_bc_women_ub <-.97

ref_bc_m <- NA
rr_bc_men_ub <-NA
rr_bc_men <- NA
rr_bc_men_lb <- NA

# Cardiovascular disease (Hamer & Chida, 2008)
ref_cvd = ref_cvd_m = ref_cvd_w <- 180                       # 3h per week of physical activity of moderate intensity
rr_cvd_lb = rr_cvd_men_lb =rr_cvd_women_lb<-.79
rr_cvd = rr_cvd_men = rr_cvd_women<- .84
rr_cvd_ub = rr_cvd_men_ub =rr_cvd_women_ub<-.90

# type 2 diabetes (Jeon et al., 2007)
ref_diab2 = ref_diab2_m =ref_diab2_w <- 168*10/11.25
rr_diab2_lb = rr_diab2_men_lb =rr_diab2_women_lb<-.75
rr_diab2 = rr_diab2_men = rr_diab2_women <- .83
rr_diab2_ub = rr_diab2_men_ub =rr_diab2_women_ub<- .91

# Depression (Pearce et al, 2022)
ref_dep = ref_dep_m = ref_dep_w <- 168                           # 168 minutes per week for walking
rr_dep_women_lb <- 0.7
rr_dep_women <-  0.75
rr_dep_women_ub <-  0.8

rr_dep_men_lb <- 0.76
rr_dep_men <- 0.8
rr_dep_men_ub <- 0.85

# mortality (Kelly et al. 2014) 
ref_mort = ref_mort_m =ref_mort_w <- 168                         # 168 minutes per week for walking
rr_mort_lb=  rr_mort_men_lb = rr_mort_women_lb<-.85         
rr_mort = rr_mort_men =  rr_mort_women<- .90
rr_mort_ub =rr_mort_men_ub =rr_mort_women_ub <-.95




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



##############################################################
#                   Good students targets                    #
##############################################################
# Part modale kilométrique après transformation x Distance journalière moyenne
  # (Rapport ADEME 2015 - "Contribution de la marche et du vélo à la décarbonation et à l'amélioration de la qualité de l'air en France")
urban_target <- 9/100 * 15 / step_length        # 1.35 km per person per day, in steps
periurban_target <- 4/100 * 22 / step_length    # 0.88 km per person per day, in steps
rural_target <- 2/100 * 32 / step_length        # 0.64 km per person per day, in steps
