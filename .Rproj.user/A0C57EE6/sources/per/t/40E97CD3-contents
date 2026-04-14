###############################################
############ EXTRACR EMP DATABASE #############
###############################################

#This code imports EMP data files from csv and puts together a dataset at the individual level, with information on:
# - indiviudal characteristics
# - total number of kilometres by bike in local mobility + split by regular vs e-bike
# - household characteristics, including their bike equipment

################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  tidyverse,
  haven)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Household data
household <- import(here("emp_2019_donnees_individuelles_anonymisees_novembre2024", "EMP", "tcm_men_public_V3.csv" ))

# Individual data
ind_kish <- import(here("emp_2019_donnees_individuelles_anonymisees_novembre2024", "EMP", "k_individu_public_V3.csv" ))

# Individual characteristics
ind <- import(here("emp_2019_donnees_individuelles_anonymisees_novembre2024", "EMP", "tcm_ind_kish_public_V3.csv" ))

# Trip data
trip <- import(here("emp_2019_donnees_individuelles_anonymisees_novembre2024", "EMP", "5. k_deploc_public_V4.csv" ))




################################################################################################################################
#                                                     3. SELECT DATA                                                           #
################################################################################################################################

household <- household %>%
  select(ident_men, 
         dep_res, 
         nuts_res,
         catcom_aa_res, 
         tuu2017_res,
         statutcom_uu_res, 
         type_uu_res)
# degré de commune de résidence

ind_kish <- ind_kish %>% 
  select(ident_men, 
         ident_ind, 
         pond_indc)
  

ind <- ind %>% 
  select(
    ident_ind,
    ident_men,
    sexe, 
    age, 
    ddipl, 
    typolog, 
    handicap,
    cs42, 
    cs24,
    typemploi
  )


trip <- trip %>% 
  select(ident_ind, 
         ident_men, 
         ident_dep,
         pond_jour,
         mdisttot_fin,                               # Trip length
         mtempsmap,                                  # Walking duration
         starts_with("mmoy"),                        # Means of transportation
         mtp,                                        # Main mean of transportation
         stautcom_uu_ori) %>%                        # Departure commune status        
  mutate(mmoy2s = na_if(mmoy2s, "vide"),
         mmoy2s = as.numeric(mmoy2s))



##############################################################
#                          WALKING                           #
##############################################################
walk <- trip %>% 
  mutate(walking_lower = mtp %in% c(1.1, 1.2, 1.3, 1.4),
         
         walking_upper = mtp %in% c(1.1, 1.2, 1.3, 1.4) |
           mmoy1s %in% c(1.1, 1.2, 1.3, 1.4) |
           mmoy2s %in% c(1.1, 1.2, 1.3, 1.4),
         
         nbkm_walking_lower = if_else(walking_lower, mdisttot_fin, 0),
         nbkm_walking_upper = if_else(walking_upper, mdisttot_fin, 0)) %>%
  group_by(ident_men, ident_ind, pond_jour) %>%
  summarise(nbkm_walking_lower = sum(nbkm_walking_lower, na.rm = TRUE),
            nbkm_walking_upper = sum(nbkm_walking_upper, na.rm = TRUE),
            .groups = "drop")



##############################################################
#                          DRIVING                           #
##############################################################

select(
  ident_men, ident_ind, pond_jour,
  mdisttot_fin, voiture_passager,
  n_trip, nb_trip_voit,
  conduit_moins10km_slt, conduit_moins10km,
  conduit_moins5km, conduit_moins5km_slt
) %>%  pivot_wider(
  id_cols = c(
    ident_men, ident_ind, pond_jour,
    nb_trip_voit),
  names_from = n_trip,
  values_from = c(mdisttot_fin, voiture_passager))




