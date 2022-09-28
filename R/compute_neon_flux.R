#' @title Compute NEON fluxes at a site

#' @author
#' John Zobitz \email{zobitz@augsburg.edu}
#' based on code developed by Edward Ayres \email{eayres@battelleecology.org}

#' @description
#' Given a site filename (from acquire_neon_data), process and compute fluxes.
#' This file takes a saved data file from acquire:
#' 1) Takes the needed components (QF and measurement flags) for soil water, temperature, co2, binding them together in a tidy data frame
#' 2) Interpolates across the measurements
#' 3) Merges air pressure data into this data frame
#' 4) Does a final QF check so we should have only timeperiods where all measurements exist
#' 5) Adds in the megapit data so we have bulk density, porosity measurements at the interpolated depth.
#' 6) Saves the data

#' @param file_name Required. Path of location for the save file from acquire_neon_data. Must end in .Rda (a string)
#'
#' @example acquire_neon_data("SJER","2020-05","2020-08","my-file.Rda")
#'
#' @import dplyr

#' @return Data frame of fluxes from the timeperiod

# changelog and author contributions / copyrights
#   John Zobitz (2021-07-22)
#     original creation




compute_neon_flux <- function(input_file_name) {

  ################
  # 1) Takes the needed components (QF and measurement flags) for soil water, temperature, co2, binding them together in a tidy data frame

  # Load up the data (this may take a while)
  load(input_file_name)

  soil_water <- measurement_merge(site_swc,"SWS","00094","VSWCMean","VSWCFinalQF") %>%
    mutate(measurement = "soil_water") %>%
    rename(value = VSWCMean,
           finalQF = VSWCFinalQF)

  temperature <- measurement_merge(site_temp,"ST","00041","soilTempMean","finalQF") %>%
    mutate(measurement = "temperature") %>%
    rename(value = soilTempMean)

  co2 <- measurement_merge(site_co2,"SCO2C","00095","soilCO2concentrationMean","finalQF") %>%
    mutate(measurement = "co2") %>%
    rename(value = soilCO2concentrationMean)


  # Bind these all up together
  site <- rbind(soil_water,temperature,co2)



  # Get the sensor positions for CO2:
  site_co2_positions <- sensor_positions(site_co2,"SCO2C","00095")

  ################


  ################
  # 2) Interpolates across the measurements

  # Filters out measurements that don't have enough QF flags
  site_filtered <- measurement_detect(site)

  if(nrow(site_filtered) > 0) {

    # Interpolate all the measurements together in one nested function
    site_interp <- depth_interpolate(site_filtered)

    ################
    # 3) Merges air pressure data into this data frame
    # Next we need to join the pressure data - this is just one measurement by time
    pressure <- site_press$BP_30min %>%   # This is where we get the different time periods
      select(startDateTime,staPresMean,staPresFinalQF) %>%
      mutate(measurement = "pressure") %>%
      rename(value = staPresMean,
             finalQF = staPresFinalQF) %>%
      filter(finalQF == 0) %>%
      select(-finalQF) %>%
      group_by(startDateTime)# %>%
#      nest()

    # Then we will want to join the times where the pressure is - for ease of use it will be at all the depths using some of the joining and pivoting skills here.

    site_depth_nest <- site_interp %>% pivot_wider(names_from = "measurement",values_from = "value") %>%
      inner_join(pressure,b=c("startDateTime")) %>%
      #unnest(cols=c("data")) %>%
      select(-measurement) %>%
      rename(pressure=value) %>%
      pivot_longer(cols=c("pressure","co2","temperature","soil_water"),
                   names_to = "measurement", values_to = "value") %>%
      group_by(zOffset) %>%
      nest()



    ################

    ################
    # 5) Adds in the megapit data so we have bulk density, porosity measurements at the interpolated depth.

    # Pull in the megapit data

    # Merge the soil properties into a single data frame
    biogeo_sample <- site_megapit$mgp_perbiogeosample %>%
      inner_join(site_megapit$mgp_perbulksample , by=c("horizonID", "pitID", "domainID", "siteID", "pitNamedLocation", "horizonName",  "laboratoryName", "labProjID", "setDate", "collectDate")) %>%
      select(c("horizonID", "pitID",
               "coarseFrag2To5","coarseFrag5To20","biogeoTopDepth", "bulkDensExclCoarseFrag","biogeoBottomDepth",
               "biogeoCenterDepth")) %>%
      mutate(across(.cols=matches("biogeo"),~-.x/100)) %>%
      drop_na()

    # Now we should go across the nested depths and have the horizon
    # Kicking it old school with the double loop (the indices are small, so that is ok.)
    for(i in seq_along(site_depth_nest$zOffset)) {
      for(j in 1:dim(biogeo_sample)[1]) {
        if(between(site_depth_nest$zOffset[i],biogeo_sample$biogeoBottomDepth[j],biogeo_sample$biogeoTopDepth[j])) {
          site_depth_nest$data[[i]] <- site_depth_nest$data[[i]] %>% mutate(biogeo_sample[j,])
        }
      }
    }

    ################

    ################
    # 6) Saves the data
    site_final_interp <- site_depth_nest %>%
      unnest(cols=c(data))

    #7) Compute the fluxes
    out_dates <- tibble(startDateTime = seq(min(site_final_interp$startDateTime),max(site_final_interp$startDateTime),by="30 min"))

    # Fill in where there is no flux measurement
    out_fluxes <- neon_site_flux(site_final_interp,site_co2_positions) %>%
      group_by(horizontalPosition) %>%
      nest() %>% # Need to nest the data so that we have the horizontal positions correct
      mutate(data = map(.x=data,.f=~right_join(.x,out_dates,by="startDateTime") %>% arrange(startDateTime))) %>%  # Join the dates to each horizontal position, arrange by date
      unnest(cols=c(data))



    return(out_fluxes)

  } else {

    out_fluxes <- NA
    return(NA)
  }


  ################



}
