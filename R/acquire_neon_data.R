#' @title Acquire NEON data for processing

#' @author
#' John Zobitz \email{zobitz@augsburg.edu}

#' @description
#' Given a site code and dates, apply the neonUtilities package to download the data from NEON API
#' @param site_name Required. NEON code for a particular site (a string)
#' @param start_date Required. Date where we end getting NEON data. Format: YYYY-MM (can't specify day).  So "2020-05" means it will grab data for the entire 5th month of 2020. (a string)
#' @param end_date Required. Date where we end getting NEON data. Format: YYYY-MM (can't specify day).  So "2020-08" means it will grab data for the entire 8th month of 2020. (a string)
#' @param data_file_name Required. Path of location for the file of environmental measurements (includes megapit data and nested data). Must end in .Rda (a string)
#' @param time_frequency Required. Will you be using 30 minute ("30_minute") or 1 minute ("1_minute") recorded data? Defaults to 30 minutes.
#' @param column_selectors Required. Types of measurements we will be computing (typically column_selectors = c("Mean","Minimum","Maximum","ExpUncert","StdErMean"))

#'
#' @example acquire_neon_data("SJER","2020-05","2020-08","my-file.Rda")
#'
#' @import neonUtilities

#' @return Nothing is returned - the file is saved to the location provided

# changelog and author contributions / copyrights
#   John Zobitz (2021-07-22)
#     original creation
#     update to fix auto download (2021-07-25)
#     2022-06-10: update to correct flags on swc

acquire_neon_data <- function(site_name,
                              start_date,
                              end_date,
                              data_file_name,
                              time_frequency = "30_minute",
                              column_selectors = c("Mean","Minimum","Maximum","ExpUncert","StdErMean")
                              ) {

  site_megapit <- neonUtilities::loadByProduct(dpID="DP1.00096.001",
                                               site=site_name,
                                               package="expanded",
                                               check.size = F)


  site_temp <- neonUtilities::loadByProduct(dpID="DP1.00041.001",
                                            site=site_name,
                                            startdate=start_date,
                                            enddate=end_date,
                                            package="expanded",
                                            check.size = F)


  site_swc <- neonUtilities::loadByProduct(dpID="DP1.00094.001",
                                           site=site_name,
                                           startdate=start_date,
                                           enddate=end_date,
                                           package="expanded",
                                           check.size = F)
  # Then correct the swc
  site_swc <- swc_correct(site_swc,site_name)



  site_press <- neonUtilities::loadByProduct(dpID="DP1.00004.001",
                                             site=site_name,
                                             startdate=start_date,
                                             enddate=end_date,
                                             package="expanded",
                                             check.size = F)

  site_co2 <- neonUtilities::loadByProduct(dpID="DP1.00095.001",
                                           site=site_name,
                                           startdate=start_date,
                                           enddate=end_date,
                                           package="expanded",
                                           check.size = F)



  # Process each site measurement
    co2 <- site_co2 |>
      pluck(paste0("SCO2C_",time_frequency)) |>
      select(domainID,siteID,horizontalPosition,verticalPosition,startDateTime,matches(str_c("soilCO2concentration",column_selectors)),finalQF) |>
      rename(soilCO2concentrationFinalQF = finalQF)


    # Determine a data frame of the different horizontal and vertical positions
    co2_positions <- site_co2 |>
      pluck(paste0("sensor_positions_","00095"))

    # Add on the positions for co2
    co2 <- determine_position(co2_positions,co2)

    # Apply monthly means
    co2_monthly_mean <- compute_monthly_mean(co2)

    temperature <- site_temp %>%
      pluck(paste0("ST_",time_frequency)) |>
      select(domainID,siteID,horizontalPosition,verticalPosition,startDateTime,matches(str_c("soilTemp",column_selectors)),finalQF)  |>
      rename(soilTempFinalQF = finalQF)

    # Determine a data frame of the different horizontal and vertical positions
    temperature_positions <- site_temp |>
      pluck(paste0("sensor_positions_","00041"))


    # Add on the positions for temperature
    temperature <- determine_position(temperature_positions,temperature)


    # Apply monthly means
    temperature_monthly_mean <- compute_monthly_mean(temperature)

    swc <- site_swc %>%
      pluck(paste0("SWS_",time_frequency)) |>
      select(domainID,siteID,horizontalPosition,verticalPosition,startDateTime,matches(str_c("VSWC",column_selectors)),VSWCFinalQF)


    # Determine a data frame of the different horizontal and vertical positions

    swc_positions <- site_swc |>
      pluck(paste0("sensor_positions_","00094"))

    # Add on the positions for swc
    swc <- determine_position(swc_positions,swc)




    # Apply monthly means
    swc_monthly_mean <- compute_monthly_mean(swc)

    time_frequency_bp <- if_else(time_frequency == "30_minute","30min","1min")

    pressure <- site_press |>
      pluck(paste0("BP_",time_frequency_bp)) |>
      select(domainID,siteID,horizontalPosition,verticalPosition,startDateTime,matches(str_c("staPres",column_selectors)),staPresFinalQF)

    pressure_positions <- site_press |>
      pluck(paste0("sensor_positions_","00004"))


    # Add on the positions for pressure
    pressure <- determine_position(pressure_positions,pressure)

    # Apply monthly means
    pressure_monthly_mean <- compute_monthly_mean(pressure)


    # Put everything in a nested data frame
    site_data <- tibble(
      data = list(co2,swc,temperature,pressure),
      monthly_mean = list(co2_monthly_mean,swc_monthly_mean,temperature_monthly_mean,pressure_monthly_mean),
      measurement=c("soilCO2concentration","VSWC","soilTemp","staPres")) |>
      mutate(data = map(.x=data,.f=~(.x |> mutate(startDateTime = lubridate::force_tz(startDateTime,tzone="UTC"))))) # Make sure the time zone stamp is in universal time


    save(site_data,site_megapit,file=data_file_name)


}
