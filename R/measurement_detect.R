#' @title Makes sure for each time, position, and depth we have at least two data points for soil temp and soil h20, 3 for soil co2

#' @author
#' John Zobitz \email{zobitz@augsburg.edu}
#' based on code developed by Edward Ayres \email{eayres@battelleecology.org}

#' @description
#' Given a merged data frame of co2, water, and temperature:
#' 1) Filters on QF measurement flags
#' 2) Filters if we have at least 2 soil h20 and temperature measurements, 3 co2 measurements at each time, horizontal position, and vertical depth
#' 3) Filters if we have at least 3 distinct measurements at each time and horizontal position
#' 4) Returns the resulting data frame.
#'
#' This as created to speed up data processing.

#' @param input_data Required. Data frame of merged soil water, temperature, and co2 across different NEON locations and depths
#'
#' @example TBD
#'
#' @import dplyr

#' @return Data frame of fluxes from the timeperiod

# changelog and author contributions / copyrights
#   John Zobitz (2022-06-11)
#     original creation



measurement_detect <- function(input_data) {

  # Filter and nest the data
  nested_data <- input_data %>%
    as_tibble() %>%
    filter(finalQF == 0,!is.na(value)) %>%  # Make sure we have the data
    group_by(horizontalPosition, startDateTime, measurement) %>%
    nest()

  # Can we interpolate data?
  have_interpolate <- nested_data %>%
    mutate(tot = map_dbl(data, nrow)) %>%  # count measurement totals
    filter((if_else(measurement == 'co2',tot>2,tot>=1))) %>%  # Filter out the measurements
    select(-tot)

  # Do we have three different measurements (co2, pressure, temp) to compute flux?
  have_measurement <- have_interpolate %>%
    unnest(cols=c(data)) %>%
    group_by(horizontalPosition, startDateTime) %>%
    nest() %>%
    mutate(tot_meas = map_dbl(.x=data, .f=~(length(unique(.x$measurement))) )) %>%
    filter(tot_meas == 3) %>%  # Keep only when we have co2, temp, and water
    select(-tot_meas) %>%
    unnest(cols=c(data))

  return( have_measurement )
}

