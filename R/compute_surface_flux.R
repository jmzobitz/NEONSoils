#' @title Compute CO2 surface flux

#' @author
#' John Zobitz \email{zobitz@augsburg.edu}

#' @description
#' Computation function. Given a measurement of the co2 and diffusive flux at different levels, return the surface flux

#' @import dplyr

#' @param input_data Required. A data frame containing zOffsets, diffusivity, and co2 (umol mol-1)

#' @return A value of the surface CO2 flux (umol m-2 s-1)

# changelog and author contributions / copyrights
#   John Zobitz (2021-07-21)
#     original creation
#   John Zobitz (2022-06-21)
#     revision to account for correct directionality when computing fluxes

compute_surface_flux <- function(input_data) {

  # Select the first three levels (closest to the surface)
  input_data_rev <- input_data %>%
    arrange(desc(zOffset)) %>%
    slice(1:3) %>%
    mutate(zOffset = abs(zOffset))  # Make these positive, just to avoid problems down the road.

  # Soil co2 grad check
  if(input_data_rev$co2mol[1] > input_data_rev$co2mol[2] | input_data_rev$co2mol[2] > input_data_rev$co2mol[3] | input_data_rev$co2mol[1] > input_data_rev$co2mol[3]){
    return(NA) #if sensors are in the wrong order return NA
  }

  # Calculate diffusivity of the first two levels as harmonic mean
  diffusive_mean <- c( (sum(input_data_rev$diffus[1:2]^(-1))/2)^(-1), (sum(input_data_rev$diffus[2:3]^(-1))/2)^(-1) )

  # Calculate diffusive flux between CO2 measurement levels (µmol m-2 s-1).
  # This is Fick's law, so Flux = -Da *(gradient)
  # Note: Because the depth is a positive number, graphs of CO2 vs depth will be increasing, so a negative diffusive flux means that the flux moves DOWN (negative)  ** sigh Calculus **
  diffusive_flux <- -diff(input_data_rev$co2mol)/diff(input_data_rev$zOffset) * diffusive_mean

  # Determine the average depth (m) between measurement levels so that this can be assigned to the flux measurement
  diffusive_flux_depth <- diff(input_data_rev$zOffset)/2 + input_data_rev$zOffset[1:2]


  # Now compute the fluxes, adding on the harmonic mean, fluxes, and depths
  flux_depth <- input_data_rev %>%
    slice(1:2) %>%
    mutate(diff_hm = diffusive_mean, # hm = harmonic mean
           flux = diffusive_flux,
           diff_depth = diffusive_flux_depth)






  # Calculate soil CO2 diffusive flux at the soil surface (i.e., flux to the atmosphere) by linear extrapolation of the fluxes to the soil surface (µmol m-2 s-1).
  surface_flux <- flux_depth$flux[1] - (((flux_depth$flux[1] - flux_depth$flux[2]) / (flux_depth$diff_depth[1] - flux_depth$diff_depth[2])) * flux_depth$diff_depth[1])

  # Now when we have a positive surface flux, that means CO2 is moving down to the soil. (or UP the depth vs concentration plot) Flux gradient method will be violated, let's just return NA

  if(surface_flux > 0) {return(NA)}
  else {
    return(-surface_flux) # Make the surface flux positive for consistency * sigh * minus signs
  }



}




