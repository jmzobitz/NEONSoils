library(tidyverse)
library(NEONSoils)

# --> YOU ADD: Define the name of the site:
sites <- c("WREF")  # Add in the 4 letter codes for the other sites WREF

# --> Create the dates vector You can shorten these to a single month if you want.
dates <- c("2021-01","2021-02")

# Combine these two together:
site_date_combine <- expand_grid(sites,dates)

# Now we go through and do the dirty work of saving and computing fluxes. Yay .... :)
for(i in seq_along(site_date_combine)) {

  # Name current month (you will need to adjust this on your computer)
  curr_site <- site_date_combine$sites[[i]]
  curr_date <- site_date_combine$dates[[i]]
  acquire_name <- paste0('data/',curr_site,"-NEON-flux-",curr_date,".Rda")
  env_name <- paste0('data/',curr_site,"-env-meas-",curr_date,".Rda")
  flux_name <- paste0('data/',curr_site,"-out-flux-",curr_date,".Rda")
  print(env_name)
  # Process
  acquire_neon_data(site_name =curr_site,
                        start_date = curr_date,
                        end_date = curr_date,
                        data_file_name = acquire_name,
                        env_file_name = env_name)

   # out_fluxes <- compute_neon_flux("data/curr-flux.Rda")

  #  save(out_fluxes,file=flux_name)

}

