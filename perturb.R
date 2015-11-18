source('seismology.R')

dir.create('perturb', showWarnings=FALSE)

speed_of_light <- 299792 # km/s

### Obtain properties of real stars varied within their uncertainties 
monte_carlo_perturbations <- function(obs_data_file, freqs_data_file,
        n_perturbations=10000) {
    freqs <- read.table(freqs_data_file, header=TRUE)
    obs_data <- read.table(obs_data_file, header=TRUE)
    noisy_freqs <- freqs
    parallelStartMulticore(max(1, detectCores()))
    do.call(rbind, with(obs_data, {
        parallelMap(function(n) {
            # Perturb observations by their uncertainties
            
            obs.DF <- data.frame(
                radius     = rnorm(1, value[name=="radius"], 
                ifelse(n==1, 0, uncertainty[name=="radius"])),
                L          = rnorm(1, value[name=="L"], 
                ifelse(n==1, 0, uncertainty[name=="L"])),
                Teff       = rnorm(1, value[name=="Teff"], 
                ifelse(n==1, 0, uncertainty[name=="Teff"])),
                Fe_H       = rnorm(1, value[name=="Fe/H"], 
                ifelse(n==1, 0, uncertainty[name=="Fe/H"])),
                nu_max     = rnorm(1, value[name=="nu_max"], 
                ifelse(n==1, 0, uncertainty[name=="nu_max"]))
            )
            
            # Correct frequencies for Doppler shift
            radial_velocity <- with(obs_data, 
                rnorm(1, value[name=="radial_velocity"], 
                ifelse(n==1, 0, uncertainty[name=="radial_velocity"])))
            doppler_beta <- radial_velocity/speed_of_light
            doppler_shift <- sqrt((1+doppler_beta)/(1-doppler_beta))
            
            # Perturb frequencies
            noisy_freqs$nu <- rnorm(nrow(freqs), freqs$nu * doppler_shift, 
                ifelse(n==1, 0, freqs$dnu * doppler_shift))
            
            # Calculate Dnu, dnus, and ratios
            seis.DF <- seismology(noisy_freqs, obs.DF$nu_max)
            merge(rbind(obs.DF), rbind(seis.DF))
        }, 1:n_perturbations)
    }))
}

# Perturb every star 10k times and save the results
star_names <- c("Tagesstern", "16CygA", "16CygB", "Sun")
stars <- list()
for (star in star_names) {
    fname <- file.path('perturb', paste0(star, "_perturb.dat"))
    #stars[[star]] <- #if (file.exists(fname)) {
        #read.table(fname, header=TRUE)
    #} else {
    perturbations <- monte_carlo_perturbations(
        file.path("..", "data", paste0(star, "-obs.dat")),
        file.path("..", "data", paste0(star, "-freqs.dat")))
    write.table(perturbations, fname, quote=FALSE, 
                sep='\t', row.names=FALSE)
    #}
}
