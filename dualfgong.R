#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(RColorBrewer)
library(magicaxis)
library(sfsmisc)

fgong_cols <- c('l', 'n', 'nu')
profile_pattern <- 'profile_.+.data$'
fgong_pattern <- 'profile_.+-freqs.dat$'

ell_cl <- brewer.pal(4, "BrBG")

exp_dir <- file.path('dualexp')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

seismology <- function(profile, profile_header, fgong, history) {
    history <- history[history$model_number==profile_header$model_number,]
    
    DF <- NULL
    DF["acoustic_cutoff"] <- history$acoustic_cutoff/(2*pi)
    DF["delta_nu_scaling"] <- history$delta_nu
    DF["nu_max"] <- history$nu_max
    DF["radius"] <- max(profile$radius)
    DF["age"] <- profile_header$star_age
    
    converted_fwhm <- (0.66*DF["nu_max"]**0.88)/(2*sqrt(2*log(2)))
    
    # fix radial modes because ADIPLS breaks sometimes
    for (l_mode in 0:3) {
        # grab the relevant l's and n's 
        ell <- fgong[fgong$l==l_mode,]
        ns <- ell$n[ell$n>0]
        # check if any n's are duplicated and if so, shift them down
        if (any(duplicated(ns))) {
            dup <- which(duplicated(ns))[1] # grab duplicated (hopef. only one)
            toshift <- ns[which(ns>0) < dup] # find the ones to shift 
            ell$n[ell$n>0][toshift] <- ns[toshift] - 1 # calculate new n vals
            fgong[fgong$l==l_mode,] <- ell # replace the old with the new 
            fgong <- fgong[!(fgong$l==l_mode & fgong$n==0),] # overwrite data
        }
    }
    
    # calculate delta_nu and d_delta_nu
    for (l_mode in -1:3) {
        strname <- "delta_nu"
        ell <- fgong[fgong$n > 0,] # obtain radial modes
        if (l_mode >= 0) { 
            strname <- paste0(strname, "_", l_mode)
            ell <- fgong[fgong$l==l_mode,]
        }
        
        gaussian_env <- dnorm(ell$nu, DF["nu_max"], converted_fwhm)
        fit <- lm(ell$nu ~ ell$n, weights=gaussian_env)
        DF[strname] <- coef(fit)[2]
        DF[paste0("d_", strname)] <- coef(summary(fit))[2, "Std. Error"]
    }
    
    return(DF)
}

parse_dir <- function(directory) {
    # parse dirname string e.g. "M=1.0;Y=0.28"
    params <- NULL
    for (var in unlist(strsplit(basename(directory), ';'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params[nameval[1]] <- as.numeric(nameval[2])
    }
    
    # obtain history
    log_dir <- file.path(directory, "LOGS")
    history <- read.table(file.path(log_dir, 'history.data'), 
        header=TRUE, skip=5)
    
    # obtain seismology information
    logs <- list.files(log_dir)
    seis <- do.call(rbind, Map(function(profile, fgong)
            seismology(read.table(profile, header=TRUE, skip=5), 
                       read.table(profile, header=TRUE, nrows=1, skip=1),
                       read.table(fgong, col.names=fgong_cols), 
                       history), 
        profile=file.path(log_dir, logs[grep(profile_pattern, logs)]), 
        fgong=file.path(log_dir, logs[grep(fgong_pattern, logs)])))
    
    return(merge(rbind(params), seis))
}

data <- do.call(rbind, Map(parse_dir, 
    directory=list.dirs(exp_dir, recursive=FALSE)))

solar.age <- data[data$age > 100000000,]

### plot Y vs DNU


