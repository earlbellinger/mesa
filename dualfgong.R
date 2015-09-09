#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(RColorBrewer)
library(magicaxis)
library(sfsmisc)
source('utils.R')

freqs_cols <- c('l', 'n', 'nu')
profile_pattern <- 'profile.+.data$'
freqs_pattern <- 'profile.+-freqs.dat$'

cl <- heat.colors(8)#brewer.pal(5, "Dark2")

exp_dir <- file.path('dualexp')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

seismology <- function(profile_header, freqs, history) {
    history <- history[history$model_number==profile_header$model_number,]
    
    DF <- NULL
    
    DF["acoustic_cutoff"] <- history$acoustic_cutoff/(2*pi)
    DF["delta_nu_scaling"] <- history$delta_nu
    DF["nu_max"] <- history$nu_max
    DF["log_g"] <- history$log_g
    DF["log_surf_z"] <- history$log_surf_z
    
    #DF["radius"] <- max(profile_data$radius)
    DF["age"] <- profile_header$star_age
    DF["mass"] <- profile_header$star_mass
    DF["Teff"] <- profile_header$Teff
    
    converted_fwhm <- (0.66*DF["nu_max"]**0.88)/(2*sqrt(2*log(2)))
    
    freqs <- freqs[freqs$nu < DF["acoustic_cutoff"],]
    # fix radial modes because ADIPLS breaks sometimes
    for (l_mode in 0:3) {
        # grab the relevant l's and n's 
        ell <- freqs[freqs$l==l_mode,]
        ns <- ell$n[ell$n>0]
        # check if any n's are duplicated and if so, shift them down
        if (any(duplicated(ns))) {
            dup <- which(duplicated(ns))[1] # grab duplicated (hopef. only one)
            toshift <- ns[which(ns>0) < dup] # find the ones to shift 
            ell$n[ell$n>0][toshift] <- ns[toshift] - 1 # calculate new n vals
            freqs[freqs$l==l_mode,] <- ell # replace the old with the new 
            freqs <- freqs[!(freqs$l==l_mode & freqs$n==0),] # overwrite data
        }
    }
    
    # calculate delta_nu and d_delta_nu
    for (l_mode in -1:3) {
        strname <- "delta_nu"
        ell <- freqs[freqs$n > 0,] # obtain radial modes
        if (l_mode >= 0) { 
            strname <- paste0(strname, "_", l_mode)
            ell <- freqs[freqs$l==l_mode,]
        }
        
        gaussian_env <- dnorm(ell$nu, DF["nu_max"], converted_fwhm)
        fit <- lm(ell$nu ~ ell$n, weights=gaussian_env)
        DF[strname] <- coef(fit)[2]
        DF[paste0("d_", strname)] <- coef(summary(fit))[2, "Std. Error"]
        if (l_mode > 0)
            DF[paste0("ratio_", l_mode)] <- DF[strname] / DF["delta_nu_0"]
    }
    
    return(DF)
}

parse_dir <- function(directory) {
    print(directory)
    
    # parse dirname string e.g. "M=1.0;Y=0.28"
    params <- NULL
    for (var in unlist(strsplit(basename(directory), ';'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params[nameval[1]] <- as.numeric(nameval[2])
    }
    
    # obtain history
    log_dir <- file.path(directory, "LOGS")
    logs <- list.files(log_dir)
    if (length(logs) <= 1) return(NA)
    history <- read.table(file.path(log_dir, 'history.data'), 
        header=TRUE, skip=5)
    
    # figure out which profiles & frequency files to use
    profile_candidates <- logs[grep(profile_pattern, logs)]
    freq_file_candidates <- logs[grep(freqs_pattern, logs)]
    profile_files <- c()
    freq_files <- c()
    for (profile_file in profile_candidates) {
        freq_name <- sub(".data", "-freqs.dat", profile_file, fixed=TRUE)
        if (freq_name %in% freq_file_candidates) {
            profile_files <- c(profile_files, profile_file)
            freq_files <- c(freq_files, freq_name)
        }
    }
    
    # obtain seismology information
    seis <- do.call(rbind, Map(function(profile_file, freqs_file)
            seismology(#read.table(profile_file, header=TRUE, skip=5), 
                       read.table(profile_file, header=TRUE, nrows=1, skip=1),
                       read.table(freqs_file, col.names=freqs_cols), 
                       history), 
        profile_file=file.path(log_dir, profile_files), 
        freqs_file=file.path(log_dir, freq_files)))
    
    return(merge(rbind(params), seis))
}

DF <- do.call(rbind, Map(parse_dir, 
    directory=list.dirs(exp_dir, recursive=FALSE)))
DF <- DF[complete.cases(DF),]


### correlogram
cairo_pdf(file.path(plot_dir, 'correlogram.pdf'),
    width=11.69, height=8.27, family=font)
plot(DF[,c(1:3, 6:12)], lower.panel=NULL, pch=3, cex=0.1)
dev.off()

### linear model
cairo_pdf(file.path(plot_dir, 'Y_lm.pdf'),
    width=plot_width, height=plot_height, family=font)
fmla <- lm(Y ~ nu_max + exp(log_g) + exp(log_surf_z) + mass + Teff + 
    delta_nu_0 + d_delta_nu_0 + delta_nu_1 + d_delta_nu_1 + 
    delta_nu_3 + delta_nu + d_delta_nu, data=DF)
boxplot(resid(fmla), ylab="Residuals", main="Y ~ observables", pch=3,
    xlab=expression(atop(Y~"~"~M+g+z+T[eff]+nu[max]+
        Delta*nu+Delta*nu[0]+Delta*nu[1]+Delta*nu[3]+"",
        delta*Delta*nu+delta*Delta*nu[0]+delta*Delta*nu[1]+delta*Delta*nu[3])))
abline(h=0, lty=2)
abline(h=0.01, lty=3)
abline(h=-0.01, lty=3)
dev.off()


fmla <- lm(Y ~ nu_max + 
    delta_nu_0 + d_delta_nu_0 + delta_nu_1 + d_delta_nu_1 + 
    delta_nu_3 + delta_nu + d_delta_nu, data=DF)


Y ~ nu_max + exp(log_g) + exp(log_surf_z) + mass + Teff
+ delta_nu_0 + d_delta_nu_0 + delta_nu_1 + d_delta_nu_1 +
delta_nu_3 + delta_nu + d_delta_nu

data <- DF[c(FALSE, TRUE),]#
#data <- DF[DF$age == 4570000000,]
masses <- as.integer(factor(data$M))

### plot Y vs DNU
cairo_pdf(file.path(plot_dir, 'Y_dnu_ZAMS.pdf'), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 3, 1), mgp=c(2, 0.25, 0))
plot(data$Y, data$delta_nu, 
    pch=20, col=cl[masses], tck=0, 
    xlab=expression("helium"~Y), 
    ylab=expression("large frequency separation"~Delta*nu~"["~mu*Hz~"]"),
    main="Large frequency separation of\nZAMS stars by helium and mass")
magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
legend("bottomleft", col=cl, pch=20, paste("M =", unique(data$M)), bty='n')
for (mass_i in masses) {
    mass <- unique(data$M)[mass_i]
    m <- data$M == mass
    lines(data$Y[m], data$delta_nu[m], col=cl[mass_i])
}
dev.off()

### plot Y vs DNUs
cairo_pdf(file.path(plot_dir, 'Y_dnus_ZAMS.pdf'), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 3, 1), mgp=c(2, 0.25, 0))
plot(data$Y, data$ratio_1, 
    ylim=c(0.9995, max(data$ratio_3)),
    pch=2, col=cl[masses], tck=0, 
    xlab=expression("helium"~Y), 
    ylab=expression("large frequency separation ratio"),
    main="Large frequency separation ratios of\nZAMS stars by helium and mass")
magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
legend("bottomleft", ncol=2, bty='n',
       col=c(cl[1:5], rep(1, 5)), 
       pch=c(rep(20, 5), NA, NA, 4:2), 
       legend=c(paste("M =", unique(data$M)), "", "", 
           expression(Delta*nu[3]/Delta*nu[0]), 
           expression(Delta*nu[2]/Delta*nu[0]), 
           expression(Delta*nu[1]/Delta*nu[0])))
for (l_mode in 2:3)
    points(data$Y, unlist(data[paste0("ratio_", l_mode)]),
           pch=l_mode+1, col=cl[masses])
for (mass_i in masses) {
    mass <- unique(data$M)[mass_i]
    m <- data$M == mass
    for (l_mode in 1:3) {
        lines(data$Y[m], unlist(data[paste0("ratio_", l_mode)])[m],
              col=cl[mass_i])
    }
}
dev.off()

