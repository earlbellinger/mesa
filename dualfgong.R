#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(RColorBrewer)
library(magicaxis)
library(sfsmisc)
library(nnet)
source('utils.R')

freqs_cols <- c('l', 'n', 'nu')
profile_pattern <- 'profile.+.data$'
freqs_pattern <- 'profile.+-freqs.dat$'

cl <- heat.colors(8)#brewer.pal(5, "Dark2")

exp_dir <- file.path('dexp-old2')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

seismology <- function(freqs, nu_max, acoustic_cutoff=Inf) {
    seis.DF <- NULL
    
    converted_fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
    
    freqs <- freqs[freqs$nu < acoustic_cutoff,]
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
    for (l_mode in c(-1, sort(unique(freqs$l)))) {
        strname <- "delta_nu"
        ell <- freqs[freqs$n > 0,] # obtain radial modes
        if (l_mode >= 0) { 
            strname <- paste0(strname, "_", l_mode)
            ell <- freqs[freqs$l==l_mode,]
        }
        
        gaussian_env <- dnorm(ell$nu, nu_max, converted_fwhm)
        fit <- lm(ell$nu ~ ell$n, weights=gaussian_env)
        seis.DF[strname] <- coef(fit)[2]
        seis.DF[paste0("d_", strname)] <- coef(summary(fit))[2, "Std. Error"]
    }
    
    return(seis.DF)
}

get_obs <- function(profile_header, freqs, history) {
    hstry_info <- history[history$model_number==profile_header$model_number,]
    if (nrow(hstry_info) == 0) return(NA)
    
    obs.DF <- NULL
    
    obs.DF["age"] <- profile_header$star_age
    obs.DF["radius"] <- profile_header$photosphere_r
    obs.DF["mass"] <- profile_header$star_mass
    obs.DF["L"] <- profile_header$photosphere_L
    obs.DF["Teff"] <- profile_header$Teff
    
    acoustic_cutoff <- hstry_info$acoustic_cutoff/(2*pi)
    #obs.DF["delta_nu_scaling"] <- hstry_info$delta_nu
    
    obs.DF["log_g"] <- hstry_info$log_g
    obs.DF["log_surf_z"] <- hstry_info$log_surf_z
    obs.DF["nu_max"] <- hstry_info$nu_max
    
    seis.DF <- seismology(freqs, obs.DF["nu_max"], acoustic_cutoff)
    
    return(merge(rbind(obs.DF), rbind(seis.DF)))
}

parse_dir <- function(directory) {
    print(directory)
    
    # parse dirname string e.g. "M=1.0_Y=0.28"
    params.DF <- NULL
    for (var in unlist(strsplit(basename(directory), '_'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params.DF[nameval[1]] <- as.numeric(nameval[2])
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
    if (length(profile_files) <= 0) return(NA)
    
    # obtain observable information
    obs.DF <- do.call(plyr:::rbind.fill, Map(function(profile_file, freqs_file) {
            print(profile_file)
            get_obs(read.table(profile_file, header=TRUE, nrows=1, skip=1),
                    read.table(freqs_file, col.names=freqs_cols), 
                    history) }, 
        profile_file=file.path(log_dir, profile_files), 
        freqs_file=file.path(log_dir, freq_files)))
    
    return(merge(rbind(params.DF), obs.DF))
}

DF <- do.call(rbind, Map(parse_dir, 
    directory=list.dirs(exp_dir, recursive=FALSE)))
DF <- DF[complete.cases(DF),]

if (0) {

### inputs
cairo_pdf(file.path(plot_dir, 'correlogram-sobol-inputs.pdf'),
    width=11.69, height=8.27, family=font)
plot(DF[,1:5], lower.panel=NULL, pch=3)
dev.off()

### correlogram
cairo_pdf(file.path(plot_dir, 'correlogram-sobol.pdf'),
    width=11.69, height=8.27, family=font)
plot(DF[,c(1:6, 8:10, 12:14)], lower.panel=NULL, pch=3, cex=0.1)
dev.off()

### helium formula
He.fmla <- Y ~ nu_max + exp(log_g) + exp(log_surf_z) + L + Teff + 
    delta_nu_0 + d_delta_nu_0 + delta_nu_1 + d_delta_nu_1 + 
    delta_nu_2 + d_delta_nu_2 + delta_nu + d_delta_nu
He.lm <- lm(He.fmla, data=DF)
He.lm.resid <- resid(He.lm)
He.nn <- nnet(He.fmla, data=DF, linout=TRUE, size=100, decay=0.1, 
    maxit=10000, MaxNWts=2000)
He.nn.resid <- He.nn$residuals

### lm resid plot
cairo_pdf(file.path(plot_dir, 'Y_lm_resid.pdf'),
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
plot(He.lm.resid ~ DF$Y, pch=3, cex=0.25, tck=0, 
    ylim=c(min(He.nn.resid, He.lm.resid), max(He.nn.resid, He.lm.resid)),
    xlab=expression("initial helium"~Y), 
    ylab=expression(hat(Y) - Y))#,
    #main="A Linear Model for Initial Helium")
magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
abline(h=0, lty=2, col="red")
abline(h=fivenum(He.lm.resid)[2], lty=3, col="red")
abline(h=fivenum(He.lm.resid)[4], lty=3, col="red")
dev.off()

### neural network
cairo_pdf(file.path(plot_dir, 'Y_nn_resid.pdf'),
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
plot(He.nn.resid ~ DF$Y, pch=3, cex=0.25, tck=0, 
    ylim=c(min(He.nn.resid, He.lm.resid), max(He.nn.resid, He.lm.resid)),
    xlab=expression("initial helium"~Y), 
    ylab=expression(hat(Y) - Y))#,
    #main="A Neural Network for Initial Helium")
magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
abline(h=0, lty=2, col="red")
abline(h=fivenum(He.nn.resid)[2], lty=3, col="red")
abline(h=fivenum(He.nn.resid)[4], lty=3, col="red")
dev.off()

### box plots
cairo_pdf(file.path(plot_dir, 'Y_boxplots.pdf'),
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
boxplot(data.frame("Linear Model"=He.lm.resid, 
                   "Neural Network"=He.nn.resid), 
    ylab="Residuals", #main="Y ~ observables", 
    pch=3, cex=0.25)
abline(h=0, lty=2, col="red")
dev.off()



### age formula
age.fmla <- log10(age) ~ nu_max + exp(log_g) + exp(log_surf_z) + L + Teff + 
    delta_nu_0 + d_delta_nu_0 + delta_nu_1 + d_delta_nu_1 + 
    delta_nu_2 + d_delta_nu_2 + delta_nu + d_delta_nu
age.nn <- nnet(age.fmla, data=DF, linout=TRUE, size=100, decay=0.1, 
    maxit=10000, MaxNWts=2000)
age.nn.resid <- age.nn$residuals

### neural network
cairo_pdf(file.path(plot_dir, 'age_nn_resid.pdf'),
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
plot(age.nn.resid ~ log10(DF$age), pch=3, cex=0.25, tck=0, 
    ylim=c(-0.5, 0.5),
    xlab=expression(log[10]~"Age"~"["~G*yr~"]"), 
    ylab=expression(log~hat(t) - log~t))#,
    #main="A Neural Network for Stellar Age")
magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
abline(h=0, lty=2, col="red")
abline(h=fivenum(age.nn.resid)[2], lty=3, col="red")
abline(h=fivenum(age.nn.resid)[4], lty=3, col="red")
dev.off()

### percent difference
cairo_pdf(file.path(plot_dir, 'age_nn_percent_diff.pdf'),
    width=plot_width, height=plot_height, family=font)
par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
percent.diff <- (DF$age-10**age.nn$fitted.values)/DF$age*100
plot(percent.diff ~ log10(DF$age), 
    pch=3, cex=0.25, ylim=c(-100, 100), tck=0,
    xlab=expression(log[10]~"Age"~"["~G*yr~"]"), 
    ylab=expression("Percent Difference"~"["~(hat(t) - t)/t%*%100~"]"))#,
    #main="A Neural Network for Stellar Age")
magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
abline(h=0, lty=2, col="red")
abline(h=fivenum(percent.diff)[2], lty=3, col="red")
abline(h=fivenum(percent.diff)[4], lty=3, col="red")
dev.off()

}

age.fmla2 <- log10(age) ~ nu_max + log_g + log_surf_z + L + Teff + 
    mass + radius + delta_nu_0 + d_delta_nu_0 + delta_nu_1 + d_delta_nu_1 + 
    delta_nu_2 + d_delta_nu_2 + delta_nu + d_delta_nu
age.nn2 <- nnet(age.fmla2, data=DF, linout=TRUE, size=100, decay=0.1, 
    maxit=10000, MaxNWts=2000)
### parse 16 Cyg data
for (cyg in c("16CygA", "16CygB")) {
    obs_data <- read.table(file.path("data", paste0(cyg, "-obs.dat")), 
        header=TRUE)
    attach(obs_data)
    freqs <- read.table(file.path("data", paste0(cyg, "-freqs.dat")), 
        header=TRUE)
    
    ages <- c()
    for (i in 1:5000) {
        obs.DF <- data.frame(
            nu_max     = value[name=="nu_max"],
            log_g      = rnorm(1, value[name=="log_g"], 
                            uncertainty[name=="log_g"]),
            log_surf_z = rnorm(1, value[name=="Fe/H"], 
                            uncertainty[name=="Fe/H"]) - 1.765,
            L          = rnorm(1, value[name=="L"], 
                            uncertainty[name=="L"]),
            Teff       = rnorm(1, value[name=="Teff"], 
                            uncertainty[name=="Teff"]),
            mass       = rnorm(1, value[name=="mass"], 
                            uncertainty[name=="mass"]),
            radius     = rnorm(1, value[name=="radius"], 
                            uncertainty[name=="radius"])
        )
        
        noisy_freqs <- freqs
        noisy_freqs$nu <- rnorm(nrow(freqs), freqs$nu, freqs$dnu)
        
        seis.DF <- seismology(noisy_freqs, 
            obs_data[obs_data$name=='nu_max', 2])
        
        ages <- c(ages, 
            10**predict(age.nn2, merge(rbind(obs.DF), rbind(seis.DF))))
    }
    
    cairo_pdf(file.path(plot_dir, paste0('age_', cyg, '.pdf')),
        width=plot_width, height=plot_height, family=font)
    par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    hist(ages, breaks=50, xlab="Age [ years ]", main="")
    abline(v=mean(ages), col='red')
    abline(v=mean(ages)+sd(ages), col='red', lty=2)
    abline(v=mean(ages)-sd(ages), col='red', lty=2)
    text(#diff(par("usr")[1:2])/2, diff(par("usr")[3:4])/2,
        #4100000000, 150,
        mean(ages)+2*sd(ages), 100,
        paste("Age:", formatC(mean(ages)/1e9), 
              "\n           +/-", formatC(sd(ages)/1e9), "Gyr"))
    dev.off()
    
    detach(obs_data)
}

