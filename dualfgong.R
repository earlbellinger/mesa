#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

options(error=traceback)

library(RColorBrewer)
library(magicaxis)
library(sfsmisc)
library(nnet)
library(parallel)
library(parallelMap)
library(stargazer)
library(plyr)
library(matrixStats)
library(mblm)

plot_width <- 9
plot_height <- 6
font <- "Palatino"

freqs_cols <- c('l', 'n', 'nu', 'inertia')
profile_pattern <- 'profile.+.data$'
freqs_pattern <- 'profile.+-freqs.dat$'

exp_dir <- file.path('diffusion')
plot_dir <- file.path('plots3')
dir.create(plot_dir, showWarnings=FALSE)

Z_div_X_solar = 0.02293

latex_labs <- c(expression(M[0]), 
                expression(Y[0]),
                expression(Z[0]),
                expression(alpha["MLT"]),
                "age",
                "Y",
                "radius",
                expression("L"),
                expression(T["eff"]),
                expression(log~g),
                expression("Fe/H"),
                expression(nu[max]),
                expression(Delta*nu[0]),
                expression(epsilon[0]),
                expression(delta*nu[0*","*2]),
                expression(Delta*nu[1]),
                expression(epsilon[1]),
                expression(delta*nu[1*","*3]),
                expression(Delta*nu[2]),
                expression(epsilon[2]),
                expression(Delta*nu[3]),
                expression(epsilon[3]))

### Obtain seismic information from frequencies 

separation <- function(first_l, first_n, second_l, second_n, df) {
  first <- df$l == first_l & df$n == first_n
  second <- df$l == second_l & df$n == second_n
  if (sum(first) == 1 && sum(second) == 1)
    return(df[first,]$nu - df[second,]$nu)
  return(NA)
}

dnu <- function(l, n, df) separation(l, n, l+2, n-1, df)
Dnu <- function(l, n, df) separation(l, n, l, n-1, df)
r02 <- function(n, df) dnu(0, n, df) / Dnu(1, n, df)
r13 <- function(n, df) dnu(1, n, df) / Dnu(0, n+1, df)

get_averages <- function(f, df, ell, freqs, l=NA) {
  sep_name <- deparse(substitute(f))
  if (sep_name == 'dnu' || sep_name == 'Dnu') {
    a <- sapply(unique(ell$n), function(n) f(l, n, freqs))
    sep_name <- paste0(sep_name, '_', l)
  } else {
    a <- sapply(unique(ell$n), function(n) f(n, freqs))
  }
  not.nan <- complete.cases(a)
  a <- a[not.nan]
  b <- ell$n[not.nan]
  df[paste0(sep_name, "_median")] <- median(a)
  fit <- mblm(a~b, repeated=1)
  df[paste0(sep_name, "_slope")] <- coef(fit)[2]
  df[paste0(sep_name, "_intercept")] <- coef(fit)[1]
  df
}

seismology <- function(freqs, nu_max, acoustic_cutoff=Inf) {
  if (nrow(freqs) == 0) return(NULL)
  
  seis.DF <- NULL
  
  converted_fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
  
  freqs <- unique(freqs[freqs$nu < acoustic_cutoff,])
  # fix radial modes because ADIPLS breaks sometimes
  for (l_deg in 0:3) {
    # grab the relevant l's and n's 
    ell <- freqs[freqs$l==l_deg,]
    ns <- ell$n[ell$n>0]
    # check if any n's are duplicated and if so, shift them down
    if (any(duplicated(ns))) {
      return(NULL) # screw it, just discard this data
      dup <- which(duplicated(ns))[1] # grab duplicated (hopef. only one)
      toshift <- ns[which(ns>0) < dup] # find the ones to shift 
      ell$n[ell$n>0][toshift] <- ns[toshift] - 1 # calculate new n vals
      freqs[freqs$l==l_deg,] <- ell # replace the old with the new 
      freqs <- freqs[!(freqs$l==l_deg & freqs$n==0),] # overwrite data
    }
  }
  
  # calculate frequency separations
  ells <- sort(unique(freqs$l))
  for (l_deg in ells) {
    ell <- freqs[freqs$n > 0 & freqs$l==l_deg,] # pressure modes
    #gaussian_env <- dnorm(ell$nu, nu_max, converted_fwhm)
    #fit <- lm(ell$nu ~ ell$n, weights=gaussian_env)
    fit <- mblm(nu ~ n, dataframe=ell, repeated=TRUE)
    seis.DF[paste0("Dnu_", l_deg)] <- coef(fit)[2]
    seis.DF[paste0("eps_", l_deg)] <- coef(fit)[1]/coef(fit)[2]
    
    if (any(freqs$l==l_deg+2)) {
      seis.DF <- get_averages(dnu, seis.DF, ell, freqs, l_deg)
      #a <- sapply(unique(ell$n), function(n) small_sep(l_deg, n, freqs))
      #not.nan <- complete.cases(a)
      #a <- a[not.nan]
      #b <- ell$n[not.nan]
      #seis.DF[paste0("dnu_", l_deg)] <- median(a)
      #fit <- mblm(a~b, repeated=1)
      #seis.DF[paste0("dnu_slope_", l_deg)] <- coef(fit)[2]
      #seis.DF[paste0("dnu_intercept_", l_deg)] <- coef(fit)[1]
    #}
    }
  }
  seis.DF <- get_averages(r02, seis.DF, ell, freqs)
  seis.DF <- get_averages(r13, seis.DF, ell, freqs)
  
#   # calculate frequency ratios
#   a <- sapply(unique(ell$n), function(n) r02(n, freqs))
#   not.nan <- complete.cases(a)
#   a <- a[not.nan]
#   b <- freqs[freqs$l==0,]$n[not.nan]
#   seis.DF["r02_median"] <- median(a)
#   fit <- mblm(a~b, repeated=1)
#   seis.DF["r02_slope"] <- coef(fit)[2]
#   seis.DF["r02_intercept"] <- coef(fit)[1]
#   
#   a <- sapply(unique(ell$n), function(n) r13(n, freqs))
#   not.nan <- complete.cases(a)
#   a <- a[not.nan]
#   b <- freqs[freqs$l==0,]$n[not.nan]
#   seis.DF["r02_median"] <- median(a)
#   fit <- mblm(a~b, repeated=1)
#   seis.DF["r13_slope"] <- coef(fit)[2]
#   seis.DF["r13_intercept"] <- coef(fit)[1]
    
#     ell2 <- freqs[freqs$n > 0 & freqs$l==l_deg+2,]
#     if (nrow(ell2) > 0) {
#       diffs <- c()
#       nus <- c()
#       for (ii in 1:nrow(ell)) {
#         row_i <- ell[ii,]
#         if (any(ell2$n == row_i$n-1)) {
#           dnu_n <- (row_i$nu - ell2$nu[ell2$n == (row_i$n-1)])[1]
#           if (dnu_n > 0 && dnu_n < 50) {
#             diffs <- c(diffs, dnu_n)
#             nus <- c(nus, row_i$nu)
#           }
#         }
#       }
#       seis.DF[paste0("dnu_", l_deg)] <- 
#         if (length(nus)>0 && length(diffs) == length(nus)) {
#           gaussian_env <- dnorm(nus, nu_max, converted_fwhm)
#           #shifted <- nus - nu_max
#           #fit <- lm(diffs ~ shifted, weights=gaussian_env)
#           #coef(fit)[2]
#           #fit <- mblm(diffs ~ shifted)
#           #coef(fit)[1]
#           #weighted.mean(diffs, gaussian_env)
#           weightedMedian(diffs, gaussian_env)
#         } else { NA }
#     }
#     #seis.DF[paste0("d_", strname)] <- coef(summary(fit))[2, "Std. Error"]
# }
  return(seis.DF)
}

### Obtain observable properties from models 
get_obs <- function(profile_file, freqs_file, ev_history) {
    print(profile_file)
    
    profile_header <- read.table(profile_file, header=TRUE, nrows=1, skip=1)
    freqs <- read.table(freqs_file, col.names=freqs_cols, fill=TRUE)
    freqs <- freqs[complete.cases(freqs),]
    hstry <- ev_history[ev_history$model_number==profile_header$model_number,]
    
    if (nrow(hstry) == 0) return(NA)
    
    obs.DF <- NULL
    
    obs.DF["age"] <- profile_header$star_age
    
    #obs.DF["mass"] <- profile_header$star_mass
    obs.DF["He"] <- (profile_header$star_mass_he3 + 
        profile_header$star_mass_he4)/profile_header$star_mass
    
    obs.DF["radius"] <- profile_header$photosphere_r
    obs.DF["L"] <- profile_header$photosphere_L
    obs.DF["Teff"] <- profile_header$Teff
    
    acoustic_cutoff <- hstry$acoustic_cutoff/(2*pi)
    #obs.DF["delta_nu_scaling"] <- hstry$delta_nu
    
    obs.DF["log_g"] <- hstry$log_g
    #obs.DF["log_surf_z"] <- hstry$log_surf_z
    
    Z <- 10**hstry$log_surf_z
    H <- hstry$surface_h1
    obs.DF["Fe_H"] <- log10(Z/H/Z_div_X_solar)
    
    obs.DF["nu_max"] <- hstry$nu_max
    nu_max <- obs.DF["nu_max"]
    
    seis.DF <- seismology(freqs, nu_max, acoustic_cutoff)
    
    return(merge(rbind(obs.DF), rbind(seis.DF)))
}

### Obtain evolutionary tracks from a MESA directory
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
    ev_history <- read.table(file.path(log_dir, 'history.data'), 
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
    if (length(profile_files) <= 2) return(NA)
    
    # obtain observable information
    obs.DF <- do.call(plyr:::rbind.fill, 
        Map(function(profile_file, freqs_file)
                get_obs(profile_file, freqs_file, ev_history), 
            profile_file=file.path(log_dir, profile_files), 
            freqs_file=file.path(log_dir, freq_files)))
    
    return(merge(rbind(params.DF), obs.DF))
}





### Obtain grid of models 
fname <- 'sobol_grid-diffusion.dat'
DF <- if (file.exists(fname)) {
    read.table(fname, header=TRUE)
} else {
    parallelStartMulticore(max(1, detectCores()-1))
    DF <- do.call(rbind, parallelMap(parse_dir, 
        directory=list.dirs(exp_dir, recursive=FALSE)))
    DF <- DF[complete.cases(DF),]
    write.table(DF, fname, quote=FALSE, sep='\t', row.names=FALSE)
    write.table(t(sapply(DF, fivenum)), "fivenums.dat", quote=FALSE, sep='\t',
        row.names=TRUE, col.names=FALSE)

    if (FALSE) {
    ### inputs
    cairo_pdf(file.path(plot_dir, 'correlogram-sobol-inputs.pdf'),
        width=11.69, height=8.27, family=font)
    plot(DF[,1:5], lower.panel=NULL, pch=3, cex=0.5)
    dev.off()
    
    ### correlogram
    png(file.path(plot_dir, 'correlogram-sobol.png'),
        width=1169, height=827, family=font)
    plot(DF[,c(1:6, 8:10, 12:14)], lower.panel=NULL, pch=3, cex=0.1)
    dev.off()
    
    age_col <- which(names(DF)=="age")
    rbPal <- colorRampPalette(c('yellow','red'))
    for (method in c("pearson", "spearman")) {
        cairo_pdf(file.path(plot_dir, paste0('correlation-', method,'.pdf')),
            width=5, height=5, family=font)
        par(mar=c(5, 5, 0.1, 0.1))
        corrs <- cor(DF, DF$age, method=method)[-5,]
        sorted <- order(corrs)
        barplot(corrs[sorted], names.arg=latex_labs[-age_col][sorted], 
            xlim=c(0.3, -0.9),#rev(c(round(min(corrs)-0.05, 1), 
                 #  round(max(corrs)+0.05, 1))),
            horiz=1, las=1, #col=heat.colors(ncol(DF)-1),
            col=rbPal(ncol(DF)-1)[as.numeric(cut(abs(corrs[sorted]), 
                breaks=ncol(DF)-1))],
            xlab=if (method=="pearson") {
                expression("Pearson product-moment correlation coefficient"~r)
            } else if (method=="kendall") {
                expression("Kendall rank correlation coefficient"~tau)
            } else if (method=="spearman") {
                expression("Spearman rank correlation coefficient"~rho)
            })
        dev.off()
    }
    
    
    cairo_pdf(file.path(plot_dir, paste0('corr-boxplot-spearman.pdf')),
            width=5, height=5, family=font)
    par(mar=c(5, 5, 1, 1))
    spearmans <- do.call(cbind, Map(function(metal) 
        cor(DF[DF$Z==metal,-1:-5], DF[DF$Z==metal,]$age, method="s"),
        metal=unique(DF$Z)))
    medians <- apply(spearmans, 1, function(x) median(x, na.rm=1))
    sorted <- order(medians)
    boxplot(t(spearmans[sorted,]), pch=3, cex=0.1, 
        names=latex_labs[-1:-5][sorted],
        horizontal=1, las=1,
        col=adjustcolor(rbPal(ncol(DF)-5)[as.numeric(cut(abs(medians[sorted]), 
                breaks=ncol(DF)-5))], alpha.f=0.5),
        xlab=expression("Spearman rank correlation coefficient"~rho))
    dev.off()
    
    for (variable in c("M", "Y")) {
        quantiles <- quantile(unique(DF[[variable]]))
        for (quant in 1:4) {
            ss <- DF[DF[[variable]] <= quantiles[quant+1] &
                     DF[[variable]] >= quantiles[quant],]
            cairo_pdf(file.path(plot_dir, 
                    paste0('corrbplot-', variable, quant, '.pdf')),
                width=5, height=5, family=font)
            par(mar=c(5, 5, 1, 1))
            spearmans <- do.call(cbind, Map(function(metal) 
                cor(ss[ss$Z==metal,-1:-5], ss[ss$Z==metal,]$age, method="s"),
                metal=unique(ss$Z)))
            medians <- apply(spearmans, 1, function(x) median(x, na.rm=1))
            sorted <- order(names(ss)[-1:-5])
            boxplot(t(spearmans[sorted,]), pch=3, cex=0.1, 
                names=latex_labs[-1:-5][sorted],
                horizontal=1, las=1,
                col=adjustcolor(
                        rbPal(ncol(ss)-5)[as.numeric(cut(abs(medians[sorted]), 
                    breaks=ncol(ss)-5))], alpha.f=0.5),
                xlab=expression("Spearman rank correlation coefficient"~rho))
            dev.off()
        }
    }
    }
}





### Obtain properties of real stars varied within their uncertainties 
monte_carlo_perturbations <- function(obs_data_file, freqs_data_file,
        n_perturbations=10000) {
    freqs <- read.table(freqs_data_file, header=TRUE)
    obs_data <- read.table(obs_data_file, header=TRUE)
    noisy_freqs <- freqs
    do.call(rbind, with(obs_data, {replicate(n_perturbations, {
        obs.DF <- data.frame(
            radius     = rnorm(1, value[name=="radius"], 
                            uncertainty[name=="radius"]),
            L          = rnorm(1, value[name=="L"], 
                            uncertainty[name=="L"]),
            Teff       = rnorm(1, value[name=="Teff"], 
                            uncertainty[name=="Teff"]),
            log_g      = rnorm(1, value[name=="log_g"], 
                            uncertainty[name=="log_g"]),
            Fe_H       = rnorm(1, value[name=="Fe/H"], 
                            uncertainty[name=="Fe/H"]),
            nu_max     = rnorm(1, value[name=="nu_max"], 
                            uncertainty[name=="nu_max"])
        )
        noisy_freqs$nu <- rnorm(nrow(freqs), freqs$nu, freqs$dnu)
        seis.DF <- seismology(noisy_freqs, obs.DF$nu_max)
        merge(rbind(obs.DF), rbind(seis.DF))
    }, simplify=FALSE)}))
}

star_names <- c("16CygA", "16CygB", "Sun")
stars <- list()
for (star in star_names) {
    fname <- file.path('perturb', paste0(star, "_perturb.dat"))
    stars[[star]] <- if (file.exists(fname)) {
        read.table(fname, header=TRUE)
    } else {
        perturbations <- monte_carlo_perturbations(
            file.path("..", "data", paste0(star, "-obs.dat")),
            file.path("..", "data", paste0(star, "-freqs.dat")))
        write.table(perturbations, fname, quote=FALSE, 
                    sep='\t', row.names=FALSE)
    }
}





if(0) {
### Train a neural network and get predictions on real stars
make_predictions <- function(target) {
    use_log <- 0
    if (target == "M") {
        description <- "initial mass"
        symbol <- "M"
    } else if (target == "age") {
        description <- "age"
        symbol <- "t"
        use_log <- 1
    } else if (target == "Y") {
        description <- "initial helium"
        symbol <- "Y"
    } else if (target == "alpha") {
        description <- "mixing length"
        symbol <- "alpha"
    }
    
    fmla <- as.formula(paste(
        ifelse(use_log, paste("log10(", target, ")"), target), 
        "~ exp(log_g) + Fe_H + L + Teff +", 
        "radius + nu_max + dnu_0 + dnu_1 +",
        "Dnu_0 + Dnu_1 + Dnu_2 + Dnu_3 +",
        "eps_0 + eps_1 + eps_2 + eps_3"))
    lm. <- lm(fmla, data=DF)
    lm.resid <- resid(lm.)
    nn. <- nnet(fmla, data=DF, linout=TRUE, size=100, decay=0.1, 
                maxit=50000, MaxNWts=2000)
    nn.resid <- nn.$residuals
    
    save(nn., file=paste0('nn.', target))
    
    # lm resid
    cairo_pdf(file.path(plot_dir, paste0(target, '_lm_resid.pdf')),
        width=plot_width, height=plot_height, family=font)
    par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    plot(lm.resid ~ DF[[target]], pch=3, cex=0.25, tck=0, 
        ylim=c(min(nn.resid, lm.resid), max(nn.resid, lm.resid)),
        xlab=as.expression(bquote(.(description) ~ .(symbol))), 
        ylab=as.expression(bquote(hat(.(symbol)) - .(symbol))))
    magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
    abline(h=0, lty=2, col="red")
    abline(h=fivenum(lm.resid)[2], lty=3, col="red")
    abline(h=fivenum(lm.resid)[4], lty=3, col="red")
    dev.off()
    
    # nn resid
    cairo_pdf(file.path(plot_dir, paste0(target, '_nn_resid.pdf')),
        width=plot_width, height=plot_height, family=font)
    par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    plot(nn.resid ~ DF[[target]], pch=3, cex=0.25, tck=0, 
        ylim=c(min(nn.resid, lm.resid), max(nn.resid, lm.resid)),
        xlab=paste(description, symbol), 
        ylab=as.expression(bquote(hat(.(symbol)) - .(symbol))))
    magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
    abline(h=0, lty=2, col="red")
    abline(h=fivenum(nn.resid)[2], lty=3, col="red")
    abline(h=fivenum(nn.resid)[4], lty=3, col="red")
    dev.off()
    
    # box plots
    cairo_pdf(file.path(plot_dir, paste0(target, '_boxplots.pdf')),
        width=plot_width, height=plot_height, family=font)
    par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    boxplot(data.frame("Linear Model"=lm.resid, 
                       "Neural Network"=nn.resid), 
        ylab="Residuals",
        pch=3, cex=0.25)
    abline(h=0, lty=2, col="red")
    dev.off()
    
    # percent difference
    cairo_pdf(file.path(plot_dir, paste0(target, '_nn_percent_diff.pdf')),
        width=plot_width, height=plot_height, family=font)
    par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    fitted_ys <- nn.$fitted.values
    ys <- DF[[target]]
    if (use_log) {
        fitted_ys <- 10**fitted_ys
        ys <- log10(ys)
    }
    percent.diff <- (DF[[target]] - fitted_ys) / DF[[target]]*100
    plot(percent.diff ~ ys,
        pch=3, cex=0.25, tck=0,
        ylim=c(min(min(percent.diff), -100), max(max(percent.diff, 100))),
        xlab=ifelse(use_log, 
             as.expression(bquote(log[10]~.(description)~.(symbol))),
             as.expression(bquote(.(description)~.(symbol)))),
        ylab=as.expression(bquote("Percent Difference"~
            "["~(hat(.(symbol)) - .(symbol))/.(symbol)%*%100~"]")))
    magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
    abline(h=0, lty=2, col="red")
    abline(h=fivenum(percent.diff)[2], lty=3, col="red")
    abline(h=fivenum(percent.diff)[4], lty=3, col="red")
    dev.off()
    
    # predict cyg
    do.call(rbind, Map(function(star_name) {
        predictions <- adply(stars[[star_name]], 1, function(x) 
            predict(nn., x), .expand=FALSE)[,2]
        if (use_log) predictions <- 10**predictions
        
        pred.summary <- data.frame(row.names=star_name)
        pred.summary[[target]] <- mean(predictions)
        pred.summary[[paste0('d_', target)]] <- sd(predictions)
        
        cairo_pdf(file.path(plot_dir, paste0(target, '_', star_name, '.pdf')),
            width=plot_width, height=plot_height, family=font)
        par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
        hist(predictions, breaks=50)
        abline(v=mean(predictions), col='red')
        abline(v=mean(predictions)+sd(predictions), col='red', lty=2)
        abline(v=mean(predictions)-sd(predictions), col='red', lty=2)
        dev.off()
        pred.summary
    }, star_name=star_names))
}

parallelStartMulticore(4)
best_guesses <- do.call(cbind, 
    parallelMap(make_predictions, target=c('M', 'Y', 'alpha', 'age')))

stargazer(t(best_guesses), type='text', summary=FALSE)









### helium formula
He.fmla <- Y ~ exp(log_g) + Fe_H + L + Teff + nu_max + 
    Dnu_0 + Dnu_1 + Dnu_2 + Dnu_3 + dnu_0 + dnu_1 +
    eps_0 + eps_1 + eps_2 + eps_3
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
    ylab=expression(hat(Y) - Y))
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
    ylab=expression(hat(Y) - Y))
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
age.fmla <- log10(age) ~ exp(log_g) + Fe_H + L + Teff + nu_max + 
    Dnu_0 + Dnu_1 + Dnu_2 + Dnu_3 + dnu_0 + dnu_1 +
    eps_0 + eps_1 + eps_2 + eps_3
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

#}

age.fmla2 <- log10(age) ~ radius + exp(log_g) + Fe_H + L + Teff + nu_max + 
    Dnu_0 + Dnu_1 + Dnu_2 + Dnu_3 + dnu_0 + dnu_1 +
    eps_0 + eps_1 + eps_2 + eps_3
age.nn2 <- nnet(age.fmla2, data=DF, linout=TRUE, size=100, decay=0.1, 
    maxit=100000, MaxNWts=2000)
nn <- age.nn2
### parse 16 Cyg data
for (cyg in c("16CygA", "16CygB")) {
    obs_data <- read.table(file.path("..", "data", paste0(cyg, "-obs.dat")), 
        header=TRUE)
    attach(obs_data)
    freqs <- read.table(file.path("..", "data", paste0(cyg, "-freqs.dat")), 
        header=TRUE)
    
    predictions <- c()
    for (i in 1:5000) {
        obs.DF <- data.frame(
            nu_max     = rnorm(1, value[name=="nu_max"],
                            uncertainty[name=="nu_max"]),
            log_g      = rnorm(1, value[name=="log_g"], 
                            uncertainty[name=="log_g"]),
            Fe_H       = rnorm(1, value[name=="Fe/H"], 
                            uncertainty[name=="Fe/H"]),
            L          = rnorm(1, value[name=="L"], 
                            uncertainty[name=="L"]),
            Teff       = rnorm(1, value[name=="Teff"], 
                            uncertainty[name=="Teff"]),
            radius     = rnorm(1, value[name=="radius"], 
                            uncertainty[name=="radius"])
        )
        
        noisy_freqs <- freqs
        noisy_freqs$nu <- rnorm(nrow(freqs), freqs$nu, freqs$dnu)
        
        seis.DF <- seismology(noisy_freqs, obs.DF$nu_max)
        
        predictions <- c(predictions, 
            10**predict(nn, merge(rbind(obs.DF), rbind(seis.DF))))
    }
    
    cairo_pdf(file.path(plot_dir, paste0('age_', cyg, '.pdf')),
        width=plot_width, height=plot_height, family=font)
    par(mar=c(4, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    hist(predictions, breaks=50, xlab="Age [ years ]", main="")
    abline(v=mean(predictions), col='red')
    abline(v=mean(predictions)+sd(predictions), col='red', lty=2)
    abline(v=mean(predictions)-sd(predictions), col='red', lty=2)
    text(#diff(par("usr")[1:2])/2, diff(par("usr")[3:4])/2,
        #4100000000, 150,
        mean(predictions)+2*sd(predictions), 100,
        paste("Age:", formatC(mean(predictions)/1e9), 
              "\n           +/-", formatC(sd(predictions)/1e9), "Gyr"))
    dev.off()
    
    detach(obs_data)
}

}
