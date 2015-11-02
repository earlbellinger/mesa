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

exp_dir <- file.path('diffusion_off')
plot_dir <- file.path('plots')
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

#dd_01= 1/8( nu_[n-1,0] - 4*nu_[n-1,1] + 6*nu_[n,0] - 4*nu[n,  1] + nu_[n+1,0] )
#dd_10=-1/8( nu_[n-1,1] - 4*nu_[n,  0] + 6*nu_[n,1] - 4*nu[n+1,0] + nu_[n+1,1] )
dd <- function(l0, l1, n, df) {
    ell.0 <- df[df$l==0 & df$n>0,]
    ell.1 <- df[df$l==1 & df$n>0,]
    n. <- df[df$n==n,]
    n.minus.one <- df[df$n==n-1,]
    n.plus.one <- df[df$n==n+1,]
    val <- if (l0 == 0 && l1 == 1) {
        ( merge(n.minus.one, ell.0)$nu -
        4*merge(n.minus.one, ell.1)$nu +
        6*merge(n., ell.0)$nu -
        4*merge(n., ell.1)$nu +
          merge(n.plus.one, ell.0)$nu )/8
    } else if (l1 == 0 && l0 == 1) {
        -( merge(n.minus.one, ell.1)$nu -
         4*merge(n., ell.0)$nu +
         6*merge(n., ell.1)$nu -
         4*merge(n.plus.one, ell.0)$nu +
           merge(n.plus.one, ell.1)$nu )/8
    } else NA
    if (length(val) == 0) NA
    else val
}

dnu <- function(l, n, df) separation(l, n, l+2, n-1, df)
Dnu <- function(l, n, df) separation(l, n, l, n-1, df)
r02 <- function(n, df) dnu(0, n, df) / Dnu(1, n, df)
r13 <- function(n, df) dnu(1, n, df) / Dnu(0, n+1, df)
r01 <- function(n, df) dd(0, 1, n, df) / Dnu(1, n, df)
r10 <- function(n, df) dd(1, 0, n, df) / Dnu(0, n+1, df)

get_averages <- function(f, df, freqs, l_degs, plot=FALSE) {
  sep_name <- deparse(substitute(f))
  a <- c()
  b <- c()
  pchs <- c()
  for (l_deg in l_degs) {
      ell <- freqs[freqs$n > 1 & freqs$l==l_deg,]
      if (sep_name == 'dnu' || sep_name == 'Dnu') {
        vals <- sapply(unique(ell$n), function(n) f(l_deg, n, freqs))
      } else {
        vals <- sapply(unique(ell$n), function(n) f(n, freqs))
      }
      not.nan <- complete.cases(vals)
      a <- c(a, vals[not.nan])
      b <- c(b, ell$n[not.nan])
      pchs = c(pchs, rep(l_deg+1, length(a)))
  }
  
  ylab <- if (sep_name=='Dnu' && length(l_degs) == 1) 
                               bquote(Delta*nu[.(l_degs)]~"["*mu*Hz*"]")
     else if (sep_name=='Dnu' && length(l_degs) >= 1)
                               bquote(Delta*nu~"["*mu*Hz*"]")
     else if (sep_name=='dnu') bquote(delta*nu[.(l_degs)*','
                                              *.(l_degs+2)]~"["*mu*Hz*"]")
     else if (sep_name=='r02') bquote(r[0*","*2]~"["*mu*Hz*"]")
     else if (sep_name=='r13') bquote(r[1*","*3]~"["*mu*Hz*"]")
     else if (sep_name=='r01') bquote(r[0*","*1]~"["*mu*Hz*"]")
     else if (sep_name=='r10') bquote(r[1*","*0]~"["*mu*Hz*"]")
  
  if (sep_name == 'Dnu' && length(l_degs) == 1)
      sep_name <- paste0(sep_name, l_deg)
  if (sep_name == 'dnu')
      sep_name <- paste0(sep_name, l_deg, l_deg+2)
  
  df[paste0(sep_name, "_median")] <- median(a)
  fit <- mblm(a~b, repeated=1)
  df[paste0(sep_name, "_slope")] <- coef(fit)[2]
  df[paste0(sep_name, "_intercept")] <- coef(fit)[1]
  
  if (plot != FALSE) {
    cairo_pdf(file.path('separation_plots-known', 
            paste0(sep_name, '-', plot, '.pdf')),
        width=4, height=2.5, family=font)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    plot(a~b, tck=0, ylab=ylab, 
         pch=if (length(l_degs)==1) 3 else pchs, 
         xlab=expression("radial order"~n))
    abline(fit, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
    if (length(l_degs)>1)
        legend("topright", pch=l_degs+1, legend=paste0("l=", l_degs), bty="n")
    dev.off()
  }
  
  df
}

seismology <- function(freqs, acoustic_cutoff=Inf, plot=FALSE) {
  if (nrow(freqs) == 0) return(NULL)
  
  seis.DF <- NULL
  
  #converted_fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
  
  freqs <- unique(freqs[freqs$nu < acoustic_cutoff,])
  for (l_deg in 0:3) # sometimes ADIPLS repeats 'n's
    if (any(duplicated(freqs[freqs$l==l_deg,]$n)))
      return(NULL) # just discard this data
  
  # calculate frequency separations
  #ells <- sort(unique(freqs$l))
  #for (l_deg in ells) {
    ##ell <- freqs[freqs$n > 0 & freqs$l==l_deg,] # pressure modes
    ##fit <- mblm(nu ~ n, dataframe=ell, repeated=TRUE)
    ##seis.DF[paste0("Dnu_", l_deg)] <- coef(fit)[2]
    ##seis.DF[paste0("eps_", l_deg)] <- coef(fit)[1]/coef(fit)[2]
    #seis.DF <- get_averages(Dnu, seis.DF, freqs, l_deg, plot)
    #if (any(freqs$l==l_deg+2))
    #  seis.DF <- get_averages(dnu, seis.DF, freqs, l_deg, plot)
  #}
  seis.DF <- get_averages(Dnu, seis.DF, freqs, 0, plot=plot)
  seis.DF <- get_averages(Dnu, seis.DF, freqs, 0:3, plot=plot)
  seis.DF <- get_averages(dnu, seis.DF, freqs, 0, plot=plot)
  seis.DF <- get_averages(dnu, seis.DF, freqs, 1, plot=plot)
  seis.DF <- get_averages(r02, seis.DF, freqs, 0, plot=plot)
  seis.DF <- get_averages(r13, seis.DF, freqs, 0, plot=plot)
  seis.DF <- get_averages(r10, seis.DF, freqs, 0, plot=plot)
  seis.DF <- get_averages(r01, seis.DF, freqs, 0, plot=plot)
  return(seis.DF)
}

### Obtain observable properties from models 
get_obs <- function(profile_file, freqs_file, ev_history) {
    print(profile_file)
    
    profile_header <- read.table(profile_file, header=TRUE, nrows=1, skip=1)
    freqs <- read.table(freqs_file, col.names=freqs_cols, fill=TRUE)
    freqs <- freqs[complete.cases(freqs),]
    
    # discard models with convective cores
    #prev_hstry <- ev_history$model_number<=profile_header$model_number
    #conv_cores <- ev_history[prev_hstry,]$mass_conv_core
    #if (any(conv_cores > 0)) return(NA)
    
    hstry <- ev_history[ev_history$model_number==profile_header$model_number,]
    if (nrow(hstry) == 0 || hstry$mass_conv_core > 0) return(NULL)
    
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
    
    seis.DF <- seismology(freqs, acoustic_cutoff,
        plot=ifelse(sample(0:5000, 1)==0, gsub("/", "-", freqs_file), FALSE))
    
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
fname <- 'grids/sobol_grid-no_diffusion-eddington-no_ccore-ratios.dat'
DF <- if (file.exists(fname)) {
    read.table(fname, header=TRUE)
} else {
    parallelStartMulticore(max(1, detectCores()))
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
