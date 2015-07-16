#### Plot HR diagrams for simulations in all subdirectories 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(Hmisc)

#simulation_dir <- 'exp_alpha'
#simulation_dir <- 'exp_Y'
for (simulation_dir in c('exp_alpha', 'exp_Y')) {

dirs <- dir(simulation_dir)
plot_dir <- paste0('plots_', simulation_dir)
dir.create(plot_dir, showWarnings=FALSE)

if (simulation_dir == 'exp_Y') {
    exp <- 'helium content'
    sun_num <- 9
    per_off <- 0.3
} else if (simulation_dir == 'exp_alpha') {
    exp <- 'mixing lengths'
    sun_num <- 7
    per_off <- 0.1
}

Teff_sun = log10(5777)
solar_age = 4.57e9

mar <- c(4, 4.5, 3, 1)
mgp <- c(3, 0.25, 0)
approx_xout <- seq(0, 1, .04)
font <- "Palatino"

freep <- as.name(sub('.+_', '', simulation_dir))
offset <- 6
cl <- heat.colors(length(dirs)+offset)[1:length(dirs)]
labels <- sapply(sub('.+_', '', dirs), 
    function(x) { as.expression(bquote(.(freep) ~ "=" ~ .(x))) })

### Plot HR diagram
print("Plotting HR diagram")
cairo_pdf(file.path(plot_dir, 'HR.pdf'), width=6, height=6, family=font)
sim_i <- 0
model_nos <- c()
for (simulation in dirs) {
    data_file <- file.path(simulation_dir, simulation, 'LOGS', 'history.data')
    if (file.exists(data_file)) { 
        sim_i <- sim_i + 1
        data <- read.table(data_file, header=TRUE, skip=5)
        HR <- data$log_L ~ data$log_Teff
        if (sim_i == 1) {
            par(bty="l", las=1, mar=mar, mgp=mgp)
            plot(HR, type='l', col=cl[sim_i],
                 xlab=expression(log~T[eff]),
                 ylab=expression(log~L/L["⊙"]), 
                 main=paste("HR diagram of Sun-like stars with varied", exp),
                 xlim=rev(c(min(data$log_Teff)-0.01, 
                            max(data$log_Teff)+0.05)),
                 ylim=c(floor(min(data$log_L)),
                        ceiling(max(data$log_L))),
                 xaxs='i', yaxs='i', tck=0.01)
            minor.tick(nx=5, ny=5, tick.ratio=-0.25)
        } else {
            lines(HR, col=cl[sim_i])
        }
    }
}
points(Teff_sun, 0, col="black")
points(Teff_sun, 0, col="black", pch=21, cex=0.1)
legend("topleft", col=cl, lty=1, bty='n', legend=labels)
dev.off()

bcl <- c(cl[1:(sun_num-1)], "black", cl[(sun_num+1):length(cl)])

### Find profile files at the solar age
print("Finding profile files at the solar age")
modelS <- read.table('modelS.dat', header=TRUE) 
profile_nos <- c()
for (simulation in dirs) {
    log_dir <- file.path(simulation_dir, simulation, 'LOGS')
    log_files <- dir(log_dir)
    log_files <- log_files[order(as.numeric(gsub('\\D', '', log_files)))]
    best_log <- NA
    best_age <- Inf
    for (log_i in grep('profile\\d+.data$', log_files)) {
        log_file <- file.path(log_dir, log_files[log_i])
        age <- read.table(log_file, skip=1, nrows=1, header=TRUE)$star_age
        age_diff <- abs(age - solar_age)
        print(c(log_file, best_age, age_diff))
        if (age_diff < best_age) {
            best_log <- log_file
            best_age <- age_diff
        }
        if (age_diff > best_age) break
    }
    profile_nos <- c(profile_nos, best_log)
}

### Plot internal profiles
print("Plotting internal profiles")
for (col_name in c("csound", "logRho", "pressure", "temperature")) {
    if (col_name=="csound") {
        modelS_y <- log10(modelS[['csound']])
        main_label <- "Sound speed"
        symbol <- bquote(c)
        ylabel <- expression(log~sound~speed~c~"["*cm/sec*"]")
    } else if (col_name=="logRho") {
        modelS_y <- log10(modelS[['rho']])
        main_label <- "Density"
        symbol <- bquote(rho)
        ylabel <- expression(log~density~rho~"["*g/cm^3*"]")
    } else if (col_name=="pressure") {
        modelS_y <- log10(modelS[['pressure']])
        main_label <- "Pressure"
        symbol <- bquote(p)
        ylabel <- expression(log~pressure~p~"["*dyn/cm^2*"]")
    } else if (col_name=="temperature") {
        modelS_y <- log10(modelS[['temperature']])
        main_label <- "Temperature"
        symbol <- bquote(T)
        ylabel <- expression(log~temperature~T~"["*K*"]")
    }
    interps <- matrix(nrow=length(approx_xout), ncol=length(dirs)) # Preallocate
    ## Find plot limits 
    y_min <- min(modelS_y)
    y_max <- max(modelS_y)
    x_max <- 1
    for (simulation in dirs) {
        data_file <- file.path(profile_nos[grep(simulation, profile_nos)])
        if (file.exists(data_file)) {
            data <- read.table(data_file, header=TRUE, skip=5)
            y <- data[[col_name]]
            if (col_name=='csound' || 
                col_name=='pressure' ||
                col_name=='temperature') {
                y <- log10(y)
            }
            if (min(y) < y_min) {
                y_min <- min(y)
            }
            if (max(y) > y_max) {
                y_max <- max(y)
            }
            if (round(max(data$radius)+.05, 1) > x_max) {
                x_max <- round(max(data$radius)+.05, 1)
            }
        }
    }
    
    ## Make profile plots 
    cairo_pdf(file.path(plot_dir, paste0('profile_', col_name, '.pdf')), 
              width=6, height=6, family="Palatino")
    sim_i <- 0
    for (simulation in dirs) {
        data_file <- file.path(profile_nos[grep(simulation, profile_nos)])
        if (file.exists(data_file)) {
            sim_i <- sim_i + 1
            data <- read.table(data_file, header=TRUE, skip=5)
            y <- data[[col_name]]
            if (col_name=='csound' || 
                col_name=='pressure' ||
                col_name=='temperature') {
                y <- log10(y)
            }
            interps[,sim_i] <- approx(data$radius, y, n=length(approx_xout), 
                                        xout=approx_xout, rule=2)$y
            if (sim_i == 1) {
                par(bty="l", las=1, mar=mar, mgp=mgp)
                plot(data$radius, y, type='l', col=cl[sim_i],
                     xlim=c(0, x_max),
                     ylim=c(y_min, y_max+0.01*y_max), 
                     xaxs='i', yaxs='i', tck=0.01,
                     xlab=expression(r/R["⊙"]), 
                     ylab=ylabel,
                     main=paste(main_label, "profile as a function of\n", 
                                exp, "at the solar age"))
                minor.tick(nx=5, ny=5, tick.ratio=-0.25)
            } else {
                lines(data$radius, y, col=cl[sim_i])
            }
        }
    }
    ## plot modelS
    lines(modelS$radius, modelS_y, col='black')
    legend("bottomleft", lty=1, bty='n',
           legend=c(labels, "Model S"), col=c(cl, "black"))
    dev.off()
    
    ## plot relative differences
    cairo_pdf(file.path(plot_dir, paste0('diff_', col_name, '.pdf')), 
              width=6, height=6, family=font)
    S_interp <- approx(modelS$radius, modelS_y, n=length(approx_xout), 
                       xout=approx_xout, rule=2)$y
    rel_diff <- apply(interps, 2, 
    #                  function(x) {abs((x-S_interp)/max(x, S_interp))})
                      function(x) {(x-S_interp)/S_interp})
    #core_max <- max(rel_diff[1,])
    #r_max <- min(which(rel_diff[2:nrow(rel_diff),]>core_max, 
    #                   arr.ind=TRUE)[,1])+1
    core_var <- var(rel_diff[1,])
    offset <- round(nrow(rel_diff) * 0.45, 0)
    r_max <- min(which(apply(rel_diff[offset:nrow(rel_diff),], 1, 
                             var)>core_var))+offset
    x_max <- max(approx_xout[1:r_max])
    for (sim_i in 1:ncol(rel_diff)) {
        if (sim_i == 1) {
            par(bty='l', las=1, mar=mar+c(0,0,0,6), mgp=mgp, xpd=TRUE)
            plot(approx_xout[1:r_max], rel_diff[1:r_max, sim_i], 
                 type='l', col=cl[sim_i],
                 xaxs='i', yaxs='i', 
                 tck=0.01,
                 xlim=c(0, x_max),
                 ylim=range(rel_diff[1:r_max,]),
                 #ylab=bquote(abs((.(symbol)[.(freep)] - .(symbol)[S])
                 #    /max*"("*.(symbol)[S]*","~.(symbol)[.(freep)]*")")),
                 ylab=bquote((log~.(symbol)[.(freep)] - log~.(symbol)[S]) 
                            / log~.(symbol)[S]),
                 xlab=expression(r/R["⊙"]),
                 main=paste("Relative", tolower(main_label),
                            "difference with\nModel S at the solar age"))
            minor.tick(nx=5, ny=5, tick.ratio=-0.25)
            lines(c(x_max,0), c(0, 0), lty=2)
        } else {
            lines(approx_xout[1:r_max], rel_diff[1:r_max, sim_i], 
                  col=bcl[sim_i])
        }
    }
    legend("right", legend=labels, col=bcl, lty=1, bty='n', inset=c(-0.35,0))
    dev.off()
}

}
