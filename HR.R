#### Plot HR diagrams for simulations in all subdirectories 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(Hmisc)

Teff_sun = log10(5777)
mar <- c(3, 3.5, 3, 1)
mgp <- c(2, 0.25, 0)
approx_n <- 80
approx_xout <- seq(0.01, 0.8, .01)
font <- "Palatino"

dirs <- dir('exp')
offset <- 4
cl <- rev(heat.colors(length(dirs)+offset))[offset:(offset+length(dirs)-1)]

### Plot HR diagram
cairo_pdf(file.path('plots', 'HR.pdf'), width=6, height=6, family=font)
counter <- 0
model_nos <- c()
for (simulation in dirs) {
    data_file <- file.path('exp', simulation, 'LOGS', 'history.data')
    if (file.exists(data_file)) { 
        counter <- counter + 1
        data <- read.table(data_file, header=TRUE, skip=5)
        HR <- data$log_L ~ data$log_Teff
        if (counter == 1) {
            par(bty="l", las=1, mar=mar, mgp=mgp)
            plot(HR, type='l', col=cl[counter],
                 xlab=expression(log~T[eff]),
                 ylab=expression(log~L/L["⊙"]), 
                 main="HR diagram of Sun-like stars with varied mixing lengths",
                 xlim=rev(c(min(data$log_Teff), 
                            max(data$log_Teff)+0.05)),
                 ylim=c(floor(min(data$log_L)),
                        ceiling(max(data$log_L))),
                 xaxs='i', yaxs='i', tck=0.01)
            minor.tick(nx=5, ny=5, tick.ratio=-0.25)
        } else {
            lines(HR, col=cl[counter])
        }
        sun_like <- Teff_sun < data$log_Teff+.005 & 0 < data$log_L+.05 & 
                    Teff_sun > data$log_Teff-.005 & 0 > data$log_L-.05
        model_nos <- c(model_nos, data$model_number[sun_like])
    }
}
points(Teff_sun, 0, col="black")
points(Teff_sun, 0, col="black", pch=21, cex=0.1)
labels <- sapply(sub('.+_', '', dirs), 
                 function(x) { as.expression(bquote(alpha ~ "=" ~ .(x))) })
legend("topleft", col=cl, lty=1, bty='n', legend=labels)
dev.off()

### Plot sound speed and density profiles
modelS <- read.table('modelS.dat', header=TRUE) 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
model_no <- Mode(model_nos) # find nearest solar-like model
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
    interps <- matrix(nrow=approx_n, ncol=length(dirs)) # Preallocate
    ## Find y limits 
    y_min <- min(modelS_y)
    y_max <- max(modelS_y)
    for (simulation in dirs) {
        profile_index <- file.path('exp', simulation, 'LOGS', 'profiles.index')
        profiles <- read.table(profile_index, skip=1)
        profile_no <- profiles[which(profiles[,1]==model_no), 3]
        data_file <- file.path('exp', simulation, 'LOGS', 
                               paste0('profile', profile_no, '.data'))
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
        }
    }
    
    ## Make plots 
    cairo_pdf(file.path('plots', paste0('profile_', col_name, '.pdf')), 
              width=6, height=6, family="Palatino")
    counter <- 0
    for (simulation in dirs) {
        profile_index <- file.path('exp', simulation, 'LOGS', 'profiles.index')
        profiles <- read.table(profile_index, skip=1)
        profile_no <- profiles[which(profiles[,1]==model_no),3]
        data_file <- file.path('exp', simulation, 'LOGS', 
                               paste0('profile', profile_no, '.data'))
        if (file.exists(data_file)) {
            counter <- counter + 1
            data <- read.table(data_file, header=TRUE, skip=5)
            y <- data[[col_name]]
            if (col_name=='csound' || 
                col_name=='pressure' ||
                col_name=='temperature') {
                y <- log10(y)
            }
            interps[,counter] <- approx(data$radius, y, 
                                        n=approx_n, xout=approx_xout)$y
            if (counter == 1) {
                par(bty="l", las=1, mar=mar, mgp=c(2, 0.25, 0))
                plot(data$radius, y, type='l', col=cl[counter],
                     xlim=c(0, round(max(data$radius)+.05, 1)),
                     ylim=c(y_min, y_max+0.01*y_max), 
                     xaxs='i', yaxs='i', tck=0.01,
                     xlab=expression(r/R["⊙"]), 
                     ylab=ylabel,
                     main=paste(main_label, 
                                "profile as a function of mixing length"))
                minor.tick(nx=5, ny=5, tick.ratio=-0.25)
            } else {
                lines(data$radius, y, col=cl[counter])
            }
        }
    }
    ## plot modelS
    lines(modelS$radius, modelS_y, col='black')
    legend("bottomleft", lty=1, bty='n',
           legend=c(labels, "Model S"), col=c(cl, "black"))
    dev.off()
    
    cairo_pdf(file.path('plots', paste0('profile_diff_', col_name, '.pdf')), 
              width=6, height=6, family=font)
    S_interp <- approx(modelS$radius, modelS_y, n=approx_n, xout=approx_xout)$y
    rel_diff <- apply(interps, 2, 
                      function(x) {abs((x-S_interp)/max(x, S_interp))})
    for (col in 1:ncol(interps)) {
        if (col == 1) {
            par(bty="l", las=1, mar=mar, mgp=mgp)
            plot(approx_xout, rel_diff[,col], type='l', col=cl[col],
                 xaxs='i', yaxs='i', tck=0.01,
                 ylim=c(0, max(rel_diff)+0.01*max(rel_diff)),
                 ylab=bquote(abs((.(symbol)[alpha] - .(symbol)[S])
                             /max*"("*.(symbol)[S]*","~.(symbol)[alpha]*")")),
                 xlab=expression(r/R["⊙"]),
                 main=paste("Relative", tolower(main_label),
                            "difference with Model S"))
            minor.tick(nx=5, ny=5, tick.ratio=-0.25)
        } else {
            lines(approx_xout, rel_diff[,col], col=cl[col])
        }
    }
    legend("top", legend=labels, col=cl, lty=1, bty='n')
    dev.off()
}
