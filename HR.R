#### Plot HR diagrams for simulations in all subdirectories 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(Hmisc)

minor.tick <- function (nx = 2, ny = 2, tick.ratio = 0.5) {
    ax <- function(w, n, tick.ratio) {
        range <- par("usr")[if (w == "x") 
            1:2
        else 3:4]
        tick.pos <- if (w == "x") 
            par("xaxp")
        else par("yaxp")
        distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
        possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
        low.candidates <- possible.minors >= range[1]
        low.minor <- if (any(low.candidates)) {
            min(possible.minors[low.candidates])
        } else {
            tick.pos[1]
        }
        possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
        hi.candidates <- possible.minors <= range[2]
        hi.minor <- if (any(hi.candidates)) {
            max(possible.minors[hi.candidates])
        } else {
            tick.pos[2]
        }
        axis(if (w == "x") 
            1
        else 2, seq(low.minor, hi.minor, by = distance.between.minor), 
            labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1) 
        ax("x", nx, tick.ratio = tick.ratio)
    if (ny > 1) 
        ax("y", ny, tick.ratio = tick.ratio)
    invisible()
}

Teff_sun = log10(5777)
solar_age = 4.57e9
font <- "Palatino"
approx_xout <- seq(0, 1, .04)
color_offset <- 6

plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

make_legend <- function(labels, ev_stage, position='left') {
    par(mar=rep(0,4))
    plot.new()
    legend(position, bty='n', inset=0, #pt.cex=1, cex=1, 
        lty=if(ev_stage=="solar-like") c(rep(1, length(labels)), 2) else 1, 
        col=if(ev_stage=="solar-like") c(cl, "black") else cl,
        legend=if(ev_stage=="solar-like") c(labels, "Model S") else labels)
}

start_dev <- function(maintext, fname) {
    cairo_pdf(file.path(subdir, 
        paste0(ev_stage, '-', simulation_dir, '_', fname, '.pdf')),
              width=6, height=6, family=font)
    layout(matrix(c(1,1,1,2,3,4,5,6,4), ncol=3, byrow=TRUE), 
           heights=c(0.1,0.45,0.45), widths=c(.425,.425,.15))
    par(mar=rep(0,4))
    plot.new()
    text(0.5, 0.5, 
         paste(maintext, ev_stage, "stars by", exp),
         cex=2, font=2)
}

for (simulation_dir in c('exp_alpha', 'exp_Y')) {
    dirs <- dir(simulation_dir)
    if (simulation_dir == 'exp_Y') {
        exp <- 'helium'
        sun_num <- grep('.28', dirs)
    } else if (simulation_dir == 'exp_alpha') {
        exp <- 'mixing length'
        sun_num <- grep('2.10', dirs)
    }
    
    freep <- as.name(sub('.+_', '', simulation_dir))
    cl <- heat.colors(length(dirs)+color_offset)[1:length(dirs)]
    cl <- c(cl[1:(sun_num-1)], "black", cl[(sun_num+1):length(cl)])
    labels <- sapply(sub('.+_', '', dirs), 
        function(x) { as.expression(bquote(.(freep) ~ "=" ~ .(x))) })
    
    ###########################################################
    ### Find profile files at different stages of evolution ###
    ###########################################################
    print("Finding profile files")
    modelS <- read.table('modelS.dat', header=TRUE) 
    saved_pro_file <- paste0(".", simulation_dir, "_profile_nos")
    if (file.exists(saved_pro_file)) {
        load(saved_pro_file)
    } else {
        profile_nos <- list()
        ## Find profiles at the pre-main sequence and solar age
        for (simulation in dirs) {
            log_dir <- file.path(simulation_dir, simulation, 'LOGS')
            log_files <- dir(log_dir)
            log_nos <- as.numeric(gsub('\\D', '', log_files))
            log_files <- log_files[order(log_nos)]
            
            best_PMS_log <- NA
            best_solar_log <- NA
            best_bump_log <- NA
            
            #found_sun <- FALSE
            lum_val <- -Inf
            
            min_lum <- Inf
            best_age_diff <- Inf
            
            for (log_i in grep('profile\\d+.data$', log_files)) {
                log_file <- file.path(log_dir, log_files[log_i])
                profile <- read.table(log_file, skip=1, nrows=1, header=TRUE)
                
                lum <- profile$photosphere_L
                if (lum < min_lum) {
                    best_PMS_log <- log_file
                    min_lum <- lum
                }
                
                age <- profile$star_age
                age_diff <- abs(age - solar_age)
                if (age_diff < best_age_diff) {
                    best_solar_log <- log_file
                    best_age_diff <- age_diff
                }
                if (age_diff > best_age_diff) found_sun <- TRUE
                
                if (as.numeric(gsub('\\D', '', log_files[log_i])) > 1000) {
                    if (lum > lum_val) {
                        lum_val <- lum
                    } else {
                        best_bump_log <- log_file
                        break
                    }
                }
            }
            
            profile_nos[["PMS"]] <- c(profile_nos[["PMS"]], best_PMS_log)
            profile_nos[["solar-like"]] <- c(profile_nos[["solar-like"]], 
                best_solar_log)
            profile_nos[["RGB bump"]] <- c(profile_nos[["RGB bump"]], 
                best_bump_log)
        }
        save(profile_nos, file=saved_pro_file)
    }
    print(profile_nos)
        
    ##############################
    ### Plot internal profiles ###
    ##############################
    print("Plotting internal profiles")
    for (ev_stage in names(profile_nos)) {
        subdir <- file.path(plot_dir, ev_stage)
        dir.create(subdir, showWarnings=FALSE)
        
        ## Start profile difference plot device
        start_dev("Differences in", "diff")
        start_dev("Profiles of", "profile")
        
        ## Iterate through the different internal profiles
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
            # Preallocate interpolation matrix
            interps <- matrix(nrow=length(approx_xout), ncol=length(dirs))
            ## Find plot limits 
            y_min <- ifelse(ev_stage=="solar-like", min(modelS_y), Inf)
            y_max <- ifelse(ev_stage=="solar-like", max(modelS_y), -Inf)
            x_max <- 1
            for (simulation in dirs) {
                sim_no <- grep(simulation, profile_nos[[ev_stage]])
                data_file <- file.path(profile_nos[[ev_stage]][sim_no])
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
            sim_i <- 0
            for (simulation in dirs) {
                sim_no <- grep(simulation, profile_nos[[ev_stage]])
                data_file <- file.path(profile_nos[[ev_stage]][sim_no])
                if (file.exists(data_file)) {
                    sim_i <- sim_i + 1
                    data <- read.table(data_file, header=TRUE, skip=5)
                    y <- data[[col_name]]
                    if (col_name=='csound' || 
                        col_name=='pressure' ||
                        col_name=='temperature') {
                        y <- log10(y)
                    }
                    interps[,sim_i] <- approx(data$radius, y, 
                        n=length(approx_xout), xout=approx_xout, rule=2)$y
                    if (sim_i == 1) {
                        par(bty="l", las=1, 
                            mar=c(3, 3.2, 1, 1), mgp=c(1.8, 0.25, 0))
                        plot(data$radius, y, type='l', col=cl[sim_i],
                             xlim=c(0, x_max),
                             ylim=c(y_min, y_max+0.01*y_max), 
                             xaxs='i', yaxs='i', tck=0.01,
                             #xlab=expression(r/R["⊙"]), 
                             xlab="",
                             ylab=ylabel)
                        title(xlab=expression(r/R["⊙"]))
                        minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                    } else {
                        lines(data$radius, y, col=cl[sim_i])
                    }
                }
            }
            if (ev_stage=="solar-like") ## plot Model S
                lines(modelS$radius, modelS_y, col='black', lty=2)
            if (col_name=="logRho") make_legend(labels, ev_stage)
            dev.set(dev.prev())
            
            ## plot relative differences
            
            ref <- if(ev_stage=="solar-like") {
                approx(modelS$radius, modelS_y, n=length(approx_xout), 
                       xout=approx_xout, rule=2)$y
            } else {
                interps[,sun_num]
            }
            rel_diff <- apply(interps, 2, function(x) {(x-ref)/ref})
            row_offset <- round(nrow(rel_diff) * 0.5, 0)
            #core_var <- var(rel_diff[1,])
            core_var <- max(apply(rel_diff[1:(row_offset-1),],1,var))
            endpts <- apply(rel_diff[row_offset:nrow(rel_diff),],1,var)>core_var
            row_max <- ifelse(any(endpts), min(which(endpts))+row_offset, 
                nrow(rel_diff))
            x_max <- max(approx_xout[1:row_max])
            for (sim_i in 1:ncol(rel_diff)) {
                if (sim_i == 1) {
                    par(bty='l', las=1, xpd=TRUE, 
                        mar=c(3, 4.5, 1, 1), mgp=c(3, 0.25, 0))
                    plot(approx_xout[1:row_max], rel_diff[1:row_max, sim_i], 
                         type='l', col=cl[sim_i],
                         xaxs='i', yaxs='i', 
                         tck=0.01,
                         xlim=c(0, x_max),
                         ylim=range(rel_diff[1:row_max,]),
                         ylab=bquote(delta * .(symbol) / .(symbol)),
                         xlab="")
                    title(xlab=expression(r/R["⊙"]), mgp=par()$mgp-c(1.3,0,0))
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                    lines(c(x_max,0), c(0, 0), lty=2)
                } else {
                    lines(approx_xout[1:row_max], rel_diff[1:row_max, sim_i], 
                          col=cl[sim_i])
                }
            }
            if (col_name=="logRho") make_legend(labels, ev_stage)
            dev.set(dev.next())
        }
        dev.off()
        dev.off()
    }
    
    #######################
    ### Plot HR diagram ###
    #######################
    print("Plotting HR diagram")
    cairo_pdf(file.path(plot_dir, paste0(simulation_dir, '_HR.pdf')), 
              width=6, height=6, family=font)
    sim_i <- 0
    for (simulation in dirs) {
        data_file <- file.path(simulation_dir, simulation, 
            'LOGS', 'history.data')
        if (file.exists(data_file)) { 
            sim_i <- sim_i + 1
            data <- read.table(data_file, header=TRUE, skip=5)
            start <- which(data$log_L == min(data$log_L))
            HR <- data$log_L[-1:-start] ~ data$log_Teff[-1:-start]
            if (sim_i == 1) {
                par(bty="l", las=1, mar=c(3, 3.2, 3, 1), mgp=c(1.8, 0.25, 0))
                plot(HR, type='l', col=cl[sim_i], 
                     xlab=expression(log~T[eff]), 
                     ylab=expression(log~L/L["⊙"]), 
                     main=paste("HR diagram of Sun-like stars by", exp),
                     xlim=rev(c(round(min(data$log_Teff), 2), 
                                round(max(data$log_Teff)+0.05, 2))),
                     ylim=c(floor(min(data$log_L)),
                            ceiling(max(data$log_L))),
                     xaxs='i', yaxs='i', tck=0.01)
                minor.tick(nx=5, ny=5, tick.ratio=-0.25)
            } else {
                lines(HR, col=cl[sim_i])
            }
            for (ev_stage in names(profile_nos)) {
                sim_num <- grep(simulation, profile_nos[[ev_stage]])
                print(sim_num)
                pro_num <- as.numeric(sub('.data', '', 
                    sub(".+profile", '', profile_nos[[ev_stage]][sim_num])))
                print(pro_num)
                print(data$log_Teff[data$model_number==pro_num])
                print(data$log_L[data$model_number==pro_num])
                points(data$log_Teff[data$model_number==pro_num],
                       data$log_L[data$model_number==pro_num], 
                       col="black", pch=21, cex=0.1)
            }
        }
    }
    points(Teff_sun, 0, col="black")
    points(Teff_sun, 0, col="black", pch=21, cex=0.1)
    legend("topleft", col=cl, lty=1, bty='n', legend=labels)
    dev.off()
}
warnings()
