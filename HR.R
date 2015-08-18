#### Plot HR diagrams for simulations in all subdirectories 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(magicaxis)

source('utils.R')

experiments <- 'exp'
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)
fgong_dir <- file.path('fgongs')
dir.create(fgong_dir, showWarnings=FALSE)

modelS <- read.table('modelS.dat', header=TRUE) 

for (experiment in list.dirs(experiments, recursive=FALSE)) {
    simulations <- list.dirs(experiment, recursive=FALSE)
    simulations <- simulations[order(as.numeric(basename(simulations)))]
    load_experiment_info(experiment, simulations)
    
    ###########################################################
    ### Find profile files at different stages of evolution ###
    ###########################################################
    print("Finding profile files")
    pro_files <- list()
    for (simulation in simulations) {
        files <- list.files(file.path(simulation, "LOGS"))
        ev_pro_files <- grep('profile_.+.data$', files)
        for (ev_pro_file in ev_pro_files) {
            ev_stage <- sub('.data','',sub('profile_','',files[ev_pro_file]))
            pro_files[[ev_stage]] <- c(pro_files[[ev_stage]], 
                file.path(simulation, 'LOGS', files[ev_pro_file]))
        }
    }
    
    #######################
    ### Plot HR diagram ###
    #######################
    print("Plotting HR diagram")
    cairo_pdf(file.path(plot_dir, paste0('HR_', basename(experiment), '.pdf')), 
              width=plot_width, height=plot_height, family=font)
    seismology <- list()
    sim_i <- 0
    for (simulation in simulations) {
        data_file <- file.path(simulation, 'LOGS', 'history.data')
        if (file.exists(data_file)) { 
            sim_i <- sim_i + 1
            data <- read.table(data_file, header=TRUE, skip=5)
            start <- which(data$log_L == min(data$log_L))
            log_L <- data$log_L[-1:-start]
            log_Teff <- data$log_Teff[-1:-start]
            HR <- log_L ~ log_Teff
            if (sim_i == 1) {
                par(bty="l", las=1, mar=c(3, 3.2, 3, 1), cex.lab=1.3, 
                    mgp=c(1.8, 0.25, 0))
                plot(HR, type='l', col=cl[sim_i], 
                     xlab=expression(log~T[eff]), 
                     ylab=expression(log~L/L["⊙"]), 
                     main=paste("HR diagram of Sun-like stars by", 
                                experiment_name),
                     xlim=rev(round(range(log_Teff)+c(0,0.05), 1)),
                     ylim=round(range(log_L)*2+c(0,0.5), 0)/2,
                     xaxs='i', yaxs='i', tck=0.01)
                magaxis(side=1:4, tcl=0.25, labels=FALSE)
            } else {
                lines(HR, col=cl[sim_i])
            }
            for (ev_stage in names(pro_files)) {
                pros <- grep(simulation, pro_files[[ev_stage]])
                for (pro_file in pro_files[[ev_stage]][pros]) {
                    model <- data$model_number == read.table(pro_file, 
                        header=TRUE, nrows=1, skip=1)$model_number
                    points(data$log_Teff[model], 
                           data$log_L[model], 
                           col="black", pch=21, cex=0.1)
                    seismology[[ev_stage]] <- c(seismology[[ev_stage]], 
                        basename(simulation), 
                        data$acoustic_cutoff[model]/(2*pi), 
                        data$delta_nu[model], data$nu_max[model])
                }
            }
        }
    }
    points(Teff_sun, 0, col="black")
    points(Teff_sun, 0, col="black", pch=21, cex=0.1)
    legend("topright", col=cl, lty=1, bty='n', legend=labls)
    dev.off()
    
    print("Plotting internal profiles")
    for (ev_stage in names(pro_files)) {
        plot_subdir <- file.path(plot_dir, ev_stage)
        dir.create(plot_subdir, showWarnings=FALSE)
        
        fgong_subdir <- file.path(fgong_dir, basename(experiment), ev_stage)
        dir.create(fgong_subdir, showWarnings=FALSE, recursive=TRUE)
        for (pro_file in pro_files[[ev_stage]]) {
            fgong_file <- paste0(pro_file, '.FGONG')
            new_path <- file.path(fgong_subdir, 
                paste0(basename(dirname(dirname(fgong_file))), '.FGONG'))
            if (file.exists(fgong_file) && !file.exists(new_path))
                file.copy(fgong_file, new_path)
        }
        write.table(matrix(seismology[[ev_stage]], 
                    nrow=length(simulations), byrow=TRUE), 
            row.names=FALSE, quote=FALSE, 
            col.names=c('name', 'acoustic_cutoff', 'delta_nu', 'nu_max'), 
            file=file.path(dirname(fgong_subdir), 
                paste0('seismology_', ev_stage, '.dat')))
        
        ## Start profile plot devices
        start_dev("Profile differences in", "diff", experiment, ev_stage)
        start_dev("Profiles of", "profile", experiment, ev_stage)
        
        ## Iterate through the different internal profiles
        for (col_name in c("csound", "density", "pressure", "temperature")) {
            if (col_name=="csound") {
                main_label <- "Sound speed"
                symbol <- bquote(c)
                ylabel <- expression(sound~speed~c~"["*cm/sec*"]")
            } else if (col_name=="density") {
                main_label <- "Density"
                symbol <- bquote(rho)
                ylabel <- expression(density~rho~"["*g/cm^3*"]")
            } else if (col_name=="pressure") {
                main_label <- "Pressure"
                symbol <- bquote(p)
                ylabel <- expression(pressure~p~"["*dyn/cm^2*"]")
            } else if (col_name=="temperature") {
                main_label <- "Temperature"
                symbol <- bquote(T)
                ylabel <- expression(temperature~T~"["*K*"]")
            }
            modelS_yt <- modelS[[col_name]]
            modelS_y <- log10(modelS_yt)
            
            ## Find plot limits 
            y_min <- ifelse(ev_stage=="solar-age", min(modelS_y), Inf)
            y_max <- ifelse(ev_stage=="solar-age", max(modelS_y), -Inf)
            x_max <- 1
            for (simulation in simulations) {
                sim_no <- grep(simulation, pro_files[[ev_stage]])
                data_file <- file.path(pro_files[[ev_stage]][sim_no])
                if (file.exists(data_file)) {
                    data <- read.table(data_file, header=TRUE, skip=5)
                    y <- log10(data[[col_name]])
                    if (min(y) < y_min) y_min <- min(y)
                    if (max(y) > y_max) y_max <- max(y)
                    if (round(max(data$radius)+.05, 1) > x_max) 
                        x_max <- round(max(data$radius)+.05, 1)
                }
            }
            
            # Preallocate interpolation matrix
            interps <- matrix(nrow=length(approx_xout), 
                              ncol=length(simulations))
            
            ##############################
            ### Plot internal profiles ###
            ##############################
            sim_i <- 0
            for (simulation in simulations) {
                sim_no <- grep(simulation, pro_files[[ev_stage]])
                data_file <- file.path(pro_files[[ev_stage]][sim_no])
                if (file.exists(data_file)) {
                    sim_i <- sim_i + 1
                    data <- read.table(data_file, header=TRUE, skip=5)
                    yt <- data[[col_name]]
                    y <- log10(yt)
                    interps[,sim_i] <- approx(data$radius, yt, 
                        n=length(approx_xout), xout=approx_xout, rule=2)$y
                    if (sim_i == 1) {
                        par(bty="l", las=1, cex.lab=1.3,
                            mar=c(3, 3.5, 1, 1), mgp=c(1.8, 0.25, 0))
                        plot(data$radius, y, type='l', col=cl[sim_i],
                             xlim=c(0, x_max),
                             ylim=c(y_min, y_max),#+0.01*y_max), 
                             xaxs='i', tck=0.01, yaxt='n', yaxs='i', 
                             xlab="",
                             ylab=ylabel)
                        title(xlab=expression(fractional~radius~r/R["⊙"]))
                        magaxis(side=1:2, unlog='y', family=font, tcl=0.25,
                            mgp=c(1.8, 0.25, 0),
                            labels=c(FALSE, TRUE))
                    } else {
                        lines(data$radius, y, col=cl[sim_i])
                    }
                }
            }
            if (ev_stage=="solar-age") ## plot Model S
                lines(modelS$radius, modelS_y, col='black', lty=2)
            if (col_name=="density") make_legend(labls, ev_stage)
            dev.set(dev.prev())
            
            ################################
            ### Plot profile differences ###
            ################################
            ref <- if(ev_stage=="solar-age")
                approx(modelS$radius, modelS_yt, n=length(approx_xout), 
                       xout=approx_xout, rule=2)$y
            else interps[,sun_num]
            rel_diff <- apply(interps, 2, function(x) {(x-ref)/ref})
            row_offset <- round(nrow(rel_diff) * 0.5, 0)
            rd <- function(x) diff(range(x))
            core_var <- apply(rel_diff[1:(row_offset-1),],1,rd)
            row_max <- if (any(core_var > 1)) {
                min(which(core_var>0.1))
            } else {
                outer_var <- apply(rel_diff[row_offset:nrow(rel_diff),],1,rd)
                endpts <- outer_var > max(core_var) & outer_var > 0.01 | 
                          outer_var > 1
                ifelse(any(endpts), min(which(endpts))+row_offset-1, 
                       nrow(rel_diff))
            }
            if (row_max <= 1 || row_max >= nrow(rel_diff)) 
                row_max <- nrow(rel_diff)
            x_max <- max(approx_xout[1:row_max])
            for (sim_i in 1:ncol(rel_diff)) {
                if (sim_i == 1) {
                    par(bty='l', las=1, xpd=TRUE, cex.lab=1.3,
                        mar=c(3, 5, 1, 1), mgp=c(2.5, 0.25, 0))
                    plot(approx_xout[1:row_max], rel_diff[1:row_max, sim_i], 
                         type='l', col=cl[sim_i],
                         xaxs='i', yaxs='i',
                         tck=0.01,
                         xlim=c(0, x_max),
                         ylim=range(rel_diff[1:row_max,]),
                         ylab=bquote("rel."~.(tolower(main_label))~
                             "diff."~delta * .(symbol) / .(symbol)),
                         xlab="")
                    title(xlab=expression("fractional radius"~r/R["⊙"]), 
                          mgp=par()$mgp-c(0.5,0,0))
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                    lines(c(x_max,0), c(0, 0), lty=2)
                } else {
                    lines(approx_xout[1:row_max], rel_diff[1:row_max, sim_i], 
                          col=cl[sim_i])
                }
            }
            if (col_name=="density") make_legend(labls, ev_stage)
            dev.set(dev.next())
        }
        dev.off()
        dev.off()
        
        ###########################################
        ### Plot Lamb/Brunt-Vaisala frequencies ###
        ###########################################
        print("Plotting Lamb frequencies")
        cairo_pdf(file.path(plot_subdir, paste0("lamb", '_', 
            basename(experiment), '_', ev_stage, '.pdf')),
            width=plot_width, height=plot_height, family=font)
        main_label <- "Frequency"
        ylabel <- expression(frequency~nu~"["*mu*Hz*"]")
        y_min <- Inf
        y_max <- -Inf
        x_max <- 1
        for (simulation in simulations) {
            sim_no <- grep(simulation, pro_files[[ev_stage]])
            data_file <- file.path(pro_files[[ev_stage]][sim_no])
            if (file.exists(data_file)) {
                data <- read.table(data_file, header=TRUE, skip=5)
                y <- data[["csound"]]**2 /
                    (data[["radius"]]*solar_radius)**2
                y <- y[is.finite(y)]
                y_s1 <- log10(10**6/(2*pi) * sqrt(1*(1+1)*y))
                y_s3 <- log10(10**6/(2*pi) * sqrt(3*(3+1)*y))
                if (min(y_s1) < y_min) y_min <- min(y_s1)
                if (max(y_s3) > y_max) y_max <- max(y_s3)
                if (round(max(data$radius)+.05, 1) > x_max) 
                    x_max <- round(max(data$radius)+.05, 1)
            }
        }
        sim_i <- 0
        for (simulation in simulations) {
            sim_no <- grep(simulation, pro_files[[ev_stage]])
            data_file <- file.path(pro_files[[ev_stage]][sim_no])
            if (file.exists(data_file)) {
                sim_i <- sim_i + 1
                data <- read.table(data_file, header=TRUE, skip=5)
                brunt_N2 <- data[["brunt_N2"]]
                stable <- which(brunt_N2 >= 0)
                indices <-  stable[which(diff(stable)!=1)[1]+1] : max(stable)
                y <- log10(10**6/(2*pi) * sqrt(brunt_N2[indices]))
                if (sim_i == 1) {
                    par(las=1, mar=c(3, 3.2, 3, 1), cex.lab=1.3, 
                        mgp=c(1.8, 0.25, 0))
                    plot(data$radius[indices], y, type='l', cex=0.25,
                         col=cl[sim_i],
                         xlim=c(0, x_max),
                         ylim=c(y_min, y_max+0.01*y_max), 
                         xaxs='i', yaxt='n', tck=0, #yaxs='i', 
                         xlab="", ylab=ylabel,
                         main=paste("Critical frequencies of", ev_stage, 
                                    "stars\nby", experiment_name))
                    title(xlab=expression(r/R["⊙"]))
                    magaxis(side=1:4, unlog='y', family=font, tcl=0.25,
                            mgp=c(1.8, 0.25, 0),
                            labels=c(FALSE, TRUE, FALSE, FALSE))
                } else {
                    lines(data$radius[indices], y, 
                          pch=20, cex=0.25, col=cl[sim_i])
                }
                for (ell in 1:3) {
                    S <- data[["csound"]]**2 / 
                        (data[["radius"]]*solar_radius)**2
                    y <- log10(10**6/(2*pi) * sqrt(ell*(ell+1)*S))
                    lines(data$radius, y, col=cl[sim_i], lty=ell+1)
                }
            }
        }
        legend("topleft", bty='n', inset=0,
               lty=1:4, legend=c("Brunt-Väisälä", 
                        expression("Lamb"~S[1]),
                        expression("Lamb"~S[2]),
                        expression("Lamb"~S[3])))
        legend("topright", bty='n', inset=0, col=cl, lty=1, legend=labls)
        dev.off()
    }
}
warnings()
