#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source('utils.R')

fgong_dir <- file.path('fgongs')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

modelS <- read.table('fgong.l5bi.d.dat')

for (experiment in list.dirs(fgong_dir, recursive=FALSE)) {
    for (ev_stage_dir in list.dirs(experiment, recursive=FALSE)) {
        ev_stage <- basename(ev_stage_dir)
        plot_subdir <- file.path(plot_dir, ev_stage)
        dir.create(plot_subdir, showWarnings=FALSE)
        
        files <- list.files(ev_stage_dir)
        print(files)
        data_files <- grep('*dat$', files)
        print(data_files)
        if (!any(data_files)) next
        
        simulations <- list.dirs(ev_stage_dir, recursive=FALSE)
        simulations <- simulations[order(as.numeric(basename(simulations)))]
        load_experiment_info(experiment, simulations)
        
        widths <- c(.44, .44, .12)
        start_dev("Frequency differences in", "freqdiffs", 
                  basename(experiment), ev_stage, width=widths)
        start_dev("Frequencies of", "freqs", basename(experiment), ev_stage,
                  widths)
        y_min <- Inf
        y_max <- -Inf
        x_max <- -Inf
        x_min <- Inf
        for (model_i in data_files) {
            data <- read.table(file.path(ev_stage_dir, files[model_i]))
            if (min(data[,3]) < y_min) y_min <- min(data[,3])
            if (max(data[,3]) > y_max) y_max <- max(data[,3])
            if (max(data[,2]) > x_max) x_max <- max(data[,2])
            if (min(data[,2]) < x_min) x_min <- min(data[,2])
        }
        
        for (l_mode in 0:3) {
            l_name <- if (l_mode == 0) "radial"
                 else if (l_mode == 1) "dipole"
                 else if (l_mode == 2) "quadrupole"
                 else "octupole"
            counter <- 1
            if (experiment_name=='mixing length') {
                fs <- rev(data_files)
                tcl <- rev(cl)
            } else {
                fs <- data_files
                tcl <- cl
            }
            
            ############################
            ### Make frequency plots ###
            ############################
            for (model_i in fs) {
                data <- read.table(file.path(ev_stage_dir, files[model_i]))
                ell <- data[data[,1]==l_mode,]
                if (!grepl('solar-like', ev_stage) && counter == sun_num)
                    ref <- ell
                relation <- ell[,3] ~ ell[,2]
                if (counter == 1) {
                    par(bty="l", las=1, mar=c(3, 3.6, 1, 1), 
                        mgp=c(2.2, 0.25, 0))
                    plot(relation, pch=20, cex=0.01, lwd=0.5,
                         xlim=c(min(0, x_min-1), x_max), 
                         ylim=c(0, round(y_max+50, -2)),
                         yaxs='i', xaxs='i',
                         col=tcl[counter], 
                         tck=0.01, 
                         xlab='',
                         ylab=bquote(.(l_name)~"frequency"~nu~"["*mu*Hz*"]"))
                                     #~omega~"["*rad/s*"]"))
                                     #"["*mu*Hz*"]"))
                    title(xlab=expression("radial order"~n), 
                          mgp=par()$mgp-c(0.4,0,0))
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    points(relation, pch=20, cex=0.01, 
                           col=tcl[counter], lwd=0.5)
                }
                counter <- counter + 1
            }
            if (grepl('solar-like', ev_stage)) {
                ell <- modelS[modelS[,1]==l_mode,]
                relation <- ell[,3] ~ ell[,2]
                points(relation, pch=1, cex=0.5, col="black", lwd=0.5)
                ref <- ell
            }
            if (l_mode==1) make_legend(labls, ev_stage, lty=FALSE)
            dev.set(dev.prev())
            
            #############################
            ### Make difference plots ###
            #############################
            n_range <- min(ref[,2]):max(ref[,2])
            diffs <- matrix(nrow=max(ref[,2])-min(ref[,2])+1, #x_max-x_min+1,
                            ncol=length(simulations))
            counter <- 1
            for (model_i in fs) {
                data <- read.table(file.path(ev_stage_dir, files[model_i]))
                ell <- data[data[,1]==l_mode,]
                x_counter <- 1
                for (x in n_range) {#x_min:x_max) {
                    diffs[x_counter, counter] <- 
                        if (x%in%ref[,2] && x%in%ell[,2]) {
                            a <- 2*pi*Mode(unlist(ell[ell[,2]==x,][3]))
                            b <- 2*pi*Mode(unlist(ref[ref[,2]==x,][3]))
                            (a-b)/b
                            } else NA
                    x_counter <- x_counter + 1
                }
                counter <- counter + 1
            }
            counter <- 1
            for (model_i in fs) {
                if (counter == 1) {
                    par(bty="l", las=1, mar=c(3, 3.6, 1, 1), 
                        mgp=c(2.2, 0.25, 0))
                    plot(n_range,#x_min:x_max, 
                         diffs[,counter], pch=20, cex=0.01,
                         xlim=c(min(0, min(n_range)-1), max(n_range)), 
                         ylim=range(diffs[complete.cases(diffs)]),
                         xaxs='i', lwd=0.5,
                         col=tcl[counter], 
                         tck=0.01, 
                         xlab='',
                         ylab=bquote({delta * omega["\u2113"==.(l_mode)] /
                                      omega["\u2113"==.(l_mode)]}))
                    title(xlab=expression("radial order"~n), 
                          mgp=par()$mgp-c(0.4,0,0))
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    points(n_range,#x_min:x_max, 
                           diffs[,counter], pch=20, cex=0.01, lwd=0.5,
                           col=tcl[counter])
                }
                counter <- counter + 1
            }
            if (grepl('solar-like', ev_stage)) {
                ell <- modelS[modelS[,1]==l_mode,]
                relation <- rep(0, length(ell[,2])) ~ ell[,2]
                points(relation, pch=1, cex=0.5, col="black", lwd=0.5)
            }
            if (l_mode==1) make_legend(labls, ev_stage, lty=FALSE)
            dev.set(dev.next())
        }
        dev.off()
        dev.off()
    }
}
