#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(RColorBrewer)
library(magicaxis)

source('utils.R')

fgong_dir <- file.path('fgongs')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

modelS <- read.table('fgong.l5bi.d.dat')

for (experiment in list.dirs(fgong_dir, recursive=FALSE)) {
    for (ev_stage_dir in list.dirs(experiment, recursive=FALSE)) {
        ev_stage <- basename(ev_stage_dir)
        #if (grepl('RGB', ev_stage)) next()
        
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
        
        seismology <- read.table(file.path(experiment, 
            paste0('seismology_', ev_stage, '.dat')), header=TRUE)
        print(seismology)
        
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('echelle_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        
        widths <- c(.44, .44, .12)
        start_dev("Frequency differences in", "freqdiffs", 
            basename(experiment), ev_stage, width=widths)
        start_dev("Frequencies of", "freqs", basename(experiment), ev_stage,
            widths)
        
        y_min <- Inf
        y_max <- -Inf
        x_max <- -Inf
        x_min <- Inf
        max_delta_nu <- -Inf
        nu_min <- Inf
        for (model_i in data_files) {
            seismo_row <- which(seismology$name == 
                    as.numeric(sub('.dat', '', files[model_i])))
            delta_nu <- seismology$delta_nu[seismo_row] 
            if (delta_nu > max_delta_nu) max_delta_nu <- delta_nu
            
            nu_max <- seismology$nu_max[seismo_row] 
            cutoff_freq <- seismology$acoustic_cutoff[seismo_row] 
            
            data <- read.table(file.path(ev_stage_dir, files[model_i]))
            data <- data[data[,3] < cutoff_freq,]
            if (nrow(data)==0) next
            if (min(data[,3]) < y_min) y_min <- min(data[,3])
            if (max(data[,3]) > y_max) y_max <- max(data[,3])
            if (max(data[,2]) > x_max) x_max <- max(data[,2])
            if (min(data[,2]) < x_min) x_min <- min(data[,2])
            
            nu <- data[(data[,3]<cutoff_freq) & (data[,2]>=0)
                     & (data[,3]>(nu_max-5*delta_nu)),]
            if (min(nu[,3]) < nu_min) nu_min <- min(nu[,3])
        }
        
        delta_nus <- list()
        deltas <- list()
        noisy_delta <- list()
        for (l_mode in 0:3) {
            l_name <- if (l_mode == 0) "radial"
                 else if (l_mode == 1) "dipole"
                 else if (l_mode == 2) "quadrupole"
                 else "octupole"
            
            ############################
            ### Make frequency plots ###
            ############################
            counter <- 1
            for (model_i in data_files) {
                seismo_row <- which(seismology$name == 
                    as.numeric(sub('.dat', '', files[model_i])))
                cutoff_freq <- seismology$acoustic_cutoff[seismo_row]
                
                data <- read.table(file.path(ev_stage_dir, files[model_i]))
                ell <- data[data[,1]==l_mode,]
                ell <- ell[ell[,3]<cutoff_freq,]
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
                         col=cl[counter], 
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
                           col=cl[counter], lwd=0.5)
                }
                counter <- counter + 1
            }
            if (grepl('solar-like', ev_stage)) {
                ell <- modelS[modelS[,1]==l_mode,]
                ell <- ell[ell[,3]<y_max,]
                relation <- ell[,3] ~ ell[,2]
                points(relation, pch=1, cex=0.5, col="black", lwd=0.5)
                ref <- ell
            }
            if (l_mode==1) make_legend(labls, ev_stage, lty=FALSE)
            dev.set(dev.prev())
            
            #############################
            ### Make difference plots ###
            #############################
            print(ref[,2])
            n_range <- if (length(ref[,2])>0) min(ref[,2]):max(ref[,2]) else c()
            diffs <- matrix(nrow=max(ref[,2])-min(ref[,2])+1, #x_max-x_min+1,
                            ncol=length(simulations))
            counter <- 1
            for (model_i in data_files) {
                seismo_row <- which(seismology$name == 
                    as.numeric(sub('.dat', '', files[model_i])))
                cutoff_freq <- seismology$acoustic_cutoff[seismo_row]
                
                data <- read.table(file.path(ev_stage_dir, files[model_i]))
                ell <- data[data[,1]==l_mode,]
                ell <- ell[ell[,3]<cutoff_freq,]
                
                x_counter <- 1
                for (x in n_range) {
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
            for (model_i in data_files) {
                if (counter == 1) {
                    par(bty="l", las=1, mar=c(3, 3.6, 1, 1), 
                        mgp=c(2.2, 0.25, 0))
                    plot(n_range,#x_min:x_max, 
                         diffs[,counter], pch=20, cex=0.01,
                         xlim=c(min(0, min(n_range)-1), max(n_range)), 
                         ylim=range(diffs[complete.cases(diffs)]),
                         xaxs='i', lwd=0.5,
                         col=cl[counter], 
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
                           col=cl[counter])
                }
                counter <- counter + 1
            }
            if (grepl('solar-like', ev_stage)) {
                ell <- modelS[modelS[,1]==l_mode,]
                relation <- rep(0, length(ell[,2])) ~ ell[,2]
                points(relation, pch=1, cex=0.5, col="black", lwd=0.5)
            }
            if (l_mode==1) make_legend(labls, ev_stage, lty=FALSE)
            dev.set(dev.prev())
            
            ##########################
            ### Make Echelle plots ###
            ##########################
            counter <- 1
            for (model_i in data_files) {
                seismo_row <- which(seismology$name == 
                    as.numeric(sub('.dat', '', files[model_i]))) 
                delta_nu <- seismology$delta_nu[seismo_row] 
                nu_max <- seismology$nu_max[seismo_row] 
                cutoff_freq <- seismology$acoustic_cutoff[seismo_row] 
                
                data <- read.table(file.path(ev_stage_dir, files[model_i]))
                ell <- data[data[,1]==l_mode,]
                ell <- ell[(ell[,3]<cutoff_freq) & (ell[,2]>=0)
                         & (ell[,3]>(nu_max-5*delta_nu)),]
                #delta_nu <- median(diff(ell[,3]))
                converted_fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
                gaussian_env <- dnorm(ell[,3], nu_max, converted_fwhm)
                fit <- lm(ell[,3] ~ ell[,2], weights=gaussian_env)
                delta_nu <- coef(fit)[2]
                l_name <- toString(l_mode)
                delta_nus[[l_name]] <- c(delta_nus[[l_name]], delta_nu)
                deltas[[l_name]] <- c(deltas[[l_name]], 
                    summary(fit)$coefficients[2,2])
                
                ### estimate std errors
                out <- rep(0, 100)
                d <- if (l_mode == 0) 0.91
                else if (l_mode == 1) 0.61
                else if (l_mode == 2) 0.87
                else if (l_mode == 3) 1.4
                for (ii in 1:100) {
                    new_fit <- lm(ell[,3] + rnorm(nrow(ell), 0, d) ~ ell[,2], 
                        weights=gaussian_env)
                    out[ii] <- summary(new_fit)$coefficients[2,2]
                    #coef(new_fit)[2]
                    #noisy_nu <- ell[,3] + rnorm(nrow(ell), 0, 1)
                    #new_ns <- 1:nrow(ell)
                    ##new_ns <- sample(nrow(ell),replace=T)
                    #new_fit <- lm(noisy_nu[new_ns] ~ ell[new_ns,2], 
                    #    weights=gaussian_env[new_ns])
                    #print(summary(new_fit)$sigma**2)
                    #out[ii] <- summary(new_fit)$coefficients[2,2]
                    ##coef(new_fit)[2]
                }
                noisy_delta[[l_name]] <- c(noisy_delta[[l_name]],
                    mean(out))
                print(c(summary(fit)$coefficients[2,2], mean(out), max(out)))
                #    sqrt((0.1 + summary(fit)$sigma**2)
                #          * summary(fit)$cov.unscaled[2,2])))
                
                relation <- ell[,3] ~ ell[,3]%%delta_nu
                if (l_mode == 0 && counter == 1) {
                    layout(matrix(c(1,1,2,3), ncol=2, byrow=TRUE), 
                           heights=c(0.14,0.86), widths=c(.85, .15))
                    par(mar=rep(0,4))
                    plot.new()
                    text(0.5, 0.5, paste("Echelle diagram of", ev_stage, 
                         "stars\nby", experiment_name), cex=2, font=2)
                    par(bty="l", las=1, mar=c(3, 3.4, 0.1, 0.1), 
                        mgp=c(2, 0.25, 0))
                    plot(relation, pch=l_mode, col=cl[counter],
                         main="", cex=0.75,
                         ylim=c(round(nu_min-500,-3), round(y_max+1000, -3)),
                         xlim=c(0, round(max_delta_nu+5, -1)),
                         xlab="", xaxs='i', yaxs='i', tck=0.01, 
                         ylab=bquote("frequency"~nu~"["*mu*Hz*"]"))
                    title(xlab=expression(nu~mod~Delta*nu), 
                          mgp=par()$mgp-c(0.4,0,0))
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    points(relation, pch=l_mode, col=cl[counter], cex=0.75)
                }
                abline(v=delta_nu, col=cl[counter], lty=l_mode+1, lwd=0.5)
                counter <- counter+1
            }
            if (l_mode==3) {
                par(mar=rep(0,4))
                plot.new()
                legend(-0.1, 0.9, bty='n', inset=0, col=cl, pch=20, 
                       legend=labls)
                legend(0.1, 0.45, bty='n', inset=0, pch=0:3, 
                       legend=c(expression("\u2113"==0), 
                                expression("\u2113"==1),
                                expression("\u2113"==2),
                                expression("\u2113"==3)))
                
                legend(-0.1, 0.25, bty='n', inset=0, lty=1:4, 
                       legend=c(expression(Delta*nu["\u2113"==0]), 
                                expression(Delta*nu["\u2113"==1]),
                                expression(Delta*nu["\u2113"==2]),
                                expression(Delta*nu["\u2113"==3])))
            }
            dev.set(dev.next())
            dev.set(dev.next())
        }
        dev.off()
        dev.off()
        dev.off()
        
        ################################
        ### Make plots of exp vs dnu ###
        ################################
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('ell_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        layout(matrix(c(1,1,2,3), ncol=2, byrow=TRUE), 
               heights=c(0.14,0.86), widths=c(.85, .15))
        par(mar=rep(0,4))
        plot.new()
        text(0.5, 0.5, paste("Large frequency separation of", ev_stage,
             "\nstars by", experiment_name), cex=2, font=2)
        par(bty="l", las=1, mar=c(3, 3.4, 0.1, 0.1), 
            mgp=c(2, 0.25, 0))
        for (l_mode in 0:3) {
            relation <- delta_nus[[toString(l_mode)]] ~ exp_vals
            if (l_mode == 0) {
                plot(relation, pch=l_mode, main="", xlab="", 
                     tck=0.01, ylim=range(delta_nus),
                     ylab=expression("large frequency separation"~Delta*nu))
                title(xlab=bquote(.(experiment_name)~.(freep)), 
                      mgp=par()$mgp-c(0.4,0,0))
                lines(relation, lty=l_mode+1)
                minor.tick(nx=5, ny=5, tick.ratio=-0.15)   
            } else {
                points(relation, pch=l_mode)
                lines(relation, lty=l_mode+1)
            }
        }
        par(mar=rep(0,4))
        plot.new()
        legend("left", bty='n', inset=0, pch=0:3, lty=1:4,
               legend=c(expression("\u2113"==0), 
                        expression("\u2113"==1),
                        expression("\u2113"==2),
                        expression("\u2113"==3)))
        dev.off()
        
        #####################################
        ### Make plots of exp vs dnu diff ###
        #####################################
        ell_cl <- brewer.pal(4, "BrBG")
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('elldiff_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights=c(0.14,0.86))
        par(mar=rep(0,4))
        plot.new()
        text(0.5, 0.5, paste("Differences in Δν of", ev_stage,
             "\nstars by degree and", experiment_name), cex=2, font=2)
        par(bty="l", las=1, mar=c(3, 3.4, 0.1, 1), 
            mgp=c(2, 0.25, 0))
        y_max <- -Inf
        y_min <- Inf
        for (l_mode in 1:3) {
            y <- delta_nus[[toString(l_mode)]] - delta_nus[["0"]]
            if (max(y) > y_max) y_max <- max(y)
            if (min(y) < y_min) y_min <- min(y)
        }
        for (l_mode in 1:3) {
            i_mid <- round(median(1:length(exp_vals)))
            l_name <- toString(l_mode)
            difference <- delta_nus[[l_name]] - delta_nus[["0"]]
            relation <- difference ~ exp_vals
            if (l_mode == 1) {
                plot(relation, pch=l_mode, main="", xlab="", yaxs='i',
                     col=ell_cl[l_mode+1],
                     tck=0.01, ylim=c(round(y_min), round(y_max+0.5)),
                     ylab=expression("large frequency separation difference"
                                     ~delta*Delta*nu["\u2113"]==
                                     Delta*nu["\u2113"]-
                                     Delta*nu[0]))
                title(xlab=bquote(.(experiment_name)~.(freep)), 
                      mgp=par()$mgp-c(0.4,0,0))
                lines(relation, lty=l_mode+1, col=ell_cl[l_mode+1])
                errbar(exp_vals[i_mid], 0, 0, 
                    noisy_delta[["0"]][i_mid]/2, 
                    add=TRUE, lty=5, pch='-')
                magaxis(side=1:4, labels=FALSE)
            } else {
                points(relation, pch=l_mode, col=ell_cl[l_mode+1])
                lines(relation, lty=l_mode+1, col=ell_cl[l_mode+1])
            }
            i_off <- 1
            if (l_mode == 2) i_off <- -1
            if (l_mode == 3) i_off <- 0
            errbar(exp_vals[i_mid+i_off], difference[i_mid+i_off], 
                 difference[i_mid+i_off] - noisy_delta[[l_name]][i_mid+i_off]/2,
                 difference[i_mid+i_off] + noisy_delta[[l_name]][i_mid+i_off]/2,
                 add=TRUE, pch=l_mode, lty=5, 
                 col=ell_cl[l_mode+1], errbar.col=ell_cl[l_mode+1])
        }
        legend("topright", bty='n', inset=0, pch=1:3, lty=2:4, col=ell_cl[2:4],
            legend=c(expression(delta*Delta*nu[1]), 
                     expression(delta*Delta*nu[2]),
                     expression(delta*Delta*nu[3])))
        dev.off()
        
    }
}
warnings()
