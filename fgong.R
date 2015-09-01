#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(RColorBrewer)
library(magicaxis)
library(sfsmisc)

source('utils.R')
ell_cl <- brewer.pal(4, "BrBG")

fgong_dir <- file.path('fgongs')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

#widths <- c(.44, .44, .12)
picker <- c(TRUE)#, FALSE, FALSE)
n_samples <- 100
col_names <- c('l', 'n', 'nu')

modelS <- read.table(file.path('data', 'fgong.l5bi.d.dat'), col.names=col_names)

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
        
        i_mid <- round(median(1:length(exp_vals)))
        
        seismology <- read.table(file.path(experiment, 
            paste0('seismology_', ev_stage, '.dat')), header=TRUE)
        print(seismology)
        
        start_dev("Gradient diagram of", "gradient", 
            basename(experiment), ev_stage)
        
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('echelle_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        
        start_dev("Frequency differences in", "angdiff", 
            basename(experiment), ev_stage)
        start_dev("Frequencies of", "quadfreq", basename(experiment), ev_stage)
        
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('freq_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        
        cache <- list()
        scaler <- list()
        nu_maxs <- list()
        #y_min <- Inf
        y_max <- -Inf
        x_max <- -Inf
        x_min <- Inf
        #max_delta_nu <- -Inf
        p_nu_min <- Inf
        p_nu_max <- -Inf
        for (model_i in data_files) {
            exp_value <- as.numeric(sub('.dat', '', files[model_i]))
            M <- if (basename(experiment) == 'M') exp_value else 1
            M <- M * solar_mass
            seismo_row <- which(seismology$name == exp_value)
            R <- seismology$radius[seismo_row] * solar_radius
            scaler[[model_i]] <- solar_scale / sqrt(M/R^3)
            
            #delta_nu <- seismology$delta_nu[seismo_row] * scaler[[model_i]]
            #if (delta_nu > max_delta_nu) max_delta_nu <- delta_nu
            
            nu_max <- seismology$nu_max[seismo_row] * scaler[[model_i]]
            cutoff_freq <- seismology$acoustic_cutoff[seismo_row] * 
                scaler[[model_i]]
            nu_maxs[[model_i]] <- nu_max
            
            data <- read.table(file.path(ev_stage_dir, files[model_i]),
                col.names=col_names)
            data <- data[data$nu < cutoff_freq,]
            data$nu <- data$nu * scaler[[model_i]]
            
            ## There is a bug in ADIPLS that causes l=1 frequencies to sometimes
            ## duplicate. The workaround is to shift all n at the duplicate
            ## down by one and remove the n=0.
            ell <- data[data$l==1,]
            ns <- ell$n[ell$n>0]
            if (any(duplicated(ns))) {
                dup <- which(duplicated(ns))[1]
                toshift <- ns[which(ns>0) < dup]
                ell$n[ell$n>0][toshift] <- ns[toshift] - 1
                data[data$l==1,] <- ell
                data <- data[!(data$l==1 & data$n==0),]
            }
            
            cache[[model_i]] <- data
            
            if (nrow(data)==0) next
            #if (min(data$nu) < y_min) y_min <- min(data$nu)
            if (max(data$nu) > y_max) y_max <- max(data$nu)
            if (max(data$n) > x_max) x_max <- max(data$n)
            if (min(data$n) < x_min) x_min <- min(data$n)
            
            p_nu <- data[data$n>0,] #& data$nu>(nu_max-5*delta_nu),]
            if (min(p_nu$nu) < p_nu_min) p_nu_min <- min(p_nu$nu)
            if (max(p_nu$nu) > p_nu_max) p_nu_max <- max(p_nu$nu)
        }
        
        delta_nus <- list()
        deltas <- list()
        dnu_err <- list()
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
                data <- cache[[model_i]]
                ell <- data[data$l==l_mode,]
                #ell <- ell[!duplicated(ell$n),]
                
                if (!grepl('solar-age', ev_stage) && counter == sun_num)
                    ref <- ell
                relation <- ell$nu ~ ell$n
                if (counter == 1 && l_mode == 0) {
                    make_layout(paste("Frequencies of", ev_stage, 
                                      "stars\nby", experiment_name))
                    set_par()
                    plot(relation, pch=l_mode, cex=0.25, lwd=0.5,
                         xlim=c(min(0, x_min-1), x_max), 
                         ylim=c(0, round(y_max+50, -2)),
                         yaxs='i', xaxs='i',
                         col=cl[counter], 
                         tck=0.01, 
                         xlab='',
                         ylab=bquote("frequency"~nu~"["*mu*Hz*"]"))
                    title(xlab=expression("radial order"~n))
                    magaxis(side=1:4, tcl=0.25, labels=FALSE)
                    legend("topleft", bty='n', inset=0, col=lgnd_cl, pch=20, 
                       legend=labls)
                    legend("bottomright", bty='n', inset=0, pch=0:3, 
                       legend=c(expression("\u2113"==0), 
                                expression("\u2113"==1),
                                expression("\u2113"==2),
                                expression("\u2113"==3)))
                } else {
                    points(relation, cex=0.25, pch=l_mode,
                           col=cl[counter], lwd=0.5)
                }
                counter <- counter + 1
            }
            dev.set(dev.prev())
            
            ######################################
            ### Make separated frequency plots ###
            ######################################
            counter <- 1
            for (model_i in data_files) {
                data <- cache[[model_i]]
                ell <- data[data$l==l_mode,]
                #ell <- ell[!duplicated(ell$n),]
                
                if (!grepl('solar-age', ev_stage) && counter == sun_num)
                    ref <- ell
                relation <- ell$nu ~ ell$n
                if (counter == 1) {
                    set_par()
                    plot(relation, pch=20, cex=0.01, lwd=0.5,
                         xlim=c(min(0, x_min-1), x_max), 
                         ylim=c(0, round(y_max+50, -2)),
                         yaxs='i', xaxs='i',
                         col=cl[counter], 
                         tck=0.01, 
                         xlab='',
                         ylab=bquote(.(l_name)~"frequency"~nu~"["*mu*Hz*"]"))
                    title(xlab=expression("radial order"~n))#, 
                          #mgp=par()$mgp-mgp_xoff)
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    points(relation, pch=20, cex=0.01, 
                           col=cl[counter], lwd=0.5)
                }
                counter <- counter + 1
            }
            if (grepl('solar-age', ev_stage)) {
                ell <- modelS[modelS$l==l_mode,]
                ell <- ell[ell$nu<y_max,]
                relation <- ell$nu ~ ell$n
                points(relation, pch=1, cex=0.5, col="black", lwd=0.5)
                ref <- ell
            }
            if (l_mode==1) make_legend(labls, ev_stage, lty=FALSE)
            dev.set(dev.prev())
            
            #############################
            ### Make difference plots ###
            #############################
            print(ref$n)
            n_range <- if (length(ref$n)>0) min(ref$n):max(ref$n) else c()
            #n_range <- n_range[n_range >= 5]
            diffs <- matrix(nrow=max(ref$n)-min(ref$n)+1, 
                            ncol=length(simulations))
            n_counter <- 1
            for (model_i in data_files) {
                data <- cache[[model_i]]
                ell <- data[data$l==l_mode,]
                #ell <- ell[!duplicated(ell$n),]
                
                x_counter <- 1
                for (x in n_range) {
                    diffs[x_counter, n_counter] <- 
                        if (x%in%ref$n && x%in%ell$n) {
                            a <- 2*pi*Mode(unlist(ell[ell$n==x,][3]))
                            b <- 2*pi*Mode(unlist(ref[ref$n==x,][3]))
                            (a-b)/b
                        } else NA
                    x_counter <- x_counter + 1
                }
                n_counter <- n_counter + 1
            }
            counter <- 1
            for (model_i in data_files) {
                if (counter == 1) {
                    set_par()
                    plot(n_range,#x_min:x_max, 
                         diffs[,counter], pch=20, cex=0.01,
                         xlim=range(n_range), 
                         ylim=range(diffs[complete.cases(diffs)]),
                         xaxs='i', lwd=0.5,  col=cl[counter], 
                         tck=0.01, xlab='',
                         ylab=if (l_mode==0) as.expression(
                                 bquote("rel. ang. freq. diff."~
                                        {delta * omega["\u2113"==.(l_mode)] /
                                                 omega["\u2113"==.(l_mode)]}))
                             else bquote({delta * omega["\u2113"==.(l_mode)] /
                                                  omega["\u2113"==.(l_mode)]}))
                    title(xlab=expression("radial order"~n))#, 
                          #mgp=par()$mgp-mgp_xoff)
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    points(n_range,
                           diffs[,counter], pch=20, cex=0.01, lwd=0.5,
                           col=cl[counter])
                }
                counter <- counter + 1
            }
            if (grepl('solar-age', ev_stage)) {
                ell <- modelS[modelS$l==l_mode,]
                relation <- rep(0, length(ell$n)) ~ ell$n
                points(relation, pch=1, cex=0.5, col="black", lwd=0.5)
            }
            if (l_mode==1) make_legend(labls, ev_stage, lty=FALSE)
            dev.set(dev.prev())
            
            ##########################
            ### Make Echelle plots ###
            ##########################
            counter <- 1
            for (model_i in data_files) {
                data <- cache[[model_i]]
                ell <- data[data$l==l_mode & data$n>0,]
                #ell <- ell[!duplicated(ell$n),]
                nu_max <- nu_maxs[[model_i]]
                
                converted_fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
                gaussian_env <- dnorm(ell$nu, nu_max, converted_fwhm)
                fit <- lm(ell$nu ~ ell$n, weights=gaussian_env)
                delta_nu <- coef(fit)[2]
                l_name <- toString(l_mode)
                delta_nus[[l_name]] <- c(delta_nus[[l_name]], delta_nu)
                deltas[[l_name]] <- c(deltas[[l_name]], 
                    summary(fit)$coefficients[2,2])
                
                ### estimate std errors
                out <- rep(0, n_samples)
                d <- if (l_mode == 0) 0.07
                else if (l_mode == 1) 0.08
                else if (l_mode == 2) 0.17
                else if (l_mode == 3) 0.38
                for (ii in 1:n_samples) {
                    new_fit <- lm(ell$nu + rnorm(nrow(ell), 0, d) ~ ell$n, 
                        weights=gaussian_env)
                    out[ii] <- summary(new_fit)$coefficients[2,2]
                }
                dnu_err[[l_name]] <- c(dnu_err[[l_name]],
                    mean(out))
                print(c(summary(fit)$coefficients[2,2], mean(out), max(out)))
                
                relation <- ell$nu ~ ell$nu%%delta_nu
                if (l_mode == 0 && counter == 1) {
                    make_layout(paste("Echelle diagram of", ev_stage, 
                                      "stars\nby", experiment_name), 
                                outside_legend=TRUE)
                    set_par()
                    plot(relation, pch=l_mode, 
                         col=adjustcolor(cl[counter], 0.5),
                         main="", cex=1000*gaussian_env,
                         ylim=c(round(p_nu_min-500,-3), 
                                round(p_nu_max+500, -3)),
                         xlim=c(0, round(max_delta_nu+5, -1)),
                         xlab="", xaxs='i', yaxs='i', tck=0.01, 
                         ylab=bquote("frequency"~nu~"["*mu*Hz*"]"))
                    title(xlab=expression(nu~mod~Delta*nu))#, 
                          #mgp=par()$mgp-mgp_xoff)
                    magaxis(side=1:4, tcl=0.25, labels=FALSE)
                    #minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    points(relation, pch=l_mode, cex=1000*gaussian_env,
                           col=adjustcolor(cl[counter], 
                               ifelse(counter==sun_num, 1, 0.5)))
                }
                abline(v=delta_nu, col=cl[counter], lty=l_mode+1, lwd=0.5)
                counter <- counter+1
            }
            if (l_mode==3) {
                par(mar=rep(0,4))
                plot.new()
                legend(-0.1, 0.9, bty='n', inset=0, col=lgnd_cl, pch=20, 
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
            dev.set(dev.prev())
            
            ###########################
            ### Make gradient plots ###
            ###########################
            alpha_cache <- list()
            new_ns_cache <- list()
            ns_cache <- list()
            log_dnus_cache <- list()
            ns_cache <- list()
            
            n_min <- Inf
            n_max <- Inf
            ii <- 1
            for (model_i in data_files) {
                data <- cache[[model_i]]
                ell <- data[data$l==l_mode & data$n>0,]
                #ell <- ell[!duplicated(ell$n),]
                ell <- ell[ell$nu>(nu_maxs[[model_i]] - 
                           5*delta_nus[[toString(l_mode)]][ii]),]
                ns <- head(ell$n, -1)
                if (max(ns) < n_max) n_max <- max(ns)
                if (min(ns) < n_max) n_min <- min(ns)
                ns_cache[[model_i]] <- ns
                ii <- ii + 1
            }
            
            grad_min <- Inf
            grad_max <- -Inf
            ii <- 1
            for (model_i in data_files) {
                data <- cache[[model_i]]
                ell <- data[data$l==l_mode & data$n>0,]
                #ell <- ell[!duplicated(ell$n),]
                ell <- ell[ell$nu>(nu_maxs[[model_i]] - 
                           5*delta_nus[[toString(l_mode)]][ii]),]
                log_dnus <- log(diff(ell$nu))
                ns <- head(ell$n, -1)
                new_ns <- seq(max(n_min, min(ns)), min(n_max, max(ns)), 0.01)
                alphas <- D1ss(ns, log_dnus, xout=new_ns)
                if (max(alphas) > grad_max) grad_max <- max(alphas)
                if (min(alphas) < grad_min) grad_min <- min(alphas)
                log_dnus_cache[[model_i]] <- log_dnus
                new_ns_cache[[model_i]] <- new_ns
                alpha_cache[[model_i]] <- alphas
                ii <- ii + 1
            }
            
            if (grepl('solar-age', ev_stage)) {
                ell <- modelS[modelS$l==l_mode,]
                ell <- ell[ell$nu<y_max & ell$n>=n_min & ell$n<=n_max,]
                #ell <- ell[!duplicated(ell$n),]
                modelS_log_dnus <- log(diff(ell$nu))
                modelS_ns <- head(ell$n, -1)
                new_modelS_ns <- seq(n_min, n_max, 0.01)
                modelS_alphas <- D1ss(modelS_ns, modelS_log_dnus, 
                    xout=new_modelS_ns)
                #if (max(modelS_alphas)>grad_max) grad_max <- max(modelS_alphas)
                #if (min(modelS_alphas)<grad_min) grad_min <- min(modelS_alphas)
            }
            
            counter <- 1
            for (model_i in data_files) {
                relation <- alpha_cache[[model_i]] ~ new_ns_cache[[model_i]]
                if (counter == 1) {
                    set_par()
                    plot(relation, col=adjustcolor(cl[counter], 0.3),
                         type='l',
                         ylim=c(grad_min, grad_max),
                         xlim=c(n_min, n_max),
                         main="", xlab="", cex=0.75, tck=0.01, 
                         ylab=if (l_mode == 0) bquote("gradient"~alpha[0]
                             =="("*d~log~Delta*nu/dn*")"["\u2113"==0])
                           else 
                             as.expression(bquote("gradient"~alpha[.(l_mode)])))
                    title(xlab=expression("radial order"~n))#, 
                          #mgp=par()$mgp-mgp_xoff)
                    minor.tick(nx=5, ny=5, tick.ratio=-0.15)
                } else {
                    lines(relation, col=adjustcolor(cl[counter], 0.3), 
                          cex=0.75, lwd=0.5)
                }
                ns <- ns_cache[[model_i]]
                log_dnus <- log_dnus_cache[[model_i]]
                points(D1ss(ns, log_dnus) ~ ns, col=cl[counter], pch=20, 
                       cex=0.5)
                counter <- counter+1
            }
            if (grepl('solar-age', ev_stage)) {
                lines(modelS_alphas ~ new_modelS_ns, lty=2, 
                      col=adjustcolor("black", 0.3))
                points(D1ss(modelS_ns, modelS_log_dnus) ~ modelS_ns)
            }
            if (l_mode==1) make_legend(labls, ev_stage)
            dev.set(dev.next())
            dev.set(dev.next())
            dev.set(dev.next())
            dev.set(dev.next())
        }
        dev.off()
        dev.off()
        dev.off()
        dev.off()
        dev.off()
        
        ################################
        ### Make plots of exp vs dnu ###
        ################################
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('dnu_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        make_layout(paste("Large separation of", ev_stage,
                    "stars\nby", experiment_name))
        set_par()
        max_err <- max(dnu_err[["0"]], dnu_err[["1"]], dnu_err[["2"]])/2
        for (l_mode in 0:2) {
            l_name <- toString(l_mode)
            relation <- delta_nus[[l_name]] ~ exp_vals
            if (l_mode == 0) {
                plot(relation, pch=l_mode, main="", xlab="", cex=0.5,
                     tck=0.01, col=ell_cl[l_mode+1],
                     ylim=range(delta_nus[["0"]], delta_nus[["1"]], 
                                delta_nus[["2"]]) + c(-max_err, max_err),
                     ylab=expression("large frequency separation"~Delta*nu))
                title(xlab=bquote(.(experiment_name)~.(freep)))#, 
                      #mgp=par()$mgp-mgp_xoff)
                lines(relation, lty=l_mode+1, col=ell_cl[l_mode+1])
                magaxis(side=1:4, tcl=0.25, labels=FALSE)
            } else {
                points(relation, pch=l_mode, col=ell_cl[l_mode+1], cex=0.5)
                lines(relation, lty=l_mode+1, col=ell_cl[l_mode+1])
            }
            y <- delta_nus[[l_name]][picker]
            yerr <- dnu_err[[l_name]][picker]/2
            errbar(exp_vals[picker], y, y - yerr, y + yerr,
                 add=TRUE, pch=l_mode, lty=5, cex=0.5,
                 col=ell_cl[l_mode+1], errbar.col=ell_cl[l_mode+1])
        }
        legend("bottom", bty='n', inset=0, pch=2:0, lty=3:1, col=ell_cl[3:1],
               legend=c(#expression("\u2113"==3), 
                        expression("\u2113"==2),
                        expression("\u2113"==1),
                        expression("\u2113"==0)))
        dev.off()
        
        #####################################
        ### Make plots of exp vs dnu diff ###
        #####################################
        cairo_pdf(file.path(plot_dir, ev_stage, 
            paste0('sepdiff_', basename(experiment), '_', ev_stage, '.pdf')), 
            width=plot_width, height=plot_height, family=font)
        make_layout(paste("Differences in Δν of", ev_stage,
                    "\nstars by degree and", experiment_name))
        set_par()
        diff_max <- -Inf
        diff_min <- Inf
        for (l_mode in 1:2) {
            y <- delta_nus[[toString(l_mode)]] - delta_nus[["0"]] + 
                max(dnu_err[[toString(l_mode)]])/2
            if (max(y) > diff_max) diff_max <- max(y)
            if (min(y) < diff_min) diff_min <- min(y)
        }
        print(c("diff min: ", diff_min))
        for (l_mode in 1:2) {
            l_name <- toString(l_mode)
            difference <- delta_nus[[l_name]] - delta_nus[["0"]]
            relation <- difference ~ exp_vals
            if (l_mode == 1) {
                plot(relation, pch=l_mode, main="", xlab="", yaxs='i',
                     col=ell_cl[l_mode+1], tck=0.01, cex=0.5,
                     ylim=c(min(0, round(diff_min-0.05, 1)), 
                            round(diff_max+0.05, 1)),
                     ylab=expression("large separation diff."
                                     ~delta*Delta*nu["\u2113"]==
                                     Delta*nu["\u2113"]-
                                     Delta*nu[0]))
                title(xlab=bquote(.(experiment_name)~.(freep)))#, 
                      #mgp=par()$mgp-mgp_xoff)
                lines(relation, lty=l_mode+1, col=ell_cl[l_mode+1])
                npoints <- length(exp_vals[picker])
                yerr <- dnu_err[["0"]][i_mid]/2
                errbar(exp_vals[picker], rep(0, npoints), -yerr, yerr, cex=0.5,
                    add=TRUE, lty=5, pch='-')
                abline(h=0)
                magaxis(side=1:4, tcl=0.25, labels=FALSE)
            } else {
                points(relation, pch=l_mode, col=ell_cl[l_mode+1], cex=0.5)
                lines(relation, lty=l_mode+1, col=ell_cl[l_mode+1])
            }
            y <- difference[picker]
            yerr <- dnu_err[[l_name]][picker]/2
            errbar(exp_vals[picker], y, y - yerr, y + yerr,
                 add=TRUE, pch=l_mode, lty=5, cex=0.5,
                 col=ell_cl[l_mode+1], errbar.col=ell_cl[l_mode+1])
        }
        legend("topleft", bty='n', inset=0, pch=2:1, lty=3:2, col=ell_cl[3:2],
            legend=c(#expression(delta*Delta*nu[3]), 
                     expression(delta*Delta*nu[2]),
                     expression(delta*Delta*nu[1])))
        dev.off()
        
    }
}
warnings()

