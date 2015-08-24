#### Plot frequency information for 16 Cyg A & B 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(stargazer)
library(RColorBrewer)
library(magicaxis)

source('utils.R')

plot_dir <- file.path('plots', 'cyg')
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

cl <- brewer.pal(4, "BrBG")#rainbow(4)[1:4]

for (cyg in c('16 Cyg A', '16 Cyg B')) {
    spaceless <- gsub(' ', '', cyg)
    data <- read.table(file.path('data', paste0(spaceless,'.dat')), header=TRUE)
    
    delta_nu <- ifelse(cyg=='16 Cyg A', 103.4, 116.97) # Lund et al 2014
    nu_max <- ifelse(cyg=='16 Cyg A', 2101, 2552) # Detection of l=4 and l=5...
    converted_fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
    gaussian_env <- dnorm(data$nu, nu_max, converted_fwhm)
    
    # start echelle device
    cairo_pdf(file.path(plot_dir, paste0('echelle_', spaceless, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    layout(matrix(c(1,2), ncol=1), heights=c(0.14,0.86))
    par(mar=rep(0,4))
    plot.new()
    text(0.5, 0.5, paste("Echelle diagram of", cyg), cex=2, font=2)
    par(las=1, mar=c(3, 3.4, 0.1, 1), mgp=c(2, 0.25, 0))
    
    # start frequency device
    cairo_pdf(file.path(plot_dir, paste0('nu_', spaceless, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    layout(matrix(c(1,2), ncol=1), heights=c(0.14,0.86))
    par(mar=rep(0,4))
    plot.new()
    text(0.5, 0.5, paste("Frequencies of", cyg), cex=2, font=2)
    par(las=1, mar=c(3, 3.4, 0.1, 1), mgp=c(2, 0.25, 0))#bty="l", 
    plot(data$nu ~ data$n, col=cl[data$l+1], pch=data$l, 
         main="", xlab="", tck=0, 
         cex=40*dnorm(data$nu, nu_max, converted_fwhm)/data$dnu,
         ylab=bquote("frequency"~nu~"["*mu*Hz*"]"))
    title(xlab=expression("radial order"~n), 
          mgp=par()$mgp-c(0.4,0,0))
    magaxis(side=1:4, labels=FALSE)
    abline(h=nu_max, lty=2)
    
    all_fit <- lm(data$nu ~ data$n, weights=gaussian_env/data$dnu)
    abline(all_fit)
    
    ells <- sapply(data$l+1, toString)
    print(anova(lm(data$nu ~ data$n, weights=gaussian_env/data$dnu),
                lm(data$nu ~ data$n + ells + data$n:ells, 
                   weights=gaussian_env/data$dnu)))
    
    ###############
    ### Echelle ###
    ###############
    delta_nus <- c()
    stderrors <- c()
    fits <- c()
    intercepts <- c()
    for (l_mode in 0:3) {
        ell <- data[data$l==l_mode,]
        #ell <- ell[ell$nu>(nu_max-5*delta_nu) & ell$nu<(nu_max+5*delta_nu),]
        attach(ell)
        gaussian_env <- dnorm(nu, nu_max, converted_fwhm)
        fit <- lm(nu ~ n, weights=gaussian_env/dnu)
        detach(ell)
        
        delta_nu <- coef(fit)[2]
        delta_nus <- c(delta_nus, delta_nu)
        stderrors <- c(stderrors, summary(fit)$coefficients[2,2])
        fits <- c(fits, fit)
        intercepts <- c(intercepts, coef(fit)[1])
        if (l_mode == 0) first_fit <- fit
        
        abline(fit, lty=l_mode+3, col=cl[l_mode+1])
        dev.set(dev.prev())
        
        abline(h=nu_max, lty=2)
        
        relation <- ell$nu ~ ell$nu%%delta_nu
        if (l_mode == 0) {
            plot(relation, pch=l_mode, col=cl[l_mode+1],
                 main="", 
                 cex=50*gaussian_env/ell$dnu,
                 ylim=range(data$nu),#c(1500, 3500),
                 xlim=c(0, 125),
                 xlab="", xaxs='i', tck=0.01, #yaxs='i', 
                 ylab=bquote("frequency"~nu~"["*mu*Hz*"]"))
            title(xlab=expression(nu~mod~Delta*nu), 
                  mgp=par()$mgp-c(0.4,0,0))
            magaxis(side=1:4, tcl=0.25, labels=FALSE)
            #abline(v=coef(all_fit)[2], lty=1, lwd=0.5)
            legend(3, par("usr")[4], bty='n', inset=0, 
                   pch=c(0:3, rep(NA, 6)), 
                   lty=c(rep(NA, 5), 2, 3:6),
                   col=c(cl, "black", "black", cl),
                   legend=c(expression("\u2113"==0), 
                            expression("\u2113"==1),
                            expression("\u2113"==2),
                            expression("\u2113"==3),
                            "",
                            expression(nu[max]),
                            expression(Delta*nu["\u2113"==0]), 
                            expression(Delta*nu["\u2113"==1]),
                            expression(Delta*nu["\u2113"==2]),
                            expression(Delta*nu["\u2113"==3])))
        } else {
            points(relation, col=cl[l_mode+1], pch=l_mode, 
                   cex=50*gaussian_env/ell$dnu)
        }
        #points(ell$nu ~ ell$nu%%coef(all_fit)[2], 
        #       col=rgb(0,0,0,0.5), pch=l_mode+3, 
        #       cex=50*gaussian_env/ell$dnu)
        abline(v=delta_nu, lty=l_mode+3, lwd=0.5, col=cl[l_mode+1])
        dev.set(dev.next())
    }
    
    legend("topleft", bty='n', inset=0,
           pch=0:3, col=cl,
           legend=c(expression("\u2113" == 0),
                    expression("\u2113" == 1),
                    expression("\u2113" == 2),
                    expression("\u2113" == 3)))
    
    dnus <- sapply(0:3, function(l) as.expression(bquote(
        Delta*nu[.(l)] == 
            .(sprintf("%.2f", round(delta_nus[l+1],2)))~"\U00b1"~
            .(sprintf("%.2f", round(stderrors[l+1],2))))))
    print(summary(all_fit)$coefficients[2,2])
    legend("bottomright", bty='n', inset=0,
           lty=c(2, 1, 3:6),
           col=c("black", "black", cl),
           legend=c(bquote(nu[max]==.(nu_max)~mu*Hz), 
           as.expression(bquote(Delta*nu == 
            .(sprintf("%.2f", round(coef(all_fit)[2],2)))~"\U00b1"~
            .(sprintf("%.2f", round(summary(all_fit)$coefficients[2,2],2))))),
           dnus))
    
    dev.off()
    dev.off()
    stargazer(data.frame(ell=0:3, dnu=delta_nus, std=stderrors),
              title=cyg, type="text", summary=FALSE, rownames=FALSE)
    
    # start difference device
    cairo_pdf(file.path(plot_dir, paste0('diff_', spaceless, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    layout(matrix(c(1,2), ncol=1), heights=c(0.14,0.86))
    par(mar=rep(0,4))
    plot.new()
    text(0.5, 0.5, paste("Residuals of", cyg), cex=2, font=2)
    par(las=1, mar=c(3, 3.4, 0.1, 1), mgp=c(2, 0.25, 0))
    #fit <- fits[1]
    df <- data.frame()
    for (l_mode in 0:3) {
        ell <- data[data$l == l_mode,]
        nus <- ell$nu - predict(first_fit, newdata=data.frame(n=ell$n))
        df <- rbind(df, data.frame(l=l_mode, n=ell$n, nu=nus, 
            nu_orig=ell$nu, dnu=ell$dnu))
    }
    plot(df$nu ~ df$n, col=cl[df$l+1], pch=df$l, 
         main="", xlab="", tck=0, 
         cex=100*dnorm(df$nu_orig, nu_max, converted_fwhm)/df$dnu,
         ylab=bquote("frequency residuals"~
             nu["\u2113,"*n]-Delta*nu["0,"*n]~"["*mu*Hz*"]"))
    title(xlab=expression("radial order"~n), 
          mgp=par()$mgp-c(0.4,0,0))
    magaxis(side=1:4, tcl=0.25, labels=FALSE)
    for (l_mode in 0:3) {
        ell <- df[df$l==l_mode,]
        w <- dnorm(ell$nu_orig, nu_max, converted_fwhm)/ell$dnu
        fit <- lm(ell$nu ~ ell$n, weights=w)
        #abline(fit, lty=2)
        abline(h=weighted.mean(ell$nu, w), lty=2)
        ci <- predict(fit, interval="confidence")
        lines(ci[,2] ~ ell$n, lty=3)
        lines(ci[,3] ~ ell$n, lty=3)
    }
    print(c(cl, rep(0, 4)))
    legend("center", bty='n', inset=0, ncol=2, 
           pch=c(3:0, rep(NA, 4)), 
           col=c(rev(cl), rep("black", 4)), 
           lty=c(rep(0, 5), 2, 3, 0),
           legend=c(expression("\u2113" == 3),
                    expression("\u2113" == 2),
                    expression("\u2113" == 1),
                    expression("\u2113" == 0),
                    "", 
                    "weighted mean", 
                    "95% confidence interval",
                    ""))
    dev.off()
}
