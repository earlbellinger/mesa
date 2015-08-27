#### Plot Lamb frequencies and large frequency separations of Model S
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(magicaxis)
library(RColorBrewer)
library(Bolstad)

source('utils.R')

plot_dir <- file.path('plots', 'modelS')
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

#cl <- brewer.pal(3, "Dark2")
cl <- brewer.pal(4, "BrBG")

solar_radius <- 6.955*10**10
modelS <- read.table(file.path('data', 'modelS.dat'), header=TRUE) 
modelS <- modelS[modelS$radius > 0,]
attach(modelS)

modelS_freqs <- read.table(file.path('data', 'fgong.l5bi.d.dat'))
modelS_freqs <- data.frame(l=modelS_freqs[,1], n=modelS_freqs[,2], 
    nu=modelS_freqs[,3])
attach(modelS_freqs)

s <- function(ell) 
    10**6/(2*pi) * sqrt(ell*(ell+1) * csound**2 / (radius*solar_radius)**2)

#############################
### Plot Lamb frequencies ###
#############################
cairo_pdf(file.path(plot_dir, paste0('lamb_modelS.pdf')), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(3, 3.5, 3, 1), mgp=c(2, 0.5, 0))
plot(radius, log10(s(1)), type='l', lty=2, 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     main="Lamb frequencies of Model S", 
     xlab=expression(fractional~radius~r/R), 
     ylab=expression(frequency~nu~"["*mu*Hz*"]"),
     xlim=c(0, 1), ylim=c(min(log10(s(1))), max(log10(s(3)))))
magaxis(side=1:4, family=font, unlog='y', tcl=0.25,
        mgp=c(1.8, 0.25, 0), labels=c(TRUE, TRUE, FALSE, FALSE))
dnu <- list()
for (ell in 1:3) {
    ess <- log10(s(ell))
    lines(radius, ess, lty=ell+1)
    n_is <- which(l==ell)
    linept <- if (ell == 3) min(n_is)
        else if (ell == 2) ceil(mean(n_is))
        else if (ell == 1) max(n_is)
    for (n_i in n_is) {
        freq <- log10(nu[n_i])
        rt <- radius[max(which(freq > ess))]
        dnu[[toString(ell)]] <- c(dnu[[toString(ell)]], 
            1/(2*sintegral(solar_radius*radius[radius>=rt], 
               1/(csound[radius>=rt]*10**6))$value))
        points(rt, freq, col=cl[ell+1], pch=ell, cex=0.5)
        if (n_i == linept) 
            segments(rt, freq, 1, freq, lty=ell+1, col=cl[ell+1])
    }
}
legend("bottomleft", bty='n', lty=3:1+1,
       legend=c(expression(Lamb~S[3]),
                expression(Lamb~S[2]),
                expression(Lamb~S[1])))
legend("topright", bty='n', col=rev(cl[1:3+1]), pch=3:1,
       legend=c(expression("\u2113"==3),
                expression("\u2113"==2),
                expression("\u2113"==1)))
dev.off()

########################################
### Plot large frequency separations ###
########################################
cairo_pdf(file.path(plot_dir, paste0('dnu_integral_modelS.pdf')), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(7, 5, 3, 1), mgp=c(5, 0.75, 0))
boxplot(dnu, horizontal=TRUE, pch=1:3, col=cl[2:4], ylim=c(136, max(dnu$`3`)),
        main="Large frequency separations of Model S", 
        xlab=expression(Delta*nu["\u2113,"*n] == 
            bgroup("(", 2*integral(frac(dr, c(r)), r[t], R), ")")^-1),
        ylab="")
title(ylab=expression("harmonic degree"~"\u2113"), mgp=par()$mgp-c(2,0,0))
abline(v=135, lty=2)
dnu0 <- 1/(2*sintegral(solar_radius*radius, 1/(csound*10**6))$value)
abline(v=dnu0, lty=3)
legend("bottomright", bty='n', lty=3:2, 
       legend=c(bquote(Delta*nu[0] == .(round(dnu0, 2)) ~ mu*Hz),
                expression(Delta*nu["⊙"] %~~% 135 ~ mu*Hz)))
dev.off()

cairo_pdf(file.path(plot_dir, paste0('dnu_pairwise_modelS.pdf')), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(7, 5, 3, 1), mgp=c(5, 0.75, 0))
pairwise_dnu = list(`0`=diff(nu[l==0 & n > 0]),
                    `1`=diff(nu[l==1 & n > 0]),
                    `2`=diff(nu[l==2 & n > 0]),
                    `3`=diff(nu[l==3 & n > 0]))
boxplot(pairwise_dnu, horizontal=TRUE, pch=0:3, 
        col=cl, #ylim=c(136, max(pairwise_dnu$`3`)),
        main="Large frequency separations of Model S", 
        xlab=expression(Delta*nu["\u2113,"*n] == 
             nu["\u2113,"*n+1] - nu["\u2113,"*n+1]),
        ylab="")
title(ylab=expression("harmonic degree"~"\u2113"), mgp=par()$mgp-c(2,0,0))
abline(v=135, lty=2)
legend("bottomright", bty='n', lty=2, 
       legend=expression(Delta*nu["⊙"] %~~% 135 ~ mu*Hz))
dev.off()

cairo_pdf(file.path(plot_dir, paste0('dnu_both_modelS.pdf')), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(7, 5, 3, 1), mgp=c(5, 0.75, 0))
boxplot(pairwise_dnu, horizontal=TRUE, pch=1, outcol=cl[1],
        col=adjustcolor(cl[1], 0.5), ylim=c(130, 180),
        main="Large frequency separations of Model S", 
        xlab=expression(Delta*nu["\u2113,"*n]), ylab="")
boxplot(c(dnu0, dnu), horizontal=TRUE, pch=3, add=TRUE, outcol=cl[3],
        col=adjustcolor(cl[3], 0.5))
title(ylab=expression("harmonic degree"~"\u2113"), mgp=par()$mgp-c(2,0,0))
abline(v=135, lty=2)
legend("bottomright", bty='n', lty=c(NA, NA, 2), pch=c(1, 3, NA),
       col=c(cl[1], cl[3], 1),
       legend=c(expression("pairwise"~Delta*nu),
                expression("integral"~Delta*nu),
                expression(Delta*nu["⊙"] %~~% 135 ~ mu*Hz)))
dev.off()

#data <- approx(radius, log(s(1)), n=500)
#plot(radius, log(s(1)), type='l')
#f <- function(r, b0, b1, b2, b3, b4, b5) {
#    b0 + b1*log(r) + b2*exp(b3*r) + b4*r
#}
#init.fit <- lm(y ~ log(x)+exp(x)+x, data=data)
#lines(radius, predict(init.fit, newdata=data.frame(x=radius)), 
#       col='blue', lty=3)
#coefs <- coef(init.fit)
#fit <- nls(y ~ f(x, b0, b1, b2, b3, b4), data=data, 
#           control=nls.control(maxiter=10000, tol=1e-6, minFactor=1e-10),
#           start=list(b0=coefs[1], b1=coefs[2], b2=coefs[3], b3=1, 
#                      b4=coefs[4]))
#lines(radius, predict(fit, newdata=data.frame(x=radius)), col='red', lty=2)

