library(magicaxis)
library(RColorBrewer)
library(akima)

font <- 'Palatino'

sim_dir <- 'deleter'
simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]
sims.DF <- do.call(rbind, Map(function (f) read.table(f, header=1), simulations))
sims.DF <- sims.DF[sims.DF$dnu02_median > 0 & sims.DF$dnu02_median < 15,]
write.table(sims.DF, file.path('grids', 'deleter.dat'), quote=FALSE, sep='\t', row.names=FALSE)
combos <- unique(sims.DF[,1:4])

# L_lim <- c(0)
# T_lim <- c()
# Dnu_lim <- c()
# dnu_lim <- c()
# Dnus <- c()
# dnus <- c()
# r02 <- c()
# r01 <- c()
# r13 <- c()
# ages <- c()
# max_R <- 0
# max_M <- 0
# for (simulation in simulations) {
    # DF <- read.table(simulation, header=1)
    # DF <- DF[DF$dnu02_median > 0 & DF$dnu02_median < 15,]
    # L_lim <- range(L_lim, DF$L)
    # T_lim <- range(T_lim, DF$Teff)
    # Dnu_lim <- range(Dnu_lim, DF$Dnu_median)
    # dnu_lim <- range(dnu_lim, DF$dnu02_median)
    # Dnus <- c(Dnus, DF$Dnu_median)
    # dnus <- c(dnus, DF$dnu02_median)
    # r02 <- c(r02, DF$r_sep02_median)
    # r01 <- c(r01, DF$r_avg01_median)
    # r13 <- c(r13, DF$r_sep13_median)
    # ages <- c(ages, DF$age)
    # if (max(DF$radius) > max_R) max_R <- max(DF$radius)
    # if (max(DF$mass) > max_M) max_M <- max(DF$mass)
# }

cairo_pdf('HR.pdf', width=4, height=3, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
#for (simulation_i in 1:length(simulations)) {
    #DF <- read.table(simulations[simulation_i], header=1)
    #DF <- DF[DF$dnu02_median > 0 & DF$dnu02_median < 15,]
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(combos[simulation_i,], sims.DF)
    HR <- DF$L ~ DF$Teff
    if (simulation_i == 1) {
        plot(HR, type='l', tcl=0,
            #ylim=L_lim, xlim=rev(T_lim),
            ylim=range(sims.DF$L, 0), xlim=rev(range(sims.DF$Teff)),
            xlab=expression(T[eff]), xaxt='n', yaxt='n',
            ylab=expression(L / L['\u0298']))
        abline(v=5777, lty=3, col='lightgray')
        abline(h=1, lty=3, col='lightgray')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0))
    } else {
        lines(HR)
    }
    points(HR, pch=1, lwd=0.1, 
        col=brewer.pal(11,"Spectral")[floor(DF$age/13.9*11)+1],
        cex=0.5*DF$radius/max(sims.DF$radius))
}
dev.off()

cairo_pdf('JCD.pdf', width=4, height=3, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
#for (simulation_i in 1:length(simulations)) {
#    DF <- read.table(simulations[simulation_i], header=1)
#    DF <- DF[DF$dnu02_median > 0 & DF$dnu02_median < 15,]
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(combos[simulation_i,], sims.DF)
    JCD <- DF$dnu02_median ~ DF$Dnu_median
    jcdcol <- brewer.pal(11,"Spectral")[floor(DF$age/13.9*11)+1]
    jcdcex <- 0.5*DF$mass/max(sims.DF$mass)#max_M
    if (simulation_i == 1) {
        plot(JCD, pch=20, lwd=0.1, 
            col=jcdcol,
            cex=jcdcex, tcl=0,
            #ylim=dnu_lim, xlim=Dnu_lim,
            ylim=range(sims.DF$dnu02_median), xlim=range(sims.DF$Dnu_median),
            xaxt='n', yaxt='n',
            xlab=expression(Delta*nu~"["*mu*Hz*"]"),
            ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"))
        #abline(h=9, lty=3, col='lightgray')
        #abline(v=136, lty=3, col='lightgray')
        points(135, 9, pch=1, cex=1)
        points(135, 9, pch=20, cex=0.1)
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0))
    } else {
        points(JCD, col=jcdcol, pch=20, cex=jcdcex)
    }
    
}
dev.off()

cairo_pdf('JCD-mesh.pdf', width=8, height=6, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
#my.matrix  <- interp(Dnus,dnus,ages)
my.matrix  <- interp(sims.DF$Dnu_median, sims.DF$dnu02_median, sims.DF$age)
#ind.mat.na <- which(is.na(c(my.matrix$z)))
#my.matrix$z[ind.mat.na] <- 0
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }
filled.contour(my.matrix,
    #color=function(x) { brewer.pal(11,"Spectral")[x] },
    color=my.heat.colors,
    key.axes={
        axis(4, cex.axis=1, tcl=0)
        mtext("Age [Gyr]", side=4, las=3, line=3, cex=2)
    },
    plot.axes={
        #abline(v=136, lty=3, col='black')
        #abline(h=9, lty=3, col='black')
        contour(my.matrix, add=TRUE, nlevels=13, labcex=1)
        points(135, 9, pch=1, cex=1)
        points(135, 9, pch=20, cex=0.1)
        magaxis(side=1:4, tcl=0.5, labels=c(1,1,0,0), cex.lab=0.5, family='Palatino')
    },
    plot.title={
        title(xlab=expression(Delta*nu~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
    })
dev.off()


cairo_pdf('ratios-mesh.pdf', width=8, height=6, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
#my.matrix  <- interp(r02,r01,ages)
my.matrix  <- interp(sims.DF$r_sep02_median, sims.DF$r_avg01_median, sims.DF$age)
#ind.mat.na <- which(is.na(c(my.matrix$z)))
#my.matrix$z[ind.mat.na] <- 0
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }
filled.contour(my.matrix,
    #color=function(x) { brewer.pal(11,"Spectral")[x] },
    color=my.heat.colors,
    key.axes={
        axis(4, cex.axis=1, tcl=0)
        mtext("Age [Gyr]", side=4, las=3, line=3, cex=2)
    },
    plot.axes={
        #abline(v=136, lty=3, col='black')
        #abline(h=9, lty=3, col='black')
        points(6.642727e-02, 2.238694e-02, pch=1, cex=1)
        points(6.642727e-02, 2.238694e-02, pch=20, cex=0.1)
        contour(my.matrix, add=TRUE, nlevels=13, labcex=1)
        magaxis(side=1:4, tcl=0.5, labels=c(1,1,0,0), cex.lab=0.5, family='Palatino')
    },
    plot.title={
        title(xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(r[0*","*1]~"["*mu*Hz*"]"), cex.lab=2, line=3)
    })
dev.off()


#corrs <- do.call(cbind, Map(function (f) { a <- read.table(f, header=1); cor(a, a$age) } , simulations))

