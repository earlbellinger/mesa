library(magicaxis)
library(RColorBrewer)
library(akima)
library(parallel)
library(parallelMap)
library(data.table)

font <- 'Palatino'

plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

sim_dir <- 'deleter'
simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

parallelStartMulticore(max(1, detectCores()))
seis.DF <- do.call(rbind, 
    Map(function (f) read.table(f, header=1), simulations))

print("Before outlier removal:")
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

get_outliers <- function(DF) {
    unique(unlist(parallelMap(
        function(ii) which(DF[,ii] %in% boxplot.stats(DF[,ii], coef=10)$out),
        1:ncol(DF))))
}
outliers <- get_outliers(seis.DF)
while (length(outliers) > 0) {
    seis.DF <- seis.DF[!(1:nrow(seis.DF)) %in% outliers,]
    singletons <- c()
    for (ii in which(!duplicated(seis.DF[,1:4])))
        if (!(ii+1 < nrow(seis.DF) && all(seis.DF[ii,1:4] == seis.DF[ii+1,1:4])
            | ii-1 > 0 && all(seis.DF[ii,1:4] == seis.DF[ii-1,1:4])))
                singletons <- c(singletons, ii)
    if (!is.null(singletons)) {
        print(paste("Removing", length(singletons), "singletons"))
        seis.DF <- seis.DF[-singletons,]
    }
    outliers <- get_outliers(seis.DF)
}

print("After outlier removal:")
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

write.table(seis.DF, file.path('grids', 'deleter.dat'), quote=FALSE, 
    sep='\t', row.names=FALSE)
combos <- unique(seis.DF[,1:4])
#ages <- sapply(1:nrow(combos), function(i) max(merge(combos[i,], seis.DF)$age))
ages <- parallelMap(function(i) max(merge(combos[i,], seis.DF)$age), 
    1:nrow(combos))
combos <- combos[order(unlist(ages)),]


png(file.path(plot_dir, 'HR.png'), width=600, height=450, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(combos[simulation_i,], seis.DF)
    HR <- DF$L ~ log10(DF$Teff)
    hrcex <- 0.1*DF$M/max(seis.DF$M)
    hrcol <- brewer.pal(10,"Spectral")[floor(DF$age/13.9*10)+1]
    if (simulation_i == 1) {
        plot(HR, #type='l', 
            tcl=0, pch=20, cex=hrcex, col=hrcol, #lwd=0.1,
            ylim=range(seis.DF$L, 0), xlim=rev(log10(range(seis.DF$Teff))),
            xlab=expression(log[10]~T[eff]), xaxt='n', yaxt='n',
            ylab=expression(log[10]~L / L['\u0298']))
        abline(v=5777, lty=3, col='black')
        abline(h=1, lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0))
    } else {
        points(HR, pch=20, col=hrcol, cex=hrcex)
        #lines(HR)
    }
}
points(1, log10(5777), pch=1, cex=1)
points(1, log10(5777), pch=20, cex=0.1)
dev.off()


png(file.path(plot_dir, 'JCD.png'), width=600, height=450, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(combos[simulation_i,], seis.DF)
    JCD <- DF$dnu02_median ~ DF$Dnu_median
    jcdcol <- brewer.pal(10,"Spectral")[floor(DF$age/13.9*10)+1]
    jcdcex <- 0.1*DF$M/max(seis.DF$M)#max_M
    if (simulation_i == 1) {
        plot(JCD, pch=20, 
            col=jcdcol,
            cex=jcdcex, tcl=0,
            ylim=range(seis.DF$dnu02_median), xlim=range(seis.DF$Dnu_median),
            xaxt='n', yaxt='n',
            xlab=expression(Delta*nu~"["*mu*Hz*"]"),
            ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"))
        abline(v=135, lty=3, col='black')
        abline(h=9, lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0))
    } else {
        points(JCD, col=jcdcol, pch=20, cex=jcdcex)
    }
}
points(135, 9, pch=1, cex=1)
points(135, 9, pch=20, cex=0.1)
dev.off()


cairo_pdf(file.path(plot_dir, 'JCD-mesh.pdf'), width=8, height=6, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
my.matrix  <- interp(seis.DF$Dnu_median, seis.DF$dnu02_median, seis.DF$age)
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }
filled.contour(my.matrix,
    color=my.heat.colors,
    key.axes={
        axis(4, cex.axis=1, tcl=0)
        mtext("Age [Gyr]", side=4, las=3, line=3, cex=2)
    },
    plot.axes={
        contour(my.matrix, add=TRUE, nlevels=13, labcex=1)
        points(135, 9, pch=1, cex=1)
        points(135, 9, pch=20, cex=0.1)
        magaxis(side=1:4, tcl=0.5, labels=c(1,1,0,0), cex.lab=0.5, 
                family='Palatino')
    },
    plot.title={
        title(xlab=expression(Delta*nu~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, 
              line=3)
    })
dev.off()


cairo_pdf(file.path(plot_dir, 'ratios-mesh.pdf'), 
    width=8, height=6, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
my.matrix  <- interp(seis.DF$r_sep02_median, seis.DF$r_avg01_median, 
    seis.DF$age)
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }
filled.contour(my.matrix,
    color=my.heat.colors,
    key.axes={
        axis(4, cex.axis=1, tcl=0)
        mtext("Age [Gyr]", side=4, las=3, line=3, cex=2)
    },
    plot.axes={
        points(6.642727e-02, 2.238694e-02, pch=1, cex=1)
        points(6.642727e-02, 2.238694e-02, pch=20, cex=0.1)
        contour(my.matrix, add=TRUE, nlevels=13, labcex=1)
        magaxis(side=1:4, tcl=0.5, labels=c(1,1,0,0), cex.lab=0.5, 
                family='Palatino')
    },
    plot.title={
        title(xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(r[0*","*1]~"["*mu*Hz*"]"), cex.lab=2, line=3)
    })
dev.off()


cairo_pdf(file.path(plot_dir, 'rcd-mesh.pdf'), width=8, height=6, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
my.matrix  <- interp(seis.DF$r_sep02_median, seis.DF$dnu02_median, seis.DF$age)
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }
filled.contour(my.matrix,
    color=my.heat.colors,
    key.axes={
        axis(4, cex.axis=1, tcl=0)
        mtext("Age [Gyr]", side=4, las=3, line=3, cex=2)
    },
    plot.axes={
        points(0.06642727, 9, pch=1, cex=1)
        points(0.06642727, 9, pch=20, cex=0.1)
        contour(my.matrix, add=TRUE, nlevels=13, labcex=1)
        magaxis(side=1:4, tcl=0.5, labels=c(1,1,0,0), cex.lab=0.5, 
                family='Palatino')
    },
    plot.title={
        title(xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(dnu[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
    })
dev.off()


cairo_pdf(file.path(plot_dir, 'r02-age.pdf'), width=4, height=3, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
plot(seis.DF$age ~ seis.DF$r_sep02_median, pch=3, tcl=0, cex=0.1,
     xaxt='n', yaxt='n',
     xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), ylab='Age [Gyr]')
magaxis(side=1:4, tcl=0.5, labels=c(1,1,0,0), family='Palatino')
dev.off()


predict(lm(age ~ poly(r_sep02_median,6)*poly(r_avg01_median,6), data=seis.DF), 
        data.frame(r_sep02_median=6.642727e-02, r_avg01_median=2.238694e-02), 
        interval="predict", level=0.68)
#corrs <- do.call(cbind, Map(function (f) { 
#    a <- read.table(f, header=1); cor(a, a$age) } , simulations))
