library(magicaxis)
library(RColorBrewer)
library(akima)
library(parallel)
library(parallelMap)
library(data.table)
library(lattice)

font <- 'Palatino'

plot_dir <- file.path('plots')
hrcdr_dir <- file.path(plot_dir, 'hrcdr')
dir.create(plot_dir, showWarnings=FALSE)
dir.create(hrcdr_dir, showWarnings=FALSE)

sim_dir <- 'deleter'
simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

plot_width = 6.97522
plot_height = 4.17309

# Load data
seis.DF <- data.table(do.call(rbind, 
    Map(function (f) read.table(f, header=1), simulations)))
setkey(seis.DF, M, Y, Z, alpha)
keys <- key(seis.DF)

# Remove outliers
print("Before outlier removal:")
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

parallelStartMulticore(max(1, detectCores()))
get_outliers <- function(DF) {
    unique(unlist(parallelMap(
        function(ii) which(DF[[ii]] %in% boxplot.stats(DF[[ii]], coef=10)$out),
        names(DF))))
}
repeat {
    outliers <- get_outliers(seis.DF)
    if (length(outliers) <= 0) break
    print(paste("Removing", length(outliers), "models"))
    seis.DF <- seis.DF[-outliers,]
    
    singletons <- c()
    for (ii in which(!duplicated(seis.DF[, keys, with=0])))
        if (!(ii+1 < nrow(seis.DF) && 
              all(seis.DF[ii, keys, with=0]==seis.DF[ii+1, keys, with=0])
            | ii-1 > 0 && 
              all(seis.DF[ii, keys, with=0]==seis.DF[ii-1, keys, with=0])))
                singletons <- c(singletons, ii)
    if (!is.null(singletons)) {
        print(paste("Removing", length(singletons), "singletons"))
        seis.DF <- seis.DF[-singletons,]
    }
}

print("After outlier removal:")
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, file.path('grids', 'deleter.dat'), quote=FALSE, 
    sep='\t', row.names=FALSE)

solar_vals <- read.table(file.path('perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- parallelMap(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos))
combos <- combos[order(unlist(ages)),]

# Make inputs diagram
cairo_pdf(file.path(plot_dir, 'inputs.pdf'), 
          width=plot_width, height=plot_width,#height=plot_height*2, 
          family=font)
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0, 0, 0, 0))
splom(combos, pch=20, cex=0.05, 
      xlab=NULL, ylab=NULL, 
      axis.text.cex=1,
      axis.text.lineheight=0.1,
      axis.line.tck=1,
      varname.cex=2,
      varnames=c(expression(M[0]), expression(Y[0]), expression(Z[0]), 
               expression(alpha["MLT"])))
#pairs(combos, upper.panel=NULL, 
#      pch=20, cex=0.25, tcl=-0.25, gap=1.5,
#      cex.labels=2, font.labels=2,
#      labels=c(expression(M[0]), expression(Y[0]), expression(Z[0]), 
#               expression(alpha["MLT"])))
dev.off()

# HR Diagram
png(file.path(plot_dir, 'HR.png'), family=font, res=400,
    width=plot_width*250, height=plot_height*250)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(seis.DF, combos[simulation_i,])
    HR <- log10(DF$L) ~ log10(DF$Teff)
    hrcex <- 0.1*DF$M/max(seis.DF$M)
    hrcol <- brewer.pal(10,"Spectral")[floor(DF$age/13.9*10)+1]
    if (simulation_i == 1) {
        plot(HR, #type='l', 
            tcl=0, pch=20, cex=hrcex, col=hrcol, #lwd=0.1,
            ylim=range(log10(seis.DF$L), 0, 1), 
            xlim=rev(log10(range(seis.DF$Teff))),
            xlab=expression(T[eff]), xaxt='n', yaxt='n',
            ylab=expression(L / L['\u0298']))
        abline(v=log10(5777), lty=3, col='black')
        abline(h=log10(1), lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                unlog='xy', mgp=c(2, 0.25, 0), prettybase=100)
    } else {
        points(HR, pch=20, col=hrcol, cex=hrcex)
        #lines(HR)
    }
}
points(log10(5777), log10(1), pch=1, cex=1)
points(log10(5777), log10(1), pch=20, cex=0.1)
dev.off()


scatter_mesh <- function(plotname, var1, var2, var3, 
        label1, label2, label3) {
    #my.matrix <- interp(seis.DF[[var1]], seis.DF[[var2]], seis.DF[[var3]])
    cairo_pdf(file.path(hrcdr_dir, paste0('mesh-', plotname, '-', var3, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1)
    filled.contour(my.matrix,
        levels=sort(unique(round(seis.DF[[third]], 1)))[c(TRUE, FALSE, FALSE)], 
        color=colorRampPalette(brewer.pal(10, "Spectral")),
        key.axes={
            axis(4, cex.axis=1.5, tcl=0, line=0)
            mtext(label3, side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(my.matrix, add=TRUE, labcex=0.5, levels=levels[[third]])
            points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.5, 0), cex.axis=1.5)
        },
        plot.title={
            title(xlab=label1, cex.lab=2, line=3)
            title(ylab=label2, cex.lab=2, line=3)
        })
    dev.off()
    
    
    varmax <- max(seis.DF[[third]])
    varmin <- min(seis.DF[[third]])
    png(file.path(hrcdr_dir, paste0(plotname, '-', var3, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        relation <- DF[[var2]] ~ DF[[var1]]
        color <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*9)+1]
        cex <- 0.01
        if (simulation_i == 1) {
            plot(relation, 
                pch=20, xaxt='n', yaxt='n',
                col=color, cex=cex, tcl=0,
                ylim=range(seis.DF[[var2]]), 
                xlim=range(seis.DF[[var1]]),
                xlab=label1, ylab=label2)
            abline(v=solar_vals[[var1]], lty=3, col='black')
            abline(h=solar_vals[[var2]], lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(relation, col=color, pch=20, cex=cex)
        }
    }
    points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
    points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
    dev.off()
}



## Color by age, mass, helium, metallicity, mix length, and core hydrogen
col.pal <- brewer.pal(10, "Spectral")
thirds <- c('age', 'M', 'Y', 'Z', 'alpha')
labels <- c('Age [Gyr]', expression(M[0]), expression(Y[0]), expression(Z[0]), 
            expression(alpha["MLT"]))
levels <- list(age=0:14)



scatter_mesh('JCD', 'Dnu_median', 'dnu02_median', 'age', 
    expression(Delta*nu~"["*mu*Hz*"]"),
    expression(delta*nu[0*","*2]~"["*mu*Hz*"]"),
    labels[which(third==thirds)])




for (third in thirds) {
    var1 <- 'Dnu_median'
    var2 <- 'dnu02_median'
    my.matrix <- interp(seis.DF[[var1]], seis.DF[[var2]], seis.DF[[third]])
    cairo_pdf(file.path(hrcdr_dir, paste0('mesh-JCD-', third, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1)
    filled.contour(my.matrix,
        levels=sort(unique(round(seis.DF[[third]], 1)))[c(TRUE, FALSE, FALSE)], 
        color=colorRampPalette(brewer.pal(10, "Spectral")),
        key.axes={
            axis(4, cex.axis=1.5, tcl=0, line=0)
            mtext(labels[which(third==thirds)], side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(my.matrix, add=TRUE, labcex=0.5, levels=levels[[third]])
            points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.5, 0), cex.axis=1.5)
        },
        plot.title={
            title(xlab=expression(Delta*nu~"["*mu*Hz*"]"), cex.lab=2, line=3)
            title(ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"), 
                  cex.lab=2, line=3)
        })
    dev.off()
    
    
    cairo_pdf(file.path(hrcdr_dir, paste0('mesh-r_sep-', third, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
    var1 <- 'r_sep02_median'
    var2 <- 'r_sep13_median'
    my.matrix <- interp(seis.DF[[var1]], seis.DF[[var2]], seis.DF[[third]])
    filled.contour(my.matrix,
        color=function(x) { rev(heat.colors(x, alpha=1)) },
        key.axes={
            axis(4, cex.axis=1.5, tcl=0)
            mtext(labels[which(third==thirds)], side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(my.matrix, add=TRUE, labcex=0.5)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.5, 0), cex.axis=1.5)
        },
        plot.title={
            title(xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
            title(ylab=expression(r[1*","*3]~"["*mu*Hz*"]"), cex.lab=2, 
                  line=3)
        })
    dev.off()
    
    
    cairo_pdf(file.path(hrcdr_dir, paste0('mesh-r_avg01-', third, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
    var1 <- 'Dnu_median'
    var2 <- 'r_avg01_median'
    my.matrix <- interp(seis.DF[[var1]], seis.DF[[var2]], seis.DF[[third]])
    filled.contour(my.matrix,
        color=function(x) { rev(heat.colors(x, alpha=1)) },
        key.axes={
            axis(4, cex.axis=1.5, tcl=0)
            mtext(labels[which(third==thirds)], side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(my.matrix, add=TRUE, labcex=0.5)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.5, 0), cex.axis=1.5)
        },
        plot.title={
            title(xlab=expression(Delta*nu~"["*mu*Hz*"]"), cex.lab=2, line=3)
            title(ylab=expression(r[0*","*1]~"["*mu*Hz*"]"), cex.lab=2, 
                  line=3)
        })
    dev.off()
    
    
    cairo_pdf(file.path(hrcdr_dir, paste0('mesh-HR-', third, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
    var1 <- 'Teff'
    var2 <- 'L'
    my.matrix <- interp(seis.DF[[var1]], seis.DF[[var2]], 
        seis.DF[[third]])
    filled.contour(my.matrix,
        xlim=rev(range(seis.DF[[var1]])),
        color=function(x) { rev(heat.colors(x, alpha=1)) },
        key.axes={
            axis(4, cex.axis=1.5, tcl=0)
            mtext(labels[which(third==thirds)], side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(my.matrix, add=TRUE, labcex=0.5)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
            points(solar_vals[[var1]], solar_vals[[var2]], pch=20, 
                cex=0.1)
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.5, 0), cex.axis=1.5)
        },
        plot.title={
            title(xlab=expression(T[eff]), cex.lab=2, line=3)
            title(ylab=expression(L / L['\u0298']), cex.lab=2, 
                  line=3)
        })
    dev.off()
}


for (third in thirds) {
    varmax <- max(seis.DF[[third]])
    varmin <- min(seis.DF[[third]])
    
    ## Ratio separations
    png(file.path(hrcdr_dir, paste0('r_sep-', third, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        var1 <- 'r_sep02_median'
        var2 <- 'r_sep13_median'
        fmla <- DF[[var2]] ~ DF[[var1]]
        cols <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*9)+1]
        cexs <- 0.01#*DF$M/max(seis.DF$M)#max_M
        if (simulation_i == 1) {
            plot(fmla, pch=20, 
                col=cols, cex=cexs, tcl=0, 
                xlim=range(seis.DF[[var1]]),
                ylim=range(seis.DF[[var2]]),
                xaxt='n', yaxt='n',
                xlab=expression(r[0*","*2]~"["*mu*Hz*"]"),
                ylab=expression(r[1*","*3]~"["*mu*Hz*"]"))
            abline(v=solar_vals[[var1]], lty=3, col='black')
            abline(h=solar_vals[[var2]], lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(fmla, col=cols, pch=20, cex=cexs)
        }
    }
    points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
    points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
    dev.off()
    
    ## Ratio averages
    png(file.path(hrcdr_dir, paste0('r_avg10-', third, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        var1 <- 'Dnu_median'
        var2 <- 'r_avg10_median'
        fmla <- DF[[var2]] ~ DF[[var1]]
        cols <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*9)+1]
        cexs <- 0.01#*DF$M/max(seis.DF$M)#max_M
        if (simulation_i == 1) {
            plot(fmla, pch=20, 
                col=cols, cex=cexs, tcl=0,
                xlim=range(seis.DF[[var1]]),
                ylim=range(seis.DF[[var2]]),
                xaxt='n', yaxt='n',
                xlab=expression(Delta*nu~"["*mu*Hz*"]"),
                ylab=expression(r[1*","*0]~"["*mu*Hz*"]"))
            abline(v=solar_vals[[var1]], lty=3, col='black')
            abline(h=solar_vals[[var2]], lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(fmla, col=cols, pch=20, cex=cexs)
        }
    }
    points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
    points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
    dev.off()
    
    ## Ratio averages
    png(file.path(hrcdr_dir, paste0('r_avg01-', third, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        var1 <- 'Dnu_median'
        var2 <- 'r_avg01_median'
        fmla <- DF[[var2]] ~ DF[[var1]]
        cols <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*9)+1]
        cexs <- 0.01#*DF$M/max(seis.DF$M)#max_M
        if (simulation_i == 1) {
            plot(fmla, pch=20, 
                col=cols, cex=cexs, tcl=0,
                xlim=range(seis.DF[[var1]]),
                ylim=range(seis.DF[[var2]]),
                xaxt='n', yaxt='n',
                xlab=expression(Delta*nu~"["*mu*Hz*"]"),
                ylab=expression(r[0*","*1]~"["*mu*Hz*"]"))
            abline(v=solar_vals[[var1]], lty=3, col='black')
            abline(h=solar_vals[[var2]], lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(fmla, col=cols, pch=20, cex=cexs)
        }
    }
    points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
    points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
    dev.off()
    
    ## HR
    png(file.path(hrcdr_dir, paste0('HR-', third, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        HR <- DF$L ~ DF$Teff
        hrcex <- 0.01
        hrcol <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*9)+1]
        if (simulation_i == 1) {
            plot(HR, #type='l', 
                tcl=0, pch=20, cex=hrcex, col=hrcol, 
                ylim=range(seis.DF$L, 0), xlim=rev(range(seis.DF$Teff)),
                xlab=expression(T[eff]), xaxt='n', yaxt='n',
                ylab=expression(L / L['\u0298']))
            abline(v=5777, lty=3, col='black')
            abline(h=1, lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(HR, pch=20, col=hrcol, cex=hrcex)
        }
    }
    points(5777, 1, pch=1, cex=1)
    points(5777, 1, pch=20, cex=0.1)
    dev.off()
    
    ## JCD
    png(file.path(hrcdr_dir, paste0('JCD-', third, '.png')), 
        family=font, res=400, width=plot_width*500, height=plot_height*500)
    par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        JCD <- DF$dnu02_median ~ DF$Dnu_median
        jcdcol <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*9)+1]
        jcdcex <- 0.01#*DF$M/max(seis.DF$M)#max_M
        if (simulation_i == 1) {
            plot(JCD, pch=20, 
                col=jcdcol,
                cex=jcdcex, tcl=0,
                ylim=range(seis.DF$dnu02_median), 
                xlim=range(seis.DF$Dnu_median),
                xaxt='n', yaxt='n',
                xlab=expression(Delta*nu~"["*mu*Hz*"]"),
                ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"))
            abline(v=solar_vals$Dnu_median, lty=3, col='black')
            abline(h=solar_vals$dnu02_median, lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(JCD, col=jcdcol, pch=20, cex=jcdcex)
        }
    }
    points(solar_vals$Dnu_median, solar_vals$dnu02_median, pch=1, cex=1)
    points(solar_vals$Dnu_median, solar_vals$dnu02_median, pch=20, cex=0.1)
    dev.off()
}




if (FALSE) {


# JCD Diagram
png(file.path(plot_dir, 'JCD.png'), family=font, res=400,
    width=plot_width*250, height=plot_height*250)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(seis.DF, combos[simulation_i,])
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
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                mgp=c(2, 0.25, 0))
    } else {
        points(JCD, col=jcdcol, pch=20, cex=jcdcex)
    }
}
points(135, 9, pch=1, cex=1)
points(135, 9, pch=20, cex=0.1)
dev.off()


# Ratio Diagram
png(file.path(plot_dir, 'RD.png'), width=400, height=300, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(seis.DF, combos[simulation_i,])
    JCD <- DF$r_sep02_median ~ DF$r_sep13_median
    jcdcol <- brewer.pal(10,"Spectral")[floor(DF$age/13.9*10)+1]
    jcdcex <- 0.1*DF$M/max(seis.DF$M)#max_M
    if (simulation_i == 1) {
        plot(JCD, pch=20, 
            col=jcdcol,
            cex=jcdcex, tcl=0,
            ylim=range(seis.DF$r_sep02_median), 
            xlim=range(seis.DF$r_sep13_median),
            xaxt='n', yaxt='n',
            xlab=expression(r[0*","*2]~"["*mu*Hz*"]"),
            ylab=expression(r[1*","*3]~"["*mu*Hz*"]"))
        abline(v=0.118065474714876, lty=3, col='black')
        abline(h=0.0664272689314098, lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                mgp=c(2, 0.25, 0))
    } else {
        points(JCD, col=jcdcol, pch=20, cex=jcdcex)
    }
}
points(0.118065474714876, 0.0664272689314098, pch=1, cex=1)
points(0.118065474714876, 0.0664272689314098, pch=20, cex=0.1)
dev.off()


# Ratio Diagram 2
png(file.path(plot_dir, 'RD2.png'), width=400, height=300, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(seis.DF, combos[simulation_i,])
    JCD <- DF$r_avg10_median ~ DF$r_avg01_median
    jcdcol <- brewer.pal(10,"Spectral")[floor(DF$age/13.9*10)+1]
    jcdcex <- 0.1*DF$M/max(seis.DF$M)#max_M
    if (simulation_i == 1) {
        plot(JCD, pch=20, 
            col=jcdcol,
            cex=jcdcex, tcl=0,
            ylim=range(seis.DF$r_avg10_median), 
            xlim=range(seis.DF$r_avg01_median),
            xaxt='n', yaxt='n',
            xlab=expression(r[0*","*1]~"["*mu*Hz*"]"),
            ylab=expression(r[1*","*0]~"["*mu*Hz*"]"))
        abline(v=0.0223869425875988, lty=3, col='black')
        abline(h=0.0224873559599195, lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                mgp=c(2, 0.25, 0))
    } else {
        points(JCD, col=jcdcol, pch=20, cex=jcdcex)
    }
}
points(0.0223869425875988, 0.0224873559599195, pch=1, cex=1)
points(0.0223869425875988, 0.0224873559599195, pch=20, cex=0.1)
dev.off()




# JCD Mesh
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
        contour(my.matrix, add=TRUE, labcex=0.5)
        points(135, 9, pch=1, cex=1)
        points(135, 9, pch=20, cex=0.1)
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                mgp=c(2, 0.25, 0), family='Palatino')
    },
    plot.title={
        title(xlab=expression(Delta*nu~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(delta*nu[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, 
              line=3)
    })
dev.off()

# Ratios Mesh
cairo_pdf(file.path(plot_dir, 'ratios-mesh.pdf'), 
    width=8, height=6, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1.3)
my.matrix  <- interp(seis.DF$r_sep02_median, seis.DF$r_sep13_median, 
    seis.DF$age)
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }
filled.contour(my.matrix,
    color=my.heat.colors,
    key.axes={
        axis(4, cex.axis=1, tcl=0)
        mtext("Age [Gyr]", side=4, las=3, line=3, cex=2)
    },
    plot.axes={
        points(0.118065474714876, 0.0664272689314098, pch=1, cex=1)
        points(0.118065474714876, 0.0664272689314098, pch=20, cex=0.1)
        contour(my.matrix, add=TRUE, labcex=0.5)
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                mgp=c(2, 0.25, 0), family='Palatino')
    },
    plot.title={
        title(xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(r[1*","*3]~"["*mu*Hz*"]"), cex.lab=2, line=3)
    })
dev.off()


# Ratio-dnu Mesh
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
        contour(my.matrix, add=TRUE, labcex=0.5)
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                mgp=c(2, 0.25, 0), family='Palatino')
    },
    plot.title={
        title(xlab=expression(r[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
        title(ylab=expression(dnu[0*","*2]~"["*mu*Hz*"]"), cex.lab=2, line=3)
    })
dev.off()


predict(lm(age ~ poly(r_sep02_median,6)*poly(r_avg01_median,6), data=seis.DF), 
        data.frame(r_sep02_median=6.642727e-02, r_avg01_median=2.238694e-02), 
        interval="predict", level=0.68)
#corrs <- do.call(cbind, Map(function (f) { 
#    a <- read.table(f, header=1); cor(a, a$age) } , simulations))

}
