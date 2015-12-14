library(magicaxis)
library(RColorBrewer)
library(akima)
library(parallel)
library(parallelMap)
library(data.table)
library(lattice)
library(lpSolve)

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
load_data <- function(f) {
    DF <- read.table(f, header=1)
    n <- 100
    if (nrow(DF) < n) return(NULL)
    N <- length(DF$Hc)
    ideal <- seq(max(DF$Hc), min(DF$Hc), length=n)
    cost.mat  <- outer(ideal, DF$Hc, function(x, y) abs(x-y))
    row.signs <- rep("==", n)
    row.rhs   <- rep(1, n)
    col.signs <- rep("<=", N)
    col.rhs   <- rep(1, N)
    sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
        col.signs, col.rhs)$solution
    DF[apply(sol, 1, which.max),]
}
seis.DF <- data.table(do.call(rbind, Map(load_data, simulations)))
setkey(seis.DF, M, Y, Z, alpha)
keys <- key(seis.DF)

# Remove outliers
print("Before outlier removal:")
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

parallelStartMulticore(max(1, detectCores()))
repeat {
    outliers <- unique(unlist(parallelMap(
        function(ii) which(
            seis.DF[[ii]] %in% boxplot.stats(seis.DF[[ii]], coef=100)$out),
        names(seis.DF))))
    if (length(outliers) <= 0) break
    print(paste("Removing", length(outliers), "models"))
    seis.DF <- seis.DF[-outliers,]
    
    combos <- unique(seis.DF[,keys, with=0])
    remove_list <- unlist(parallelMap(function(i) nrow(merge(seis.DF, combos[i,])), 
        1:nrow(combos))) < 20
    
    if (sum(remove_list) > 0) {
        print(paste("Removing", sum(remove_list), "singletons"))
        seis.DF <- seis.DF[!combos[remove_list]]
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

# Plot histograms
d <- melt(seis.DF[,1:8, with=0])
ggplot(d,aes(x = value)) +
    geom_histogram() + 
    facet_wrap(~variable,scales = "free_x", nrow=2)

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(parallelMap(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]

col.pal <- colorRampPalette(brewer.pal(11, "Spectral"))(1000)

thirds <- c('age', 'M', 'Y', 'Z', 'alpha', 'He', 'Hc')
labels <- c(
    'Age [Gyr]', 
    expression(M/M["\u0298"]), 
    expression(Y[0]), 
    expression(Z[0]),
    expression(alpha["MLT"]), 
    expression(X(He)), 
    expression(H[c])
)
levels <- list(
    age=0:14, 
    M=seq(0.7, 1.3, 0.1),
    Y=seq(0.22, 0.34, 0.01),
    Z=log10(seq(10**1e-04, 10**0.04, length=10)),
    alpha=seq(1.5, 2.5, 0.1),
    He=seq(0.22, 0.45, 0.02),
    Hc=seq(0, 0.78, 0.05)
)
color_levels <- list(
    age=seq(0, 13.8, 0.25),
    M=seq(0.7, 1.3, 0.025),
    Y=seq(0.22, 0.34, 0.0025),
    Z=log10(seq(10**1e-04, 10**0.04, length=20)),
    alpha=seq(1.5, 2.5, 0.025),
    He=seq(0.22, 0.45, 0.005),
    Hc=seq(0, 0.78, 0.015)
)


# Make inputs diagram
#cairo_pdf(file.path(plot_dir, 'inputs.pdf'), 
#          width=plot_width+0.25*plot_width, height=plot_width, 
#          family=font)
png(file.path(plot_dir, 'inputs.png'), res=400, 
    width=150*plot_width, height=150*plot_width, 
    family=font)
par(mar=c(0, 0, 0, 0), mgp=c(0, 0, 0), oma=c(0, 0, 0, 0))
H <- 1-combos$Y-combos$Z
varmax <- max(H)
varmin <- min(H)
cols <- col.pal[floor((H-varmin) / (varmax-varmin) * (length(col.pal)-1))+1]
splom(combos, cex=0.001, pch=3,
      col=cols,
      #col.pal[floor(ages/max(ages)*length(col.pal))],
      xlab=NULL, ylab=NULL, 
      axis.text.cex=0.25,
      axis.text.lineheight=0.0001,
      axis.line.tck=0.25,
      xaxs='n', yaxs='n',
      varname.cex=0.5,
      varnames=c(expression(M[0]), expression(Y[0]), expression(Z[0]), 
               expression(alpha["MLT"])))
dev.off()

png(file.path(plot_dir, 'inputs-legend.png'), res=400, 
    width=150*plot_width/8, height=150*plot_width, 
    family=font)
par(mar=c(0, 0, 0, 0), mgp=c(0, 0, 0), oma=c(0, 0, 0, 0))
color.legend(par()$usr[2], par()$usr[1], par()$usr[4], par()$usr[3], 
             signif(quantile(seq(varmin, varmax, length=1000), 
                    c(0.05, 0.275, 0.5, 0.725, 0.95)), 3), 
             col.pal[1:length(col.pal)], gradient='y', align='rb')
mtext(expression(H_0), 4, line=4.5, cex=1.3)
dev.off()


# HR scatter
third <- 'M'
varmax <- round(max(seis.DF[[third]]), 2)
varmin <- round(min(seis.DF[[third]]), 2)
png(file.path(plot_dir, 'HR-M-linear.png'), 
    family=font, res=400, width=plot_width*250, height=plot_height*250)
par(mar=c(3, 4, 1, 6), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(seis.DF, combos[simulation_i,])
    relation <- log10(DF$L) ~ DF$Teff
    color <- col.pal[
        floor((DF[[third]]-varmin)/(varmax-varmin)*length(col.pal))+1]
    cex <- 0.01
    if (simulation_i == 1) {
        plot(relation, 
            pch=1, axes=FALSE,
            col=color, cex=cex, tcl=0,
            ylim=range(log10(seis.DF$L)),#, 0, 1),
            xlim=rev(range(seis.DF$Teff)), 
            xlab=expression(T["eff"]~"["*K*"]"), 
            ylab=expression(L / L['\u0298']))
        abline(v=5777, lty=3, col='black')
        abline(h=0, lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                majorn=c(4, 3, 4, 3),
                unlog='y', mgp=c(2, 0.25, 0))
    } else {
        points(relation, col=color, pch=20, cex=cex)
    }
}
points(5777, 0, pch=1, cex=1)
points(5777, 0, pch=20, cex=0.1)
var1range <- diff(par()$usr)[1]
color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
             par()$usr[2]+0.10*var1range, par()$usr[4], 
    signif(quantile(seq(varmin, varmax, length=1000), 
        c(0.05, 0.275, 0.5, 0.725, 0.95)), 3), 
    col.pal[1:length(col.pal)], gradient='y', align='rb')
mtext(expression(M/M['\u0298']), 4, line=4.5, cex=1.3)
dev.off()




# HR mesh
mesh <- interp(seis.DF$Teff, log10(seis.DF$L), seis.DF$M,
    xo=seq(min(seis.DF$Teff), max(seis.DF$Teff), length=100),
    yo=seq(log10(min(seis.DF$L)), log10(max(seis.DF$L)), length=100))
cairo_pdf(file.path(plot_dir, 'mesh-HR-M-linear.pdf'), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1)
filled.contour(mesh,
    ylim=range(log10(seis.DF$L)),
    xlim=rev(range(seis.DF$Teff)), 
    xaxs='r', yaxs='r',
    levels=color_levels[['M']], 
    color=colorRampPalette(brewer.pal(11, "Spectral")),
    key.axes={
        axis(4, cex.axis=1.5, tcl=0, line=0)
        mtext(expression(M/M['\u0298']), side=4, las=3, line=4, cex=2)
    },
    plot.axes={
        contour(mesh, add=TRUE, labcex=1, levels=levels[['M']],
           method="simple")
        points(5777, 0, pch=1, cex=1)
        points(5777, 0, pch=20, cex=0.1)
        abline(v=5777, lty=3, col=adjustcolor('black', alpha.f=0.25))
        abline(h=0, lty=3, col=adjustcolor('black', alpha.f=0.25))
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                majorn=c(4, 3, 4, 3), cex.axis=1.5,
                unlog='y')#, mgp=c(2, 0.25, 0))
    },
    plot.title={
        title(xlab=expression(T["eff"]~"["*K*"]"), cex.lab=2, line=3)
        title(ylab=expression(L / L['\u0298']), cex.lab=2, line=3)
    })
dev.off()




## Color by age, mass, Y0, X(He), metallicity, mix length, and core hydrogen
scatter_mesh <- function(plotname, var1, var2, var3, label1, label2, label3) { 
    # scatter 
    varmax <- max(seis.DF[[third]])
    varmin <- min(seis.DF[[third]])
    png(file.path(hrcdr_dir, paste0(plotname, '-', var3, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 6), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        relation <- DF[[var2]] ~ DF[[var1]]
        color <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*length(col.pal))+1]
        cex <- 0.01
        if (simulation_i == 1) {
            plot(relation, 
                pch=20, axes=FALSE,
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
    var1max <- max(seis.DF[[var1]])
    var1range <- diff(par()$usr)[1]
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4], 
        signif(quantile(seq(varmin, varmax, length=1000), 
            c(0.05, 0.275, 0.5, 0.725, 0.95)), 2), 
        col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(label3, 4, line=4.5, cex=1.3)
    dev.off()
    
       
    # mesh
    mesh <- interp(seis.DF[[var1]], seis.DF[[var2]], seis.DF[[var3]])
    cairo_pdf(file.path(hrcdr_dir, 
            paste0('mesh-', plotname, '-', var3, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1)
    filled.contour(mesh,
        levels=color_levels[[third]], 
        color=colorRampPalette(brewer.pal(11, "Spectral")),
        key.axes={
            axis(4, cex.axis=1.5, tcl=0, line=0)
            mtext(label3, side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(mesh, add=TRUE, labcex=0.5, levels=levels[[third]])
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
}

for (third in thirds) {
    scatter_mesh('JCD', 'Dnu_median', 'dnu02_median', third, 
        expression(Delta*nu~"["*mu*Hz*"]"), 
        expression(delta*nu[0*","*2]~"["*mu*Hz*"]"), 
        labels[which(third==thirds)])
    
    #scatter_mesh('HR', 'Teff', 'L', third, 
    #    expression(T[eff]~"["*K*"]"), 
    #    expression(L/L['\u0298']), 
    #    labels[which(third==thirds)])
}


if (FALSE) {

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

}
