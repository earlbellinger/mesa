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

solar_vals <- read.table(file.path('perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

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

labs <- expression(M, Y[0], Z[0], alpha["MLT"], tau, "mass", R, 
    H[c], "X(He)", log~g, L, T["eff"], "Fe"/"H", 
    "<"*Delta*nu*">", "<"*d*Delta*nu/d*nu*">", #"<"*Delta*nu^b*">", 
    "<"*Delta*nu[0]*">", "<"*d*Delta*nu[0]/d*nu*">", #"<"*Delta*nu[0]^b*">", 
    "<"*delta*nu[0*","*2]*">", "<"*d*delta*nu[0*","*2]/d*nu*">", 
        #"<"*delta*nu[0*","*2]^b*">", 
    "<"*r[0*","*2]*">", "<"*d*r[0*","*2]/d*nu*">", #"<"*r[0*","*2]^b*">", 
    "<"*r[0*","*1]*">", "<"*d*r[0*","*1]/d*nu*">", #"<"*r[0*","*1]^b*">", 
    "<"*delta*nu[1*","*3]*">", "<"*d*delta*nu[1*","*3]/d*nu*">", 
        #"<"*delta*nu[1*","*3]^b*">", 
    "<"*r[1*","*3]*">", "<"*d*r[1*","*3]/d*nu*">",
    "<"*r[1*","*0]*">", "<"*d*r[1*","*0]/d*nu*">"#, "<"*r[0*","*1]^b*">", 
)

latex_labs <- c("M", "$Y_0$", "$Z_0$", "$\\alpha_{\\text{\"MLT\"}}$", 
    "$\\tau$", "mass", "R", "$H_c$", "X(He)", "$\\log g$", "L", 
    "$T_{\text{\"eff\"}}$", "Fe/H", 
    
    "$\\langle\\Delta\\nu\\rangle$", 
    "$\\langle\\frac{d\\Delta\\nu}{d\nu}\\rangle$", 
    
    "$\\langle\\Delta\\nu_0\\rangle$", 
    "$\\langle\\frac{d\\Delta\\nu_0}{d\nu}\\rangle$", 
    
    "$\\langle\\delta\\nu_{02}\\rangle$", 
    "$\\langle\\frac{d\\delta\\nu_{02}}{d\nu}\\rangle$", 
    
    "$\\langle r_{02}\\rangle$", 
    "$\\langle\\frac{dr_{02}}{d\nu}\\rangle$", 
    
    "$\\langle r_{01}\\rangle$", 
    "$\\langle\\frac{dr_{01}}{d\nu}\\rangle$", 
    
    "$\\langle\\delta\\nu_{13}\\rangle$", 
    "$\\langle\\frac{d\\delta\\nu_{13}}{d\nu}\\rangle$", 
    
    "$\\langle r_{13}\\rangle$", 
    "$\\langle\\frac{dr_{13}}{d\nu}\\rangle$", 
    
    "$\\langle r_{10}\\rangle$", 
    "$\\langle\\frac{dr_{10}}{d\nu}\\rangle$"
)

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
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

exclude <- which(grepl("mass|Dnu_", names(seis.DF)))
seis.DF <- seis.DF[,-exclude, with=0]

# Save data
write.table(seis.DF, file.path('grids', 'deleter.dat'), quote=FALSE, 
    sep='\t', row.names=FALSE)

## Plot histograms
#tmp <- data.frame(seis.DF[,1:8, with=0])
#colnames(tmp) <- labs[-exclude][1:8]
#d <- melt(tmp)
#ggplot(d, aes(x = value)) +#, y=..density..)) +
#    geom_histogram(aes(y=..ncount..), fill="#c0392b", alpha=0.75) + 
#    fte_theme() + 
#    geom_density(aes(y = ..scaled..)) +#col=2) + 
#    scale_y_continuous(labels=comma) + 
#    geom_hline(yintercept=0, size=0.4, color="black") +
#    facet_wrap(~variable,scales = "free_x", nrow=4) +
#    ggtitle(labs[-exclude][1:8])

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(parallelMap(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]


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

#png(file.path(plot_dir, 'inputs-legend.png'), res=400, 
#    width=150*plot_width/8, height=150*plot_width, 
#    family=font)
#par(mar=c(0, 0, 0, 0), mgp=c(0, 0, 0), oma=c(0, 0, 0, 0))
#color.legend(par()$usr[2], par()$usr[1], par()$usr[4], par()$usr[3], 
#             signif(quantile(seq(varmin, varmax, length=1000), 
#                    c(0.05, 0.275, 0.5, 0.725, 0.95)), 3), 
#             col.pal[1:length(col.pal)], gradient='y', align='rb')
#mtext(expression(H_0), 4, line=4.5, cex=1.3)
#dev.off()


# HR scatter
third <- 'M'
varmax <- round(max(seis.DF[[third]]), 2)
varmin <- round(min(seis.DF[[third]]), 2)
png(file.path(plot_dir, 'HR-M.png'), 
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
            xlab=expression(T["eff"]/K), 
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
        c(0, 0.25, 0.5, 0.75, 1)), 3), 
    col.pal[1:length(col.pal)], gradient='y', align='rb')
mtext(expression(M/M['\u0298']), 4, line=4.5, cex=1.3)
dev.off()




# HR mesh
mesh <- with(seis.DF, interp(Teff, log10(L), M,
    xo=seq(min(Teff), max(Teff), length=40),
    yo=seq(log10(min(L)), log10(max(L)), length=40)))
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
        color <- col.pal[floor((DF[[third]]-varmin)/(varmax-varmin)*
            length(col.pal))+1]
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
