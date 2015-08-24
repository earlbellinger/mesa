#### Utility file for HR.R and fgong.R
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(Hmisc)

#options(warn=2)

minor.tick <- function (nx = 2, ny = 2, tick.ratio = 0.5) {
    ax <- function(w, n, tick.ratio) {
        range <- par("usr")[if (w == "x") 
            1:2
        else 3:4]
        tick.pos <- if (w == "x") 
            par("xaxp")
        else par("yaxp")
        distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
        possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
        low.candidates <- possible.minors >= range[1]
        low.minor <- if (any(low.candidates)) {
            min(possible.minors[low.candidates])
        } else {
            tick.pos[1]
        }
        possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
        hi.candidates <- possible.minors <= range[2]
        hi.minor <- if (any(hi.candidates)) {
            max(possible.minors[hi.candidates])
        } else {
            tick.pos[2]
        }
        axis(if (w == "x") 
            1
        else 2, seq(low.minor, hi.minor, by = distance.between.minor), 
            labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1) 
        ax("x", nx, tick.ratio = tick.ratio)
    if (ny > 1) 
        ax("y", ny, tick.ratio = tick.ratio)
    invisible()
}

Teff_sun = log10(5777)
solar_age = 4.57e9
solar_radius = 6.955*10**10
solar_mass = 1.9891*10**30
solar_scale = sqrt(solar_mass/solar_radius^3)
font = "Palatino"
approx_xout = seq(0, 1, .04)
color_offset = 6
plot_width = 6
plot_height = 6
layout_width_three = c(.425,.425,.15)
layout_width_two = c(.85, .15)
cex_lab = 1.3
title_cex = 1.3
par_mar = c(3.5, 4.5, 0.1, 1)
par_mgp = c(2.5, 0.25, 0)
mgp_xoff = c(0.5,0,0)
layout_heights = c(0.14,0.86)

set_par <- function() 
    par(bty="l", las=1, mar=par_mar, cex.lab=cex_lab, mgp=par_mgp)

make_layout <- function(title_text, outside_legend=FALSE) {
    if (outside_legend)
        layout(matrix(c(1,1,2,3), ncol=2, byrow=TRUE), 
               heights=layout_heights, 
               widths=layout_width_two)
    else
        layout(matrix(c(1,2), ncol=1, byrow=TRUE), 
               heights=layout_heights)
    make_title(title_text, 
        title.cex=ifelse(outside_legend, 1.2*title_cex, title_cex))
}

make_title <- function(title_text, title.cex=title_cex) {
    par(mar=rep(0,4))
    plot.new()
    text(0.5, 0.5, title_text, cex=title.cex, font=2)
}

make_legend <- function(labls, ev_stage, position='left', lty=TRUE) {
    par(mar=rep(0,4))
    plot.new()
    cls <- if(ev_stage=="solar-age") c(lgnd_cl, "black") else lgnd_cl
    lgnd <- if(ev_stage=="solar-age") c(labls, "Model S") else labls
    if (lty) {
        legend(position, bty='n', inset=0, col=cls, legend=lgnd,
            lty=if(ev_stage=="solar-age") 
                c(rep(1, length(labls)), 2) else 1)
    } else {
        legend(position, bty='n', inset=0, col=cls, legend=lgnd,
            pch=if(ev_stage=="solar-age") 
                c(rep(20, length(labls)), 1) else 20)
    }
}

start_dev <- function(maintext, plot_name, experiment, ev_stage, 
                      width=layout_width_three) {
    cairo_pdf(file.path(plot_subdir, 
        paste0(plot_name, '_', basename(experiment), '_', ev_stage, '.pdf')),
              width=plot_width, height=plot_height, family=font)
    layout(matrix(c(1,1,1,2,3,4,5,6,4), ncol=3, byrow=TRUE), 
           heights=c(0.14,0.43,0.43), widths=width)
    #layout(matrix(c(1,2,3,4,5,3), ncol=3, byrow=TRUE), 
    #       heights=c(0.5,0.5), widths=width)
    make_title(paste(maintext, sub("_", " ", ev_stage), 
                     "stars\nby", experiment_name),
               title.cex=title_cex*1.5)
}

load_experiment_info <- function(experiment, simulations) {
    print(experiment)
    print(simulations)
    if (grepl('Y', experiment)) {
        experiment_name <<- 'helium'
        sun_num <<- grep('.280', simulations)[1]
    } else if (grepl('alpha', experiment)) {
        experiment_name <<- 'mixing length'
        sun_num <<- grep('2.10', simulations)[1]
    } else if (grepl('M', experiment)) {
        experiment_name <<- 'mass'
        sun_num <<- grep('1.000', simulations)[1]
    }
    
    freep <<- as.name(basename(experiment))
    cl <<- heat.colors(length(simulations)+color_offset)[1:length(simulations)]
    cl <<- c(cl[1:(sun_num-1)], "black", cl[(sun_num+1):length(cl)])
    lgnd_cl <<- cl
    exp_vals <<- as.numeric(basename(simulations))
    labls <<- sapply(simulations,
        function(x) { as.expression(bquote(.(freep) ~ "=" ~ .(basename(x)))) })
    
    while (length(lgnd_cl) > 11) {
        lgnd_cl <<- lgnd_cl[c(TRUE, FALSE)]
        labls <<- labls[c(TRUE, FALSE)]
    }
    
    print(freep)
    print(labls)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
