#### Utility file for HR.R and fgong.R
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(Hmisc)

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
font <- "Palatino"
approx_xout <- seq(0, 1, .04)
color_offset <- 6

make_legend <- function(labls, ev_stage, position='left', lty=TRUE) {
    par(mar=rep(0,4))
    plot.new()
    cls <- if(ev_stage=="solar-like") c(cl, "black") else cl
    lgnd <- if(ev_stage=="solar-like") c(labls, "Model S") else labls
    if (lty) {
        legend(position, bty='n', inset=0, col=cls, legend=lgnd,
            lty=if(ev_stage=="solar-like") 
                c(rep(1, length(labls)), 2) else 1)
    } else {
        legend(position, bty='n', inset=0, col=cls, legend=lgnd,
            pch=if(ev_stage=="solar-like") 
                c(rep(20, length(labls)), 1) else 20)
    }
}

start_dev <- function(maintext, fname, exp, ev_stage, width=c(.425,.425,.15)) {
    cairo_pdf(file.path(plot_subdir, 
        paste0(ev_stage, '-', exp, '_', fname, '.pdf')),
              width=6, height=6, family=font)
    layout(matrix(c(1,1,1,2,3,4,5,6,4), ncol=3, byrow=TRUE), 
           heights=c(0.1,0.45,0.45), widths=width)
    par(mar=rep(0,4))
    plot.new()
    text(0.5, 0.5, 
         paste(maintext, sub("_", " ", ev_stage), "stars by", experiment),
         cex=2, font=2)
}

load_experiment_info <- function(exp, dirs) {
    print(exp)
    print(dirs)
    if (grepl('exp_Y', exp)) {
        experiment <<- 'helium'
        sun_num <<- grep('.28', dirs)[1]
    } else if (grepl('exp_alpha', exp)) {
        experiment <<- 'mixing length'
        sun_num <<- grep('2.10', dirs)[1]
    }
    
    freep <<- as.name(sub('.+_', '', basename(exp)))
    cl <<- heat.colors(length(dirs)+color_offset)[1:length(dirs)]
    cl <<- c(cl[1:(sun_num-1)], "black", cl[(sun_num+1):length(cl)])
    labls <<- sapply(sub('.+_', '', dirs), 
        function(x) { as.expression(bquote(.(freep) ~ "=" ~ .(x))) })
}

add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, function(x) 
        rgb(x[1], x[2], x[3], alpha=alpha))
}
