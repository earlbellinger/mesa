#### Plot frequency information for stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(RColorBrewer)
library(magicaxis)
library(sfsmisc)

ell_cl <- brewer.pal(4, "BrBG")

exp_dir <- file.path('dualexp')
plot_dir <- file.path('plots')
dir.create(plot_dir, showWarnings=FALSE)

for (experiment in list.dirs(exp_dir, recursive=FALSE)) {
    for (var in unlist(strsplit(basename(experiment), ';'))) { 
        nameval <- unlist(strsplit(var, "=")) # obtain M and Y
        assign(nameval[1], as.numeric(nameval[2]))
    }
    
    
}
