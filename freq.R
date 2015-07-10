#### Plot l=0, n=1 modes for simulated Suns  
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

dirs <- dir('exp')
n_max <- 5
l0 <- matrix(nrow=length(dirs), ncol=n_max)
mix_lengths = c()
for (exp_i in 1:length(dirs)) {
    exp <- dirs[exp_i]
    mix_lengths <- c(mix_lengths, as.numeric(exp))
    data <- read.table(file.path('exp', exp, paste0('exp', exp, '.dat')))
    l0[exp_i,] = data[1:n_max,3]
}

offset <- 5
cl <- rev(heat.colors(n_max+offset))[offset:(offset+n_max)]
#cl <- heat.colors(n_max+offset)

cairo_pdf('freq.pdf', family="CM Roman")
for (n_i in 1:n_max) {
    fmla <- l0[,n_i] ~ mix_lengths
    if (n_i == 1) {
        par(las=1, cex.main=1.25, cex.lab=1.25)
        plot(fmla, pch=3, col=cl[1],
             xlab="Mixing Length (α)", ylab="Frequency (Hz)",
             #xlim=c(floor(min(mix_lengths)), ceiling(max(mix_lengths))),
             ylim=c(floor(min(l0)), ceiling(max(l0))),
             main=paste0("Radial frequencies of Sun-like ",
                         "stars with varied mixing lengths"))
    } else {
        points(fmla, pch=3, col=cl[n_i])
    }
    abline(lm(fmla), col=cl[n_i], lty=3)
}
legend("topleft", col=cl, pch=3, lty=3, bty='n',
       legend=c('ℓ=0, n=1',
                'ℓ=0, n=2',
                'ℓ=0, n=3',
                'ℓ=0, n=4',
                'ℓ=0, n=5'))
dev.off()
