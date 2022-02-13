set.seed(145477)

noy_meir <- function(r=0.75, b=0.1, h=0.75, p=2, F=1, t_start=1, t_end=10000, x0=6, sigma=0, noise_on=0)
{
    F <- if(length(F)==1){rep(F,t_end)}else{F}

    ricker <- function(N, r, b, h, p, F, noise_add=0)
    {
        out <- N*exp(r-b*N)-F*(N^p/(N^p+h^p)) + noise_add
        
        return(out)
    }
    
    L <- matrix(NaN, nrow = t_end+1-t_start, ncol = 1)
    L[1] <- x0
    
    for(i in 1:(NROW(L)-1))
    {
        L_noise_add <- if(noise_on!=0){sigma*L[i]*rnorm(1)}else{0}
        L[i+1] <- ricker(N=L[i], noise_add=L_noise_add, r=r, b=b, h=h, p=p, F=F[i])
    }
    
    return(L)
}

F_bifurcation <- seq(1,2,by=1/10000)
ts_noy_meir <- noy_meir(F=F_bifurcation, sigma=0.01, noise_on=1)

### plot rho vs window size
pdf("Variable_vs_time.pdf", width = 4, height = 3)

par(mar=c(4,4.5,2,1)+0.1)

plot(ts_noy_meir, lty=1, lwd=2, type="l", col="blue", las=0, cex.axis=1.5, cex.lab=1.5, main="", xlab="Time", ylab="Variable, N")

dev.off()
