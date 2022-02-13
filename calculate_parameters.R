### Load SSR ###
library(rEDM)

###set up parameter
raw_time_series <- ts_noy_meir[1:8200] #from noy_meir.R
E <- seq(1,12,by=1)
tau <- seq(1,12,by=1)
theta <- 0
window_size <- seq(100, 600, by=50)
step_size <- 50

### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time

### start algorithm
matrix_results <- data.frame()

for(i in 1:length(window_size))
{
    for(j in 1:length(E))
    {
        for(k in 1:length(tau))
        {
            for(l in 1:length(theta))
            {
                matrix_rho <- c()
                m <- 0
                
                while(m <= length(raw_time_series) - window_size[i] - step_size)
                {
                    raw_ts_part <- raw_time_series[(m + 1):(m + window_size[i])]
                    time_series <- (raw_ts_part - mean(raw_ts_part, na.rm=TRUE))/sd(raw_ts_part, na.rm=TRUE)
                    
                    smap <- s_map(time_series, E=E[j], tau=tau[k], theta=theta[l], silent=TRUE)
                    matrix_rho <- cbind(matrix_rho, smap$rho)
                    
                    m <- m + step_size
                }
                
                matrix_results <- rbind(matrix_results, data.frame(
                w = window_size[i],
                E = E[j],
                tau = tau[k],
                theta = theta[l],
                rho = mean(matrix_rho)
                ))
                
            }
        }
    }
}

results <- matrix_results

### sort results
matrix_rho <- c()
for(i in 1:length(window_size))
{
    best_w <- order(-results[results$w == window_size[i],]$rho)[1]
    matrix_rho <- cbind(matrix_rho, results[results$w == window_size[i],]$rho[best_w])
}


### plot rho vs window size
pdf(sprintf("rho_vs_window.pdf",plot_name), width = 4, height = 3)

par(mar=c(4,4.5,2,1)+0.1)

plot(window_size, matrix_rho, lty=1, lwd=2, type="l", col="blue", las=0, cex.axis=1.5, cex.lab=1.5, xlim=c(min(window_size), max(window_size)), main="", xlab="Window size", ylab=expression("Predictability, "~rho))

dev.off()
