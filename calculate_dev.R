### Load SSR ###
library(rEDM)

###set up parameter
raw_time_series <- ts_noy_meir[1:8200] #from noy_meir.R
E <- 2
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 250
step_size <- 50

### start algorithm
window_indices <- seq(window_size, NROW(raw_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

for(j in window_indices)
{
    index <- index + 1
    rolling_window <- raw_time_series[(j-window_size+1):j]
    
    norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
    
    # calculate best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta, silent=TRUE)
    best <- order(-smap$rho)[1]
    theta_best <- smap[best,]$theta
    
    # calculate eigenvalues for best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
    
    smap_co <- smap[[1]]$smap_coefficients
    
    matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
    
    for(k in 1:NROW(smap_co))
    {
        if(!is.na(smap_co[k,1]))
        {
            M <- rbind(smap_co[k, 1:E], cbind(diag(E - 1), rep(0, E - 1)))
            M_eigen <- eigen(M)$values
            lambda1 <- M_eigen[order(abs(M_eigen))[E]]
            
            matrix_eigen[k,1] <- abs(lambda1)
            matrix_eigen[k,2] <- Re(lambda1)
            matrix_eigen[k,3] <- Im(lambda1)
        }
    }
    
    # save results
    matrix_result[index,1] <- j
    matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
    matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
    matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

### plot dev vs time and in complex plain
pdf("dev_results.pdf", width = 7, height = 3)

par(mfrow=c(1,2),oma=c(0,0,0,0))

## dev vs time
par(mar=c(4,4.5,2,1)+0.1)

plot(matrix_result[,1], matrix_result[,2], ylim=c(0.25,1.1), lty=1, lwd=0.5, type="p", cex=1, pch=20, col="blue", las=0, cex.axis=1.5, cex.lab=1.5, main="", xlab="Time", ylab="| DEV |")

points(c(0,10000),c(1,1), type="l", col="gray", lwd=3, lty=2)

## complex plain
par(mar=c(4,4.5,2,1)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.5, cex.axis=1.5, type="l", xlab="Re(DEV)", ylab="Im(DEV)", xlim=c(0,1), ylim=c(-0.5,0.5))

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(matrix_result[,3],matrix_result[,4],type="p", cex=1, pch=20, col="blue")

dev.off()
