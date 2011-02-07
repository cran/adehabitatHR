kerneloverlap <- function(xy, method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
                          percent=95, conditional=FALSE, ...)
{
    ## Verifications
    method <- match.arg(method)
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    if (percent >100)
        stop("percent should not be >100")

    ## UD estimation
    x <- kernelUD(xy, same4all=TRUE, ...)
    vol <- getvolumeUD(x)
    gp <- gridparameters(vol[[1]])

    ## Matrix of results
    res <- matrix(0, ncol=length(x), nrow=length(x))

    ## loop for each animal
    for (i in 1:length(x)) {
        for (j in 1:i) {

            if (method=="HR") {
                vi <- as.data.frame(vol[[i]])[,1]
                vj <- as.data.frame(vol[[j]])[,1]
                vi[vi<=percent] <- 1
                vi[vi>percent] <- 0
                vj[vj<=percent] <- 1
                vj[vj>percent] <- 0
                vk <- vi*vj
                res[i,j] <- sum(vk)/sum(vi)
                res[j,i] <- sum(vk)/sum(vj)
            }

            if (method=="PHR") {
                vi <- as.data.frame(x[[i]])[,1]
                vj <- as.data.frame(x[[j]])[,1]
                ai <- as.data.frame(vol[[i]])[,1]
                aj <- as.data.frame(vol[[j]])[,1]
                ai[ai<=percent] <- 1
                ai[ai>percent] <- 0
                aj[aj<=percent] <- 1
                aj[aj>percent] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- sum(vi*aj)*(gp[1,2]^2)
                    res[i,j] <- sum(vj*ai)*(gp[1,2]^2)
                } else {
                    res[j,i] <- sum(vi*aj)*(gp[1,2]^2)
                    res[i,j] <- sum(vj*ai)*(gp[1,2]^2)
                }
            }



            if (method=="VI") {
                vi <- c(as.data.frame(x[[i]])[,1])
                vj <- c(as.data.frame(x[[j]])[,1])
                ai <- as.data.frame(vol[[i]])[,1]
                aj <- as.data.frame(vol[[j]])[,1]
                ai[ai<=percent] <- 1
                ai[ai>percent] <- 0
                aj[aj<=percent] <- 1
                aj[aj>percent] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(gp[1,2]^2)
                } else {
                    res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(gp[1,2]^2)
                }
            }

            if (method=="BA") {
                vi <- c(as.data.frame(x[[i]])[,1])
                vj <- c(as.data.frame(x[[j]])[,1])
                ai <- as.data.frame(vol[[i]])[,1]
                aj <- as.data.frame(vol[[j]])[,1]
                ai[ai<=percent] <- 1
                ai[ai>percent] <- 0
                aj[aj<=percent] <- 1
                aj[aj>percent] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(gp[1,2]^2)
                } else {
                    res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(gp[1,2]^2)
                }
            }

            if (method=="UDOI") {
                vi <- c(as.data.frame(x[[i]])[,1])
                vj <- c(as.data.frame(x[[j]])[,1])
                ai <- as.data.frame(vol[[i]])[,1]
                aj <- as.data.frame(vol[[j]])[,1]
                ai[ai<=percent] <- 1
                ai[ai>percent] <- 0
                aj[aj<=percent] <- 1
                aj[aj>percent] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    ak <- sum(ai*aj)*(gp[1,2]^2)
                    res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(gp[1,2]^2)
                } else {
                    ak <- sum(ai*aj)*(gp[1,2]^2)
                    res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(gp[1,2]^2)
                }
            }

            if (method=="HD") {
                vi <- c(as.data.frame(x[[i]])[,1])
                vj <- c(as.data.frame(x[[j]])[,1])
                ai <- as.data.frame(vol[[i]])[,1]
                aj <- as.data.frame(vol[[j]])[,1]
                ai[ai<=percent] <- 1
                ai[ai>percent] <- 0
                aj[aj<=percent] <- 1
                aj[aj>percent] <- 0
                if (conditional) {
                    vi <- vi*ai
                    vj <- vj*aj
                    res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(gp[1,2]^2)))
                } else {
                    res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(gp[1,2]^2)))
                }
            }
        }
    }
    rownames(res) <- names(x)
    colnames(res) <- names(x)
    return(res)
}
