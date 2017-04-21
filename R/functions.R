meanlist <- function(x) rowMeans(sapply(x,I))

varlist  <- function(x) diag(var((t(sapply(x, I)))))

mycumsum <- function(x) {
    out <- NULL
    for(i in 1:length(x)) out <- c(out,sum(x[1:i],na.rm=T))
    out
}

bound <- function(x){
    xx <- x
    xx[x > 0.999] <- 0.999
    xx[x < 0.001] <- 0.001
    as.numeric(xx)
}


#' Compute log odds ratio
#'
#' @param theta is the vector with multinomial counterfactual distribution, P(Y <= k)
#'
#' @return point estimate
compBeta <- function(theta){
    mean(log(theta[1,]/(1-theta[1,]) * (1-theta[2,])/theta[2,]))
}

#' Inverse Probability Weighted Estimator
#'
#' @param data is a data frame with n rows with variables A (treatment variable), D (missingness indicator), and Y (outcome)
#' @param gDn1 is a vector of sicze n containing Pr(D=1 | A=1, W)
#' @param gDn0 is a vector containing Pr(D=1 | A=0, W)
#' @param gAn1 is a ectoe containing the propensity scores Pr(A=1|W)
#'
#' @return point estimate
iptw <- function(data, gDn1, gDn0, gAn1){

    A <- data$A
    D <- data$D
    Y <- data$Y
    Y[D==0] <- 9999
    if(length(unique(D))==1) K <- length(unique(Y)) else K <- length(unique(Y))-1
    gAn0 <- 1 - gAn1
    if(length(unique(gAn1))==1) gAn0 <- gAn1

    # Compute the IPTW
    numer <- sapply(1:(K-1), function(k)colSums(cbind(D*A/(gAn1*gDn1) * as.numeric(Y<=k), D*(1-A)/(gAn0*gDn0) * as.numeric(Y<=k))))
    denom <- sapply(1:(K-1), function(k)colSums(cbind(D*A/(gAn1*gDn1), D*(1-A)/(gAn0*gDn0))))
    return(numer/denom)
}

#' Augemented Inverse Probability Weighted Estimator of the Counterfactual Multinomial Distributions
#'
#' @param data is a data frame with n rows with variables A (treatment variable), D (missingness indicator), and Y (outcome)
#' @param Qkn1 n x (k-1) matrix containing P(Y=k|A=1,Y>=k,W)
#' @param Qkn0 Qkn0 same as above for A=0
#' @param gDn1 is a vector of sicze n containing Pr(D=1 | A=1, W)
#' @param gDn0 is a vector containing Pr(D=1 | A=0, W)
#' @param gAn1 is a ectoe containing the propensity scores Pr(A=1|W)
#'
#' @return point estimate
aiptw <- function(data, Qkn1, Qkn0, gDn1, gDn0, gAn1){

    A <- data$A
    D <- data$D
    Y <- data$Y
    Y[D==0] <- 9999
    if(length(unique(D))==1) K <- length(unique(Y)) else K <- length(unique(Y))-1
    n <- length(A)
    id <- rep(1:n, rep(K-1, n))

    if(K>2) {
        pkn1 <- 1 - t(apply(1-Qkn1, 1, cumprod))
        pkn0 <- 1 - t(apply(1-Qkn0, 1, cumprod))
    }
    if(K==2) {
        pkn1 <- data.frame(Qkn1)
        pkn0 <- data.frame(Qkn0)
    }

    gAn0 <- 1 - gAn1
    return(sapply(1:(K-1), function(k)colMeans(cbind(D*A/(gAn1*gDn1) * (as.numeric(Y<=k) - pkn1[,k]) +
                                                     pkn1[,k], D*(1-A)/(gAn0*gDn0) * (as.numeric(Y<=k) - pkn0[,k]) + pkn0[,k]))))
}

#' TMLE of the Counterfactual Multinomial Distributions
#'
#' @param data is a data frame with n rows with variables A (treatment variable), D (missingness indicator), and Y (outcome)
#' @param Qkn1 n x (k-1) matrix containing P(Y=k|A=1,Y>=k,W)
#' @param Qkn0 Qkn0 same as above for A=0
#' @param gDn1 is a vector of sicze n containing Pr(D=1 | A=1, W)
#' @param gDn0 is a vector containing Pr(D=1 | A=0, W)
#' @param gAn1 is a ectoe containing the propensity scores Pr(A=1|W)
#'
#' @return point estimate
estimaTheta <- function(data, Qkn1, Qkn0, gDn1, gDn0, gAn1){

    A <- data$A
    D <- data$D
    Y <- data$Y
    Y[D==0] <- 9999
    if(length(unique(D))==1) K <- length(unique(Y)) else K <- length(unique(Y))-1
    n <- length(Y)

    id <- rep(1:n, rep(K-1, n))
    Al <- A[id]
    Dl <- D[id]
    kl <- as.factor(rep(1:(K-1), n))
    Yl <- as.numeric(Y[id] == kl)
    Rl <- unlist(lapply(tapply(Yl, id, cumsum),
                        function(x)as.numeric(cumsum(x)<=1)))

    gAn0 <- 1 - gAn1

    Qkn <- A*Qkn1 + (1-A)*Qkn0

    Qkn1 <- do.call(c, as.data.frame(t(Qkn1)))
    Qkn0 <- do.call(c, as.data.frame(t(Qkn0)))
    Qkn  <- do.call(c, as.data.frame(t(Qkn)))

    Z1 <- Z0 <- matrix(0, ncol = K-1, nrow = n*(K-1))

    H1 <- 1/(gAn1[id] * gDn1[id])
    H0 <- 1/(gAn0[id] * gDn0[id])

    m    <- length(id)
    crit <- TRUE
    iter <- 1
    psin <- Dn <- list()

    compTheta <- function(Qkn1, Qkn0, K){

        tmp1 <- tapply(1-Qkn1, id, cumprod, simplify = FALSE)
        tmp0 <- tapply(1-Qkn0, id, cumprod, simplify = FALSE)

        if(K>2) {
            theta1 <- 1 - meanlist(tmp1)
            theta0 <- 1 - meanlist(tmp0)
        }
        if(K==2) {
            theta1 <- 1 - mean(sapply(tmp1,I))
            theta0 <- 1 - mean(sapply(tmp0,I))
        }

        return(rbind(theta1=theta1,theta0=theta0))

    }

    while(crit && iter <= 100){

        ind <- outer(as.numeric(kl), 1:(K-1), '<=')

        tmp1 <- tapply(1-Qkn1, id, cumprod, simplify = FALSE)
        tmp0 <- tapply(1-Qkn0, id, cumprod, simplify = FALSE)

        prodk1 <- do.call('rbind', tmp1[id])
        prodk0 <- do.call('rbind', tmp0[id])

        prodj1 <- unlist(tmp1)
        prodj0 <- unlist(tmp0)

        Z1 <- ind * prodk1 / prodj1 * H1
        Z0 <- ind * prodk0 / prodj0 * H0

        Z   <- Dl * cbind(Al * Z1, (1 - Al) * Z0)

        eps <- coef(glm(Yl ~ 0 + offset(qlogis(Qkn)) + Z, family = binomial(), subset = Rl == 1 & Dl == 1))

        Qkn1old <- Qkn1
        Qkn0old <- Qkn0

        Qkn1 <- bound(plogis(qlogis(Qkn1) + as.vector(cbind(Z1, 0 * Z0) %*% eps)))
        Qkn0 <- bound(plogis(qlogis(Qkn0) + as.vector(cbind(0 * Z1, Z0) %*% eps)))

        Qkn  <- Al*Qkn1 + (1-Al)*Qkn0

        iter <-  iter + 1
        crit <- any(colMeans((cbind(Qkn1old - Qkn1, Qkn0old - Qkn0))^2) > 1e-4/n)

    }

    tmp1 <- matrix(Qkn1, nrow = m, ncol = K-1)
    tmp0 <- matrix(Qkn0, nrow = m, ncol = K-1)
    DnY  <- Reduce('+', split(as.data.frame(Rl * Z * (Yl - cbind(tmp1, tmp0))), kl))

    tmp1 <- do.call(rbind, tapply(1-Qkn1, id, cumprod, simplify = FALSE))
    tmp0 <- do.call(rbind, tapply(1-Qkn0, id, cumprod, simplify = FALSE))

    theta <- compTheta(Qkn1, Qkn0,K)

    Dnl <- list()
    for(k in 1:(K-1)){
        Dnl[[k]]     <- matrix(0, nrow=n, ncol=2)
        Dnl[[k]][,1] <- DnY[, k]       - tmp1[,k]
        Dnl[[k]][,2] <- DnY[, k + K-1] - tmp0[,k]
    }

    sdn <- sqrt(sapply(Dnl, function(x)diag(var(x)))/n)


	# theta matrix 2 x (k-1), row 1 A = 1, row 2 A = 0
	# sdn matrix 2 x (k-1), standard errors
    return(list(tmle = theta, sdn = sdn))

}

#' Enhanded efficiency TMLE of the Log-odds ratio
#'
#' @param data is a data frame with n rows with variables A (treatment variable), D (missingness indicator), and Y (outcome)
#' @param Qkn1 n x (k-1) matrix containing P(Y=k|A=1,Y>=k,W)
#' @param Qkn0 Qkn0 same as above for A=0
#' @param gDn1 is a vector of sicze n containing Pr(D=1 | A=1, W)
#' @param gDn0 is a vector containing Pr(D=1 | A=0, W)
#' @param gAn1 is a ectoe containing the propensity scores Pr(A=1|W)
#'
#' @return point estimates and standard error
estimaBeta <- function(data, Qkn1, Qkn0, gDn1, gDn0, gAn1){

    A <- data$A
    D <- data$D
    Y <- data$Y
    Y[D==0] <- 9999
    if(length(unique(D))==1) K <- length(unique(Y)) else K <- length(unique(Y))-1
    n <- length(Y)

    id <- rep(1:n, rep(K-1, n))
    Al <- A[id]
    Dl <- D[id]
    kl <- as.factor(rep(1:(K-1), n))
    Yl <- as.numeric(Y[id] == kl)
    Rl <- unlist(lapply(tapply(Yl, id, cumsum),
                        function(x)as.numeric(cumsum(x)<=1)))

    gAn0 <- 1 - gAn1

    Qkn <- A*Qkn1 + (1-A)*Qkn0
    gDn <- A*gDn1 + (1-A)*gDn0
    gAn <- A*gAn1 + (1-A)*gAn0

    Qkn1 <- do.call(c, as.data.frame(t(Qkn1)))
    Qkn0 <- do.call(c, as.data.frame(t(Qkn0)))
    Qkn  <- do.call(c, as.data.frame(t(Qkn)))

    m    <- length(id)
    crit <- TRUE
    iter <- 1
    psin <- Dn <- list()

    compTheta <- function(Qkn1, Qkn0,K){

        tmp1 <- tapply(1-Qkn1, id, cumprod, simplify = FALSE)
        tmp0 <- tapply(1-Qkn0, id, cumprod, simplify = FALSE)

        if(K>2) {
            theta1 <- 1 - meanlist(tmp1)
            theta0 <- 1 - meanlist(tmp0)
        }
        if(K==2) {
            theta1 <- 1 - mean(sapply(tmp1,I))
            theta0 <- 1 - mean(sapply(tmp0,I))
        }

        return(rbind(theta1=theta1,theta0=theta0))

    }

    while(crit && iter <= 100){

        ind <- outer(as.numeric(kl), 1:(K-1), '<=')

        tmp1 <- tapply(1-Qkn1, id, cumprod, simplify = FALSE)
        tmp0 <- tapply(1-Qkn0, id, cumprod, simplify = FALSE)

        prodk1 <- do.call('rbind', tmp1[id])
        prodk0 <- do.call('rbind', tmp0[id])

        prodj1 <- unlist(tmp1)
        prodj0 <- unlist(tmp0)

        theta <- compTheta(Qkn1, Qkn0,K)

        Z1 <- ind * prodk1 / (prodj1 * gAn1[id] * gDn1[id])
        Z0 <- ind * prodk0 / (prodj0 * gAn0[id] * gDn0[id])

        Z   <- Dl * cbind(Al * Z1, (1 - Al) * Z0)

        M1 <- do.call(rbind, tapply(1-Qkn1, id, cumprod, simplify = FALSE))
        M0 <- do.call(rbind, tapply(1-Qkn0, id, cumprod, simplify = FALSE))


        H1 <- (M1 - theta[1,]) / (gAn1 * gDn1)
        H0 <- (M0 - theta[2,]) / (gAn0 * gDn0)
        H1 <- colSums(t(H1)/(theta[1,] * (1-theta[1,])))
        H0 <- colSums(t(H0)/(theta[2,] * (1-theta[2,])))
        H  <- A * H1 - (1 - A) * H0

        M1 <- (M1 - theta[1,]) / gAn1
        M0 <- (M0 - theta[2,]) / gAn0
        M1 <- colSums(t(M1)/(theta[1,] * (1-theta[1,])))
        M0 <- colSums(t(M0)/(theta[2,] * (1-theta[2,])))
        M  <- M1 + M0

        eps   <- coef(glm(Yl ~ 0 + offset(qlogis(Qkn)) + Z, family = binomial(), subset = Rl == 1 & Dl == 1,
                           start = rep(0, ncol(Z))))
        gamma <- coef(glm(D  ~ 0 + offset(qlogis(gDn)) + H, family = binomial(), start = 0))
        nu    <- coef(glm(A ~ 0 + offset(qlogis(gAn1)) + M, family = binomial(), start = 0))

        Qkn1old <- Qkn1
        Qkn0old <- Qkn0
        gDn1old <- gDn1
        gDn0old <- gDn0
        gAn1old <- gAn1

        eps[is.na(eps)] <- 0
        gamma[is.na(gamma)] <- 0
        nu[is.na(nu)] <- 0
        Qkn1 <- bound(plogis(qlogis(Qkn1) + as.vector(cbind(Z1, 0 * Z0) %*% eps)))
        Qkn0 <- bound(plogis(qlogis(Qkn0) + as.vector(cbind(0 * Z1, Z0) %*% eps)))

        gDn1 <- bound(plogis(qlogis(gDn1) + H1 * gamma))
        gDn0 <- bound(plogis(qlogis(gDn0) - H0 * gamma))
        gAn1 <- bound(plogis(qlogis(gAn1) + M * nu))
        gAn0 <- 1 - gAn1

        Qkn <- Al*Qkn1 + (1-Al)*Qkn0
        gDn <- A*gDn1 + (1-A)*gDn0
        gAn <- A*gAn1 + (1-A)*gAn0

        iter <-  iter + 1
        crit <- any(colMeans((cbind(Qkn1old - Qkn1, Qkn0old - Qkn0, gDn1old - gDn1,
                                    gDn0old - gDn0, gAn1old - gAn1))^2) > 1e-4/n)

    }

    Z1 <- colSums(t(Z1)/(theta[1,] * (1-theta[1,])))
    Z0 <- colSums(t(Z0)/(theta[2,] * (1-theta[2,])))
    Z  <- Dl * (Al * Z1 - (1 - Al) * Z0)
    DnY  <- tapply(Rl * Z * (Yl - Qkn), id, mean)

    tmp1 <- colMeans(do.call(cbind, tapply(1-Qkn1, id, cumprod, simplify = FALSE))/(theta[1,] * (1-theta[1,])))
    tmp0 <- colMeans(do.call(cbind, tapply(1-Qkn0, id, cumprod, simplify = FALSE))/(theta[2,] * (1-theta[2,])))

    Dnl <- DnY - tmp1 + tmp0
    theta <- compTheta(Qkn1, Qkn0,K)
    beta <- compBeta(theta)

    sdn <- sqrt(var(Dnl)/n)

    return(list(tmle = beta, theta = theta, sdn = sdn))
}
