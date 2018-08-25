extendedkalmanf = function(y, A, B, Q, C, R, mu0, Sigma0, alpha) {

    library(numDeriv)

    dy = nrow(y)
    T = ncol(y)
    dx = length(mu0)
    I = diag(dx)

    ## p = PREDICTOR / f = KALMAN FILTER / s = SMOOTHED KALMAN FILTER##

    ##INITIALISE##
    mu.p = matrix(0, nrow = dx, ncol = T)
    Sigma.p = array(0, c(dx, dx, T))
    ci.p = matrix(0, nrow = 2, ncol = T)
    mu.f = matrix(0, nrow = dx, ncol = T)
    Sigma.f = array(0, c(dx, dx, T))
    ci.f = matrix(0, nrow = 2, ncol = T)
    mu.s = matrix(0, nrow = dx, ncol = T)
    Sigma.s = array(0, c(dx, dx, T))
    ci.s = matrix(0, nrow = 2, ncol = T)

    ##TIME 1##
    #PREDICTION#
    dA = grad(A, mu0)
    mu.p[, 1] = dA %*% mu0
    Sigma.p[,, 1] = dA %*% Sigma0 %*% t(dA) + B %*% Q %*% t(B)

    #UPDATE#
    dC = grad(C, mu.p[, 1])
    nu = y[, 1] - dC %*% mu.p[, 1]
    S = dC %*% Sigma.p[,, 1] + R
    K = Sigma.p[,, 1] %*% t(dC) %*% solve(S)
    mu.f[, 1] = mu.p[, 1] + K %*% nu
    Sigma.f[,, 1] = (I - K %*% dC) %*% Sigma.p[,, 1]

    ##RECURSION TIME 2:T##
    for (t in (2:T)) {
        #PREDICTION#
        dA = grad(A, mu.f[, t - 1])
        mu.p[, t] = dA %*% mu.f[, t - 1]
        Sigma.p[,, t] = dA %*% Sigma.f[,, t - 1] %*% t(dA) + B %*% Q %*% t(B)

        #UPDATE#
        dC = grad(C, mu.p[, t])
        nu = y[, t] - dC %*% mu.p[, t]
        S = dC %*% Sigma.p[,, t] + R
        K = Sigma.p[,, t] %*% t(dC) %*% solve(S)
        mu.f[, t] = mu.p[, t] + K %*% nu
        Sigma.f[,, t] = (I - K %*% dC) %*% Sigma.p[,, t]
    }

    ##SMOOTHING/BACKWARD RECURSION##
    mu.s[, T] = mu.f[, T]
    Sigma.s[,, T] = Sigma.f[,, T]
    for (t in (T - 1):1) {
        dA = grad(A, mu.s[, t + 1])
        J = Sigma.f[,, t] %*% t(dA) %*% solve(Sigma.p[,, t + 1])
        mu.s[, t] = mu.f[, t] + J %*% (mu.s[, t + 1] - mu.p[, t + 1])
        Sigma.s[,, t] = Sigma.f[,, t] + J %*% (Sigma.s[,, t + 1] - Sigma.p[,, t + 1]) %*% t(J)
    }

    ##CREDIBLE INTERVAL##
    #UPPER = 1/LOWER = 2#
    for (t in (1:T)) {
        ci.p[1, t] = qnorm(alpha / 2, mean = mu.p[, t], sd = Sigma.p[1, 1, t])
        ci.p[2, t] = qnorm(1 - alpha / 2, mean = mu.p[, t], sd = Sigma.p[1, 1, t])
        ci.f[1, t] = qnorm(alpha / 2, mean = mu.f[, t], sd = Sigma.f[1, 1, t])
        ci.f[2, t] = qnorm(1 - alpha / 2, mean = mu.f[, t], sd = Sigma.f[1, 1, t])
        ci.s[1, t] = qnorm(alpha / 2, mean = mu.s[, t], sd = Sigma.s[1, 1, t])
        ci.s[2, t] = qnorm(1 - alpha / 2, mean = mu.s[, t], sd = Sigma.s[1, 1, t])
    }

    ##PLOT##
    png(filename = "ekalmanfilter.png", width = 768, height = 1280)
    par(mfrow = c(3, 1))
    #p#
    matplot((1:T),
            t(rbind(mu.p, ci.p, y)),
            type = "l",
            col = c(1, 2, 2, 3),
            lty = c(1, 2, 2, 3),
            main = "Predictor",
            xlab = "Time",
            ylab = "mu.p",
            cex = 8
            )
    #f#
    matplot((1:T),
            t(rbind(mu.f, ci.f, y)),
            type = "l",
            col = c(1, 2, 2, 3),
            lty = c(1, 2, 2, 3),
            main = "Kalman Filter",
            xlab = "Time",
            ylab = "mu.f",
            cex = 8
            )
    #s#
    matplot((1:T),
            t(rbind(mu.s, ci.s, y)),
            type = "l",
            col = c(1, 2, 2, 3),
            lty = c(1, 2, 2, 3),
            main = "Smoothed Kalman Filter",
            xlab = "Time",
            ylab = "mu.s",
            cex = 8
            )
    dev.off()


}