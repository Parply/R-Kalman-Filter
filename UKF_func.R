unscentedkalmanf = function(lambda, kappa, y, A, B, Q, C, R, mu0, Sigma0, alpha) {

    dx = length(mu0)
    I = diag(dx)
    dy = nrow(y)
    T = ncol(y)
    N = NROW(mu0)

    ## p = PREDICTOR / f = KALMAN FILTER / s = SMOOTHED KALMAN FILTER##

    ##INITIALISE##
    Wc0 = kappa + lambda / (lambda + T)
    Wc = cbind(matrix(Wc0, nrow = dx, ncol = 1), matrix(1 / (2 * (lambda + dx)), nrow = dx, ncol = (2 * dx)))
    Wm0 = lambda / (lambda + T)
    Wm = cbind(matrix(Wm0, nrow = dx, ncol = 1), matrix(1 / (2 * (lambda + dx)), nrow = dx, ncol = (2 * dx)))

    ## p = PREDICTOR / f = KALMAN FILTER / s = SMOOTHED KALMAN FILTER##

    mu.p = matrix(0, nrow = dx, ncol = T)
    Yk = array(0, c(dx, (2 * dx + 1), T))
    Yk.p = matrix(0, nrow = dx, ncol = T)
    Xk.p = array(0, c(dx, (2 * dx + 1), T))
    Sigma.p = array(0, c(dx, dx, T))
    ci.p = matrix(0, nrow = 2, ncol = T)
    mu.f = matrix(0, nrow = dx, ncol = T)
    Sigma.f = array(0, c(dx, dx, T))
    ci.f = matrix(0, nrow = 2, ncol = T)
    mu.s = matrix(0, nrow = dx, ncol = T)
    Sigma.s = array(0, c(dx, dx, T))
    ci.s = matrix(0, nrow = 2, ncol = T)

    ##TIME 1##
    #CALCULATE SIGMAPOINTS#
    Qc = chol(Q)
    Sigmapoints = matrix(mu0, nrow = N, ncol = (2 * dx + 1))
    for (i in dx) {
        Sigmapoints[, i + 1] = Sigmapoints[, i + 1] + sqrt(dx + lambda) %*% Qc[, i]
        Sigmapoints[, i + dx + 1] = Sigmapoints[, i + dx + 1] + sqrt(dx + lambda) %*% Qc[, i]
    }

    #PREDICTION#
    Xk.p[,, 1] = A(1) %*% Sigmapoints
    mu.p[, 1] = sum(Wm[, 1:(2 * dx + 1)] * Xk.p[, 1, 1:(2 * dx + 1)])
    Sigma.p[,, 1] = sum(Wc[, 1:(2 * dx + 1)] * (Xk.p[, 1:(2 * dx + 1), 1] - mu.p[, 1]) %*% t(Xk.p[, 1:(2 * dx + 1), 1] - mu.p[, 1])) + Q
    Yk[,, 1] = C(1) %*% Xk.p[,, 1]
    Yk.p[, 1] = sum(Wm[, 1:(2 * dx + 1)] * Yk[, 1:(2 * dx + 1), 1])

    #UPDATE#
    Pyy = sum(Wc[, 1:(2 * dx + 1)] * (Yk[, 1:(2 * dx + 1), 1] - Yk.p[, 1]) %*% t(Yk[, 1:(2 * dx + 1), 1] - Yk.p[, 1])) + R
    Pxy = sum(Wc[, 1:(2 * dx + 1)] * (Xk.p[, 1:(2 * dx + 1), 1] - mu.p[, 1]) %*% t(Yk[, 1:(2 * dx + 1), 1] - Yk.p[, 1]))
    K = Pxy %*% solve(Pyy)
    mu.f[, 1] = mu.p[, 1] + K %*% (y[, 1] - Yk.p[, 1])
    Sigma.f[,, 1] = Sigma.p[,, 1] - K %*% Pyy %*% t(K)

    ##RECURSION TIME 2:T##
    for (t in (2:T)) {
        #CALCULATE SIGMAPOINTS#
        Sigmapoints = matrix(mu.f[, t - 1], nrow = N, ncol = (2 * dx + 1))
        for (i in dx) {
            Sigmapoints[, i + 1] = Sigmapoints[, i + 1] + sqrt(dx + lambda) %*% Qc[, i]
            Sigmapoints[, i + dx + 1] = Sigmapoints[, i + dx + 1] + sqrt(dx + lambda) %*% Qc[, i]
        }

        #PREDICTION#
        Xk.p[,, t] = A(t) %*% Sigmapoints
        mu.p[, t] = sum(Wm[, 1:(2 * dx + 1)] * Xk.p[, 1:(2 * dx + 1), t])
        Sigma.p[,, t] = sum(Wc[, 1:(2 * dx + 1)] * (Xk.p[, 1:(2 * dx + 1), t] - mu.p[, t]) %*% t(Xk.p[, 1:(2 * dx + 1), t] - mu.p[, t])) + Q
        Yk[,, t] = C(t) %*% Xk.p[,, t]
        Yk.p[, t] = sum(Wm[, 1:(2 * dx + 1)] * Yk[, 1:(2 * dx + 1), t])

        #UPDATE#
        Pyy = sum(Wc[, 1:(2 * dx + 1)] * (Yk[, 1:(2 * dx + 1), t] - Yk.p[, t]) %*% t(Yk[, 1:(2 * dx + 1), t] - Yk.p[, t])) + R
        Pxy = sum(Wc[, 1:(2 * dx + 1)] * (Xk.p[, 1:(2 * dx + 1), t] - mu.p[, t]) %*% t(Yk[, 1:(2 * dx + 1), t] - Yk.p[, t]))
        K = Pxy %*% solve(Pyy)
        mu.f[, t] = mu.p[, t] + K %*% (y[, t] - Yk.p[, t])
        Sigma.f[,, t] = Sigma.p[,, t] - K %*% Pyy %*% t(K)

    }

    ##SMOOTHING/BACKWARD RECURSION##
    mu.s[, T] = mu.f[, T]
    Sigma.s[,, T] = Sigma.f[,, T]
    for (t in (T - 1):1) {
        J = Sigma.f[,, t] %*% t(F(t)) %*% solve(Sigma.p[,, t + 1])
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
    png(filename = "ukalmanfilter.png", width = 768, height = 1280)
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
