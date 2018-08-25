# kalman function
kalman = function(y, A, B, Q, C, R, mu0, Sigma0) {
    dy = nrow(y)
    T = ncol(y)
    dx = length(mu0)
    I = diag(dx)

    # make placeholder whcih is zero matrix as below;
    mu.p = matrix(0, nrow = dx, ncol = T)
    Sigma.p = array(0, c(dx, dx, T))
    mu.f = matrix(0, nrow = dx, ncol = T)
    Sigma.f = array(0, c(dx, dx, T))

    # predict at time 1
    mu.p[, 1] = A %*% mu0
    Sigma.p[,, 1] = A %*% Sigma0 %*% t(A) + B %*% Q %*% t(B)
    nu = y[, 1] - H %*% mu.p[, 1]
    S = H %*% Sigma.p[,, 1] %*% t(C) + R
    K = Sigma.p[,, 1] %*% t(C) %*% solve(S)
    mu.f[, 1] = mu.p[, 1] + K %*% nu
    Sigma.f[,, 1] = (I - K %*% C) %*% Sigma.p[,, 1]

    # for loop for time 2:T
    for (t in (2:T)) {
        # Prediction
        mu.p[, t] = F %*% mu.f[, t - 1]
        Sigma.p[,, t] = A %*% Sigma.f[,, t - 1] %*% t(A) + B %*% Q %*% t(B)
        # Update
        nu = y[, t] - C %*% mu.p[, t]
        S = C %*% Sigma.p[,, t] %*% t(C) + R
        K = Sigma.p[,, t] %*% t(C) %*% solve(S)
        mu.f[, t] = mu.p[, t] + K %*% nu
        Sigma.f[,, t] = (I - K %*% H) %*% Sigma.p[,, t]
    }
    return(list(mu.f = mu.f, Sigma.f = Sigma.f))
}

