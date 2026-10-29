#!/usr/bin/env python3
# mpmath references for the example of stan-dev/math#3469, used in the
# comment of stan-dev/math#3438: the value at theta = 0,
# ll(0) - log det B(0) / 2, and the Laplace approximation at the mode.
# Each quantity has two routes; the script prints their difference.
#
# Model: y_i ~ Poisson(exp(mu + theta_i)), theta ~ N(0, alpha^2 J + 1e-8 I).
import mpmath as mp

mp.mp.dps = 80
Y = [3, 1, 4, 1, 5]
N = len(Y)
JITTER = mp.mpf("1e-8")


def cov(alpha):
    a2 = mp.mpf(alpha) ** 2
    return mp.matrix([[a2 + (JITTER if i == j else 0) for j in range(N)]
                      for i in range(N)])


def loglik(mu, theta):
    return mp.fsum(y * (mu + t) - mp.exp(mu + t) - mp.loggamma(y + 1)
                   for y, t in zip(Y, theta))


def log_det_b(mu, theta, alpha):
    """log det(I + W^1/2 K W^1/2) with mp.det (route 1)."""
    w = [mp.exp(mu + t) for t in theta]
    k = cov(alpha)
    b = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            b[i, j] = (1 if i == j else 0) + mp.sqrt(w[i] * w[j]) * k[i, j]
    return mp.log(mp.det(b))


def log_det_b0_closed(mu, alpha):
    """log det B(0) in closed form (route 2): W = e^mu I, so the
    eigenvalues of W K are e^mu * 1e-8 (4 times) and
    e^mu * (N alpha^2 + 1e-8)."""
    w = mp.exp(mu)
    return (N - 1) * mp.log1p(w * JITTER) + mp.log1p(
        w * (N * mp.mpf(alpha) ** 2 + JITTER))


def mode_newton(mu, alpha):
    """Mode of log p(theta | y) by Newton from theta = 0 with step halving
    (route 1)."""
    kinv = cov(alpha) ** -1
    theta = mp.matrix([0] * N)

    def obj(th):
        return loglik(mu, th) - (th.T * kinv * th)[0] / 2

    for _ in range(500):
        g = mp.matrix([y - mp.exp(mu + t) for y, t in zip(Y, theta)]) \
            - kinv * theta
        h = -kinv - mp.diag([mp.exp(mu + t) for t in theta])
        step = -(h ** -1) * g
        f0 = obj(theta)
        lam = mp.mpf(1)
        while obj(theta + lam * step) < f0 and lam > mp.mpf("1e-30"):
            lam /= 2
        theta = theta + lam * step
        if mp.norm(lam * step) < mp.mpf("1e-60"):
            break
    return theta


def mode_findroot(mu, alpha, start):
    """Mode as the root of the gradient with mp.findroot (route 2)."""
    kinv = cov(alpha) ** -1

    def grad(*th):
        th = mp.matrix(th)
        g = mp.matrix([y - mp.exp(mu + t) for y, t in zip(Y, th)]) - kinv * th
        return [g[i] for i in range(N)]

    return mp.matrix(mp.findroot(grad, [start[i] for i in range(N)]))


def laplace(mu, alpha, theta):
    kinv = cov(alpha) ** -1
    return (loglik(mu, theta) - (theta.T * kinv * theta)[0] / 2
            - log_det_b(mu, theta, alpha) / 2)


def main():
    zero = mp.matrix([0] * N)
    print("mu alpha | value at theta = 0 | Laplace approximation "
          "| route differences")
    for mu in (-30, 1):
        for alpha in ("1e3", "1e5", "1e7", "1e9"):
            alpha = mp.mpf(alpha)
            v0 = loglik(mu, zero) - log_det_b(mu, zero, alpha) / 2
            v0c = loglik(mu, zero) - log_det_b0_closed(mu, alpha) / 2
            th1 = mode_newton(mu, alpha)
            th2 = mode_findroot(mu, alpha, th1 + mp.matrix([1e-6] * N))
            lap1 = laplace(mu, alpha, th1)
            lap2 = laplace(mu, alpha, th2)
            print("%4d %6s | %s | %s | %s %s" % (
                mu, mp.nstr(alpha, 1), mp.nstr(v0, 12), mp.nstr(lap1, 12),
                mp.nstr(abs(v0 - v0c), 3), mp.nstr(abs(lap1 - lap2), 3)))


if __name__ == "__main__":
    main()
