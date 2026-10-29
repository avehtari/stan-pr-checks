#!/usr/bin/env python3
# mpmath reference for example D of stan-dev/math#3438 (test
# laplace_marginal_tol.solver_3_odd_permutation_saddle_point_throws):
# the maximum of log p(theta | y) near (-0.35, 0.003) and the Laplace
# approximation of the log marginal there.  Two routes for the mode.
import mpmath as mp

mp.mp.dps = 50
K = mp.matrix([[1, -1.5], [-1.5, 4]])
S = mp.mpf("0.3")
Y = [-S * mp.sqrt(3), S * mp.sqrt(3)]
Z = mp.mpf(0)
SZ = mp.mpf("0.1")
KINV = K ** -1


def loglik(t):
    v = mp.fsum(-mp.log(mp.pi * S) - mp.log1p(((y - t[0]) / S) ** 2)
                for y in Y)
    return v - mp.log(2 * mp.pi) / 2 - mp.log(SZ) - ((Z - t[1]) / SZ) ** 2 / 2


def grad_ll(t):
    g0 = mp.fsum(2 * (y - t[0]) / (S ** 2 + (y - t[0]) ** 2) for y in Y)
    return mp.matrix([g0, (Z - t[1]) / SZ ** 2])


def w_diag(t):
    w0 = mp.fsum(2 * (S ** 2 - (y - t[0]) ** 2) / (S ** 2 + (y - t[0]) ** 2) ** 2
                 for y in Y)
    return [w0, 1 / SZ ** 2]


def post_grad(t):
    return grad_ll(t) - KINV * t


def mode_newton(t):
    for _ in range(200):
        h = -KINV - mp.diag(w_diag(t))
        step = -(h ** -1) * post_grad(t)
        t = t + step
        if mp.norm(step) < mp.mpf("1e-45"):
            break
    return t


def mode_findroot(t):
    return mp.matrix(mp.findroot(
        lambda a, b: list(post_grad(mp.matrix([a, b]))), [t[0], t[1]]))


def laplace(t):
    b = mp.eye(2) + K * mp.diag(w_diag(t))
    return loglik(t) - (t.T * KINV * t)[0] / 2 - mp.log(mp.det(b)) / 2


def main():
    t1 = mode_newton(mp.matrix([-0.35, 0.003]))
    t2 = mode_findroot(mp.matrix([-0.3, 0.0]))
    print("mode (Newton):    %s %s" % (mp.nstr(t1[0], 20), mp.nstr(t1[1], 20)))
    print("mode (findroot):  %s %s" % (mp.nstr(t2[0], 20), mp.nstr(t2[1], 20)))
    print("mode difference:  %s" % mp.nstr(mp.norm(t1 - t2), 3))
    # The analytic gradient and W against numerical differentiation.
    f = lambda a, b: loglik(mp.matrix([a, b]))
    g_num = [mp.diff(f, (t1[0], t1[1]), (1, 0)),
             mp.diff(f, (t1[0], t1[1]), (0, 1))]
    w_num = [-mp.diff(f, (t1[0], t1[1]), (2, 0)),
             -mp.diff(f, (t1[0], t1[1]), (0, 2))]
    g, w = grad_ll(t1), w_diag(t1)
    print("mp.diff check:    gradient %s, W %s" % (
        mp.nstr(max(abs(g[i] - g_num[i]) for i in range(2)), 3),
        mp.nstr(max(abs(w[i] - w_num[i]) for i in range(2)), 3)))
    h = KINV + mp.diag(w_diag(t1))
    print("eig(K^-1 + W):    %s" % [mp.nstr(e, 6) for e in mp.eigsy(h)[0]])
    print("Laplace value:    %s" % mp.nstr(laplace(t1), 20))
    print("route difference: %s" % mp.nstr(abs(laplace(t1) - laplace(t2)), 3))
    z = mp.matrix([0, 0])
    b0 = mp.eye(2) + K * mp.diag(w_diag(z))
    print("at (0, 0): det(B) %s, eig(K^-1 + W) %s" % (
        mp.nstr(mp.det(b0), 8),
        [mp.nstr(e, 6) for e in mp.eigsy(KINV + mp.diag(w_diag(z)))[0]]))


if __name__ == "__main__":
    main()
