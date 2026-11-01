#!/usr/bin/env python3
# <Purpose in one line.>  Score probe output against two independent
# mpmath routes and print the relative error of each column, with the
# agreement of the routes (STAN_MATH_DISTRIBUTIONS.md section 3).
#
# Usage: python3 refs/score.py < out/probe_<tree>.txt > out/score_<tree>.txt
#
# Rules learned on earlier PRs (STEVIEB_BOT_REPLIES.md sections 5 to 5c):
# - two routes that do not share code (a series and an integral, a sum
#   and a closed form, the lower and the upper tail);
# - every arithmetic step inside the working precision (no float()
#   before the end);
# - give the relative error; where the true value is near 0, print the
#   computed value instead;
# - one script over all regions, and print the region of each point.
import sys
import mpmath as mp

mp.mp.dps = 50


def route_a(x, a):
    # Replace: first route, for example a direct formula
    return mp.loggamma(x + a) - mp.loggamma(x)


def route_b(x, a):
    # Replace: an independent route, for example a sum or an integral
    return mp.log(mp.rf(x, a))


def refs(route, x, a):
    return [route(x, a),
            mp.diff(lambda v: route(v, a), x),
            mp.diff(lambda v: route(x, v), a)]


def rel(s, r):
    s = mp.mpf(s)
    return abs(s - r) / abs(r) if r != 0 else abs(s)


print(f"{'x':>10s} {'a':>5s} {'value':>9s} {'d/dx':>9s} {'d/da':>9s}"
      f" {'A vs B':>8s} status")
for line in sys.stdin:
    f = line.split()
    if len(f) != 6:
        continue
    x, a = mp.mpf(f[0]), mp.mpf(f[1])
    stan = f[2:5]
    status = f[5]
    ra = refs(route_a, x, a)
    rb = refs(route_b, x, a)
    agree = max(rel(p, q) for p, q in zip(ra, rb))
    errs = ["nan" if s == "nan" else mp.nstr(rel(s, r), 2) for s, r in zip(stan, ra)]
    print(f"{mp.nstr(x, 3):>10s} {mp.nstr(a, 3):>5s} "
          + " ".join(f"{e:>9s}" for e in errs)
          + f" {mp.nstr(agree, 2):>8s} {status}")
