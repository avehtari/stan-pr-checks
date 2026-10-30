"""mpmath references for the "valid ranges" table of stan-dev/math#3406.

Every value and d/dy comes from two routes that differ in erfcx(x), x >= 0:
route A is erfc(x) * exp(x^2) (asymptotic series above 1e6), route B is
the upper incomplete gamma Gamma(1/2, x^2) * exp(x^2) / sqrt(pi) (Lentz
continued fraction above 25).  Parameters: mu = 0, sigma = 1,
exp_mod_normal lambda = 1, skew_normal xi = 0, omega = 1, alpha = 1.
"""
import mpmath as mp

DPS = 40
mp.mp.dps = DPS


def erfcx_a(x):
    """erfcx(x), x >= 0: erfc * exp, or the asymptotic series."""
    assert x >= 0
    if x <= 1e6:
        with mp.extradps(20):
            return mp.erfc(x) * mp.exp(x * x)
    # sum (-1)^n (2n-1)!! / (2 x^2)^n; terms fall by 1e-12 or more
    t = mp.mpf(1)
    s = mp.mpf(1)
    n = 0
    u = 1 / (2 * x * x)
    while True:
        n += 1
        t *= -(2 * n - 1) * u
        s += t
        if abs(t) < mp.eps * abs(s) * mp.mpf(10) ** -5:
            break
    return s / (x * mp.sqrt(mp.pi))


def erfcx_b(x):
    """erfcx(x), x >= 0: incomplete gamma, or Lentz continued fraction."""
    assert x >= 0
    if x <= 25:
        with mp.extradps(20):
            return (mp.gammainc(mp.mpf(1) / 2, x * x, mp.inf)
                    * mp.exp(x * x) / mp.sqrt(mp.pi))
    # erfc(x) exp(x^2) sqrt(pi) = 1 / (x + (1/2) / (x + 1 / (x + (3/2) / ...)))
    tiny = mp.mpf(10) ** (-mp.mp.dps * 2)
    f = x
    c = x
    d = mp.mpf(0)
    k = 1
    while True:
        a = mp.mpf(k) / 2
        d = x + a * d
        d = tiny if d == 0 else d
        c = x + a / c
        c = tiny if c == 0 else c
        d = 1 / d
        delta = c * d
        f *= delta
        # c * d carries one rounding: stop at a few eps
        if abs(delta - 1) < 8 * mp.eps or k > 100000:
            break
        k += 1
    return 1 / (f * mp.sqrt(mp.pi))


class Ref:
    """All references for one erfcx route."""

    def __init__(self, erfcx):
        self.E = erfcx

    def q(self, u):  # 1 - Phi(u), u > 0
        return self.E(u / mp.sqrt(2)) / 2 * mp.exp(-u * u / 2)

    def phi_cdf(self, u):
        if u <= 0:
            return self.E(-u / mp.sqrt(2)) / 2 * mp.exp(-u * u / 2)
        return 1 - self.q(u)

    def log_phi_cdf(self, u):
        if u <= 0:
            return mp.log(self.E(-u / mp.sqrt(2)) / 2) - u * u / 2
        return mp.log1p(-self.q(u))

    def npdf(self, u):
        return mp.exp(-u * u / 2) / mp.sqrt(2 * mp.pi)

    def mills_inv(self, u):  # phi(u) / Phi(u)
        if u <= 0:
            return mp.sqrt(2 / mp.pi) / self.E(-u / mp.sqrt(2))
        return self.npdf(u) / self.phi_cdf(u)

    # ── exp_mod_normal, v = lambda sigma = 1 ──
    def emn_parts(self, z):
        """Return (F, G = 1 - F, f) as mpf.  z <= 0 uses the erfcx form
        F = exp(-z^2/2) (E(a) - E(b)) / 2 at raised precision."""
        v = mp.mpf(1)
        if z <= 0:
            extra = int(mp.log10(abs(z) + 2)) + 10
            with mp.extradps(extra):
                a = -z / mp.sqrt(2)
                b = (v - z) / mp.sqrt(2)
                ea, eb = self.E(a), self.E(b)
                w = mp.exp(-z * z / 2) / 2
                F = w * (ea - eb)
                f = w * eb
                G = 1 - F
                return +F, +G, +f, +(eb / (ea - eb))
        f = mp.exp(v * v / 2 - v * z) * self.phi_cdf(z - v)
        G = self.q(z) + f
        F = 1 - G
        return F, G, f, f / F

    def value_grad(self, name, y):
        y = mp.mpf(y)
        if name in ("std_normal_lcdf", "normal_lcdf"):
            return self.log_phi_cdf(y), self.mills_inv(y)
        if name in ("std_normal_lccdf", "normal_lccdf"):
            return self.log_phi_cdf(-y), -self.mills_inv(-y)
        if name in ("std_normal_cdf", "normal_cdf", "Phi"):
            return self.phi_cdf(y), self.npdf(y)
        if name.startswith("lognormal"):
            if y <= 0:
                return None
            z = mp.log(y)
            if name == "lognormal_lcdf":
                return self.log_phi_cdf(z), self.mills_inv(z) / y
            if name == "lognormal_lccdf":
                return self.log_phi_cdf(-z), -self.mills_inv(-z) / y
            return self.phi_cdf(z), self.npdf(z) / y
        if name == "exp_mod_normal_lpdf":
            u = y - 1
            return mp.mpf(1) / 2 - y + self.log_phi_cdf(u), -1 + self.mills_inv(u)
        if name.startswith("exp_mod_normal"):
            F, G, f, f_over_F = self.emn_parts(y)
            if name == "exp_mod_normal_lcdf":
                if y <= 0:
                    return mp.log(F), f_over_F
                return mp.log1p(-G), f_over_F
            if name == "exp_mod_normal_lccdf":
                if y <= 0:
                    return mp.log1p(-F), -f / G
                return mp.log(G), -f / G
            return F, f
        if name == "skew_normal_lpdf":
            return (mp.log(2) - mp.log(2 * mp.pi) / 2 - y * y / 2
                    + self.log_phi_cdf(y), -y + self.mills_inv(y))
        raise ValueError(name)


REF_A = Ref(erfcx_a)
REF_B = Ref(erfcx_b)


def rel_diff(a, b):
    if a == b:
        return mp.mpf(0)
    return abs(a - b) / max(abs(a), abs(b))


def reference(name, y):
    """Return (value, grad, agreement) with agreement the larger relative
    difference of the two routes, or None outside the domain."""
    ra = REF_A.value_grad(name, y)
    rb = REF_B.value_grad(name, y)
    if ra is None:
        return None
    agree = max(rel_diff(ra[0], rb[0]), rel_diff(ra[1], rb[1]))
    return ra[0], ra[1], agree
