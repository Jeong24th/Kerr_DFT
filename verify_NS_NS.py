#!/usr/bin/env python3
"""
verify_NS_NS.py
===============

Independent SymPy verification that the full rotating closed-string vacuum
{g_{mu nu}, B_{mu nu}, phi} of [Rotating_NSNS_PRL.tex, boxed eqs (4)-(6)] satisfies
the four-dimensional NS-NS vacuum equations of motion at qq != 0 generic
parameters [PRL eq:EDFE, lines 95-101]:

  E1_{mu nu} := R_{mu nu} + 2 nabla_mu partial_nu phi
                - (1/4) H_{mu rho sigma} H_nu^{rho sigma}                 = 0
  E2^{mu nu} := nabla_rho (e^{-2 phi} H^{rho mu nu})                       = 0
  E3        := R + 4 Box phi - 4 (partial phi)^2 - (1/12) H_{lambda mu nu} H^{lambda mu nu}  = 0

We work in Boyer-Lindquist coordinates and evaluate at three generic test
points (matching the Mathematica notebook test point and two independent
ones).  All residuals should be 0 to machine precision.

Run:  python verify_NS_NS.py
"""

import sympy as sp


# ============================================================
# Coordinates and parameters
# ============================================================
t, rb, th, ph = sp.symbols("t rb theta phi", real=True)
coords = (t, rb, th, ph)

m, j, q, zeta = sp.symbols("m j q zeta", positive=True, real=True)
mj = sp.sqrt(m**2 - j**2)


# ============================================================
# Auxiliary BL functions  (PRL aux, DeltaSigmaBL, Bomega)
# ============================================================
Delta_ = rb**2 - 2*m*rb + j**2
Sigma_ = rb**2 + j**2 * sp.cos(th)**2
Lambda_ = Delta_ - j**2 * sp.sin(th)**2
Xi_ = (rb**2 + j**2) * Sigma_ + 2*m*j**2*rb * sp.sin(th)**2
Kperp_ = (1 + mj**2 * sp.sin(th)**2 / Delta_) ** (-q**2 / mj**2)

# Dilaton (BL form, eq dilatonBL)
F_ = (rb - m - mj) / (rb - m + mj)
e2phi = F_**(q/mj) * sp.cos(zeta)**2 + F_**(-q/mj) * sp.sin(zeta)**2
phi_ = sp.Rational(1, 2) * sp.log(e2phi)


# ============================================================
# String-frame metric g^{str} = e^{2 phi} g^E  (eq metricBL)
# ============================================================
g = sp.zeros(4, 4)
g[0, 0] = -e2phi * Lambda_ / Sigma_
g[0, 3] = -2*m*j*rb * e2phi * sp.sin(th)**2 / Sigma_
g[3, 0] = g[0, 3]
g[3, 3] = e2phi * Xi_ * sp.sin(th)**2 / Sigma_
g[1, 1] = e2phi * Sigma_ * Kperp_ / Delta_
g[2, 2] = e2phi * Sigma_ * Kperp_

# Block-wise inverse.
# det g^{str}_{(t,phi)} = e^{4 phi} det g^E_{(t,phi)} = -e^{4 phi} Delta sin^2 vartheta
det_tphi = -e2phi**2 * Delta_ * sp.sin(th)**2
ginv = sp.zeros(4, 4)
ginv[0, 0] = g[3, 3] / det_tphi
ginv[0, 3] = -g[0, 3] / det_tphi
ginv[3, 0] = ginv[0, 3]
ginv[3, 3] = g[0, 0] / det_tphi
ginv[1, 1] = 1 / g[1, 1]
ginv[2, 2] = 1 / g[2, 2]

# Full metric determinant (string frame).
# det g^{str} = det g_{(t,phi)} * det g_{(r,th)} (block diagonal)
# Note: this is the algebraic det; for sqrt(-det g) use the absolute value.
det_g = det_tphi * g[1, 1] * g[2, 2]
sqrt_mg = sp.sqrt(-det_g)


# ============================================================
# Kalb-Ramond B-field and H = dB
# B_{t phi} = h cos vartheta with h = -2 q sin 2 zeta;  others zero
# H_{mu nu rho} = partial_mu B_{nu rho} + cyclic
# ============================================================
h_const = -2 * q * sp.sin(2 * zeta)
B = sp.zeros(4, 4)
B[0, 3] = h_const * sp.cos(th)
B[3, 0] = -B[0, 3]

H = [[[sp.S(0)] * 4 for _ in range(4)] for _ in range(4)]
for a in range(4):
    for b in range(4):
        for c in range(4):
            H[a][b][c] = (
                sp.diff(B[b, c], coords[a])
                + sp.diff(B[c, a], coords[b])
                + sp.diff(B[a, b], coords[c])
            )


# ============================================================
# Christoffel symbols  Gamma^a_{bc}
# ============================================================
print("Computing Christoffel symbols...")
Gamma = [[[sp.S(0)] * 4 for _ in range(4)] for _ in range(4)]
for a in range(4):
    for b in range(4):
        for c in range(4):
            s = sp.S(0)
            for d in range(4):
                if ginv[a, d] == 0:
                    continue
                s += ginv[a, d] * (
                    sp.diff(g[d, c], coords[b])
                    + sp.diff(g[d, b], coords[c])
                    - sp.diff(g[b, c], coords[d])
                )
            Gamma[a][b][c] = sp.Rational(1, 2) * s


# ============================================================
# Ricci tensor
# ============================================================
print("Computing Ricci tensor...")
Ric = sp.zeros(4, 4)
for mu in range(4):
    for nu in range(mu, 4):
        s = sp.S(0)
        for a in range(4):
            s += sp.diff(Gamma[a][nu][mu], coords[a])
            s -= sp.diff(Gamma[a][a][mu], coords[nu])
            for b in range(4):
                s += Gamma[a][a][b] * Gamma[b][nu][mu]
                s -= Gamma[a][nu][b] * Gamma[b][a][mu]
        Ric[mu, nu] = s
        if mu != nu:
            Ric[nu, mu] = s


# ============================================================
# Build NS-NS equation residuals
# ============================================================
print("Building E1 (graviton)...")
# nabla_mu partial_nu phi  = partial_mu partial_nu phi - Gamma^lambda_{mu nu} partial_lambda phi
dphi = [sp.diff(phi_, x) for x in coords]
NablaDphi = sp.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        s = sp.diff(dphi[nu], coords[mu])
        for lam in range(4):
            s -= Gamma[lam][mu][nu] * dphi[lam]
        NablaDphi[mu, nu] = s

# H_{mu rho sigma} H_nu^{rho sigma} = H_{mu rho sigma} g^{rho a} g^{sigma b} H_{nu a b}
HHterm = sp.zeros(4, 4)
for mu in range(4):
    for nu in range(mu, 4):
        s = sp.S(0)
        for rho in range(4):
            for sig in range(4):
                for a in range(4):
                    for b in range(4):
                        if ginv[rho, a] == 0 or ginv[sig, b] == 0:
                            continue
                        s += H[mu][rho][sig] * ginv[rho, a] * ginv[sig, b] * H[nu][a][b]
        HHterm[mu, nu] = s
        if mu != nu:
            HHterm[nu, mu] = s

E1 = sp.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        E1[mu, nu] = Ric[mu, nu] + 2 * NablaDphi[mu, nu] - sp.Rational(1, 4) * HHterm[mu, nu]

print("Building E2 (B-field)...")
# E2^{mu nu} = (1/sqrt(-g)) partial_rho (e^{-2 phi} sqrt(-g) H^{rho mu nu})
# with H^{rho mu nu} = g^{rho a} g^{mu b} g^{nu c} H_{a b c}.
# Use this form because nabla_rho (totally antisym tensor) is a partial divergence.
Hup = [[[sp.S(0)] * 4 for _ in range(4)] for _ in range(4)]
for a in range(4):
    for b in range(4):
        for c in range(4):
            s = sp.S(0)
            for p in range(4):
                for r in range(4):
                    for q_ in range(4):
                        if ginv[a, p] == 0 or ginv[b, r] == 0 or ginv[c, q_] == 0:
                            continue
                        s += ginv[a, p] * ginv[b, r] * ginv[c, q_] * H[p][r][q_]
            Hup[a][b][c] = s

E2 = [[sp.S(0)] * 4 for _ in range(4)]
for mu in range(4):
    for nu in range(mu + 1, 4):
        s = sp.S(0)
        for rho in range(4):
            s += sp.diff(sp.exp(-2 * phi_) * sqrt_mg * Hup[rho][mu][nu], coords[rho])
        E2[mu][nu] = s / sqrt_mg
        # antisymmetric, only need mu < nu

print("Building E3 (dilaton trace)...")
# Box phi = (1/sqrt(-g)) partial_mu (sqrt(-g) g^{mu nu} partial_nu phi)
BoxPhi = sp.S(0)
for mu in range(4):
    s_inner = sp.S(0)
    for nu in range(4):
        if ginv[mu, nu] == 0:
            continue
        s_inner += ginv[mu, nu] * dphi[nu]
    BoxPhi += sp.diff(sqrt_mg * s_inner, coords[mu])
BoxPhi = BoxPhi / sqrt_mg

# (partial phi)^2 = g^{mu nu} partial_mu phi partial_nu phi
dphi_sq = sp.S(0)
for mu in range(4):
    for nu in range(4):
        if ginv[mu, nu] == 0:
            continue
        dphi_sq += ginv[mu, nu] * dphi[mu] * dphi[nu]

# Ricci scalar
R_scalar = sp.S(0)
for mu in range(4):
    for nu in range(4):
        if ginv[mu, nu] == 0:
            continue
        R_scalar += ginv[mu, nu] * Ric[mu, nu]

# H^2 = H_{a b c} H^{a b c}
H2 = sp.S(0)
for a in range(4):
    for b in range(4):
        for c in range(4):
            H2 += H[a][b][c] * Hup[a][b][c]

E3 = R_scalar + 4 * BoxPhi - 4 * dphi_sq - sp.Rational(1, 12) * H2


# ============================================================
# Evaluate residuals at test points
# ============================================================
def evaluate(label, subs_dict, digits=25):
    print(f"\n--- {label} ---")
    print(f"  Substitution: " + ", ".join(f"{k}={v}" for k, v in subs_dict.items()))

    # E1: 10 independent symmetric components
    print("  Eq1 (graviton) -- R_{mu nu} + 2 nabla_mu d_nu phi - (1/4) H H = 0")
    max_E1 = 0
    for mu in range(4):
        for nu in range(mu, 4):
            v = sp.N(E1[mu, nu].subs(subs_dict), digits)
            mag = float(sp.Abs(v)) if v.is_number else float("inf")
            tag = "OK  " if mag < 1e-15 else "WARN"
            print(f"    [{tag}] E1[{mu},{nu}] = {v}")
            max_E1 = max(max_E1, mag)

    # E2: 6 independent antisymmetric components
    print("  Eq2 (B-field) -- nabla_rho (e^{-2 phi} H^{rho mu nu}) = 0")
    max_E2 = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            v = sp.N(E2[mu][nu].subs(subs_dict), digits)
            mag = float(sp.Abs(v)) if v.is_number else float("inf")
            tag = "OK  " if mag < 1e-15 else "WARN"
            print(f"    [{tag}] E2[{mu},{nu}] = {v}")
            max_E2 = max(max_E2, mag)

    # E3: scalar dilaton equation
    print("  Eq3 (dilaton) -- R + 4 Box phi - 4 (d phi)^2 - (1/12) H^2 = 0")
    v = sp.N(E3.subs(subs_dict), digits)
    mag = float(sp.Abs(v)) if v.is_number else float("inf")
    tag = "OK  " if mag < 1e-15 else "WARN"
    print(f"    [{tag}] E3 = {v}")
    max_E3 = mag

    return max_E1, max_E2, max_E3


testpoints = [
    {
        "label": "test 1 (MMA notebook point): m=5, j=3, q=2, zeta=pi/7, rb=10, theta=pi/3",
        "subs": {m: 5, j: 3, q: 2, zeta: sp.pi/7, rb: 10, th: sp.pi/3},
    },
    {
        "label": "test 2: m=1, j=4/5, q=1/3, zeta=pi/5, rb=3, theta=pi/4",
        "subs": {m: 1, j: sp.Rational(4, 5), q: sp.Rational(1, 3), zeta: sp.pi/5,
                 rb: 3, th: sp.pi/4},
    },
    {
        "label": "test 3 (large q): m=1, j=1/2, q=2, zeta=pi/8, rb=5, theta=pi/6",
        "subs": {m: 1, j: sp.Rational(1, 2), q: 2, zeta: sp.pi/8, rb: 5, th: sp.pi/6},
    },
]

global_max = 0
for tp in testpoints:
    e1, e2, e3 = evaluate(tp["label"], tp["subs"])
    global_max = max(global_max, e1, e2, e3)

print("\n" + "=" * 60)
print(f"Max residual across all eqs and test points: {global_max:.2e}")
if global_max < 1e-15:
    print("PASS: NS-NS vacuum field equations satisfied at all test points.")
else:
    print("FAIL: Residual exceeds tolerance.")
print("=" * 60)
