#!/usr/bin/env python3
"""
verify_einstein_eq.py
=====================

Independent SymPy verification that the Bogush-Galtsov Einstein-frame seed
[Rotating_NSNS_PRL.tex eq:BG_E_metric, line 415] satisfies the Einstein-scalar
field equation

    R^E_{mu nu} - 2 partial_mu phi_seed partial_nu phi_seed = 0,

with seed dilaton phi_seed = (q / m_j) ln G  (PRL line 418), where
G(r) = (2r - m_j)/(2r + m_j) in qi coordinates, equivalently
phi_seed = (q / (2 m_j)) ln F in BL with F = (rb - rb_+)/(rb - rb_-).

We work in Boyer-Lindquist coordinates (cleaner) and substitute generic
numerical test points to confirm every component of the residual is zero
to machine precision.

Run:  python verify_einstein_eq.py
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
# Auxiliary BL functions (PRL eqs aux, DeltaSigmaBL, Bomega)
# ============================================================
Delta_ = rb**2 - 2*m*rb + j**2
Sigma_ = rb**2 + j**2 * sp.cos(th)**2
Lambda_ = Delta_ - j**2 * sp.sin(th)**2          # = Sigma - 2 m rb
Xi_ = (rb**2 + j**2) * Sigma_ + 2*m*j**2*rb * sp.sin(th)**2
omega_ = -2*m*j*rb * sp.sin(th)**2 / Lambda_
Kperp_ = (1 + mj**2 * sp.sin(th)**2 / Delta_) ** (-q**2 / mj**2)


# ============================================================
# Bogush-Galtsov Einstein-frame seed metric (BL coordinates)
# ds^2_E = -(Lambda/Sigma)(dt - omega dphi)^2
#        + (Sigma K_perp / Delta) drb^2 + Sigma K_perp dvartheta^2
#        + (Sigma Delta sin^2 vartheta / Lambda) dphi^2
# After expanding (dt - omega dphi)^2 and using the identity
# Lambda Xi + 4 m^2 j^2 rb^2 sin^2 vartheta = Sigma^2 Delta,
# the (phi, phi) component collapses to Xi sin^2 vartheta / Sigma.
# ============================================================
g = sp.zeros(4, 4)
g[0, 0] = -Lambda_ / Sigma_
g[0, 3] = -2*m*j*rb * sp.sin(th)**2 / Sigma_
g[3, 0] = g[0, 3]
g[3, 3] = Xi_ * sp.sin(th)**2 / Sigma_
g[1, 1] = Sigma_ * Kperp_ / Delta_
g[2, 2] = Sigma_ * Kperp_


# ============================================================
# Block-wise inverse (analytical, avoids slow .inv() on 4x4 symbolic)
# det g^E_{(t,phi)} = - Delta sin^2 vartheta  (uses identity eq:detid)
# (rb, vartheta) block is diagonal
# ============================================================
det_tphi = -Delta_ * sp.sin(th)**2
ginv = sp.zeros(4, 4)
ginv[0, 0] = g[3, 3] / det_tphi
ginv[0, 3] = -g[0, 3] / det_tphi
ginv[3, 0] = ginv[0, 3]
ginv[3, 3] = g[0, 0] / det_tphi
ginv[1, 1] = 1 / g[1, 1]
ginv[2, 2] = 1 / g[2, 2]


# ============================================================
# Seed dilaton: phi_seed = (q / (2 m_j)) ln F with F = (rb-m-m_j)/(rb-m+m_j)
# Equivalently phi_seed = (q / m_j) ln G in qi, since G^2 = F.
# ============================================================
F_ = (rb - m - mj) / (rb - m + mj)
phi_seed = (q / (2 * mj)) * sp.log(F_)


# ============================================================
# Christoffel symbols Gamma^a_{bc} = (1/2) g^{ad}(d_b g_{dc} + d_c g_{db} - d_d g_{bc})
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

# Symmetry test: Gamma[a][b][c] == Gamma[a][c][b]  (assured by construction)


# ============================================================
# Ricci tensor R_{mu nu} = d_a Gamma^a_{nu mu} - d_nu Gamma^a_{a mu}
#                       + Gamma^a_{a b} Gamma^b_{nu mu}
#                       - Gamma^a_{nu b} Gamma^b_{a mu}
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
# Scalar source 2 partial_mu phi_seed partial_nu phi_seed
# ============================================================
dphi = [sp.diff(phi_seed, x) for x in coords]
T_scalar = sp.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        T_scalar[mu, nu] = 2 * dphi[mu] * dphi[nu]


# ============================================================
# Numerical evaluation at generic test points
# ============================================================
def evaluate(label, subs_dict, digits=25):
    print(f"\n--- {label} ---")
    print(f"  Substitution: " + ", ".join(f"{k}={v}" for k, v in subs_dict.items()))
    max_resid = 0
    for mu in range(4):
        for nu in range(mu, 4):
            val = (Ric[mu, nu] - T_scalar[mu, nu]).subs(subs_dict)
            try:
                v = sp.N(val, digits)
            except Exception:
                v = val
            mag = float(sp.Abs(v)) if v.is_number else float("inf")
            tag = "OK " if mag < 1e-20 else "WARN"
            print(f"  [{tag}] R^E[{mu},{nu}] - 2 dphi_seed[{mu}] dphi_seed[{nu}] = {v}")
            max_resid = max(max_resid, mag)
    return max_resid


# Test point 1 (matches Mathematica notebook)
testpt1 = {
    m: sp.Rational(5, 1),
    j: sp.Rational(3, 1),
    q: sp.Rational(2, 1),
    zeta: sp.pi / 7,
    rb: sp.Rational(10, 1),
    th: sp.pi / 3,
}

# Test point 2 (different angles/values; non-special radius)
testpt2 = {
    m: sp.Rational(1, 1),
    j: sp.Rational(17, 25),
    q: sp.Rational(3, 7),
    zeta: sp.pi / 5,
    rb: sp.Rational(13, 5),
    th: sp.pi / 4,
}

# Test point 3 (Schwarzschild-like, j -> 0, q small)
testpt3 = {
    m: sp.Rational(1, 1),
    j: sp.Rational(1, 100),  # near static
    q: sp.Rational(1, 10),
    zeta: sp.pi / 3,
    rb: sp.Rational(7, 2),
    th: sp.pi / 4,
}


r1 = evaluate("Test point 1 (matches MMA notebook: m=5, j=3, q=2, zeta=pi/7, rb=10, theta=pi/3)", testpt1)
r2 = evaluate("Test point 2 (m=1, j=17/25, q=3/7, zeta=pi/5, rb=13/5, theta=pi/4)", testpt2)
r3 = evaluate("Test point 3 (near-static: m=1, j=1/100, q=1/10, zeta=pi/3, rb=7/2, theta=pi/4)", testpt3)

print("\n" + "=" * 60)
print(f"Max residual across all test points: {max(r1, r2, r3):.2e}")
if max(r1, r2, r3) < 1e-20:
    print("PASS: Einstein-scalar field equation satisfied at all test points.")
else:
    print("FAIL: Residual exceeds tolerance.")
print("=" * 60)
