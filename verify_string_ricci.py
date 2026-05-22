#!/usr/bin/env python3
"""
verify_string_ricci.py
======================

Independent SymPy verification of the static-limit string-frame Ricci scalar
formula in PRL [eq:Ricci, line 459]:

    R|_{j=0} = -(2 e^{-2 phi} + 3 e^{-6 phi} sin^2 2 zeta) * R^E|_{j=0}
    R^E|_{j=0} = (2 q^2 / (rb^2 Delta)) * (1 + m^2 sin^2 vartheta / Delta)^{q^2 / m^2}

We work in Boyer-Lindquist coords at j=0 (so m_j = m, Delta = rb^2 - 2 m rb,
Sigma = rb^2, Lambda = Delta) and:
  1. Build the full string-frame metric g^{str}_{mu nu} = e^{2 phi} g^E_{mu nu}.
  2. Compute its Ricci scalar R.
  3. Build the predicted RHS using the claimed formula.
  4. Check  R - RHS = 0  at several test points.

This is the formula that gives the rotating-shifted divergence exponent
[eq:rPlusDivStatic line 467], so its correctness underlies the polar-repulsion
threshold and the curvature-singularity claim at r = r_+.

Run:  python verify_string_ricci.py
"""

import sympy as sp


# ============================================================
# Coordinates and parameters (static slice j = 0)
# ============================================================
t, rb, th, ph = sp.symbols("t rb theta phi", real=True)
coords = (t, rb, th, ph)
m, q, zeta = sp.symbols("m q zeta", positive=True, real=True)

# At j=0:  m_j -> m,  Delta -> rb^2 - 2 m rb,  Sigma -> rb^2,
# Lambda -> Delta,  Xi = (rb^2)*Sigma + 0 = rb^4,  omega = 0
Delta_ = rb**2 - 2*m*rb
Sigma_ = rb**2
Lambda_ = Delta_                    # = Delta - 0
Xi_ = rb**2 * Sigma_                # = rb^4
Kperp_ = (1 + m**2 * sp.sin(th)**2 / Delta_) ** (-q**2 / m**2)

# Dilaton e^{2 phi} (BL form, j=0 so m_j=m, exponent q/m)
F_ = (rb - 2*m) / rb                # = (rb - m - m_j) / (rb - m + m_j) at j=0
e2phi = F_**(q/m) * sp.cos(zeta)**2 + F_**(-q/m) * sp.sin(zeta)**2
phi_ = sp.Rational(1, 2) * sp.log(e2phi)


# ============================================================
# String-frame metric g^{str} = e^{2 phi} g^E   (j=0)
# ============================================================
# Einstein-frame g^E (j=0)
gE = sp.zeros(4, 4)
gE[0, 0] = -Lambda_ / Sigma_                       # = -(1 - 2m/rb)
gE[1, 1] = Sigma_ * Kperp_ / Delta_                # = K_perp / (1 - 2m/rb)
gE[2, 2] = Sigma_ * Kperp_                         # = rb^2 K_perp
gE[3, 3] = Xi_ * sp.sin(th)**2 / Sigma_            # = rb^2 sin^2 vartheta
# (omega = 0 at j=0 so no t-phi off-diagonal)

g = sp.zeros(4, 4)
for a in range(4):
    g[a, a] = e2phi * gE[a, a]

# Inverse: diagonal at j=0
ginv = sp.zeros(4, 4)
for a in range(4):
    ginv[a, a] = 1 / g[a, a]


# ============================================================
# Christoffel symbols
# ============================================================
print("Computing Christoffel symbols (j=0 slice)...")
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
# Ricci tensor and Ricci scalar
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

print("Computing Ricci scalar...")
R_scalar = sp.S(0)
for mu in range(4):
    for nu in range(4):
        if ginv[mu, nu] != 0:
            R_scalar += ginv[mu, nu] * Ric[mu, nu]


# ============================================================
# Predicted RHS from PRL eq:Ricci
# ============================================================
RE_predicted = (2 * q**2 / (rb**2 * Delta_)) * (
    1 + m**2 * sp.sin(th)**2 / Delta_
) ** (q**2 / m**2)
R_predicted = -(2 * sp.exp(-2*phi_) + 3 * sp.exp(-6*phi_) * sp.sin(2*zeta)**2) * RE_predicted


# ============================================================
# Numerical evaluation at test points
# ============================================================
def evaluate(label, subs_dict, digits=40):
    print(f"\n--- {label} ---")
    print(f"  Substitution: " + ", ".join(f"{k}={v}" for k, v in subs_dict.items()))
    val_R = sp.N(R_scalar.subs(subs_dict), digits)
    val_P = sp.N(R_predicted.subs(subs_dict), digits)
    diff = val_R - val_P
    rel = sp.Abs(diff) / sp.Max(sp.Abs(val_R), sp.Abs(val_P), sp.S(1))
    print(f"  R (computed):  {val_R}")
    print(f"  R (predicted): {val_P}")
    print(f"  Absolute diff: {sp.N(sp.Abs(diff), 8)}")
    print(f"  Relative diff: {sp.N(rel, 8)}")
    return float(rel)


testpoints = [
    {
        "label": "test 1: m=1, q=1/2, zeta=pi/7, rb=5/2, theta=pi/3",
        "subs": {m: 1, q: sp.Rational(1, 2), zeta: sp.pi/7, rb: sp.Rational(5, 2), th: sp.pi/3},
    },
    {
        "label": "test 2: m=2, q=1/3, zeta=pi/5, rb=7, theta=pi/4",
        "subs": {m: 2, q: sp.Rational(1, 3), zeta: sp.pi/5, rb: 7, th: sp.pi/4},
    },
    {
        "label": "test 3: m=1, q=1/4, zeta=pi/9, rb=10, theta=pi/6",
        "subs": {m: 1, q: sp.Rational(1, 4), zeta: sp.pi/9, rb: 10, th: sp.pi/6},
    },
    {
        "label": "test 4: m=1, q=1/2, zeta=0 (FJNW branch), rb=3, theta=pi/2",
        "subs": {m: 1, q: sp.Rational(1, 2), zeta: sp.S(0), rb: 3, th: sp.pi/2},
    },
]

max_rel = 0
for tp in testpoints:
    r = evaluate(tp["label"], tp["subs"])
    max_rel = max(max_rel, r)

print("\n" + "=" * 60)
print(f"Max relative error across all test points: {max_rel:.2e}")
if max_rel < 1e-30:
    print("PASS: eq:Ricci formula matches direct Ricci-scalar computation.")
else:
    print("FAIL: Formula does not match.")
print("=" * 60)
