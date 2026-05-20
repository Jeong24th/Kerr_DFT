# ============================================================
# Rotating_closed_string_vacuum.mpl
# This is the metric file of Rotating_closed_string_vacuum
# Coordinates: x^mu = (t, rb, vartheta, varphi)
# Parameters: m, j, q, zeta
# ============================================================

# Spacetime dimension
Ndim_ := 4: 

# Define the coordinates
x1_ := t:
x2_ := rb:
x3_ := vartheta:
x4_ := varphi:

# Choose the Signature convention: - + + + (which corresponds to 2)
sig_ := 2:

# ------------------------------------------------------------
# Auxiliary functions
# ------------------------------------------------------------

mj_ := sqrt(m^2 - j^2):

Sigma_ := rb^2 + j^2*cos(vartheta)^2:

Delta_ := rb^2 - 2*m*rb + j^2:

F_ := (rb - m - mj_)/(rb - m + mj_):

E2phi_ := F_^(q/mj_)*cos(zeta)^2
         + F_^(-q/mj_)*sin(zeta)^2:

omega_ := -2*m*j*rb*sin(vartheta)^2/(Delta_ - j^2*sin(vartheta)^2):

# ------------------------------------------------------------
# Metric components
#
# Coordinate order:
#   1 = t
#   2 = rb
#   3 = vartheta
#   4 = varphi
#
# ds^2 = g11 dt^2 + 2 g14 dt dvarphi + g44 dvarphi^2
#      + g22 dr^2 + g33 dvartheta^2
# ------------------------------------------------------------

g11_ := -E2phi_*(Delta_ - j^2*sin(vartheta)^2)/Sigma_:

g14_ := -2*m*j*rb*E2phi_*sin(vartheta)^2/Sigma_:

g22_ := E2phi_*Sigma_*((1 + mj_^2*sin(vartheta)^2/Delta_)^(-q^2/mj_^2))/Delta_:

g33_ := E2phi_*Sigma_*((1 + mj_^2*sin(vartheta)^2/Delta_)^(-q^2/mj_^2)):

g44_ := E2phi_*((rb^2 + j^2)*Sigma_ + 2*m*j^2*rb*sin(vartheta)^2)*sin(vartheta)^2/Sigma_:

Info_ := `Exact rotating closed-string vacuum metric (JKP, 2026)`: