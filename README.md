# Rotating Closed String Vacuum

This repository contains two independent verification tracks — a Maple/GRTensorIII worksheet/workbook and a Mathematica notebook — for the exact rotating closed-string NS–NS vacuum solution of the Einstein Double Field Equations in double field theory.

**Coordinate convention.** Both tracks perform the verification natively in Boyer–Lindquist coordinates $(t,\bar r,\vartheta,\varphi)$, matching the solution as written in the Supplemental Material of the accompanying paper. The quasi-isotropic chart $(t,r,\vartheta,\varphi)$ used in the body of the paper follows by the explicit pull-back $\bar r(r)=r+m+(m^{2}-j^{2})/(4r)$; the Mathematica notebook additionally checks this pull-back map.

The Maple/GRTensorIII track supplies a GRTensorIII-compatible metric file together with Maple worksheet/workbook files that perform symbolic tensor calculations — line element, determinant, inverse metric, Ricci tensor, and Ricci scalar. In the vacuum check, the Ricci tensor and Ricci scalar simplify to zero. The Mathematica notebook performs an independent symbolic verification of the dilaton, Kalb–Ramond three-form, and string-frame metric directly against the NS–NS field equations in Boyer–Lindquist coordinates.

## Repository contents

| File | Description |
| --- | --- |
| `Rotating_closed_string_vacuum.mpl` | GRTensorIII metric definition file (Boyer–Lindquist chart). Defines coordinates, signature, parameters, auxiliary functions, and nonzero metric components. |
| `Rotating_Closed_String_Vacuum.mw` | Maple worksheet containing the symbolic GRTensorIII calculation and saved outputs. |
| `Rotating_Closed_String_Vacuum.maple` | Maple workbook version of the calculation file. |
| `analytic_proof_of__qi coordinate.nb` | Mathematica notebook: independent verification track. Checks the boxed dilaton, $H_{(3)}$, and string-frame metric directly against the NS–NS field equations in Boyer–Lindquist coordinates, then verifies the pull-back to quasi-isotropic coordinates. |

## Requirements

The files were prepared for use with:

- [Maple](https://www.maplesoft.com/products/Maple/) (worksheet metadata indicates Maple 2026.1)
- [GRTensorIII](https://github.com/grtensor/grtensor)
- [Wolfram Mathematica](https://www.wolfram.com/mathematica/) (the notebook was prepared with Wolfram 14.3)

Other recent Maple/GRTensorIII/Mathematica versions may also work, but symbolic-simplification behaviour and notebook rendering can vary by version.

## Metric file summary

The metric file uses Boyer–Lindquist coordinates

```text
x^mu = (t, rb, vartheta, varphi)
```

where rb = \bar{r} and with signature convention

```text
(-, +, +, +)
```

and symbolic parameters

```text
m, j, q, zeta
```

The auxiliary quantities defined in `Rotating_closed_string_vacuum.mpl` are

```text
mj     = sqrt(m^2 - j^2)
Sigma  = rb^2 + j^2*cos(vartheta)^2
Delta  = rb^2 - 2*m*rb + j^2
F      = (rb - m - mj)/(rb - m + mj)
E2phi  = F^(q/mj)*cos(zeta)^2 + F^(-q/mj)*sin(zeta)^2
omega  = -2*m*j*rb*sin(vartheta)^2/(Delta - j^2*sin(vartheta)^2)
```

The nonzero metric components are defined in the coordinate order

```text
1 = t, 2 = rb, 3 = vartheta, 4 = varphi
```

with line element

```text
ds^2 = g11 dt^2 + 2 g14 dt dvarphi + g44 dvarphi^2
     + g22 drb^2 + g33 dvartheta^2.
```

where

```text
g11 = -E2phi*(Delta - j^2*sin(vartheta)^2)/Sigma
g14 = -2*m*j*rb*E2phi*sin(vartheta)^2/Sigma
g22 = E2phi*Sigma*(1 + mj^2*sin(vartheta)^2/Delta)^(-q^2/mj^2)/Delta
g33 = E2phi*Sigma*(1 + mj^2*sin(vartheta)^2/Delta)^(-q^2/mj^2)
g44 = E2phi*((rb^2 + j^2)*Sigma + 2*m*j^2*rb*sin(vartheta)^2)*sin(vartheta)^2/Sigma
```

The metric file also sets

```text
Info = `rotating closed-string vacuum metric (JKP, 2026)`
```

## How to run the calculation

### Maple/GRTensorIII track

1. Install Maple and GRTensorIII.
2. Clone or download this repository.
3. Open `Rotating_Closed_String_Vacuum.mw` in Maple, or open the workbook `Rotating_Closed_String_Vacuum.maple`.
4. Make sure the working directory points to the repository folder containing `Rotating_closed_string_vacuum.mpl`.
5. Execute the worksheet/workbook from the beginning.

Alternatively, load the metric file directly from a Maple/GRTensorIII session by placing `Rotating_closed_string_vacuum.mpl` in the working directory or in a directory known to GRTensorIII, then loading the metric under the name

```text
Rotating_closed_string_vacuum
```

Depending on local GRTensorIII conventions, the metric name may be treated case-sensitively by the filesystem, so keep the file name unchanged.

### Mathematica track

1. Install Wolfram Mathematica.
2. Open `analytic_proof_of__qi coordinate.nb` and evaluate the cells from top to bottom.
3. The notebook verifies the boxed dilaton, $H_{(3)}$, and string-frame metric in Boyer–Lindquist coordinates against the NS–NS field equations, and then confirms the pull-back map to the quasi-isotropic chart.

## Expected output

**Maple/GRTensorIII.** The worksheet/workbook computes tensorial quantities for the metric. In particular, the saved output includes calculations of:

- the line element `ds`
- the determinant of the metric `detg`
- the inverse metric `g(up,up)`
- the Ricci tensor `R(dn,dn)`
- the Ricci scalar `Ricciscalar`

For the vacuum check, the Ricci tensor and Ricci scalar should simplify to zero under the assumptions and simplifications used in the worksheet/workbook.

**Mathematica.** The notebook evaluates each NS–NS field equation symbolically in Boyer–Lindquist coordinates and prints simplified residuals; under the stated assumptions all residuals reduce to zero.

## Reproducibility notes

Symbolic tensor calculations can be sensitive to assumptions, simplification order, and software versions. For best reproducibility:

- use the provided worksheet/workbook/notebook as the reference execution path;
- keep the metric and notebook file names unchanged;
- run each calculation from a clean Maple or Mathematica session;
- make sure GRTensorIII can find `Rotating_closed_string_vacuum.mpl`;
- record the Maple/GRTensorIII and Mathematica versions used in any derived work;
- cite the specific GitHub release corresponding to the version used.

## Citation

If you use these files, please cite the accompanying paper and this repository.

```bibtex
@software{Rotating_closed_string_vacuum_JKP,
  author       = {Hun Jang and Minkyoo Kim and Jeong-Hyuck Park},
  title        = {Exact Rotating Closed-String Vacuum: Maple/GRTensorIII and Mathematica Verification},
  year         = {2026},
  publisher    = {GitHub},
  version      = {v1.0.0},
  url          = {https://github.com/Jeong24th/Kerr_DFT}
}
```

## License

Distributed under the [MIT License](LICENSE).

## Contact

For questions, please contact:

- Hun Jang: hun.jang@nyu.edu
- Minkyoo Kim: minkyookim@snu.ac.kr
- Jeong-Hyuck Park: park@sogang.ac.kr
