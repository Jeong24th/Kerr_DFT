# Rotating Closed String Vacuum

This repository contains three independent verification tracks — a Maple/GRTensorIII worksheet/workbook, a pair of Mathematica notebooks (analytic and numerical), and a Python/SymPy suite — for the exact rotating closed-string NS–NS vacuum solution of the Einstein Double Field Equations in double field theory. A Wolfram Language script (`figures.wl`) reproduces the five figures of the paper.

**Coordinate convention.** All three tracks perform the verification natively in Boyer–Lindquist coordinates $(t,\bar r,\vartheta,\varphi)$, matching the solution as written in the Supplemental Material of the accompanying paper. The quasi-isotropic chart $(t,r,\vartheta,\varphi)$ used in the body of the paper follows by the explicit pull-back $\bar r(r)=r+m+(m^{2}-j^{2})/(4r)$; the Mathematica notebook additionally checks this pull-back map.

**Notation.** The paper defines $\Sigma_{\pm}(\bar r,\vartheta)\equiv\bar r^{2}\pm j^{2}\cos^{2}\vartheta$. The verification code retains the legacy single-symbol $\Sigma$ for the paper's $\Sigma_{+}$ (the standard Kerr combination $\bar r^{2}+j^{2}\cos^{2}\vartheta$); the minus-variant $\Sigma_{-}$ does not enter the field-equation residuals computed here.

The Maple/GRTensorIII track supplies a GRTensorIII-compatible metric file together with Maple worksheet/workbook files that perform symbolic tensor calculations — line element, determinant, inverse metric, Ricci tensor, and Ricci scalar. For the NS–NS vacuum check, the field equation residuals — the modified Einstein tensor $R_{\mu\nu}-2\nabla_{\mu}\nabla_{\nu}\phi+\tfrac{1}{4}H_{\mu\alpha\beta}H_{\nu}{}^{\alpha\beta}$, the dilaton equation, and the H-flux Bianchi/EOM — simplify to zero. (The bare Ricci tensor and scalar of the string-frame metric are nonzero; they are sourced by the dilaton and H-flux.) The Mathematica notebook performs an independent symbolic verification of the dilaton, Kalb–Ramond three-form, and string-frame metric directly against the NS–NS field equations in Boyer–Lindquist coordinates.

## Repository contents

| File | Description |
| --- | --- |
| `Rotating_closed_string_vacuum.mpl` | GRTensorIII metric definition file (Boyer–Lindquist chart). Defines coordinates, signature, parameters, auxiliary functions, and nonzero metric components. |
| `Rotating_Closed_String_Vacuum.mw` | Maple worksheet containing the symbolic GRTensorIII calculation and saved outputs. |
| `Rotating_Closed_String_Vacuum.maple` | Maple workbook version of the calculation file. |
| `analytic_proof_of__qi coordinate.nb` | Mathematica notebook (analytic track). Checks the boxed dilaton, $H_{(3)}$, and string-frame metric directly against the NS–NS field equations in Boyer–Lindquist coordinates, then verifies the pull-back to quasi-isotropic coordinates. |
| `Numerical_check_for_qi coordinate.nb` | Mathematica notebook (numerical track). Independent numerical verification of the same NS–NS field equations on the quasi-isotropic chart, complementing the analytic notebook. |
| `figures.wl` | Wolfram Language script reproducing the five figures of the paper (`fig_polar3d.pdf`, `fig_polar_radial.pdf`, `fig_extremal.pdf`, `fig_nonsphere.pdf`, `fig_polar2d.pdf`) from the closed-form solution. Run with `wolframscript -file figures.wl`. |
| `verify_einstein_eq.py` | Third independent track (SymPy/Python). Verifies the Bogush–Galtsov Einstein-frame seed satisfies $R^E_{\mu\nu} - 2\partial_\mu\phi_{\rm seed}\partial_\nu\phi_{\rm seed} = 0$ at multiple generic numerical test points. |
| `verify_string_ricci.py` | Verifies the static-limit ($j=0$) string-frame Ricci-scalar formula $R\|_{j=0} = -(2e^{-2\phi}+3e^{-6\phi}\sin^{2}2\zeta)\,R^E\|_{j=0}$ from the SM. |
| `verify_NS_NS.py` | Full NS–NS vacuum field equation verification (graviton + B-field + dilaton) at $q \neq 0$ generic rotating parameters. |

## Requirements

The files were prepared for use with:

- [Maple](https://www.maplesoft.com/products/Maple/) (worksheet metadata indicates Maple 2026.1)
- [GRTensorIII](https://github.com/grtensor/grtensor)
- [Wolfram Mathematica](https://www.wolfram.com/mathematica/) (the notebook was prepared with Wolfram 14.3)
- [Python 3](https://www.python.org/) (3.10+) with [SymPy](https://www.sympy.org/) (1.13+) for the `.py` tracks.

Other recent Maple/GRTensorIII/Mathematica versions may also work, but symbolic-simplification behaviour and notebook rendering can vary by version. The Python tracks were tested with Python 3.12 / SymPy 1.14 and require only the SymPy standard library.

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
4. For an independent numerical confirmation, open `Numerical_check_for_qi coordinate.nb` and evaluate the cells from top to bottom; the notebook performs the same checks numerically on the quasi-isotropic chart.

### Figure reproduction

1. Install Wolfram Mathematica (or a Wolfram Engine with `wolframscript`).
2. From the repository directory, run:

   ```bash
   wolframscript -file figures.wl
   ```

3. The script produces the five PDF figures of the paper directly from the closed-form solution. Output filenames are listed in the header of `figures.wl`.

### Python/SymPy track

1. Install Python 3 and SymPy (`pip install sympy`).
2. From the repository directory, run any of:

   ```bash
   python verify_einstein_eq.py
   python verify_string_ricci.py
   python verify_NS_NS.py
   ```

3. Each script builds the metric (and B-field / dilaton where relevant) symbolically, computes Christoffel / Ricci / NS–NS residuals, then evaluates at three independent generic test points and reports residuals. Expected output: all residuals zero to ~140 decimal digits (machine precision) and a final `PASS` line.

## Expected output

**Maple/GRTensorIII.** The worksheet/workbook computes tensorial quantities for the metric. In particular, the saved output includes calculations of:

- the line element `ds`
- the determinant of the metric `detg`
- the inverse metric `g(up,up)`
- the Ricci tensor `R(dn,dn)`
- the Ricci scalar `Ricciscalar`

For the NS–NS vacuum check, the worksheet/workbook then combines these with the dilaton and Kalb–Ramond data to form the NS–NS field equation residuals — the modified Einstein tensor $R_{\mu\nu}-2\nabla_{\mu}\nabla_{\nu}\phi+\tfrac{1}{4}H_{\mu\alpha\beta}H_{\nu}{}^{\alpha\beta}$, the dilaton equation, and the H-flux Bianchi/EOM — which simplify to zero. The bare Ricci tensor and scalar of the string-frame metric are nonzero; they are sourced by the dilaton and H-flux.

**Mathematica.** The notebook evaluates each NS–NS field equation symbolically in Boyer–Lindquist coordinates and prints simplified residuals; under the stated assumptions all residuals reduce to zero.

**Python/SymPy.** Each `verify_*.py` script prints, at multiple independent test points,

- `verify_einstein_eq.py` — the 10 independent components of $R^E_{\mu\nu}-2\partial_\mu\phi_{\rm seed}\partial_\nu\phi_{\rm seed}$, all evaluating to zero;
- `verify_string_ricci.py` — the directly computed string-frame Ricci scalar at $j=0$ and the PRL-formula prediction, with zero absolute and relative difference;
- `verify_NS_NS.py` — the 10 graviton-equation components, 6 B-field-equation components, and the dilaton trace, all zero.

A final `PASS` line summarises the maximum residual.

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
