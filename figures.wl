(* ============================================================================
   figures.wl  --  Mathematica reproduction of figures.py
   "An Exact Rotating Closed-String Vacuum" (quasi-isotropic coordinates)

   Produces (PDF, matching the matplotlib originals):
     fig_polar3d.pdf      -- Fig 2 LEFT : 3D cutaway (ergosphere + r+ shell + repulsion)
     fig_polar_radial.pdf -- Fig 2 RIGHT: axial radial profile of a^r for 4 values of |q|
     fig_effective_potential.pdf -- Appendix: axial timelike effective potential W(r)
     fig_extremal.pdf     -- SM1 : displaced Kerr bound on the extremal branch
     fig_nonsphere.pdf    -- SM2 : static-limit non-sphericity ratio (heatmap)
     fig_polar2d.pdf      -- SM3 : sign(F) two-tone (r,theta) map + F=0 contour

   Radial coordinate is the quasi-isotropic r (paper's native r, "no bar");
   the metric / F-function are evaluated in their natural Boyer-Lindquist (bar r)
   form via  bar r(r) = r + m + m_j^2/(4 r).
   Run:  wolframscript -file figures.wl
   ============================================================================ *)

(* ----------------------------- palette ----------------------------------- *)
steelblue  = RGBColor[70/255, 130/255, 180/255];
crimson    = RGBColor[220/255, 20/255, 60/255];
gold       = RGBColor[1, 215/255, 0];
darkorange = RGBColor[1, 140/255, 0];
firebrick  = RGBColor[178/255, 34/255, 34/255];
ergoTxt    = RGBColor[44/255, 95/255, 141/255];

(* ----------------------- closed-form helpers ----------------------------- *)
Frad[r_, m_, mj_] := (r - m - mj)/(r - m + mj);   (* BL auxiliary F(r) *)

(* dilaton phi(r), BL form -- left ungated so D[] can differentiate it *)
phiFunc[r_, m_, mj_, q_, zeta_] :=
  (1/2) Log[Frad[r, m, mj]^(q/mj) Cos[zeta]^2 + Frad[r, m, mj]^(-q/mj) Sin[zeta]^2];

(* analytic d phi/dr : differentiate ONCE symbolically, then substitute *)
dphidrExpr = D[phiFunc[rr, mm, mjj, qq, zz], rr];
dphidr[r_?NumericQ, m_?NumericQ, mj_?NumericQ, q_?NumericQ, zeta_?NumericQ] :=
  dphidrExpr /. {rr -> r, mm -> m, mjj -> mj, qq -> q, zz -> zeta};

(* radial-acceleration coefficient F(r,theta) -- sign sets repulsion *)
FaccelScalar[r_?NumericQ, theta_?NumericQ, m_?NumericQ, j_?NumericQ,
             q_?NumericQ, zeta_?NumericQ] :=
  Module[{mj = Sqrt[m^2 - j^2], rho2, Del, phip},
    rho2 = r^2 + j^2 Cos[theta]^2;
    Del  = (r - m)^2 - mj^2;
    phip = dphidr[r, m, mj, q, zeta];
    m (r^2 - j^2 Cos[theta]^2) + phip rho2 (Del - j^2 Sin[theta]^2)];

(* proper radial acceleration a^{hat r}, incl. off-axis K_perp^{-1/2} warp *)
aHatR[r_?NumericQ, theta_?NumericQ, m_?NumericQ, j_?NumericQ,
      q_?NumericQ, zeta_?NumericQ] :=
  Module[{mj = Sqrt[m^2 - j^2], rho2, Del, Lam, Kperp, phiv, Fv},
    rho2 = r^2 + j^2 Cos[theta]^2;
    Del  = (r - m)^2 - mj^2;
    Lam  = Del - j^2 Sin[theta]^2;                 (* > 0 in stationary region *)
    If[Del <= 0 || Lam <= 0, Return[0.0]];
    Kperp = (1 + mj^2 Sin[theta]^2/Del)^(-q^2/mj^2);
    phiv  = phiFunc[r, m, mj, q, zeta];
    Fv    = FaccelScalar[r, theta, m, j, q, zeta];
    Exp[-phiv] Sqrt[Del]/(Sqrt[Kperp] rho2^1.5 Lam) Fv];

(* quasi-isotropic  <->  Boyer-Lindquist *)
rOfRbar[rbar_, m_, mj_] := rbar + m + mj^2/(4 rbar);          (* qiso -> BL *)
rbarOfR[r_, m_, mj_]    := Module[{d = (r - m)^2 - mj^2},     (* BL -> qiso *)
                             0.5 ((r - m) + Sqrt[Max[d, 0]])];

(* theta-axis tick spec reused by the two heatmaps *)
thetaTicks = {{0, "0"}, {Pi/4, "\[Pi]/4"}, {Pi/2, "\[Pi]/2"},
              {3 Pi/4, "3\[Pi]/4"}, {Pi, "\[Pi]"}};

(* ===================== Fig 2 RIGHT : axial radial profile ================= *)
makePolarRadial[file_] := Module[
  {m = 1.0, j = 0.85, zeta = 0.5, mj, rbarPlus, rbars, qs, cols,
   curves = {}, marks = {}, q, av, ip, plt},
  mj = Sqrt[m^2 - j^2]; rbarPlus = mj/2;
  rbars = Subdivide[rbarPlus + 0.003, 4.0 m, 800];
  qs    = {0.5, 1.0, 1.5, 2.0};
  cols  = RGBColor /@ {"#1f77b4", "#2ca02c", "#e67e22", "#7b1fa2"};
  Do[
    q  = qs[[k]] mj;
    av = Table[m aHatR[rOfRbar[rb, m, mj], 0.0, m, j, q, zeta], {rb, rbars}];
    AppendTo[curves, Transpose[{rbars/mj, av}]];
    If[q > mj,                                   (* mark |a| peak (most negative) *)
      ip = First@Ordering[av, 1];
      If[2 <= ip <= Length[rbars] - 1,
         AppendTo[marks, {rbars[[ip]]/mj, av[[ip]]}]]],
    {k, Length[qs]}];
  plt = ListLinePlot[curves,
    PlotStyle  -> (Directive[#, AbsoluteThickness[1.6]] & /@ cols),
    PlotRange  -> {{0.4, 3.5}, {-0.5, 1.5}},
    Frame -> True, Axes -> False, AspectRatio -> 3.0/3.5,
    FrameLabel -> {Row[{"r/", Subscript["m", "j"]}],
                   Row[{"m ", Superscript["a", "r\:0302"], "(r, \[CurlyTheta]=0)"}]},
    FrameTicks -> {{Automatic, None},
                   {{{0.5, "1/2"}, {1, "1"}, {2, "2"}, {3, "3"}}, None}},
    GridLines  -> {{{0.5, Directive[Black, Dotted]}}, {{0, Black}}},
    PlotLegends -> Placed[
       LineLegend[cols, (Row[{"|q|/", Subscript["m", "j"], " = ", #}] & /@ qs),
         LegendMarkerSize -> 18, LegendLayout -> "Column"], {0.74, 0.76}],
    Epilog -> {darkorange, AbsolutePointSize[6], Point /@ marks,
               Black, Text[Subscript["r", "+"], {0.5, -0.42}]},
    PlotLabel -> Style["Axial acceleration; dots at |a| peak", 10],
    ImageSize -> 380];
  Export[file, plt]; Print["Wrote ", file]; plt];

(* ===================== Appendix : axial effective potential =============== *)
makeEffectivePotential[file_] := Module[
  {m = 1.0, j = 0.85, zeta = 0.5, mj, xs, qs, cols, wAtX, curves, plt},
  mj = Sqrt[m^2 - j^2];
  xs = 0.5 + Exp[Subdivide[Log[2.0 10^-4], Log[5.5], 900]];
  qs = {0.5, 1.0, 1.5};
  cols = RGBColor /@ {"#1f77b4", "#2ca02c", "#e67e22"};
  wAtX[x_?NumericQ, qRatio_?NumericQ] := Module[{rq, rb, del, sig, q},
    rq = x mj; q = qRatio mj;
    rb = rOfRbar[rq, m, mj];
    del = (rq - mj^2/(4 rq))^2;
    sig = rb^2 + j^2;
    Exp[2 phiFunc[rb, m, mj, q, zeta]] del/sig];
  curves = Table[Table[{x, wAtX[x, qr]}, {x, xs}], {qr, qs}];
  plt = ListLinePlot[curves,
    PlotStyle -> (Directive[#, AbsoluteThickness[1.7]] & /@ cols),
    ScalingFunctions -> {None, "Log10"},
    PlotRange -> {{0.48, 6.0}, {10^-3, 10^2}},
    Frame -> True, Axes -> False, AspectRatio -> 2.8/3.7,
    FrameLabel -> {Row[{"r/", Subscript["m", "j"]}],
                   Row[{"W(r) = -", Subscript["g", "tt"], "|\[CurlyTheta]=0"}]},
    FrameTicks -> {{Automatic, None},
                   {{{0.5, "1/2"}, {1, "1"}, {2, "2"}, {4, "4"}, {6, "6"}}, None}},
    GridLines -> {{{0.5, Directive[Black, Dotted]}},
                  {{1, Directive[GrayLevel[0.45], Dashed]}}},
    PlotLegends -> Placed[
      LineLegend[cols, (Row[{"|q|/", Subscript["m", "j"], " = ", #}] & /@ qs),
        LegendMarkerSize -> 18, LegendLayout -> "Column"], {0.75, 0.76}],
    Epilog -> {Black, Text[Subscript["r", "+"], {0.5, 1.5 10^-3}],
               GrayLevel[0.35], Text["W \[Rule] 1", {5.7, 1.08}]},
    PlotLabel -> Style["Axial timelike effective potential", 10],
    ImageSize -> 400];
  Export[file, plt]; Print["Wrote ", file]; plt];

(* ===================== SM3 : sign(F) two-tone polar map =================== *)
makePolar2D[file_] := Module[
  {m = 1.0, j = 0.85, zeta = 0.5, mj, rbarPlus, q, Fv, dens, cont, plt},
  mj = Sqrt[m^2 - j^2]; rbarPlus = mj/2; q = 1.5 mj;
  Fv[rb_?NumericQ, th_?NumericQ] := FaccelScalar[rOfRbar[rb, m, mj], th, m, j, q, zeta];
  dens = DensityPlot[Fv[rb, th],
     {rb, rbarPlus + 0.01, 3.5 m}, {th, 0.001, Pi - 0.001},
     ColorFunctionScaling -> False,
     ColorFunction -> (If[# < 0, RGBColor["#3b6db5"], RGBColor["#c83838"]] &),
     PlotPoints -> 90, MaxRecursion -> 1, Frame -> True,
     FrameLabel -> {"r/m", "\[CurlyTheta]"},
     FrameTicks -> {{thetaTicks, None}, {Automatic, None}},
     AspectRatio -> 3.4/5, PlotRangePadding -> None];
  cont = ContourPlot[Fv[rb, th] == 0,
     {rb, rbarPlus + 0.01, 3.5 m}, {th, 0.001, Pi - 0.001},
     ContourStyle -> Directive[Black, AbsoluteThickness[1.6]], PlotPoints -> 90];
  plt = Legended[
     Show[dens, cont,
       Epilog -> {Black, Dotted, AbsoluteThickness[0.8],
         Line[{{rbarPlus, 0}, {rbarPlus, Pi}}],
         Text[Subscript["r", "+"], {rbarPlus, Pi}, {0, -1.4}],
         Style[Text["repulsive",  {0.55, 0.28}], 9, Bold, Black],
         Style[Text["repulsive",  {0.55, Pi - 0.28}], 9, Bold, Black],
         Style[Text["attractive", {2.6, Pi/2}], 9, Bold, White]}],
     SwatchLegend[{RGBColor["#3b6db5"], RGBColor["#c83838"]},
       {"\[ScriptCapitalF] < 0 (repulsive)", "\[ScriptCapitalF] > 0 (attractive)"}]];
  Export[file, plt]; Print["Wrote ", file]; plt];

(* ===================== SM2 : non-sphericity ratio ========================= *)
makeNonsphere[file_] := Module[
  {m = 1.0, q = 0.5, ratio, dens, cont, plt},
  ratio[rb_?NumericQ, th_?NumericQ] := Module[{R, d},
     R = rb + m + m^2/(4 rb);                     (* j=0 : r(bar r) *)
     d = Max[(R - m)^2 - m^2, 1.*^-12];
     (1 + m^2 Sin[th]^2/d)^(-q^2/m^2)];
  dens = DensityPlot[ratio[rb, th],
     {rb, 0.55 m, 3.0 m}, {th, 0.02, Pi - 0.02},
     ColorFunction -> "Viridis", PlotPoints -> 90, MaxRecursion -> 1,
     Frame -> True, FrameLabel -> {"r/m", "\[CurlyTheta]"},
     FrameTicks -> {{thetaTicks, None}, {Automatic, None}},
     PlotLegends -> Automatic, AspectRatio -> 3.4/5,
     PlotLabel -> Style["Static-limit non-sphericity ratio at q/m = 0.5", 10]];
  cont = ContourPlot[ratio[rb, th],
     {rb, 0.55 m, 3.0 m}, {th, 0.02, Pi - 0.02},
     Contours -> {0.5, 0.7, 0.85, 0.95, 0.99}, ContourShading -> None,
     ContourStyle -> Directive[White, Opacity[0.7], AbsoluteThickness[0.6]],
     PlotPoints -> 90];
  plt = Show[dens, cont];
  Export[file, plt]; Print["Wrote ", file]; plt];

(* ===================== SM1 : displaced Kerr bound ========================= *)
makeExtremal[file_] := Module[
  {m = 1.0, spin, fillNeg, fillPos, curve, plt},
  spin[x_] := 1/(1 + x)^2;                         (* x = q cos2zeta / m *)
  fillNeg = Plot[spin[x], {x, -1, 0}, Filling -> Axis,
     FillingStyle -> Directive[firebrick, Opacity[0.18]],
     PlotStyle -> Directive[Opacity[0]]];
  fillPos = Plot[spin[x], {x, 0, 1}, Filling -> Axis,
     FillingStyle -> Directive[steelblue, Opacity[0.18]],
     PlotStyle -> Directive[Opacity[0]]];
  curve = Plot[spin[x], {x, -1, 1},
     PlotStyle -> Directive[Black, AbsoluteThickness[1.4]]];
  plt = Show[fillNeg, fillPos, curve,
     PlotRange -> {{-1, 1}, {0, 30}}, Frame -> True, Axes -> False,
     AspectRatio -> 2.8/4.0,
     FrameLabel -> {Row[{Style["q", Italic], " cos 2\[Zeta] / m"}],
                    Row[{"J/", Superscript["(MG)", "2"]}]},
     GridLines -> {{{0, Directive[Gray, Dashed]}}, {{1, Directive[Gray, Dashed]}}},
     PlotLabel -> Style["Extremal branch: displaced Kerr bound", 10],
     Epilog -> {
        Text[Style[Row[{"J/", Superscript["(MG)", "2"], " > 1"}], 9, firebrick],  {-0.5, 20}],
        Text[Style[Row[{"J/", Superscript["(MG)", "2"], " < 1"}], 9, steelblue], { 0.5,  5}]},
     ImageSize -> 420];
  Export[file, plt]; Print["Wrote ", file]; plt];

(* ===================== Fig 2 LEFT : 3D cutaway ============================ *)
makePolar3D[file_] := Module[
  {m = 1.0, j = 0.85, zeta = 0.5, mj, q, rbarPlus, rbarErgoEq,
   rbarScan, aScan, rbarPeak, ext, ergo, shell, repCond, rep, annot, plt},
  mj = Sqrt[m^2 - j^2]; q = 1.5 mj; rbarPlus = mj/2;
  rbarErgoEq = rbarOfR[2.0 m, m, mj];              (* equator ergo, qiso *)
  rbarScan = Subdivide[rbarPlus + 0.003, 3.5 m, 3000];
  aScan    = Table[aHatR[rOfRbar[rb, m, mj], 0.0, m, j, q, zeta], {rb, rbarScan}];
  rbarPeak = rbarScan[[First@Ordering[aScan, 1]]]; (* axial repulsion peak *)
  ext = 1.25 m;

  (* (1) ergosurface : Kerr-shaped r_ergo(theta) mapped to qiso *)
  ergo = ParametricPlot3D[
    With[{rb = rbarOfR[m + Sqrt[m^2 - j^2 Cos[t]^2], m, mj]},
      {rb Sin[t] Cos[p], rb Sin[t] Sin[p], rb Cos[t]}],
    {t, 0, Pi}, {p, Pi/2, 3 Pi/2},
    PlotStyle -> Directive[steelblue, Opacity[0.22]],
    Mesh -> None, PlotPoints -> {50, 40}];

  (* (2) outer shell at constant qiso r+ = mj/2 *)
  shell = ParametricPlot3D[
    {rbarPlus Sin[t] Cos[p], rbarPlus Sin[t] Sin[p], rbarPlus Cos[t]},
    {t, 0, Pi}, {p, Pi/2, 3 Pi/2},
    PlotStyle -> Directive[crimson, Opacity[1]],
    Mesh -> None, PlotPoints -> {40, 30}];

  (* (3) polar repulsion solid : F<0, cut away on x<=0 half *)
  repCond[x_, y_, z_] := Module[{rb = Sqrt[x^2 + y^2 + z^2], th},
     If[rb <= rbarPlus, Return[False]];
     th = ArcCos[Clip[z/rb, {-1, 1}]];
     FaccelScalar[rOfRbar[rb, m, mj], th, m, j, q, zeta] < 0];
  rep = RegionPlot3D[repCond[x, y, z],
     {x, -ext, 0}, {y, -ext, ext}, {z, -ext, ext},
     PlotStyle -> Directive[gold, Opacity[0.7]],
     Mesh -> None, PlotPoints -> 30, MaxRecursion -> 1, BoundaryStyle -> None];

  (* annotations *)
  annot = Graphics3D[{
     Black, AbsoluteThickness[2],
     Line[{{0, 0, -1.2 m}, {0, 0, 1.2 m}}],
     Cone[{{0, 0, 1.2 m},  {0, 0, 1.32 m}},  0.04 m],
     Cone[{{0, 0, -1.2 m}, {0, 0, -1.32 m}}, 0.04 m],
     Text[Style["axis of rotation", 12, Black], {0, 0, 1.45 m}],
     darkorange, Sphere[{0, 0, rbarPeak}, 0.03 m], Sphere[{0, 0, -rbarPeak}, 0.03 m],
     steelblue, Dotted, AbsoluteThickness[1],
     Line[Table[{rbarErgoEq Cos[a], rbarErgoEq Sin[a], 0}, {a, 0, 2 Pi, 2 Pi/120}]],
     Text[Style["repulsion (\[ScriptCapitalF]<0)", 11, darkorange, Bold], {0, -0.95 m, 0.95 m}],
     Text[Style["repulsion peak", 11, darkorange, Bold], {0, 0.95 m, 1.05 m}],
     Text[Style["ergosphere", 11, ergoTxt, Bold], {0, -1.1 m, -0.7 m}],
     Text[Style[Row[{"outer shell (", Subscript["r", "+"], ")"}], 11, crimson, Bold],
          {0, 0.85 m, -0.65 m}]}];

  plt = Show[ergo, shell, rep, annot,
     Boxed -> False, Axes -> False,
     PlotRange -> {{-ext, ext}, {-ext, ext}, {-ext, ext}},
     ViewPoint -> {2.6, 0.7, 0.9}, ViewProjection -> "Orthographic",
     ViewVertical -> {0, 0, 1}, Lighting -> "Neutral", ImageSize -> 560];
  Export[file, plt]; Print["Wrote ", file]; plt];

(* ------------------------------- main ------------------------------------ *)
makePolar3D["fig_polar3d.pdf"];
makePolarRadial["fig_polar_radial.pdf"];
makeEffectivePotential["fig_effective_potential.pdf"];
makeExtremal["fig_extremal.pdf"];
makeNonsphere["fig_nonsphere.pdf"];
makePolar2D["fig_polar2d.pdf"];
Print["--- all figures generated ---"];
