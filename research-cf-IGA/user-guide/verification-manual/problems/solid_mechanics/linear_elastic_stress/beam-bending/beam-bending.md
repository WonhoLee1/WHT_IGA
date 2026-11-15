# Cantilever beam bending – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/linear_elastic_stress/beam-bending/beam-bending.html

---

Verification Manual
Solid mechanics
Linear elastic stress
Cantilever beam bending
Cantilever beam bending
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
\[
\newcommand{\gravityLineLoad}{784.8}
\newcommand{\pressureLineLoad}{500}
\]
1
Objective
Verify Coreform IGA for Abaqus support for the common Abaqus load and boundary-condition types on a single problem whose solution is known in closed form: a slender cantilever beam bent by gravity, body force, surface pressure (uniform and total-force), surface traction (including the follower and resultant variants), and a concentrated end moment. Each load case is checked against the corresponding Euler-Bernoulli beam solution, with pure-Abaqus body-fitted mesh-refinement sweeps retained as independent comparison points. The quantity of interest is the tip deflection; the deflection curve and the bending- and shear-stress distributions are also compared against the closed forms along line probes. The test artifact records the unit system as SI: m, kg, s, N, Pa; all dimensional inputs and results therefore use the
SI consistent unit system
.
2
Geometry
A square-section prismatic beam of width
\(b\)
, height
\(h\)
, and length
\(L\)
is shown in
Figure 1
, with its dimensions provided in
Table 1
. The cross-section is centered on the
\(x\)
/
\(y\)
origin and the beam axis runs along
\(z \in [0, L]\)
. Verification quantities are sampled along three line probes: the beam axis (
\(u_y\)
), the top fiber at
\(y = h/2\)
(
\(\sigma_{zz}\)
), and the midspan section line at
\(z = L/2\)
(
\(\sigma_{zz}\)
and
\(\tau_{zy}\)
).
Figure 1: Beam dimensions and coordinate system.
Table 1: Geometric parameters.
Quantity
Symbol
Value
Section width
\(b\)
\(0.1\ \text{m}\)
Section height
\(h\)
\(0.1\ \text{m}\)
Length
\(L\)
\(1\ \text{m}\)
Second moment of area
\(I_x = bh^3/12\)
\(8.33 \times 10^{-6}\ \text{m}^4\)
3
Material
The beam is a homogeneous, isotropic, linear-elastic material (
Table 2
). Poisson’s ratio is deliberately zero so the 3D solid matches the Euler-Bernoulli kinematics, which carry no Poisson effect.
Table 2: Material properties.
Property
Symbol
Value
Young’s modulus
\(E\)
\(200 \times 10^{9}\ Pa\)
Poisson’s ratio
\(\nu\)
0
Mass density
\(\rho\)
\(8.00 \times 10^{3}\ kg/m³\)
4
Loading & boundary conditions
The idealized beam is clamped at
\(z = 0\)
, as shown in
Figure 2
. In the three-dimensional solid models, this condition is imposed by pinning the
\(z = 0\)
face (all translations fixed). Each load case below is applied in its own study; all other model inputs are identical across cases.
Figure 2: Idealized clamped boundary condition at
\(z = 0\)
.
4.1
Gravity
A gravity load of magnitude
\(g\)
acts in
\(-y\)
on the whole body, as shown in
Figure 3
. The beam therefore carries its self-weight as the distributed line load
\(w = \rho g b h\)
.
Figure 3: Gravity loading on a beam of density
\(\rho\)
.
4.2
Body force
A uniform body force per unit volume of magnitude
\(F=\rho g\)
acts in
\(-y\)
, as shown in
Figure 4
. This load is statically identical to the gravity case but exercises the body-force load type rather than the gravity type.
Figure 4: Uniform body-force loading in the
\(-y\)
direction.
4.3
Pressure
A normal pressure is applied to the top face (
\(y = +h/2\)
), as shown in
Figure 5
. Two variants are exercised:
uniform
, where the magnitude is the pressure itself, and
total force
, where the magnitude is the total force distributed over the face. The magnitudes are chosen so both reduce to the same line load
\(w\)
.
Figure 5: Top-face pressure loading represented by the equivalent line load
\(w\)
.
4.4
Surface traction
A uniform shear traction on the end face (
\(z = L\)
) acts in
\(-y\)
, giving the resultant tip load
\(P = t\,bh\)
shown in
Figure 6
. Two additional coverage variants exercise the traction options:
follower
(
follower=ON
), for which the traction direction follows the deforming surface, and
resultant
(
resultant=ON
), for which Abaqus maintains a constant resultant by integrating the specified traction over the reference area. The
resultant
option does not reinterpret the entered traction magnitude as a total force; its magnitude retains pressure units. Because this is a small-displacement analysis, all three traction variants produce the same solution.
Figure 6: End-face shear traction represented by the resultant tip load
\(P\)
.
4.5
End moment
A concentrated moment
\(M_0\)
about
\(x\)
is applied at the beam tip, as shown in
Figure 7
. In the three-dimensional solid models, it is applied to a reference point at the centroid of the end face and transmitted through a kinematic coupling to the whole face.
Figure 7: Concentrated end moment
\(M_0\)
about the
\(x\)
axis.
Table 3: Applied load magnitudes.
Load case
Test case
Magnitude
Gravity
gravity
\(9.81\ m/s²\)
Body force
body_force
\(78.5 \times 10^{3}\ N/m³\)
Pressure (uniform)
pressure_uniform
\(5.00 \times 10^{3}\ Pa\)
Pressure (total force)
pressure_total_force
\(500\ N\)
Surface traction
traction
\(500\ Pa\)
End moment
moment
\(5\ N·m\)
5
Mesh & discretization
Coreform IGA discretizes the beam with
immersed, trimmed U-splines
: a structured rectilinear background grid with
\(p\)
padding elements per side is trimmed to the beam. The background element sizes are
\((0.1, 0.1, 0.25)/2^{k}\)
in
\((x, y, z)\)
for refinement levels
\(k = 0,\dots,4\)
, at the evidence-backed degrees
\(p = 1\)
,
\(p = 2\)
,
\(p = 3\)
with maximum continuity
\(C^{p-1}\)
. Each primary load case runs on this same mesh sweep; the follower and resultant traction variants are single-level coverage cases. The meshes are therefore shown once — the tabs below render each refinement level from the gravity study’s IGA results databases (trimmed Bezier element edges in carnation, CAD edges in black).
\(h_z\)
=0.25
\(h_z\)
=0.125
\(h_z\)
=0.0625
\(h_z\)
=0.03125
\(h_z\)
=0.015625
6
Reference & accepted solutions
The
accepted
solutions are the Euler-Bernoulli cantilever closed forms
[1]
. The distributed cases (gravity, body force, both pressures) reduce to a uniform line load
\(w\)
:
\[
V(z) = w\,(L - z), \qquad
M(z) = -\tfrac{w}{2}\,(L - z)^2, \qquad
u_y(z) = -\frac{w z^2}{24 E I_x}\left(6L^2 - 4Lz + z^2\right).
\]
The traction cases reduce to a tip load
\(P\)
:
\[
V(z) = P, \qquad
M(z) = -P\,(L - z), \qquad
u_y(z) = -\frac{P z^2\,(3L - z)}{6 E I_x},
\]
and the end moment
\(M_0\)
gives constant moment and zero shear:
\[
V(z) = 0, \qquad
M(z) = -M_0, \qquad
u_y(z) = -\frac{M_0 z^2}{2 E I_x}.
\]
Stresses are recovered from the section resultants with
\(Q(y) = \tfrac{b}{2}\left(\tfrac{h^2}{4} - y^2\right)\)
:
\[
\sigma_{zz}(y, z) = -\frac{M(z)\, y}{I_x}, \qquad
\tau_{zy}(y, z) = \frac{V(z)\, Q(y)}{I_x\, b}.
\]
Euler-Bernoulli is itself an approximation of the 3D solid — it neglects shear deformation, and the pinned end face is stiffer than the idealized support — so exact agreement is not expected, and the test gates each study on the tip deflection matching within 5% on the asymptotic refinement levels (relaxed to 10% for the slowest-converging studies; the pre-asymptotic coarse
\(p = 1\)
meshes, which lock in bending, still solve and report but are not gated). To separate that modeling gap from discretization error, each load case also includes independent pure-Abaqus body-fitted refinement sweeps using C3D8 and C3D20 elements. The reference value reported in the accepted-solution and tip-deflection tables is specifically the finest-level C3D20 result. The contour galleries retain the full C3D8 and C3D20 sweeps so their discretization behavior remains visible. Coreform IGA and the body-fitted models discretize the same elasticity problem, so they should agree more closely with each other than either does with Euler-Bernoulli.
Table 4: Equivalent loads and tip deflections: Euler-Bernoulli closed form and the finest body-fitted C3D20 result.
Load case
Equivalent load
Euler-Bernoulli
\(u_y(L)\)
Finest C3D20
\(u_y(L)\)
Gravity
\(w\)
=
\(785\)
N/m
\(-58.9 \times 10^{-6}\)
m
\(-59.3 \times 10^{-6}\)
m
Body force
\(w\)
=
\(785\)
N/m
\(-58.9 \times 10^{-6}\)
m
\(-59.3 \times 10^{-6}\)
m
Pressure (uniform)
\(w\)
=
\(500\)
N/m
\(-37.5 \times 10^{-6}\)
m
\(-37.8 \times 10^{-6}\)
m
Pressure (total force)
\(w\)
=
\(500\)
N/m
\(-37.5 \times 10^{-6}\)
m
\(-37.8 \times 10^{-6}\)
m
Surface traction
\(P\)
=
\(5\)
N
\(-1.00 \times 10^{-6}\)
m
\(-1.01 \times 10^{-6}\)
m
End moment
\(M_0\)
=
\(5\)
N·m
\(-1.50 \times 10^{-6}\)
m
\(-1.50 \times 10^{-6}\)
m
Download data (CSV)
7
Results
Each load case below reports the tip deflection against the Euler-Bernoulli value and the finest C3D20 reference result (rows are the
\(h\)
-sweep; one tab per degree, and the degree tabs switch together across load cases), the four verification profiles against Euler-Bernoulli (one figure per probe quantity, a subfigure per degree), full-span line-probe absolute
\(L^2\)
self-convergence against the same-degree
\(h_4\)
profile, and the axial-stress contours for every test case. The contour plots show the Coreform IGA sweep and the full C3D8/C3D20 Abaqus sweeps as separate synchronized tab groups on a shared color scale capped at the Euler-Bernoulli peak.
7.1
Gravity
Euler-Bernoulli tip deflection:
\(-58.9 \times 10^{-6}\)
m; finest C3D20 Abaqus reference:
\(-59.3 \times 10^{-6}\)
m.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-17.8 \times 10^{-6}\)
\(69.7\%\)
\(70\%\)
0.125
\(-31.3 \times 10^{-6}\)
\(46.8\%\)
\(47.2\%\)
0.0625
\(-48.2 \times 10^{-6}\)
\(18.2\%\)
\(18.8\%\)
0.03125
\(-56.4 \times 10^{-6}\)
\(4.14\%\)
\(4.9\%\)
0.015625
\(-58.8 \times 10^{-6}\)
\(0.119\%\)
\(0.909\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-58.9 \times 10^{-6}\)
\(0.117\%\)
\(0.675\%\)
0.125
\(-59.6 \times 10^{-6}\)
\(1.34\%\)
\(0.538\%\)
0.0625
\(-59.6 \times 10^{-6}\)
\(1.29\%\)
\(0.486\%\)
0.03125
\(-59.6 \times 10^{-6}\)
\(1.24\%\)
\(0.437\%\)
0.015625
\(-59.6 \times 10^{-6}\)
\(1.21\%\)
\(0.412\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-62.0 \times 10^{-6}\)
\(5.32\%\)
\(4.49\%\)
0.125
\(-59.8 \times 10^{-6}\)
\(1.53\%\)
\(0.728\%\)
0.0625
\(-59.6 \times 10^{-6}\)
\(1.33\%\)
\(0.532\%\)
0.03125
\(-59.6 \times 10^{-6}\)
\(1.26\%\)
\(0.455\%\)
0.015625
\(-59.6 \times 10^{-6}\)
\(1.22\%\)
\(0.418\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 8: Gravity: deflection along the axis — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 9: Gravity: bending stress along the top fiber — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 10: Gravity: bending stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 11: Gravity: shear stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a) Axis deflection
\(u_y\)
(b) Top-fiber stress
\(\sigma_{zz}\)
(c) Midspan stress
\(\sigma_{zz}\)
(d) Midspan shear
\(\tau_{zy}\)
Figure 12: Gravity: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h4 profile; the h4 point is omitted because its self-error is zero.
7.1.1
Contour plots
Coreform IGA for Abaqus
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
h0p3
h1p3
h2p3
h3p3
h4p3
Abaqus reference
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
7.2
Body force
Euler-Bernoulli tip deflection:
\(-58.9 \times 10^{-6}\)
m; finest C3D20 Abaqus reference:
\(-59.3 \times 10^{-6}\)
m.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-17.8 \times 10^{-6}\)
\(69.7\%\)
\(70\%\)
0.125
\(-31.3 \times 10^{-6}\)
\(46.8\%\)
\(47.2\%\)
0.0625
\(-48.2 \times 10^{-6}\)
\(18.2\%\)
\(18.8\%\)
0.03125
\(-56.4 \times 10^{-6}\)
\(4.14\%\)
\(4.9\%\)
0.015625
\(-58.8 \times 10^{-6}\)
\(0.119\%\)
\(0.909\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-58.9 \times 10^{-6}\)
\(0.117\%\)
\(0.675\%\)
0.125
\(-59.6 \times 10^{-6}\)
\(1.34\%\)
\(0.538\%\)
0.0625
\(-59.6 \times 10^{-6}\)
\(1.29\%\)
\(0.486\%\)
0.03125
\(-59.6 \times 10^{-6}\)
\(1.24\%\)
\(0.437\%\)
0.015625
\(-59.6 \times 10^{-6}\)
\(1.21\%\)
\(0.412\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-62.0 \times 10^{-6}\)
\(5.32\%\)
\(4.49\%\)
0.125
\(-59.8 \times 10^{-6}\)
\(1.53\%\)
\(0.728\%\)
0.0625
\(-59.6 \times 10^{-6}\)
\(1.33\%\)
\(0.532\%\)
0.03125
\(-59.6 \times 10^{-6}\)
\(1.26\%\)
\(0.455\%\)
0.015625
\(-59.6 \times 10^{-6}\)
\(1.22\%\)
\(0.418\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 13: Body force: deflection along the axis — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 14: Body force: bending stress along the top fiber — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 15: Body force: bending stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 16: Body force: shear stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a) Axis deflection
\(u_y\)
(b) Top-fiber stress
\(\sigma_{zz}\)
(c) Midspan stress
\(\sigma_{zz}\)
(d) Midspan shear
\(\tau_{zy}\)
Figure 17: Body force: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h4 profile; the h4 point is omitted because its self-error is zero.
7.2.1
Contour plots
Coreform IGA for Abaqus
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
h0p3
h1p3
h2p3
h3p3
h4p3
Abaqus reference
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
7.3
Pressure (uniform)
Euler-Bernoulli tip deflection:
\(-37.5 \times 10^{-6}\)
m; finest C3D20 Abaqus reference:
\(-37.8 \times 10^{-6}\)
m.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-11.4 \times 10^{-6}\)
\(69.7\%\)
\(70\%\)
0.125
\(-20.0 \times 10^{-6}\)
\(46.8\%\)
\(47.2\%\)
0.0625
\(-30.7 \times 10^{-6}\)
\(18.2\%\)
\(18.8\%\)
0.03125
\(-35.9 \times 10^{-6}\)
\(4.14\%\)
\(4.9\%\)
0.015625
\(-37.5 \times 10^{-6}\)
\(0.12\%\)
\(0.909\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-37.5 \times 10^{-6}\)
\(0.116\%\)
\(0.676\%\)
0.125
\(-38.0 \times 10^{-6}\)
\(1.34\%\)
\(0.538\%\)
0.0625
\(-38.0 \times 10^{-6}\)
\(1.29\%\)
\(0.485\%\)
0.03125
\(-38.0 \times 10^{-6}\)
\(1.24\%\)
\(0.437\%\)
0.015625
\(-38.0 \times 10^{-6}\)
\(1.21\%\)
\(0.412\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-39.5 \times 10^{-6}\)
\(5.32\%\)
\(4.49\%\)
0.125
\(-38.1 \times 10^{-6}\)
\(1.53\%\)
\(0.727\%\)
0.0625
\(-38.0 \times 10^{-6}\)
\(1.33\%\)
\(0.531\%\)
0.03125
\(-38.0 \times 10^{-6}\)
\(1.26\%\)
\(0.454\%\)
0.015625
\(-38.0 \times 10^{-6}\)
\(1.22\%\)
\(0.418\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 18: Uniform pressure: deflection along the axis — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 19: Uniform pressure: bending stress along the top fiber — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 20: Uniform pressure: bending stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 21: Uniform pressure: shear stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a) Axis deflection
\(u_y\)
(b) Top-fiber stress
\(\sigma_{zz}\)
(c) Midspan stress
\(\sigma_{zz}\)
(d) Midspan shear
\(\tau_{zy}\)
Figure 22: Uniform pressure: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h4 profile; the h4 point is omitted because its self-error is zero.
7.3.1
Contour plots
Coreform IGA for Abaqus
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
h0p3
h1p3
h2p3
h3p3
h4p3
Abaqus reference
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
7.4
Pressure (total force)
Euler-Bernoulli tip deflection:
\(-37.5 \times 10^{-6}\)
m; finest C3D20 Abaqus reference:
\(-37.8 \times 10^{-6}\)
m.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-11.4 \times 10^{-6}\)
\(69.7\%\)
\(70\%\)
0.125
\(-20.0 \times 10^{-6}\)
\(46.8\%\)
\(47.2\%\)
0.0625
\(-30.7 \times 10^{-6}\)
\(18.2\%\)
\(18.8\%\)
0.03125
\(-35.9 \times 10^{-6}\)
\(4.14\%\)
\(4.9\%\)
0.015625
\(-37.5 \times 10^{-6}\)
\(0.12\%\)
\(0.909\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-37.5 \times 10^{-6}\)
\(0.116\%\)
\(0.676\%\)
0.125
\(-38.0 \times 10^{-6}\)
\(1.34\%\)
\(0.538\%\)
0.0625
\(-38.0 \times 10^{-6}\)
\(1.29\%\)
\(0.485\%\)
0.03125
\(-38.0 \times 10^{-6}\)
\(1.24\%\)
\(0.437\%\)
0.015625
\(-38.0 \times 10^{-6}\)
\(1.21\%\)
\(0.412\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-39.5 \times 10^{-6}\)
\(5.32\%\)
\(4.49\%\)
0.125
\(-38.1 \times 10^{-6}\)
\(1.53\%\)
\(0.727\%\)
0.0625
\(-38.0 \times 10^{-6}\)
\(1.33\%\)
\(0.531\%\)
0.03125
\(-38.0 \times 10^{-6}\)
\(1.26\%\)
\(0.454\%\)
0.015625
\(-38.0 \times 10^{-6}\)
\(1.22\%\)
\(0.418\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 23: Total-force pressure: deflection along the axis — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 24: Total-force pressure: bending stress along the top fiber — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 25: Total-force pressure: bending stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 26: Total-force pressure: shear stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a) Axis deflection
\(u_y\)
(b) Top-fiber stress
\(\sigma_{zz}\)
(c) Midspan stress
\(\sigma_{zz}\)
(d) Midspan shear
\(\tau_{zy}\)
Figure 27: Total-force pressure: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h4 profile; the h4 point is omitted because its self-error is zero.
7.4.1
Contour plots
Coreform IGA for Abaqus
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
h0p3
h1p3
h2p3
h3p3
h4p3
Abaqus reference
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
7.5
Surface traction
Euler-Bernoulli tip deflection:
\(-1.00 \times 10^{-6}\)
m; finest C3D20 Abaqus reference:
\(-1.01 \times 10^{-6}\)
m.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-287 \times 10^{-9}\)
\(71.3\%\)
\(71.5\%\)
0.125
\(-509 \times 10^{-9}\)
\(49.1\%\)
\(49.4\%\)
0.0625
\(-810 \times 10^{-9}\)
\(19\%\)
\(19.4\%\)
0.03125
\(-955 \times 10^{-9}\)
\(4.52\%\)
\(5.08\%\)
0.015625
\(-996 \times 10^{-9}\)
\(0.431\%\)
\(1.02\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-1.00 \times 10^{-6}\)
\(0.0964\%\)
\(0.496\%\)
0.125
\(-1.01 \times 10^{-6}\)
\(1.01\%\)
\(0.409\%\)
0.0625
\(-1.01 \times 10^{-6}\)
\(0.964\%\)
\(0.366\%\)
0.03125
\(-1.01 \times 10^{-6}\)
\(0.926\%\)
\(0.328\%\)
0.015625
\(-1.01 \times 10^{-6}\)
\(0.907\%\)
\(0.309\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-1.04 \times 10^{-6}\)
\(4.05\%\)
\(3.44\%\)
0.125
\(-1.01 \times 10^{-6}\)
\(1.15\%\)
\(0.551\%\)
0.0625
\(-1.01 \times 10^{-6}\)
\(0.998\%\)
\(0.4\%\)
0.03125
\(-1.01 \times 10^{-6}\)
\(0.939\%\)
\(0.341\%\)
0.015625
\(-1.01 \times 10^{-6}\)
\(0.911\%\)
\(0.314\%\)
Download data (CSV)
The follower and resultant variants run as single-level coverage cases; in this small-displacement analysis they must reproduce the base traction solution (
Table 5
).
Table 5: Traction variants at the coarsest level.
Test case
Variant
\(u_y\)
tip (
\(k=0\)
,
\(p=2\)
)
rel. err
traction
base
\(-1.00 \times 10^{-6}\)
\(0.0964\%\)
traction_follower
follower
\(-1.00 \times 10^{-6}\)
\(0.0964\%\)
traction_resultant
resultant
\(-1.00 \times 10^{-6}\)
\(0.0964\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 28: Surface traction: deflection along the axis — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 29: Surface traction: bending stress along the top fiber — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 30: Surface traction: bending stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 31: Surface traction: shear stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a) Axis deflection
\(u_y\)
(b) Top-fiber stress
\(\sigma_{zz}\)
(c) Midspan stress
\(\sigma_{zz}\)
(d) Midspan shear
\(\tau_{zy}\)
Figure 32: Surface traction: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h4 profile; the h4 point is omitted because its self-error is zero.
7.5.1
Contour plots
Coreform IGA for Abaqus
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
h0p3
h1p3
h2p3
h3p3
h4p3
Abaqus reference
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
7.6
End moment
Euler-Bernoulli tip deflection:
\(-1.50 \times 10^{-6}\)
m; finest C3D20 Abaqus reference:
\(-1.50 \times 10^{-6}\)
m.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-411 \times 10^{-9}\)
\(72.6\%\)
\(72.6\%\)
0.125
\(-729 \times 10^{-9}\)
\(51.4\%\)
\(51.4\%\)
0.0625
\(-1.20 \times 10^{-6}\)
\(20\%\)
\(20\%\)
0.03125
\(-1.42 \times 10^{-6}\)
\(5.23\%\)
\(5.23\%\)
0.015625
\(-1.48 \times 10^{-6}\)
\(1.12\%\)
\(1.12\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-1.50 \times 10^{-6}\)
\(0.246\%\)
\(0.246\%\)
0.125
\(-1.50 \times 10^{-6}\)
\(0.329\%\)
\(0.329\%\)
0.0625
\(-1.50 \times 10^{-6}\)
\(0.25\%\)
\(0.25\%\)
0.03125
\(-1.50 \times 10^{-6}\)
\(0.22\%\)
\(0.22\%\)
0.015625
\(-1.50 \times 10^{-6}\)
\(0.207\%\)
\(0.207\%\)
\(h_z\)
\(u_y\)
tip [m]
vs Euler-Bernoulli
vs finest C3D20 reference
0.25
\(-1.54 \times 10^{-6}\)
\(2.41\%\)
\(2.41\%\)
0.125
\(-1.51 \times 10^{-6}\)
\(0.374\%\)
\(0.374\%\)
0.0625
\(-1.50 \times 10^{-6}\)
\(0.269\%\)
\(0.269\%\)
0.03125
\(-1.50 \times 10^{-6}\)
\(0.228\%\)
\(0.228\%\)
0.015625
\(-1.50 \times 10^{-6}\)
\(0.21\%\)
\(0.21\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 33: End moment: deflection along the axis — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 34: End moment: bending stress along the top fiber — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 35: End moment: bending stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 36: End moment: shear stress through the midspan section — coarse→fine sweep per degree vs Euler-Bernoulli.
(a) Axis deflection
\(u_y\)
(b) Top-fiber stress
\(\sigma_{zz}\)
(c) Midspan stress
\(\sigma_{zz}\)
(d) Midspan shear
\(\tau_{zy}\)
Figure 37: End moment: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h4 profile; the h4 point is omitted because its self-error is zero.
7.6.1
Contour plots
Coreform IGA for Abaqus
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
h0p3
h1p3
h2p3
h3p3
h4p3
Abaqus reference
h0p1
h1p1
h2p1
h3p1
h4p1
h0p2
h1p2
h2p2
h3p2
h4p2
8
Performance
Performance is reported separately for the two stages of each immersed IGA analysis:
trimming
(building the trimmed U-spline basis on the immersing grid) and the
Abaqus solve
. Each load case uses the same mesh sweep, but timings are reported per load case because each study was executed independently.
8.1
Gravity
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(63\)
\(5\)
\(2\)
\(1\)
\(1.2\)
\(0.125\)
\(275\)
\(81\)
\(74\)
\(2\)
\(2.7\)
\(0.0625\)
\(931\)
\(425\)
\(290\)
\(4\)
\(7.0\)
\(0.03125\)
\(4{,}235\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(26.6\)
\(0.015625\)
\(24{,}187\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(117.7\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(225\)
\(5\)
\(2\)
\(1\)
\(5.6\)
\(0.125\)
\(637\)
\(81\)
\(74\)
\(2\)
\(4.0\)
\(0.0625\)
\(1{,}701\)
\(425\)
\(290\)
\(4\)
\(9.6\)
\(0.03125\)
\(6{,}253\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(35.5\)
\(0.015625\)
\(30{,}429\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(128.6\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(539\)
\(5\)
\(2\)
\(1\)
\(32.8\)
\(0.125\)
\(1{,}215\)
\(81\)
\(74\)
\(2\)
\(22.5\)
\(0.0625\)
\(2{,}783\)
\(425\)
\(290\)
\(4\)
\(28.9\)
\(0.03125\)
\(8{,}775\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(55.4\)
\(0.015625\)
\(37{,}559\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(175.6\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(69\)
\(76\)
\(228\)
\(1\)
\(1\)
\(0.125\)
\(327\)
\(321\)
\(963\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}191\)
\(3{,}573\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(5{,}379\)
\(16{,}137\)
\(8\)
\(2\)
\(0.015625\)
\(33{,}283\)
\(28{,}923\)
\(86{,}769\)
\(16\)
\(4\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(51\)
\(107\)
\(321\)
\(1\)
\(1\)
\(0.125\)
\(327\)
\(436\)
\(1{,}308\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}474\)
\(4{,}422\)
\(4\)
\(2\)
\(0.03125\)
\(6{,}339\)
\(6{,}214\)
\(18{,}642\)
\(8\)
\(2\)
\(0.015625\)
\(33{,}283\)
\(31{,}726\)
\(95{,}178\)
\(16\)
\(9\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(67\)
\(184\)
\(552\)
\(1\)
\(2\)
\(0.125\)
\(457\)
\(667\)
\(2{,}001\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}745\)
\(2{,}029\)
\(6{,}087\)
\(4\)
\(2\)
\(0.03125\)
\(7{,}609\)
\(7{,}825\)
\(23{,}475\)
\(8\)
\(7\)
\(0.015625\)
\(37{,}865\)
\(37{,}081\)
\(111{,}243\)
\(16\)
\(36\)
Download data (CSV)
8.2
Body force
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(63\)
\(5\)
\(2\)
\(1\)
\(3.4\)
\(0.125\)
\(275\)
\(81\)
\(74\)
\(2\)
\(3.1\)
\(0.0625\)
\(931\)
\(425\)
\(290\)
\(4\)
\(9.0\)
\(0.03125\)
\(4{,}235\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(27.1\)
\(0.015625\)
\(24{,}187\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(140.9\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(225\)
\(5\)
\(2\)
\(1\)
\(5.4\)
\(0.125\)
\(637\)
\(81\)
\(74\)
\(2\)
\(4.2\)
\(0.0625\)
\(1{,}701\)
\(425\)
\(290\)
\(4\)
\(10.3\)
\(0.03125\)
\(6{,}253\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(34.1\)
\(0.015625\)
\(30{,}429\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(126.1\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(539\)
\(5\)
\(2\)
\(1\)
\(22.0\)
\(0.125\)
\(1{,}215\)
\(81\)
\(74\)
\(2\)
\(22.4\)
\(0.0625\)
\(2{,}783\)
\(425\)
\(290\)
\(4\)
\(28.6\)
\(0.03125\)
\(8{,}775\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(55.4\)
\(0.015625\)
\(37{,}559\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(169.5\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(69\)
\(76\)
\(228\)
\(1\)
\(1\)
\(0.125\)
\(327\)
\(321\)
\(963\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}191\)
\(3{,}573\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(5{,}379\)
\(16{,}137\)
\(8\)
\(1\)
\(0.015625\)
\(33{,}283\)
\(28{,}923\)
\(86{,}769\)
\(16\)
\(6\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(51\)
\(107\)
\(321\)
\(1\)
\(3\)
\(0.125\)
\(327\)
\(436\)
\(1{,}308\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}474\)
\(4{,}422\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(6{,}214\)
\(18{,}642\)
\(8\)
\(3\)
\(0.015625\)
\(33{,}283\)
\(31{,}726\)
\(95{,}178\)
\(16\)
\(10\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(67\)
\(184\)
\(552\)
\(1\)
\(1\)
\(0.125\)
\(457\)
\(667\)
\(2{,}001\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}745\)
\(2{,}029\)
\(6{,}087\)
\(4\)
\(3\)
\(0.03125\)
\(7{,}609\)
\(7{,}825\)
\(23{,}475\)
\(8\)
\(7\)
\(0.015625\)
\(37{,}865\)
\(37{,}081\)
\(111{,}243\)
\(16\)
\(36\)
Download data (CSV)
8.3
Pressure (uniform)
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(63\)
\(5\)
\(2\)
\(1\)
\(1.2\)
\(0.125\)
\(275\)
\(81\)
\(74\)
\(2\)
\(3.0\)
\(0.0625\)
\(931\)
\(425\)
\(290\)
\(4\)
\(8.2\)
\(0.03125\)
\(4{,}235\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(28.0\)
\(0.015625\)
\(24{,}187\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(120.8\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(225\)
\(5\)
\(2\)
\(1\)
\(5.3\)
\(0.125\)
\(637\)
\(81\)
\(74\)
\(2\)
\(4.3\)
\(0.0625\)
\(1{,}701\)
\(425\)
\(290\)
\(4\)
\(9.6\)
\(0.03125\)
\(6{,}253\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(35.2\)
\(0.015625\)
\(30{,}429\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(177.8\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(539\)
\(5\)
\(2\)
\(1\)
\(26.3\)
\(0.125\)
\(1{,}215\)
\(81\)
\(74\)
\(2\)
\(22.4\)
\(0.0625\)
\(2{,}783\)
\(425\)
\(290\)
\(4\)
\(29.0\)
\(0.03125\)
\(8{,}775\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(56.9\)
\(0.015625\)
\(37{,}559\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(179.8\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(69\)
\(76\)
\(228\)
\(1\)
\(2\)
\(0.125\)
\(327\)
\(321\)
\(963\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}363\)
\(1{,}191\)
\(3{,}573\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(5{,}379\)
\(16{,}137\)
\(8\)
\(2\)
\(0.015625\)
\(33{,}283\)
\(28{,}923\)
\(86{,}769\)
\(16\)
\(4\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(51\)
\(107\)
\(321\)
\(1\)
\(1\)
\(0.125\)
\(327\)
\(436\)
\(1{,}308\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}474\)
\(4{,}422\)
\(4\)
\(2\)
\(0.03125\)
\(6{,}339\)
\(6{,}214\)
\(18{,}642\)
\(8\)
\(2\)
\(0.015625\)
\(33{,}283\)
\(31{,}726\)
\(95{,}178\)
\(16\)
\(11\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(67\)
\(184\)
\(552\)
\(1\)
\(1\)
\(0.125\)
\(457\)
\(667\)
\(2{,}001\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}745\)
\(2{,}029\)
\(6{,}087\)
\(4\)
\(3\)
\(0.03125\)
\(7{,}609\)
\(7{,}825\)
\(23{,}475\)
\(8\)
\(8\)
\(0.015625\)
\(37{,}865\)
\(37{,}081\)
\(111{,}243\)
\(16\)
\(34\)
Download data (CSV)
8.4
Pressure (total force)
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(63\)
\(5\)
\(2\)
\(1\)
\(3.0\)
\(0.125\)
\(275\)
\(81\)
\(74\)
\(2\)
\(2.6\)
\(0.0625\)
\(931\)
\(425\)
\(290\)
\(4\)
\(8.6\)
\(0.03125\)
\(4{,}235\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(32.6\)
\(0.015625\)
\(24{,}187\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(124.3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(225\)
\(5\)
\(2\)
\(1\)
\(4.0\)
\(0.125\)
\(637\)
\(81\)
\(74\)
\(2\)
\(7.4\)
\(0.0625\)
\(1{,}701\)
\(425\)
\(290\)
\(4\)
\(9.4\)
\(0.03125\)
\(6{,}253\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(28.5\)
\(0.015625\)
\(30{,}429\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(218.3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(539\)
\(5\)
\(2\)
\(1\)
\(21.5\)
\(0.125\)
\(1{,}215\)
\(81\)
\(74\)
\(2\)
\(22.5\)
\(0.0625\)
\(2{,}783\)
\(425\)
\(290\)
\(4\)
\(28.9\)
\(0.03125\)
\(8{,}775\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(57.7\)
\(0.015625\)
\(37{,}559\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(151.2\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(69\)
\(76\)
\(228\)
\(1\)
\(3\)
\(0.125\)
\(327\)
\(321\)
\(963\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}363\)
\(1{,}191\)
\(3{,}573\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(5{,}379\)
\(16{,}137\)
\(8\)
\(1\)
\(0.015625\)
\(33{,}283\)
\(28{,}923\)
\(86{,}769\)
\(16\)
\(3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(51\)
\(107\)
\(321\)
\(1\)
\(1\)
\(0.125\)
\(327\)
\(436\)
\(1{,}308\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}363\)
\(1{,}474\)
\(4{,}422\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(6{,}214\)
\(18{,}642\)
\(8\)
\(2\)
\(0.015625\)
\(33{,}283\)
\(31{,}726\)
\(95{,}178\)
\(16\)
\(12\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(67\)
\(184\)
\(552\)
\(1\)
\(1\)
\(0.125\)
\(457\)
\(667\)
\(2{,}001\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}745\)
\(2{,}029\)
\(6{,}087\)
\(4\)
\(2\)
\(0.03125\)
\(7{,}609\)
\(7{,}825\)
\(23{,}475\)
\(8\)
\(8\)
\(0.015625\)
\(37{,}865\)
\(37{,}081\)
\(111{,}243\)
\(16\)
\(33\)
Download data (CSV)
8.5
Surface traction
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(63\)
\(5\)
\(2\)
\(1\)
\(3.3\)
\(0.125\)
\(275\)
\(81\)
\(74\)
\(2\)
\(3.0\)
\(0.0625\)
\(931\)
\(425\)
\(290\)
\(4\)
\(8.8\)
\(0.03125\)
\(4{,}235\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(27.1\)
\(0.015625\)
\(24{,}187\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(119.5\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(225\)
\(5\)
\(2\)
\(1\)
\(2.4\)
\(0.125\)
\(637\)
\(81\)
\(74\)
\(2\)
\(4.1\)
\(0.0625\)
\(1{,}701\)
\(425\)
\(290\)
\(4\)
\(8.5\)
\(0.03125\)
\(6{,}253\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(34.0\)
\(0.015625\)
\(30{,}429\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(226.5\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(539\)
\(5\)
\(2\)
\(1\)
\(20.9\)
\(0.125\)
\(1{,}215\)
\(81\)
\(74\)
\(2\)
\(22.2\)
\(0.0625\)
\(2{,}783\)
\(425\)
\(290\)
\(4\)
\(28.6\)
\(0.03125\)
\(8{,}775\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(57.8\)
\(0.015625\)
\(37{,}559\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(181.4\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(69\)
\(76\)
\(228\)
\(1\)
\(2\)
\(0.125\)
\(327\)
\(321\)
\(963\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}191\)
\(3{,}573\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(5{,}379\)
\(16{,}137\)
\(8\)
\(1\)
\(0.015625\)
\(33{,}283\)
\(28{,}923\)
\(86{,}769\)
\(16\)
\(4\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(51\)
\(107\)
\(321\)
\(1\)
\(0\)
\(0.125\)
\(327\)
\(436\)
\(1{,}308\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}363\)
\(1{,}474\)
\(4{,}422\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(6{,}214\)
\(18{,}642\)
\(8\)
\(2\)
\(0.015625\)
\(33{,}283\)
\(31{,}726\)
\(95{,}178\)
\(16\)
\(13\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(67\)
\(184\)
\(552\)
\(1\)
\(2\)
\(0.125\)
\(457\)
\(667\)
\(2{,}001\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}745\)
\(2{,}029\)
\(6{,}087\)
\(4\)
\(2\)
\(0.03125\)
\(7{,}609\)
\(7{,}825\)
\(23{,}475\)
\(8\)
\(8\)
\(0.015625\)
\(37{,}865\)
\(37{,}081\)
\(111{,}243\)
\(16\)
\(33\)
Download data (CSV)
8.6
End moment
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(63\)
\(5\)
\(2\)
\(1\)
\(3.5\)
\(0.125\)
\(275\)
\(81\)
\(74\)
\(2\)
\(2.9\)
\(0.0625\)
\(931\)
\(425\)
\(290\)
\(4\)
\(8.6\)
\(0.03125\)
\(4{,}235\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(33.8\)
\(0.015625\)
\(24{,}187\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(117.9\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(225\)
\(5\)
\(2\)
\(1\)
\(2.6\)
\(0.125\)
\(637\)
\(81\)
\(74\)
\(2\)
\(4.0\)
\(0.0625\)
\(1{,}701\)
\(425\)
\(290\)
\(4\)
\(8.7\)
\(0.03125\)
\(6{,}253\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(29.1\)
\(0.015625\)
\(30{,}429\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(214.6\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.25\)
\(539\)
\(5\)
\(2\)
\(1\)
\(22.9\)
\(0.125\)
\(1{,}215\)
\(81\)
\(74\)
\(2\)
\(22.6\)
\(0.0625\)
\(2{,}783\)
\(425\)
\(290\)
\(4\)
\(27.0\)
\(0.03125\)
\(8{,}775\)
\(2{,}673\)
\(1{,}154\)
\(8\)
\(57.4\)
\(0.015625\)
\(37{,}559\)
\(18{,}785\)
\(4{,}610\)
\(16\)
\(180.8\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(69\)
\(77\)
\(261\)
\(1\)
\(3\)
\(0.125\)
\(327\)
\(322\)
\(1{,}044\)
\(2\)
\(1\)
\(0.0625\)
\(1{,}363\)
\(1{,}192\)
\(3{,}822\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(5{,}380\)
\(17{,}010\)
\(8\)
\(1\)
\(0.015625\)
\(33{,}283\)
\(28{,}924\)
\(90{,}042\)
\(16\)
\(5\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(51\)
\(108\)
\(354\)
\(1\)
\(1\)
\(0.125\)
\(327\)
\(437\)
\(1{,}389\)
\(2\)
\(2\)
\(0.0625\)
\(1{,}363\)
\(1{,}475\)
\(4{,}671\)
\(4\)
\(1\)
\(0.03125\)
\(6{,}339\)
\(6{,}215\)
\(19{,}515\)
\(8\)
\(3\)
\(0.015625\)
\(33{,}283\)
\(31{,}727\)
\(98{,}451\)
\(16\)
\(12\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.25\)
\(67\)
\(185\)
\(597\)
\(1\)
\(2\)
\(0.125\)
\(457\)
\(668\)
\(2{,}154\)
\(2\)
\(3\)
\(0.0625\)
\(1{,}745\)
\(2{,}030\)
\(6{,}516\)
\(4\)
\(2\)
\(0.03125\)
\(7{,}609\)
\(7{,}826\)
\(24{,}888\)
\(8\)
\(8\)
\(0.015625\)
\(37{,}865\)
\(37{,}082\)
\(116{,}352\)
\(16\)
\(34\)
Download data (CSV)
9
Discussion
Each primary load-case sweep converges under mesh refinement. On the finer meshes, the tip deflection satisfies the verification tolerance and all four sampled profiles approach the closed-form solutions.
The remaining difference from Euler–Bernoulli theory does not vanish with refinement because it arises from differences between the models, not from the discretization. Euler–Bernoulli theory neglects transverse shear deformation, whereas the solid model includes it; the pinned solid end face is also stiffer than the idealized beam support. For this reason, line-probe convergence is measured against the finest same-degree Coreform IGA profile,
\(h_4\)
, rather than treating either the beam theory or the Abaqus reference solution as exact. The finest C3D20 Abaqus result is retained as an independent finite-element comparison for the tip quantities, while the body-fitted sweeps provide qualitative context in the contour galleries.
The load-equivalence checks compare cases intended to represent the same physical loading. Gravity and body force both produce a line load of
\(\gravityLineLoad\,\mathrm{N\,m^{-1}}\)
. Uniform pressure and total-force pressure both produce a line load of
\(\pressureLineLoad\,\mathrm{N\,m^{-1}}\)
. Agreement within each pair verifies the corresponding Coreform IGA load translation. The pairs need not agree with one another because their line loads differ.
For this small-displacement analysis, the follower and resultant traction variants also reproduce the base traction result.
References
[1]
S. P. Timoshenko,
Strength of materials, part i: Elementary theory and problems
, 2nd ed. New York: D. Van Nostrand Company, 1940.
Appendix
9.1
Capabilities exercised
The table is generated from current passing simulation and study artifacts. Each row must satisfy the linked capability card’s evidence contract; declaring a capability in a test does not, by itself, create evidence.
Table 6: Capabilities exercised by the current verification artifacts.
Capability
Evidence level
Linear static analysis
Numerically validated
Static stress analysis
Numerically validated
Boundary conditions and constraints
Numerically validated
Pinned boundary condition
Numerically validated
Kinematic coupling constraints
Numerically validated
Threaded Abaqus/Standard execution
Numerically validated
Parallel Coreform IGA mesh generation
Numerically validated
Parallel execution
Numerically validated
Body-force loading
Numerically validated
Gravity loading
Numerically validated
Moment loading
Numerically validated
Pressure loading
Numerically validated
Pressure specified by total force
Numerically validated
Uniform pressure loading
Numerically validated
Surface traction loading
Numerically validated
Follower surface traction
Exercised
Resultant surface traction
Exercised
Mechanical loads
Numerically validated
Material density
Numerically validated
Linear elastic material behavior
Numerically validated
Materials and section assignments
Numerically validated
Meshing
Numerically validated
Linear C0 spline basis
Exercised
Quadratic C1 spline basis
Numerically validated
Cubic C2 spline basis
Numerically validated
Immersed rectilinear meshes
Numerically validated
Spline basis degree and continuity
Numerically validated
Explicit field-output selection
Numerically validated
Field output requests
Numerically validated
Stress field output (S)
Exercised
Displacement field output (U)
Numerically validated
Individual field-output variables
Numerically validated
Line probes
Numerically validated
Coreform IGA result probes
Numerically validated
Dependent part instances
Numerically validated
9.2
Download artifacts
Download all verification artifacts (ZIP)
9.2.1
Geometry CAD files
No separate CAD file is required for this problem; the geometry is fully defined in the Abaqus/CAE journal files.
9.2.2
Abaqus/CAE journal files
152 generated Abaqus/CAE journals are available.
Show 152 individual journals
abaqus_h0p1_body_force–beam_bending_h0_model.jnl
abaqus_h0p1_gravity–beam_bending_h0_model.jnl
abaqus_h0p1_moment–beam_bending_h0_model.jnl
abaqus_h0p1_pressure_total_force–beam_bending_h0_model.jnl
abaqus_h0p1_pressure_uniform–beam_bending_h0_model.jnl
abaqus_h0p1_traction–beam_bending_h0_model.jnl
abaqus_h0p2_body_force–beam_bending_h0_model.jnl
abaqus_h0p2_gravity–beam_bending_h0_model.jnl
abaqus_h0p2_moment–beam_bending_h0_model.jnl
abaqus_h0p2_pressure_total_force–beam_bending_h0_model.jnl
abaqus_h0p2_pressure_uniform–beam_bending_h0_model.jnl
abaqus_h0p2_traction–beam_bending_h0_model.jnl
abaqus_h1p1_body_force–beam_bending_h1_model.jnl
abaqus_h1p1_gravity–beam_bending_h1_model.jnl
abaqus_h1p1_moment–beam_bending_h1_model.jnl
abaqus_h1p1_pressure_total_force–beam_bending_h1_model.jnl
abaqus_h1p1_pressure_uniform–beam_bending_h1_model.jnl
abaqus_h1p1_traction–beam_bending_h1_model.jnl
abaqus_h1p2_body_force–beam_bending_h1_model.jnl
abaqus_h1p2_gravity–beam_bending_h1_model.jnl
abaqus_h1p2_moment–beam_bending_h1_model.jnl
abaqus_h1p2_pressure_total_force–beam_bending_h1_model.jnl
abaqus_h1p2_pressure_uniform–beam_bending_h1_model.jnl
abaqus_h1p2_traction–beam_bending_h1_model.jnl
abaqus_h2p1_body_force–beam_bending_h2_model.jnl
abaqus_h2p1_gravity–beam_bending_h2_model.jnl
abaqus_h2p1_moment–beam_bending_h2_model.jnl
abaqus_h2p1_pressure_total_force–beam_bending_h2_model.jnl
abaqus_h2p1_pressure_uniform–beam_bending_h2_model.jnl
abaqus_h2p1_traction–beam_bending_h2_model.jnl
abaqus_h2p2_body_force–beam_bending_h2_model.jnl
abaqus_h2p2_gravity–beam_bending_h2_model.jnl
abaqus_h2p2_moment–beam_bending_h2_model.jnl
abaqus_h2p2_pressure_total_force–beam_bending_h2_model.jnl
abaqus_h2p2_pressure_uniform–beam_bending_h2_model.jnl
abaqus_h2p2_traction–beam_bending_h2_model.jnl
abaqus_h3p1_body_force–beam_bending_h3_model.jnl
abaqus_h3p1_gravity–beam_bending_h3_model.jnl
abaqus_h3p1_moment–beam_bending_h3_model.jnl
abaqus_h3p1_pressure_total_force–beam_bending_h3_model.jnl
abaqus_h3p1_pressure_uniform–beam_bending_h3_model.jnl
abaqus_h3p1_traction–beam_bending_h3_model.jnl
abaqus_h3p2_body_force–beam_bending_h3_model.jnl
abaqus_h3p2_gravity–beam_bending_h3_model.jnl
abaqus_h3p2_moment–beam_bending_h3_model.jnl
abaqus_h3p2_pressure_total_force–beam_bending_h3_model.jnl
abaqus_h3p2_pressure_uniform–beam_bending_h3_model.jnl
abaqus_h3p2_traction–beam_bending_h3_model.jnl
abaqus_h4p1_body_force–beam_bending_h4_model.jnl
abaqus_h4p1_gravity–beam_bending_h4_model.jnl
abaqus_h4p1_moment–beam_bending_h4_model.jnl
abaqus_h4p1_pressure_total_force–beam_bending_h4_model.jnl
abaqus_h4p1_pressure_uniform–beam_bending_h4_model.jnl
abaqus_h4p1_traction–beam_bending_h4_model.jnl
abaqus_h4p2_body_force–beam_bending_h4_model.jnl
abaqus_h4p2_gravity–beam_bending_h4_model.jnl
abaqus_h4p2_moment–beam_bending_h4_model.jnl
abaqus_h4p2_pressure_total_force–beam_bending_h4_model.jnl
abaqus_h4p2_pressure_uniform–beam_bending_h4_model.jnl
abaqus_h4p2_traction–beam_bending_h4_model.jnl
immersed_h0p1_body_force–beam_bending_h0_model.jnl
immersed_h0p1_gravity–beam_bending_h0_model.jnl
immersed_h0p1_moment–beam_bending_h0_model.jnl
immersed_h0p1_pressure_total_force–beam_bending_h0_model.jnl
immersed_h0p1_pressure_uniform–beam_bending_h0_model.jnl
immersed_h0p1_traction–beam_bending_h0_model.jnl
immersed_h0p2_body_force–beam_bending_h0_model.jnl
immersed_h0p2_gravity–beam_bending_h0_model.jnl
immersed_h0p2_moment–beam_bending_h0_model.jnl
immersed_h0p2_pressure_total_force–beam_bending_h0_model.jnl
immersed_h0p2_pressure_uniform–beam_bending_h0_model.jnl
immersed_h0p2_traction–beam_bending_h0_model.jnl
immersed_h0p2_traction_follower–beam_bending_h0_model.jnl
immersed_h0p2_traction_resultant–beam_bending_h0_model.jnl
immersed_h0p3_body_force–beam_bending_h0_model.jnl
immersed_h0p3_gravity–beam_bending_h0_model.jnl
immersed_h0p3_moment–beam_bending_h0_model.jnl
immersed_h0p3_pressure_total_force–beam_bending_h0_model.jnl
immersed_h0p3_pressure_uniform–beam_bending_h0_model.jnl
immersed_h0p3_traction–beam_bending_h0_model.jnl
immersed_h1p1_body_force–beam_bending_h1_model.jnl
immersed_h1p1_gravity–beam_bending_h1_model.jnl
immersed_h1p1_moment–beam_bending_h1_model.jnl
immersed_h1p1_pressure_total_force–beam_bending_h1_model.jnl
immersed_h1p1_pressure_uniform–beam_bending_h1_model.jnl
immersed_h1p1_traction–beam_bending_h1_model.jnl
immersed_h1p2_body_force–beam_bending_h1_model.jnl
immersed_h1p2_gravity–beam_bending_h1_model.jnl
immersed_h1p2_moment–beam_bending_h1_model.jnl
immersed_h1p2_pressure_total_force–beam_bending_h1_model.jnl
immersed_h1p2_pressure_uniform–beam_bending_h1_model.jnl
immersed_h1p2_traction–beam_bending_h1_model.jnl
immersed_h1p3_body_force–beam_bending_h1_model.jnl
immersed_h1p3_gravity–beam_bending_h1_model.jnl
immersed_h1p3_moment–beam_bending_h1_model.jnl
immersed_h1p3_pressure_total_force–beam_bending_h1_model.jnl
immersed_h1p3_pressure_uniform–beam_bending_h1_model.jnl
immersed_h1p3_traction–beam_bending_h1_model.jnl
immersed_h2p1_body_force–beam_bending_h2_model.jnl
immersed_h2p1_gravity–beam_bending_h2_model.jnl
immersed_h2p1_moment–beam_bending_h2_model.jnl
immersed_h2p1_pressure_total_force–beam_bending_h2_model.jnl
immersed_h2p1_pressure_uniform–beam_bending_h2_model.jnl
immersed_h2p1_traction–beam_bending_h2_model.jnl
immersed_h2p2_body_force–beam_bending_h2_model.jnl
immersed_h2p2_gravity–beam_bending_h2_model.jnl
immersed_h2p2_moment–beam_bending_h2_model.jnl
immersed_h2p2_pressure_total_force–beam_bending_h2_model.jnl
immersed_h2p2_pressure_uniform–beam_bending_h2_model.jnl
immersed_h2p2_traction–beam_bending_h2_model.jnl
immersed_h2p3_body_force–beam_bending_h2_model.jnl
immersed_h2p3_gravity–beam_bending_h2_model.jnl
immersed_h2p3_moment–beam_bending_h2_model.jnl
immersed_h2p3_pressure_total_force–beam_bending_h2_model.jnl
immersed_h2p3_pressure_uniform–beam_bending_h2_model.jnl
immersed_h2p3_traction–beam_bending_h2_model.jnl
immersed_h3p1_body_force–beam_bending_h3_model.jnl
immersed_h3p1_gravity–beam_bending_h3_model.jnl
immersed_h3p1_moment–beam_bending_h3_model.jnl
immersed_h3p1_pressure_total_force–beam_bending_h3_model.jnl
immersed_h3p1_pressure_uniform–beam_bending_h3_model.jnl
immersed_h3p1_traction–beam_bending_h3_model.jnl
immersed_h3p2_body_force–beam_bending_h3_model.jnl
immersed_h3p2_gravity–beam_bending_h3_model.jnl
immersed_h3p2_moment–beam_bending_h3_model.jnl
immersed_h3p2_pressure_total_force–beam_bending_h3_model.jnl
immersed_h3p2_pressure_uniform–beam_bending_h3_model.jnl
immersed_h3p2_traction–beam_bending_h3_model.jnl
immersed_h3p3_body_force–beam_bending_h3_model.jnl
immersed_h3p3_gravity–beam_bending_h3_model.jnl
immersed_h3p3_moment–beam_bending_h3_model.jnl
immersed_h3p3_pressure_total_force–beam_bending_h3_model.jnl
immersed_h3p3_pressure_uniform–beam_bending_h3_model.jnl
immersed_h3p3_traction–beam_bending_h3_model.jnl
immersed_h4p1_body_force–beam_bending_h4_model.jnl
immersed_h4p1_gravity–beam_bending_h4_model.jnl
immersed_h4p1_moment–beam_bending_h4_model.jnl
immersed_h4p1_pressure_total_force–beam_bending_h4_model.jnl
immersed_h4p1_pressure_uniform–beam_bending_h4_model.jnl
immersed_h4p1_traction–beam_bending_h4_model.jnl
immersed_h4p2_body_force–beam_bending_h4_model.jnl
immersed_h4p2_gravity–beam_bending_h4_model.jnl
immersed_h4p2_moment–beam_bending_h4_model.jnl
immersed_h4p2_pressure_total_force–beam_bending_h4_model.jnl
immersed_h4p2_pressure_uniform–beam_bending_h4_model.jnl
immersed_h4p2_traction–beam_bending_h4_model.jnl
immersed_h4p3_body_force–beam_bending_h4_model.jnl
immersed_h4p3_gravity–beam_bending_h4_model.jnl
immersed_h4p3_moment–beam_bending_h4_model.jnl
immersed_h4p3_pressure_total_force–beam_bending_h4_model.jnl
immersed_h4p3_pressure_uniform–beam_bending_h4_model.jnl
immersed_h4p3_traction–beam_bending_h4_model.jnl
9.3
Supported Abaqus keywords
The decks below are reproduced from the verification tests’ Abaqus input files with comments stripped and instance-specific names (sets, materials, files) redacted to
<placeholders>
.
Note
Gravity - beam_bending_h0.inp outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Material
, name
=
<material>
*Density
8000.,
*Elastic
2e+11,0.
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Dload
<element-set>
, U994, 0
<element-set>
, U995, -9.8100000000000005
<element-set>
, U996, 0
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
Note
Body force - beam_bending_h0.inp outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Material
, name
=
<material>
*Density
8000.,
*Elastic
2e+11,0.
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Dload
<element-set>
, U997, 0.
<element-set>
, U998, -78480.
<element-set>
, U999, 0.
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
Note
Pressure (uniform) - beam_bending_h0.inp outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Material
, name
=
<material>
*Density
8000.,
*Elastic
2e+11,0.
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
Note
Pressure (total force) - beam_bending_h0.inp outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Material
, name
=
<material>
*Density
8000.,
*Elastic
2e+11,0.
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
Note
Surface traction - beam_bending_h0.inp outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Material
, name
=
<material>
*Density
8000.,
*Elastic
2e+11,0.
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
Note
End moment - beam_bending_h0.inp outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Node
1,           0.,           0.,           1.
*Nset
, nset
=
<node-set>
1,
*Coupling
, constraint name
=moment_ref_pt_to_zmax_coupling
, ref node
=moment_ref_pt
, surface
=
<surface>
*Kinematic
*Material
, name
=
<material>
*Density
8000.,
*Elastic
2e+11,0.
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Cload
<node-set>
, 4, 5.
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
9.4
Provenance
The records below identify the execution environment and verification basis for every CTest artifact reported on this page.
Note
Test_abaqus_verification_beam_bending_abaqus_h0p1_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:14:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:23:16
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p1_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:45:52
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:47:24
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p1_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:59:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:13:50
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p1_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:22:01
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p1_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:16:37
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p1_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p1_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:48:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p2_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:17:54
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:23:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p2_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:46:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:47:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p2_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:48
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p2_body_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:58
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:16:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:42:08
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:19:24
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:45:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:34:40
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p1_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:14:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:19:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p1_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:45:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:47:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p1_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:00:08
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:13:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p1_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:59
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p1_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:16:03
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p1_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:52:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p1_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:42:08
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p2_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:17:53
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:23:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p2_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:46:51
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:47:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p2_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:33
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:23
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p2_gravity
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:53
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:52
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:39:44
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:23:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:46:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:12:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:01
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:31:51
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p1_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:15:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:23:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p1_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:47:06
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:47:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p1_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:00:19
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p1_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:44
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p1_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:16:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p1_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p1_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:37:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p2_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:20:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:19:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p2_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:46:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:46:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p2_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p2_moment
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:50
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:15:30
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:45:44
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:18:52
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:45:30
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:17:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:33:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p1_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:14:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p1_pressure_total_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:23:02
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p1_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:46:53
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p1_pressure_total_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:47:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p1_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:59:46
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p1_pressure_total_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:19
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p1_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:59
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p1_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:10:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p1_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:47:30
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p1_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:41:48
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p2_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:20:34
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:16:07
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p2_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:47:05
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:23:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p2_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:00
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:37
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p2_pressure_total_force
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:15:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:46:50
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:46:07
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p3_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:18:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p3_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:46:05
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p3_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:53
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p3_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:40
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p3_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:35:09
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p1_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:14:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p1_pressure_uniform
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:19:39
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p1_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:47:06
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p1_pressure_uniform
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:46:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p1_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:00:07
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p1_pressure_uniform
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p1_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p1_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:10:07
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p1_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:47:06
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p1_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:39:34
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p2_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:20:33
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:22:47
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p2_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:47:02
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:38:49
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p2_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:59:25
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:09
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p2_pressure_uniform
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:16:03
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:52:08
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:47:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p3_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:18:28
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p3_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:37:39
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p3_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:54
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p3_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:12:28
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p3_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:31:06
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p1_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:15:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:23:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p1_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:46:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:46:03
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p1_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:59:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:09
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p1_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:39
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p1_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:11:52
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p1_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p1_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:36:51
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h0p2_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:20:43
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:19:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h1p2_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 19:46:39
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:46:03
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h2p2_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:20:33
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:13:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h3p2_traction
Parameter
Value
Evidence type
Execution test
CTest result
Passed (Release)
Verification result
Not assessed — no structured numerical criteria
Test date
04-Nov-2026 20:21:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:16:25
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_abaqus_h4p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:47:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:46:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:19:19
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h1p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:45:54
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h2p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h3p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:11:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h4p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:37:22
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_traction_follower
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:19:27
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_beam_bending_immersed_h0p2_traction_resultant
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:23:38
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_study_beam_bending_abaqus_p1_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p1_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:52:05 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p1_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p1_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p1_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:57 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p1_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p2_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p2_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p2_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p2_pressure_total_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p2_pressure_uniform
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:52:09 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_abaqus_p2_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p1_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p1_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:42:09 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p1_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:37:19 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p1_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:41:49 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p1_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:39:34 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p1_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:36:51 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p2_body_force
Parameter
Value
Evidence type
Numerical verification with record-only diagnostics
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed; 1 record-only diagnostic
Test date
04-Nov-2026 20:42:08 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p2_gravity
Parameter
Value
Evidence type
Numerical verification with record-only diagnostics
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed; 1 record-only diagnostic
Test date
04-Nov-2026 20:39:45 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:51:58 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p2_pressure_total_force
Parameter
Value
Evidence type
Numerical verification with record-only diagnostics
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed; 1 record-only diagnostic
Test date
04-Nov-2026 20:51:57 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p2_pressure_uniform
Parameter
Value
Evidence type
Numerical verification with record-only diagnostics
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed; 1 record-only diagnostic
Test date
04-Nov-2026 20:51:56 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p2_traction
Parameter
Value
Evidence type
Numerical verification with record-only diagnostics
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed; 1 record-only diagnostic
Test date
04-Nov-2026 20:51:57 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:34:41 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:31:52 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:33:36 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p3_pressure_total_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:35:10 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p3_pressure_uniform
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:31:06 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_beam_bending_immersed_p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:37:23 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
