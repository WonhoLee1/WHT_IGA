# Cantilever bimetallic beam bending – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/linear_elastic_stress/bimetallic-beam-bending/bimetallic-beam-bending.html

---

Verification Manual
Solid mechanics
Linear elastic stress
Cantilever bimetallic beam bending
Cantilever bimetallic beam bending
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
\[
\newcommand{\modularRatio}{4}
\newcommand{\neutralAxis}{1.5}
\]
1
Objective
Verify Coreform IGA for Abaqus on a
two-material
problem: a bimetallic cantilever whose cross-section stacks a stiff layer on a flexible layer, bent by gravity, body force, surface pressure, surface traction, and a concentrated end moment. Each load case is checked against the transformed-section Euler-Bernoulli solution, whose bending stress is discontinuous at the material interface — the feature that distinguishes this problem from its single-material counterpart. The quantity of interest is the tip deflection; the deflection curve and the bending- and shear-stress distributions through the two-material section are also compared against the closed forms along line probes. Two interface workflows are exercised: the two layers merged into a single part sharing interface nodes, and separate parts joined by a tie constraint. All numerical inputs use the
millimeter–tonne–second (MMTS) consistent unit system
.
2
Geometry
A prismatic two-layer beam of width
\(b\)
and length
\(L\)
is shown in
Figure 1
, with its dimensions provided in
Table 1
. A flexible layer of height
\(h_f\)
lies below the interface plane, and a stiff layer of height
\(h_s\)
lies above it. The interface lies on the
\(y = 0\)
plane, the cross-section is centered on the
\(x\)
origin, and the beam axis runs along
\(z \in [0, L]\)
. Verification quantities are sampled along three line probes: the interface line at the section center (
\(u_y\)
), and the section lines through the two-material stack at midspan and at the clamped end (
\(\sigma_{zz}\)
and
\(\tau_{zy}\)
).
Figure 1: Bimetallic beam geometry: the stiff layer of height
\(h_s\)
lies above
\(y = 0\)
, and the flexible layer of height
\(h_f\)
lies below it; both layers have width
\(b\)
and length
\(L\)
.
Table 1: Geometric parameters.
Quantity
Symbol
Value
Unit
Section width
\(b\)
10
\(\mathrm{mm}\)
Flexible-layer height
\(h_f\)
5
\(\mathrm{mm}\)
Stiff-layer height
\(h_s\)
5
\(\mathrm{mm}\)
Length
\(L\)
100
\(\mathrm{mm}\)
Transformed second moment
\(I_t\)
\(1.52 \times 10^{3}\)
\(\mathrm{mm}^4\)
3
Material
Each layer is a homogeneous, isotropic, linear-elastic material (
Table 2
). Representative copper and tungsten moduli and densities give the properties a familiar engineering scale and provide a useful stiffness and density contrast; they are not intended as material allowables. The stiff layer is four times stiffer than the flexible layer, giving modular ratio
\(n = E_s/E_f = \modularRatio\)
. Poisson’s ratio is intentionally idealized to zero in both layers so the 3D solid remains compatible with the Euler-Bernoulli accepted solution, which carries no Poisson effect.
Table 2: Material properties per layer.
Layer
Material
Young’s modulus
\(E\)
[
\(\mathrm{N}/\mathrm{mm}^2\)
]
Poisson’s ratio
\(\nu\)
Mass density
\(\rho\)
[
\(\mathrm{tonne}/\mathrm{mm}^3\)
]
Flexible (
\(-h_f \le y < 0\)
)
copper
\(100 \times 10^{3}\)
0
\(8.90 \times 10^{-9}\)
Stiff (
\(0 \le y \le h_s\)
)
tungsten
\(400 \times 10^{3}\)
0
\(19.2 \times 10^{-9}\)
4
Loading & boundary conditions
The idealized composite beam is clamped at
\(z = 0\)
, as shown in
Figure 2
. In the three-dimensional solid models, this condition is imposed by pinning the
\(z = 0\)
face (all translations fixed). Each load case below is applied in its own study. The gravity and body-force cases use identical meshes, material definitions, boundary conditions, and perfectly tied layer interfaces; only the load definition differs.
Figure 2: Idealized clamped boundary condition at
\(z = 0\)
.
4.1
Gravity
A gravity load of magnitude
\(g = 9810\ \mathrm{mm}/\mathrm{s}^2\)
acts in
\(-y\)
on the whole body, as shown in
Figure 3
. The flexible and stiff layers contribute their own weights, so the equivalent distributed line load is
\(w = g\,b\,(\rho_f h_f + \rho_s h_s)\)
.
Figure 3: Gravity loading on flexible and stiff layers with densities
\(\rho_f\)
and
\(\rho_s\)
.
4.2
Body force
An explicitly equivalent body force per unit volume acts in
\(-y\)
, as shown in
Figure 4
. Its magnitude is defined separately in each layer as
\[
q_f = \rho_f g,
\qquad
q_s = \rho_s g.
\]
The gravity translator supplies
\(g\)
, and Coreform IGA multiplies it by the material density at each quadrature point; the explicit body-force case supplies the resulting
\(q_i\)
directly. Consequently, the two load cases produce the same force density pointwise, rather than merely the same section resultant. The two layers remain separate load targets and are joined by a perfect tie, which represents the ideal bond assumed by transformed-section beam theory.
Figure 4: Layer-specific body-force loading in the
\(-y\)
direction, with
\(q_f=\rho_f g\)
in copper and
\(q_s=\rho_s g\)
in tungsten.
4.3
Pressure
A uniform normal pressure is applied to the stiff layer’s top face (
\(y = +h_s\)
), as shown in
Figure 5
. It reduces to the equivalent distributed line load
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
\(P = t\,b\,(h_f + h_s)\)
shown in
Figure 6
. This load case is also run with the tie-constraint interface workflow.
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
Unit
Gravity
gravity
\(9.81 \times 10^{3}\)
\(\mathrm{mm}/\mathrm{s}^2\)
Body force
body_force
Copper:
\(87.3 \times 10^{-6}\)
Tungsten:
\(188 \times 10^{-6}\)
\(\mathrm{N}/\mathrm{mm}^3\)
Pressure
pressure
\(1.00 \times 10^{-3}\)
\(\mathrm{N}/\mathrm{mm}^2\)
Surface traction
traction
\(10.0 \times 10^{-3}\)
\(\mathrm{N}/\mathrm{mm}^2\)
End moment
moment
\(5\)
\(\mathrm{N}\mathbin{\cdot}\mathrm{mm}\)
5
Mesh & discretization
Coreform IGA discretizes the beam with
immersed, trimmed U-splines
: a structured rectilinear background grid with
\(p\)
padding elements per side is trimmed to the beam. The background element sizes are
\((5, 5, 10)/2^{k}\)
in
\((x, y, z)\)
for refinement levels
\(k = 0,\dots,3\)
, at degrees
\(p = 1\)
,
\(2\)
, and
\(3\)
with maximum continuity
\(C^{p-1}\)
within each layer. The material interface is resolved by the geometry rather than the grid: in the merged workflow the two layers form one part whose interface faces share nodes, and in the tie workflow the layers remain separate parts coupled by a tie constraint. Every load case runs on this same mesh sweep, so the meshes are shown once — the tabs below render each refinement level from the gravity study’s IGA results databases.
\(h_z\)
=10
\(h_z\)
=5
\(h_z\)
=2.5
\(h_z\)
=1.25
6
Reference & accepted solutions
The
accepted
solutions are the transformed-section Euler-Bernoulli closed forms
[1]
: the stiff layer is scaled by the modular ratio
\(n = E_s/E_f\)
into an equivalent single-material section, giving the neutral-axis height and the transformed second moment of area
\[
\bar{y} = \frac{\hat{y}_s\, n A_s + \hat{y}_f\, A_f}{n A_s + A_f},
\qquad
I_t = \frac{b h_f^3}{12} + A_f\,(\hat{y}_f - \bar{y})^2 + n\,\frac{b h_s^3}{12} + n A_s\,(\hat{y}_s - \bar{y})^2 ,
\]
where
\(A_f = b h_f\)
and
\(A_s = b h_s\)
are the layer areas and
\(\hat{y}_f\)
,
\(\hat{y}_s\)
their centroid heights. For this section the neutral axis sits at
\(y = \neutralAxis\)
— inside the stiff layer, above the material interface. The section resultants
\(V(z)\)
,
\(M(z)\)
and the deflection
\(u_y(z)\)
then follow the standard cantilever closed forms with flexural rigidity
\(E_f I_t\)
; the distributed cases (gravity, body force, pressure) reduce to a uniform line load
\(w\)
, the traction case to a tip load
\(P\)
, and the end moment to a constant moment
\(M_0\)
. The stresses recover the two-material behavior:
\[
\sigma_{zz}(y, z) = -\frac{E(y)}{E_f}\,\frac{M(z)\,(y - \bar{y})}{I_t},
\qquad
\tau_{zy}(y, z) = \frac{V(z)\, Q(y)}{I_t\, b},
\]
so the bending stress is
discontinuous
at the material interface — it jumps by the factor
\(n\)
— while the shear stress, built from the transformed first moment
\(Q(y)\)
, remains continuous. The transformed-section model is itself an approximation of the 3D solid (it neglects shear deformation, and the pinned end face is stiffer than the idealized support), so exact agreement is not expected; the acceptance criterion is the tip deflection matching within 5%.
Table 4: Equivalent loads and transformed-section tip deflections.
Load case
Equivalent load
Transformed-section
\(u_y(L)\)
[
\(\mathrm{mm}\)
]
Gravity
\(w\)
=
\(13.8 \times 10^{-3}\)
\(\;\mathrm{N}/\mathrm{mm}\)
\(-1.13 \times 10^{-3}\)
Body force
\(w\)
=
\(13.8 \times 10^{-3}\)
\(\;\mathrm{N}/\mathrm{mm}\)
\(-1.13 \times 10^{-3}\)
Pressure
\(w\)
=
\(10.0 \times 10^{-3}\)
\(\;\mathrm{N}/\mathrm{mm}\)
\(-822 \times 10^{-6}\)
Surface traction
\(P\)
=
\(1\)
\(\;\mathrm{N}\)
\(-2.19 \times 10^{-3}\)
End moment
\(M_0\)
=
\(5\)
\(\;\mathrm{N}\mathbin{\cdot}\mathrm{mm}\)
\(-164 \times 10^{-6}\)
Download data (CSV)
7
Results
Each load case below reports the tip deflection against the transformed-section value (rows are the
\(h\)
-sweep; one tab per degree, and the degree tabs switch together across load cases), the four verification profiles against the closed forms (one figure per probe quantity, a subfigure per degree), full-span line-probe absolute
\(L^2\)
self-convergence against the same-degree
\(h_3\)
profile, and the axial-stress contours for every test case on a shared color scale capped at the transformed-section peak. The midspan and clamped-section stress profiles cross the material interface, where the closed-form bending stress jumps by the modular ratio. The trilinear
\(p=1\)
cases can lock in bending on the coarser meshes, so those levels are retained as record-only diagnostics; the finest
\(h_3\)
result carries the same 5% tip-deflection acceptance criterion as the higher-degree cases.
7.1
Gravity
Transformed-section tip deflection:
\(-1.13 \times 10^{-3}\)
.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-500 \times 10^{-6}\)
\(55.8\%\)
5
\(-880 \times 10^{-6}\)
\(22.3\%\)
2.5
\(-1.08 \times 10^{-3}\)
\(4.27\%\)
1.25
\(-1.13 \times 10^{-3}\)
\(0.42\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-1.14 \times 10^{-3}\)
\(0.66\%\)
5
\(-1.14 \times 10^{-3}\)
\(0.689\%\)
2.5
\(-1.14 \times 10^{-3}\)
\(0.72\%\)
1.25
\(-1.14 \times 10^{-3}\)
\(0.72\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-1.14 \times 10^{-3}\)
\(0.717\%\)
5
\(-1.14 \times 10^{-3}\)
\(0.687\%\)
2.5
\(-1.14 \times 10^{-3}\)
\(0.589\%\)
1.25
\(-1.14 \times 10^{-3}\)
\(0.72\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 8: Gravity: deflection along the interface line — coarse→fine sweep per degree vs the transformed section.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 9: Gravity: bending stress through the midspan section — the closed form jumps by the modular ratio at the material interface (
\(y = 0\)
).
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 10: Gravity: shear stress through the midspan section — continuous across the interface, with its maximum at the neutral axis.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 11: Gravity: bending stress through the clamped section, where the moment (and the interface jump) is largest.
(a) Interface deflection
\(u_y\)
(b) Midspan stress
\(\sigma_{zz}\)
(c) Midspan shear
\(\tau_{zy}\)
(d) Clamped-section stress
\(\sigma_{zz}\)
Figure 12: Gravity: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h3 profile; the h3 point is omitted because its self-error is zero.
7.1.1
Contour plots
h0p1
h1p1
h2p1
h3p1
h0p2
h1p2
h2p2
h3p2
h0p3
h1p3
h2p3
h3p3
7.2
Body force
Transformed-section tip deflection:
\(-1.13 \times 10^{-3}\)
.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-500 \times 10^{-6}\)
\(55.8\%\)
5
\(-880 \times 10^{-6}\)
\(22.3\%\)
2.5
\(-1.08 \times 10^{-3}\)
\(4.27\%\)
1.25
\(-1.13 \times 10^{-3}\)
\(0.42\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-1.14 \times 10^{-3}\)
\(0.66\%\)
5
\(-1.14 \times 10^{-3}\)
\(0.689\%\)
2.5
\(-1.14 \times 10^{-3}\)
\(0.72\%\)
1.25
\(-1.14 \times 10^{-3}\)
\(0.72\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-1.14 \times 10^{-3}\)
\(0.717\%\)
5
\(-1.14 \times 10^{-3}\)
\(0.687\%\)
2.5
\(-1.14 \times 10^{-3}\)
\(0.589\%\)
1.25
\(-1.14 \times 10^{-3}\)
\(0.72\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 13: Body force: deflection along the interface line — coarse→fine sweep per degree vs the transformed section.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 14: Body force: bending stress through the midspan section — the closed form jumps by the modular ratio at the material interface (
\(y = 0\)
).
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 15: Body force: shear stress through the midspan section — continuous across the interface, with its maximum at the neutral axis.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 16: Body force: bending stress through the clamped section, where the moment (and the interface jump) is largest.
(a) Interface deflection
\(u_y\)
(b) Midspan stress
\(\sigma_{zz}\)
(c) Midspan shear
\(\tau_{zy}\)
(d) Clamped-section stress
\(\sigma_{zz}\)
Figure 17: Body force: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h3 profile; the h3 point is omitted because its self-error is zero.
7.2.1
Contour plots
h0p1
h1p1
h2p1
h3p1
h0p2
h1p2
h2p2
h3p2
h0p3
h1p3
h2p3
h3p3
7.3
Pressure
Transformed-section tip deflection:
\(-822 \times 10^{-6}\)
.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-342 \times 10^{-6}\)
\(58.4\%\)
5
\(-637 \times 10^{-6}\)
\(22.4\%\)
2.5
\(-787 \times 10^{-6}\)
\(4.27\%\)
1.25
\(-818 \times 10^{-6}\)
\(0.417\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-827 \times 10^{-6}\)
\(0.646\%\)
5
\(-828 \times 10^{-6}\)
\(0.716\%\)
2.5
\(-828 \times 10^{-6}\)
\(0.72\%\)
1.25
\(-828 \times 10^{-6}\)
\(0.72\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-828 \times 10^{-6}\)
\(0.72\%\)
5
\(-828 \times 10^{-6}\)
\(0.72\%\)
2.5
\(-828 \times 10^{-6}\)
\(0.721\%\)
1.25
\(-828 \times 10^{-6}\)
\(0.72\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 18: Pressure: deflection along the interface line — coarse→fine sweep per degree vs the transformed section.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 19: Pressure: bending stress through the midspan section — the closed form jumps by the modular ratio at the material interface (
\(y = 0\)
).
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 20: Pressure: shear stress through the midspan section — continuous across the interface, with its maximum at the neutral axis.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 21: Pressure: bending stress through the clamped section, where the moment (and the interface jump) is largest.
(a) Interface deflection
\(u_y\)
(b) Midspan stress
\(\sigma_{zz}\)
(c) Midspan shear
\(\tau_{zy}\)
(d) Clamped-section stress
\(\sigma_{zz}\)
Figure 22: Pressure: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h3 profile; the h3 point is omitted because its self-error is zero.
7.3.1
Contour plots
h0p1
h1p1
h2p1
h3p1
h0p2
h1p2
h2p2
h3p2
h0p3
h1p3
h2p3
h3p3
7.4
Surface traction
Transformed-section tip deflection:
\(-2.19 \times 10^{-3}\)
.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-874 \times 10^{-6}\)
\(60.1\%\)
5
\(-1.69 \times 10^{-3}\)
\(23\%\)
2.5
\(-2.09 \times 10^{-3}\)
\(4.51\%\)
1.25
\(-2.18 \times 10^{-3}\)
\(0.606\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-2.20 \times 10^{-3}\)
\(0.485\%\)
5
\(-2.20 \times 10^{-3}\)
\(0.537\%\)
2.5
\(-2.20 \times 10^{-3}\)
\(0.539\%\)
1.25
\(-2.20 \times 10^{-3}\)
\(0.539\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-2.20 \times 10^{-3}\)
\(0.54\%\)
5
\(-2.20 \times 10^{-3}\)
\(0.54\%\)
2.5
\(-2.20 \times 10^{-3}\)
\(0.54\%\)
1.25
\(-2.20 \times 10^{-3}\)
\(0.54\%\)
Download data (CSV)
The traction studies also run with the
tie-constraint
interface workflow: the two layers remain separate parts and are joined by a tie constraint instead of sharing interface nodes. Both workflows discretize the same problem, so their tip deflections must agree (
Table 5
).
Table 5: Interface workflows compared on the traction case.
Degree
Interface workflow
\(u_y\)
tip (finest level)
vs transformed section
p=1
Merged part (shared interface nodes)
\(-2.18 \times 10^{-3}\)
\(0.606\%\)
p=1
Tie constraint between separate parts
\(-2.18 \times 10^{-3}\)
\(0.609\%\)
p=2
Merged part (shared interface nodes)
\(-2.20 \times 10^{-3}\)
\(0.539\%\)
p=2
Tie constraint between separate parts
\(-2.20 \times 10^{-3}\)
\(0.539\%\)
p=3
Merged part (shared interface nodes)
\(-2.20 \times 10^{-3}\)
\(0.54\%\)
p=3
Tie constraint between separate parts
\(-2.20 \times 10^{-3}\)
\(0.54\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 23: Surface traction: deflection along the interface line — coarse→fine sweep per degree vs the transformed section.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 24: Surface traction: bending stress through the midspan section — the closed form jumps by the modular ratio at the material interface (
\(y = 0\)
).
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 25: Surface traction: shear stress through the midspan section — continuous across the interface, with its maximum at the neutral axis.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 26: Surface traction: bending stress through the clamped section, where the moment (and the interface jump) is largest.
(a) Interface deflection
\(u_y\)
(b) Midspan stress
\(\sigma_{zz}\)
(c) Midspan shear
\(\tau_{zy}\)
(d) Clamped-section stress
\(\sigma_{zz}\)
Figure 27: Surface traction: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h3 profile; the h3 point is omitted because its self-error is zero.
7.4.1
Contour plots
h0p1
h1p1
h2p1
h3p1
h0p2
h1p2
h2p2
h3p2
h0p3
h1p3
h2p3
h3p3
7.5
End moment
Transformed-section tip deflection:
\(-164 \times 10^{-6}\)
.
p=1
p=2
p=3
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-62.5 \times 10^{-6}\)
\(62\%\)
5
\(-125 \times 10^{-6}\)
\(23.8\%\)
2.5
\(-156 \times 10^{-6}\)
\(5.06\%\)
1.25
\(-163 \times 10^{-6}\)
\(1.14\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-164 \times 10^{-6}\)
\(0.00291\%\)
5
\(-164 \times 10^{-6}\)
\(0.00157\%\)
2.5
\(-164 \times 10^{-6}\)
\(0.0014\%\)
1.25
\(-164 \times 10^{-6}\)
\(0.00122\%\)
\(h_z\)
\(u_y\)
tip
vs transformed section
10
\(-164 \times 10^{-6}\)
\(0.00136\%\)
5
\(-164 \times 10^{-6}\)
\(0.00152\%\)
2.5
\(-164 \times 10^{-6}\)
\(0.00155\%\)
1.25
\(-164 \times 10^{-6}\)
\(0.00134\%\)
Download data (CSV)
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 28: End moment: deflection along the interface line — coarse→fine sweep per degree vs the transformed section.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 29: End moment: bending stress through the midspan section — the closed form jumps by the modular ratio at the material interface (
\(y = 0\)
).
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 30: End moment: shear stress through the midspan section — the closed-form shear is identically zero for a constant moment.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 31: End moment: bending stress through the clamped section — for a constant moment the interface jump is the same at every section.
(a) Interface deflection
\(u_y\)
(b) Midspan stress
\(\sigma_{zz}\)
(c) Midspan shear
\(\tau_{zy}\)
(d) Clamped-section stress
\(\sigma_{zz}\)
Figure 32: End moment: full-span line-probe absolute L2 self-convergence. Each degree is compared against its own h3 profile; the h3 point is omitted because its self-error is zero.
7.5.1
Contour plots
h0p1
h1p1
h2p1
h3p1
h0p2
h1p2
h2p2
h3p2
h0p3
h1p3
h2p3
h3p3
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
\(10\)
\(195\)
\(66\)
\(24\)
\(1\)
\(5.3\)
\(5\)
\(805\)
\(630\)
\(258\)
\(2\)
\(17.2\)
\(2.5\)
\(3{,}311\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(48.9\)
\(1.25\)
\(17{,}347\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(204.0\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(525\)
\(66\)
\(24\)
\(1\)
\(7.7\)
\(5\)
\(1{,}575\)
\(630\)
\(258\)
\(2\)
\(20.6\)
\(2.5\)
\(5{,}265\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(56.8\)
\(1.25\)
\(23{,}205\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(206.9\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}071\)
\(66\)
\(24\)
\(1\)
\(115.0\)
\(5\)
\(2{,}673\)
\(630\)
\(258\)
\(2\)
\(58.7\)
\(2.5\)
\(7{,}755\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(100.0\)
\(1.25\)
\(30{,}015\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(258.6\)
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
\(10\)
\(5{,}194\)
\(2{,}998\)
\(8{,}994\)
\(1\)
\(3\)
\(5\)
\(6{,}382\)
\(4{,}184\)
\(12{,}552\)
\(2\)
\(2\)
\(2.5\)
\(10{,}854\)
\(8{,}914\)
\(26{,}742\)
\(4\)
\(1\)
\(1.25\)
\(45{,}354\)
\(40{,}272\)
\(120{,}816\)
\(8\)
\(7\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}154\)
\(3{,}270\)
\(9{,}810\)
\(1\)
\(4\)
\(5\)
\(6{,}382\)
\(4{,}738\)
\(14{,}214\)
\(2\)
\(2\)
\(2.5\)
\(10{,}854\)
\(10{,}496\)
\(31{,}488\)
\(4\)
\(5\)
\(1.25\)
\(44{,}970\)
\(45{,}248\)
\(135{,}744\)
\(8\)
\(18\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(4{,}810\)
\(3{,}328\)
\(9{,}984\)
\(1\)
\(4\)
\(5\)
\(6{,}734\)
\(5{,}616\)
\(16{,}848\)
\(2\)
\(8\)
\(2.5\)
\(13{,}512\)
\(13{,}679\)
\(41{,}037\)
\(4\)
\(17\)
\(1.25\)
\(59{,}674\)
\(58{,}396\)
\(175{,}188\)
\(8\)
\(70\)
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
\(10\)
\(195\)
\(66\)
\(24\)
\(1\)
\(5.1\)
\(5\)
\(805\)
\(630\)
\(258\)
\(2\)
\(17.8\)
\(2.5\)
\(3{,}311\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(56.5\)
\(1.25\)
\(17{,}347\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(198.1\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(525\)
\(66\)
\(24\)
\(1\)
\(7.8\)
\(5\)
\(1{,}575\)
\(630\)
\(258\)
\(2\)
\(19.9\)
\(2.5\)
\(5{,}265\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(56.7\)
\(1.25\)
\(23{,}205\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(204.3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}071\)
\(66\)
\(24\)
\(1\)
\(69.4\)
\(5\)
\(2{,}673\)
\(630\)
\(258\)
\(2\)
\(57.7\)
\(2.5\)
\(7{,}755\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(99.9\)
\(1.25\)
\(30{,}015\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(260.7\)
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
\(10\)
\(5{,}194\)
\(2{,}998\)
\(8{,}994\)
\(1\)
\(3\)
\(5\)
\(6{,}382\)
\(4{,}184\)
\(12{,}552\)
\(2\)
\(2\)
\(2.5\)
\(10{,}854\)
\(8{,}914\)
\(26{,}742\)
\(4\)
\(2\)
\(1.25\)
\(45{,}354\)
\(40{,}272\)
\(120{,}816\)
\(8\)
\(7\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}154\)
\(3{,}270\)
\(9{,}810\)
\(1\)
\(2\)
\(5\)
\(6{,}382\)
\(4{,}738\)
\(14{,}214\)
\(2\)
\(2\)
\(2.5\)
\(10{,}854\)
\(10{,}496\)
\(31{,}488\)
\(4\)
\(5\)
\(1.25\)
\(44{,}970\)
\(45{,}248\)
\(135{,}744\)
\(8\)
\(17\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(4{,}810\)
\(3{,}328\)
\(9{,}984\)
\(1\)
\(3\)
\(5\)
\(6{,}734\)
\(5{,}616\)
\(16{,}848\)
\(2\)
\(8\)
\(2.5\)
\(13{,}512\)
\(13{,}679\)
\(41{,}037\)
\(4\)
\(18\)
\(1.25\)
\(59{,}674\)
\(58{,}396\)
\(175{,}188\)
\(8\)
\(72\)
Download data (CSV)
8.3
Pressure
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
\(10\)
\(325\)
\(132\)
\(66\)
\(1\)
\(7.4\)
\(5\)
\(1{,}127\)
\(630\)
\(258\)
\(2\)
\(17.2\)
\(2.5\)
\(5{,}203\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(47.0\)
\(1.25\)
\(29{,}963\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(199.1\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(66\)
\(1\)
\(8.5\)
\(5\)
\(2{,}025\)
\(630\)
\(258\)
\(2\)
\(18.5\)
\(2.5\)
\(7{,}605\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(59.1\)
\(1.25\)
\(37{,}485\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(208.8\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(66\)
\(1\)
\(54.5\)
\(5\)
\(3{,}267\)
\(630\)
\(258\)
\(2\)
\(38.6\)
\(2.5\)
\(10{,}575\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(83.9\)
\(1.25\)
\(46{,}023\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(251.6\)
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
\(10\)
\(7{,}884\)
\(4{,}537\)
\(13{,}611\)
\(1\)
\(3\)
\(5\)
\(9{,}100\)
\(5{,}676\)
\(17{,}028\)
\(2\)
\(3\)
\(2.5\)
\(14{,}180\)
\(10{,}732\)
\(32{,}196\)
\(4\)
\(2\)
\(1.25\)
\(54{,}824\)
\(45{,}258\)
\(135{,}774\)
\(8\)
\(6\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}650\)
\(4{,}651\)
\(13{,}953\)
\(1\)
\(2\)
\(5\)
\(9{,}100\)
\(6{,}230\)
\(18{,}690\)
\(2\)
\(3\)
\(2.5\)
\(14{,}180\)
\(12{,}314\)
\(36{,}942\)
\(4\)
\(4\)
\(1.25\)
\(54{,}248\)
\(50{,}135\)
\(150{,}405\)
\(8\)
\(17\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}762\)
\(5{,}033\)
\(15{,}099\)
\(1\)
\(9\)
\(5\)
\(9{,}628\)
\(7{,}200\)
\(21{,}600\)
\(2\)
\(8\)
\(2.5\)
\(17{,}586\)
\(15{,}889\)
\(47{,}667\)
\(4\)
\(18\)
\(1.25\)
\(75{,}482\)
\(66{,}641\)
\(199{,}923\)
\(8\)
\(69\)
Download data (CSV)
8.4
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
\(10\)
\(325\)
\(132\)
\(66\)
\(1\)
\(7.3\)
\(5\)
\(1{,}127\)
\(630\)
\(258\)
\(2\)
\(17.1\)
\(2.5\)
\(5{,}203\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(49.7\)
\(1.25\)
\(29{,}963\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(208.2\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(66\)
\(1\)
\(8.7\)
\(5\)
\(2{,}025\)
\(630\)
\(258\)
\(2\)
\(19.2\)
\(2.5\)
\(7{,}605\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(58.4\)
\(1.25\)
\(37{,}485\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(206.9\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(66\)
\(1\)
\(47.0\)
\(5\)
\(3{,}267\)
\(630\)
\(258\)
\(2\)
\(39.5\)
\(2.5\)
\(10{,}575\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(84.8\)
\(1.25\)
\(46{,}023\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(257.2\)
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
\(10\)
\(5{,}612\)
\(3{,}336\)
\(10{,}008\)
\(1\)
\(1\)
\(5\)
\(6{,}698\)
\(4{,}412\)
\(13{,}236\)
\(2\)
\(1\)
\(2.5\)
\(11{,}366\)
\(9{,}254\)
\(27{,}762\)
\(4\)
\(1\)
\(1.25\)
\(46{,}982\)
\(41{,}250\)
\(123{,}750\)
\(8\)
\(7\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}464\)
\(3{,}490\)
\(10{,}470\)
\(1\)
\(3\)
\(5\)
\(6{,}698\)
\(4{,}966\)
\(14{,}898\)
\(2\)
\(2\)
\(2.5\)
\(11{,}366\)
\(10{,}836\)
\(32{,}508\)
\(4\)
\(5\)
\(1.25\)
\(47{,}018\)
\(46{,}436\)
\(139{,}308\)
\(8\)
\(17\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}560\)
\(3{,}862\)
\(11{,}586\)
\(1\)
\(9\)
\(5\)
\(7{,}050\)
\(5{,}844\)
\(17{,}532\)
\(2\)
\(6\)
\(2.5\)
\(14{,}118\)
\(14{,}112\)
\(42{,}336\)
\(4\)
\(17\)
\(1.25\)
\(62{,}946\)
\(60{,}282\)
\(180{,}846\)
\(8\)
\(69\)
Download data (CSV)
8.5
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
\(10\)
\(325\)
\(132\)
\(66\)
\(1\)
\(7.3\)
\(5\)
\(1{,}127\)
\(630\)
\(258\)
\(2\)
\(17.5\)
\(2.5\)
\(5{,}203\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(47.1\)
\(1.25\)
\(29{,}963\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(208.2\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(66\)
\(1\)
\(9.0\)
\(5\)
\(2{,}025\)
\(630\)
\(258\)
\(2\)
\(19.1\)
\(2.5\)
\(7{,}605\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(57.8\)
\(1.25\)
\(37{,}485\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(216.5\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(66\)
\(1\)
\(55.9\)
\(5\)
\(3{,}267\)
\(630\)
\(258\)
\(2\)
\(39.4\)
\(2.5\)
\(10{,}575\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(81.6\)
\(1.25\)
\(46{,}023\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(284.6\)
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
\(10\)
\(5{,}612\)
\(3{,}337\)
\(10{,}479\)
\(1\)
\(2\)
\(5\)
\(6{,}698\)
\(4{,}413\)
\(13{,}842\)
\(2\)
\(2\)
\(2.5\)
\(11{,}366\)
\(9{,}255\)
\(28{,}686\)
\(4\)
\(3\)
\(1.25\)
\(46{,}982\)
\(41{,}251\)
\(126{,}492\)
\(8\)
\(7\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}464\)
\(3{,}491\)
\(10{,}932\)
\(1\)
\(4\)
\(5\)
\(6{,}698\)
\(4{,}967\)
\(15{,}504\)
\(2\)
\(2\)
\(2.5\)
\(11{,}366\)
\(10{,}837\)
\(33{,}432\)
\(4\)
\(6\)
\(1.25\)
\(47{,}018\)
\(46{,}437\)
\(142{,}680\)
\(8\)
\(17\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}560\)
\(3{,}863\)
\(12{,}096\)
\(1\)
\(10\)
\(5\)
\(7{,}050\)
\(5{,}845\)
\(18{,}138\)
\(2\)
\(9\)
\(2.5\)
\(14{,}118\)
\(14{,}113\)
\(43{,}884\)
\(4\)
\(23\)
\(1.25\)
\(62{,}946\)
\(60{,}283\)
\(186{,}210\)
\(8\)
\(69\)
Download data (CSV)
9
Discussion
Every load case converges under mesh refinement, with the tip deflection inside the acceptance tolerance on the finer meshes and the profile comparisons tracking the transformed-section closed forms — including the bending-stress discontinuity at the material interface, which the two-material discretization must capture without smearing. Against the transformed-section model the error plateaus: the residual is the beam theory’s own modeling error (no shear deformation, and the pinned end face is stiffer than the idealized support), not a discretization error, so the line-probe convergence figures use self-convergence against the same-degree finest profiles instead of treating the closed form as truth. That the merged-part and tie-constraint interface workflows land on the same solution verifies that the tie constraint transmits the interface tractions consistently with the shared-node discretization. The gravity and body-force cases have identical pointwise force-density fields and are required to produce probe results that agree to within solver roundoff, directly verifying the translation of
GRAV
into density-weighted IGA body-force components. The pressure case belongs to the same uniformly distributed-load closed-form family but has its own independently specified magnitude.
References
[1]
J. M. Gere and S. P. Timoshenko,
Mechanics of materials
, 4th ed. Boston: PWS Publishing Company, 1997.
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
Surface tie constraints
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
Uniform pressure loading
Numerically validated
Surface traction loading
Numerically validated
Mechanical loads
Numerically validated
Cell-based material section assignment
Numerically validated
Material density
Numerically validated
Multi-cell material interfaces
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
Boolean merge with retained boundaries
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
72 generated Abaqus/CAE journals are available.
Show 72 individual journals
immersed_h0p1_body_force–bimetallic_beam_bending_h0_model.jnl
immersed_h0p1_gravity–bimetallic_beam_bending_h0_model.jnl
immersed_h0p1_moment–bimetallic_beam_bending_h0_model.jnl
immersed_h0p1_pressure–bimetallic_beam_bending_h0_model.jnl
immersed_h0p1_traction–bimetallic_beam_bending_h0_model.jnl
immersed_h0p1_traction_tie–bimetallic_beam_bending_h0_model.jnl
immersed_h0p2_body_force–bimetallic_beam_bending_h0_model.jnl
immersed_h0p2_gravity–bimetallic_beam_bending_h0_model.jnl
immersed_h0p2_moment–bimetallic_beam_bending_h0_model.jnl
immersed_h0p2_pressure–bimetallic_beam_bending_h0_model.jnl
immersed_h0p2_traction–bimetallic_beam_bending_h0_model.jnl
immersed_h0p2_traction_tie–bimetallic_beam_bending_h0_model.jnl
immersed_h0p3_body_force–bimetallic_beam_bending_h0_model.jnl
immersed_h0p3_gravity–bimetallic_beam_bending_h0_model.jnl
immersed_h0p3_moment–bimetallic_beam_bending_h0_model.jnl
immersed_h0p3_pressure–bimetallic_beam_bending_h0_model.jnl
immersed_h0p3_traction–bimetallic_beam_bending_h0_model.jnl
immersed_h0p3_traction_tie–bimetallic_beam_bending_h0_model.jnl
immersed_h1p1_body_force–bimetallic_beam_bending_h1_model.jnl
immersed_h1p1_gravity–bimetallic_beam_bending_h1_model.jnl
immersed_h1p1_moment–bimetallic_beam_bending_h1_model.jnl
immersed_h1p1_pressure–bimetallic_beam_bending_h1_model.jnl
immersed_h1p1_traction–bimetallic_beam_bending_h1_model.jnl
immersed_h1p1_traction_tie–bimetallic_beam_bending_h1_model.jnl
immersed_h1p2_body_force–bimetallic_beam_bending_h1_model.jnl
immersed_h1p2_gravity–bimetallic_beam_bending_h1_model.jnl
immersed_h1p2_moment–bimetallic_beam_bending_h1_model.jnl
immersed_h1p2_pressure–bimetallic_beam_bending_h1_model.jnl
immersed_h1p2_traction–bimetallic_beam_bending_h1_model.jnl
immersed_h1p2_traction_tie–bimetallic_beam_bending_h1_model.jnl
immersed_h1p3_body_force–bimetallic_beam_bending_h1_model.jnl
immersed_h1p3_gravity–bimetallic_beam_bending_h1_model.jnl
immersed_h1p3_moment–bimetallic_beam_bending_h1_model.jnl
immersed_h1p3_pressure–bimetallic_beam_bending_h1_model.jnl
immersed_h1p3_traction–bimetallic_beam_bending_h1_model.jnl
immersed_h1p3_traction_tie–bimetallic_beam_bending_h1_model.jnl
immersed_h2p1_body_force–bimetallic_beam_bending_h2_model.jnl
immersed_h2p1_gravity–bimetallic_beam_bending_h2_model.jnl
immersed_h2p1_moment–bimetallic_beam_bending_h2_model.jnl
immersed_h2p1_pressure–bimetallic_beam_bending_h2_model.jnl
immersed_h2p1_traction–bimetallic_beam_bending_h2_model.jnl
immersed_h2p1_traction_tie–bimetallic_beam_bending_h2_model.jnl
immersed_h2p2_body_force–bimetallic_beam_bending_h2_model.jnl
immersed_h2p2_gravity–bimetallic_beam_bending_h2_model.jnl
immersed_h2p2_moment–bimetallic_beam_bending_h2_model.jnl
immersed_h2p2_pressure–bimetallic_beam_bending_h2_model.jnl
immersed_h2p2_traction–bimetallic_beam_bending_h2_model.jnl
immersed_h2p2_traction_tie–bimetallic_beam_bending_h2_model.jnl
immersed_h2p3_body_force–bimetallic_beam_bending_h2_model.jnl
immersed_h2p3_gravity–bimetallic_beam_bending_h2_model.jnl
immersed_h2p3_moment–bimetallic_beam_bending_h2_model.jnl
immersed_h2p3_pressure–bimetallic_beam_bending_h2_model.jnl
immersed_h2p3_traction–bimetallic_beam_bending_h2_model.jnl
immersed_h2p3_traction_tie–bimetallic_beam_bending_h2_model.jnl
immersed_h3p1_body_force–bimetallic_beam_bending_h3_model.jnl
immersed_h3p1_gravity–bimetallic_beam_bending_h3_model.jnl
immersed_h3p1_moment–bimetallic_beam_bending_h3_model.jnl
immersed_h3p1_pressure–bimetallic_beam_bending_h3_model.jnl
immersed_h3p1_traction–bimetallic_beam_bending_h3_model.jnl
immersed_h3p1_traction_tie–bimetallic_beam_bending_h3_model.jnl
immersed_h3p2_body_force–bimetallic_beam_bending_h3_model.jnl
immersed_h3p2_gravity–bimetallic_beam_bending_h3_model.jnl
immersed_h3p2_moment–bimetallic_beam_bending_h3_model.jnl
immersed_h3p2_pressure–bimetallic_beam_bending_h3_model.jnl
immersed_h3p2_traction–bimetallic_beam_bending_h3_model.jnl
immersed_h3p2_traction_tie–bimetallic_beam_bending_h3_model.jnl
immersed_h3p3_body_force–bimetallic_beam_bending_h3_model.jnl
immersed_h3p3_gravity–bimetallic_beam_bending_h3_model.jnl
immersed_h3p3_moment–bimetallic_beam_bending_h3_model.jnl
immersed_h3p3_pressure–bimetallic_beam_bending_h3_model.jnl
immersed_h3p3_traction–bimetallic_beam_bending_h3_model.jnl
immersed_h3p3_traction_tie–bimetallic_beam_bending_h3_model.jnl
9.3
Supported Abaqus keywords
The decks below are reproduced from the verification tests’ Abaqus input files with comments stripped and instance-specific names (sets, materials, files) redacted to
<placeholders>
.
Note
Gravity - bimetallic_beam_bending_h0.inp outline
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
0.,        -2.5,          0.,          1.,        -2.5,          0.
*System
0.,         2.5,          0.,          1.,         2.5,          0.
*System
*Tie
, name
=
<name>
, adjust
=yes
, type
=SURFACE TO SURFACE
flexible_interface, stiff_interface
*Material
, name
=
<material>
*Density
8.9e-09,
*Elastic
100000.,0.
*User Output Variables
1,
*Material
, name
=
<material>
*Density
1.92e-08,
*Elastic
400000.,0.
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
, U995, -9810
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
Body force - bimetallic_beam_bending_h0.inp outline
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
0.,        -2.5,          0.,          1.,        -2.5,          0.
*System
0.,         2.5,          0.,          1.,         2.5,          0.
*System
*Tie
, name
=
<name>
, adjust
=yes
, type
=SURFACE TO SURFACE
flexible_interface, stiff_interface
*Material
, name
=
<material>
*Density
8.9e-09,
*Elastic
100000.,0.
*User Output Variables
1,
*Material
, name
=
<material>
*Density
1.92e-08,
*Elastic
400000.,0.
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
, U998, -8.7309e-05
<element-set>
, U999, 0.
*Dload
<element-set>
, U997, 0.
<element-set>
, U998, -0.000188352
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
Pressure - bimetallic_beam_bending_h0.inp outline
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
8.9e-09,
*Elastic
100000.,0.
*User Output Variables
1,
*Material
, name
=
<material>
*Density
1.92e-08,
*Elastic
400000.,0.
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
Surface traction - bimetallic_beam_bending_h0.inp outline
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
8.9e-09,
*Elastic
100000.,0.
*User Output Variables
1,
*Material
, name
=
<material>
*Density
1.92e-08,
*Elastic
400000.,0.
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
End moment - bimetallic_beam_bending_h0.inp outline
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
1,           0.,           0.,         100.
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
8.9e-09,
*Elastic
100000.,0.
*User Output Variables
1,
*Material
, name
=
<material>
*Density
1.92e-08,
*Elastic
400000.,0.
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
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:47
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:44:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:12:47
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p1_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:07:43
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:39:44
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:43:28
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:11:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p2_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:07:00
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:35:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:42:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:55:46
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p3_body_force
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:49:22
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:45:06
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:13:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p1_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:10:03
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:39:47
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:43:46
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:12:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p2_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:02:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:35:19
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:42:22
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:59:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p3_gravity
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:50:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:44:43
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:12:33
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p1_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:04:11
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:38:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:44:11
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:12:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p2_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:59:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:22:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:39:43
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:53:54
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p3_moment
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:39:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p1_pressure
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p1_pressure
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:44:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p1_pressure
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:12:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p1_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:10:29
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p2_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:37:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p2_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:43:54
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p2_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:12:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p2_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:59:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p3_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:22:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p3_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:41:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p3_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:59:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p3_pressure
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:47:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:44:37
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:11:03
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p1_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:05:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p1_traction_tie
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p1_traction_tie
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:44:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p1_traction_tie
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:11:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p1_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:06:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:39:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:43:44
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:11:38
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p2_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:07:07
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p2_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:39:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p2_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:43:51
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p2_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:04:47
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p2_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:03:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:24:44
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:42:48
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:53:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p3_traction
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:43:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h0p3_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:32:23
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h1p3_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:42:09
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h2p3_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:57:52
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_bimetallic_beam_bending_immersed_h3p3_traction_tie
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:51:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_study_bimetallic_beam_bending_gravity_body_force_equivalence
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:13 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p1_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p1_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p1_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p1_pressure
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p1_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p1_traction_tie
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p2_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p2_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p2_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p2_pressure
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p2_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p2_traction_tie
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p3_body_force
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p3_gravity
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p3_moment
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p3_pressure
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p3_traction
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_bimetallic_beam_bending_immersed_p3_traction_tie
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
