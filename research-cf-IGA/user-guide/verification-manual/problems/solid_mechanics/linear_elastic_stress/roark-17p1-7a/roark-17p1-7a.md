# Stress concentration in a thin plate with hole – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/linear_elastic_stress/roark-17p1-7a/roark-17p1-7a.html

---

Verification Manual
Solid mechanics
Linear elastic stress
Stress concentration in a thin plate with hole
Stress concentration in a thin plate with hole
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
\[
\newcommand{\peakRoark}{4.32 \times 10^{-6}}
\newcommand{\peakCifa}{4.34 \times 10^{-6}}
\]
1
Objective
Verify the ability of Coreform IGA for Abaqus to predict stress concentrations on a classic linear-elastic stress-concentration problem: a uniaxially-loaded coupon with a central circular hole. The quantity of interest is the peak stress at the hole edge, whose
accepted
value is the stress-concentration factor tabulated by Roark
[1]
. We also compare the transverse stress and displacement profiles along the line probe against a reference body-fitted Abaqus solution. Finally, we verify
local refinement
: a coarse background refined around the hole to increasing levels, compared with uniform refinement at equal accuracy.
2
Geometry
A thin rectangular coupon with a central circular hole as shown in (
Figure 1
) and parameters provided in
Table 1
. Verification quantities are sampled along the shortest span from the hole edge to the free edge.
Figure 1: Coupon geometry: a plate of length
\(L\)
, width
\(D\)
, thickness
\(t\)
, with a central hole of radius
\(r\)
. The dashed red line indicates the location of a line probe that samples our quantities of interest. The annotated point
\(A\)
indicates the location of expected maximum principal stress, used to evaluate the peak stress in our analysis.
Table 1: Geometric parameters.
Quantity
Symbol
Value
Coupon length
\(L\)
\(200\)
Coupon width
\(D\)
\(40\)
Hole radius
\(r\)
\(10\)
Thickness
\(t\)
\(1\)
3
Material
The coupon is a homogeneous, isotropic, linear-elastic material with properties given in (
Table 2
).
Table 2: Material properties.
Property
Symbol
Value
Young’s modulus
\(E\)
1
Poisson’s ratio
\(\nu\)
0
Mass density
\(\rho\)
0
4
Loading & boundary conditions
The coupon is loaded in uniaxial tension by a uniform pressure on the right face (
Table 3
), with the left face pinned as depicted in
Figure 2
. The load is chosen to be small, compared to the Young’s modulus, so that the response is assumed linear.
(a) Idealized boundary conditions.
(b) Boundary conditions as implemented in Abaqus and Coreform IGA.
Figure 2: Boundary conditions on the coupon: (a) the idealized problem statement with uniform-pressure applied on both ends, and (b) the equivalent conditions as actually implemented in our finite element models.
Table 3: Applied loading.
Quantity
Symbol
Value
Applied end pressure
\(P\)
\(1.00 \times 10^{-6}\)
5
Mesh & discretization
Coreform IGA discretizes the coupon with
immersed, trimmed U-splines
: a structured background spline is trimmed to the coupon-with-hole geometry rather than body-fitting a mesh to it. We run an
\(h/p\)
refinement study — background element size
\(h = 6/2^{k}\)
for
\(k = 0,\dots,4\)
at polynomial degrees
\(p = 1\)
,
\(2\)
, and
\(3\)
. A separate local-refinement sweep (
Section 7.6
) keeps the
\(h = 3\)
,
\(p = 2\)
background and refines only around the hole. Per-mesh element, node, and degree-of-freedom counts are reported in
Section 8
.
The tabs below show the trimmed immersed mesh at each background element size
\(h\)
(the mesh is shared by all degrees
\(p\)
), rendered from the run’s IGA results database: trimmed Bezier element edges (carnation) on the CAD surface, CAD edges in black. The
Reference
tab shows the body-fitted Abaqus mesh, rendered from the VTU the reference test exports. The first group views the full coupon; the second is a detail of the nominal section, from the top of the hole to the top of the coupon.
Full coupon
Reference
h=6
h=3
h=1.5
h=0.75
h=0.375
Nominal-section detail
Reference
h=6
h=3
h=1.5
h=0.75
h=0.375
6
Reference & accepted solutions
Roark’s
Formulas for Stress and Strain
, Table 17.1 case 7a, gives the
accepted
stress-concentration factor for the net section,
\[
K_t = \sum_{i=0}^{3} C_i \left(\tfrac{2r}{D}\right)^i ,
\]
with the polynomial coefficients
\(C_i\)
listed in
Table 4
. The nominal stress on the net section is
\[
\sigma_{\text{nom}} = \frac{P\,t D}{t\,(D - 2r)} ,
\]
so the accepted peak max-principal stress at the hole edge is
\[
\sigma_{\max} = K_t\,\sigma_{\text{nom}} .
\]
Its evaluated value is reported alongside the reference solve in
Table 5
.
Table 4: Roark
\(K_t\)
polynomial coefficients.
Coefficient
Value
\(C_0\)
\(3.00\)
\(C_1\)
\(-3.13\)
\(C_2\)
\(3.66\)
\(C_3\)
\(-1.53\)
The
reference
solution is a normal Abaqus analysis of the same coupon using a dense, body-fitted mesh of quadratic elements. Its transverse stress and displacement profiles along the line probe are the comparison for the field quantities, and its peak max-principal stress at the hole edge provides an independent numerical estimate to sit alongside Roark’s accepted value (
Table 5
). Being itself a numerical approximation rather than the exact solution, the reference is not treated as ground truth — a distinction that matters for the convergence study (
Section 7.5
).
Table 5: Accepted and reference peak max-principal stress at the hole edge.
Source
peak
\(\sigma_{\max}\)
Accepted (Roark
\(K_t\)
)
\(4.32 \times 10^{-6}\)
Reference (body-fitted Abaqus)
\(4.35 \times 10^{-6}\)
Download data (CSV)
7
Results
7.1
Quantities of interest
For each mesh the peak stress is compared against both the accepted Roark value and the reference solve; the transverse stress and displacement profiles are compared (relative
\(L^2\)
) against the reference.
h=6, p=1
h=3, p=1
h=1.5, p=1
h=0.75, p=1
h=0.375, p=1
h=6, p=2
h=3, p=2
h=1.5, p=2
h=0.75, p=2
h=0.375, p=2
h=6, p=3
h=3, p=3
h=1.5, p=3
h=0.75, p=3
h=0.375, p=3
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(3.34 \times 10^{-6}\)
\(22.7\%\)
\(23.2\%\)
transverse stress rel.
\(L^2\)
—
—
\(11.5\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(2.15\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(35.1\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.19 \times 10^{-6}\)
\(2.98\%\)
\(3.66\%\)
transverse stress rel.
\(L^2\)
—
—
\(3.14\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.575\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(9.67\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.25 \times 10^{-6}\)
\(1.57\%\)
\(2.26\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.844\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.171\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(3.02\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.32 \times 10^{-6}\)
\(0.0822\%\)
\(0.619\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.286\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0425\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.688\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.34 \times 10^{-6}\)
\(0.542\%\)
\(0.163\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.0976\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0162\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.195\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(3.26 \times 10^{-6}\)
\(24.6\%\)
\(25.1\%\)
transverse stress rel.
\(L^2\)
—
—
\(14.1\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.599\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(9.42\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.07 \times 10^{-6}\)
\(5.66\%\)
\(6.32\%\)
transverse stress rel.
\(L^2\)
—
—
\(2.98\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0417\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.529\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.31 \times 10^{-6}\)
\(0.117\%\)
\(0.817\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.548\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0139\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.0603\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.33 \times 10^{-6}\)
\(0.228\%\)
\(0.474\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.216\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0114\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.0287\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.35 \times 10^{-6}\)
\(0.826\%\)
\(0.12\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.0984\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0109\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.0163\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.00 \times 10^{-6}\)
\(7.38\%\)
\(8.03\%\)
transverse stress rel.
\(L^2\)
—
—
\(3.32\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.143\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(3.72\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.40 \times 10^{-6}\)
\(1.82\%\)
\(1.11\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.537\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0108\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.192\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.36 \times 10^{-6}\)
\(0.903\%\)
\(0.196\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.0763\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0128\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.0547\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.35 \times 10^{-6}\)
\(0.738\%\)
\(0.032\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.0673\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0111\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.0238\%\)
Quantity
Coreform IGA
vs Roark
vs reference
peak
\(\sigma_{\max}\)
\(4.34 \times 10^{-6}\)
\(0.618\%\)
\(0.0873\%\)
transverse stress rel.
\(L^2\)
—
—
\(0.0742\%\)
\(u_x\)
profile rel.
\(L^2\)
—
—
\(0.0109\%\)
\(u_y\)
profile rel.
\(L^2\)
—
—
\(0.0217\%\)
Download data (CSV)
7.2
Stress field
Figure 3
compares the maximum-principal-stress field between the immersed IGA solution on the finest
\(p = 2\)
mesh (
Figure 3 (a)
) and the body-fitted Abaqus reference solve (
Figure 3 (b)
). Both are rendered in ParaView on a shared color scale — the Abaqus fields come from the VTU that the reference test case exports from its ODB — so the contours compare directly.
(a) Immersed IGA solution (case
immersed_h4p2
,
\(h = 0.375\)
,
\(p = 2\)
).
(b) Body-fitted Abaqus reference solution.
Figure 3: Maximum principal stress on a shared color scale.
The tabs below extend the comparison to every case in the study, on the same shared color scale — flip between the
Reference
tab and any IGA case to compare contours directly.
Full coupon
Reference
h=6, p=1
h=3, p=1
h=1.5, p=1
h=0.75, p=1
h=0.375, p=1
h=6, p=2
h=3, p=2
h=1.5, p=2
h=0.75, p=2
h=0.375, p=2
h=6, p=3
h=3, p=3
h=1.5, p=3
h=0.75, p=3
h=0.375, p=3
Nominal-section detail
Reference
h=6, p=1
h=3, p=1
h=1.5, p=1
h=0.75, p=1
h=0.375, p=1
h=6, p=2
h=3, p=2
h=1.5, p=2
h=0.75, p=2
h=0.375, p=2
h=6, p=3
h=3, p=3
h=1.5, p=3
h=0.75, p=3
h=0.375, p=3
7.3
Transverse stress profile
Figure 4
shows the transverse max-principal-stress profile along the line probe: the reference Abaqus solve against the coarse→fine IGA sweep, per degree. The IGA curves settle onto a consistent profile as
\(h\)
decreases.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 4: Transverse max-principal stress: reference vs the coarse→fine IGA sweep.
7.4
Displacement profiles
The transverse displacement components show the same
\(h\)
-refinement behavior —
Figure 5
for
\(u_x\)
and
Figure 6
for
\(u_y\)
.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 5: Transverse
\(u_x\)
: reference vs the coarse→fine IGA sweep.
(a)
\(p = 1\)
(b)
\(p = 2\)
(c)
\(p = 3\)
Figure 6: Transverse
\(u_y\)
: reference vs the coarse→fine IGA sweep.
7.5
Mesh convergence
Because the reference solve is itself an approximation — the IGA discretization may converge to a slightly different (and possibly more accurate) solution — we assess convergence against the
finest IGA solution
(self-convergence) rather than the reference.
Figure 7
shows the transverse-stress
\(L^2\)
error toward that finest solution, with fitted rates of 1.86 (
\(p=1\)
), 2.23 (
\(p=2\)
), and 2.38 (
\(p=3\)
).
Figure 7: Self-convergence of the transverse-stress profile toward the finest IGA solution (not the reference solve), per degree.
7.6
Local refinement
The local-refinement sweep keeps the
\(p = 2\)
background at
\(h = 3\)
and refines the cells near the hole’s cylindrical face, named by the
hole
set, to level
\(L = 1\)
,
\(2\)
, and
\(3\)
. Each case is a single refinement at the deepest level with two-to-one balancing, so the intermediate levels grade in automatically. The refinement offset is left unset, so it defaults to the degree times the element size at that level,
\(p\,h/2^L\)
. The finest refined cells match the uniform
\(h = 1.5\)
,
\(0.75\)
, and
\(0.375\)
meshes.
The tabs show the maximum principal stress near the hole at each level, with the refined elements drawn over it, on the same color scale as
Figure 3
.
L=1
L=2
L=3
Figure 8
plots the peak-stress error against the reference solve as a function of degrees of freedom, for the uniform
\(p = 2\)
sweep and the local-refinement sweep. At
\(L\)
= 3 the refined mesh recovers the reference peak to 0.11% with 57,300 degrees of freedom; the uniform mesh with the same finest cell (
\(h\)
= 0.375) reaches 0.12% with 516,858. The transverse stress profile is a different matter: the line probe runs from the hole edge across the ligament, mostly outside the refined band, so its accuracy is bounded by the
\(h = 3\)
background rather than by the refinement level.
Figure 8: Peak max-principal stress error against the reference solve, by degrees of freedom: uniform
\(p = 2\)
refinement and local refinement around the hole.
The tabs below pair each level with the uniform mesh that has the same finest cell. Peak errors are against the reference solve; wall times are at the processor count each case runs with.
Accuracy
Cost
\(L\)
Finest
\(h\)
DOFs (local)
DOFs (uniform)
Peak error (local)
Peak error (uniform)
\(1\)
\(1.5\)
\(14{,}280\)
\(37{,}146\)
\(1.72\%\)
\(0.817\%\)
\(2\)
\(0.75\)
\(24{,}264\)
\(138{,}207\)
\(0.61\%\)
\(0.474\%\)
\(3\)
\(0.375\)
\(57{,}300\)
\(516{,}858\)
\(0.105\%\)
\(0.12\%\)
\(L\)
CPUs
Trim [s] (local)
Trim [s] (uniform)
Solve [s] (local)
Solve [s] (uniform)
\(1\)
\(4\)
\(24.9\)
\(24.9\)
\(2\)
\(4\)
\(2\)
\(8\)
\(53.2\)
\(58.0\)
\(3\)
\(9\)
\(3\)
\(16\)
\(236.5\)
\(150.0\)
\(8\)
\(26\)
Download data (CSV)
8
Performance
Performance is reported separately for the two stages of an immersed IGA analysis:
trimming
(building the trimmed U-spline basis on the immersing grid) and the
Abaqus solve
.
8.1
Coreform IGA Mesh
For each case the immersing background grid, the active Bezier elements kept in the analysis basis, and the cells actually cut by the part boundary are read from the trim database; the counts are nested subsets (background
\(\supseteq\)
active
\(\supseteq\)
trimmed), so they decrease left to right.
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h\)
Total
Active
Trimmed
CPUs
Wall-time [s]
\(6\)
\(2{,}145\)
\(240\)
\(96\)
\(1\)
\(7.0\)
\(3\)
\(6{,}745\)
\(984\)
\(199\)
\(2\)
\(12.1\)
\(1.5\)
\(21{,}545\)
\(3{,}532\)
\(402\)
\(4\)
\(23.1\)
\(0.75\)
\(79{,}945\)
\(14{,}180\)
\(794\)
\(8\)
\(54.4\)
\(0.375\)
\(299{,}145\)
\(55{,}108\)
\(1{,}557\)
\(16\)
\(171.0\)
\(h\)
Total
Active
Trimmed
CPUs
Wall-time [s]
\(6\)
\(2{,}145\)
\(240\)
\(96\)
\(1\)
\(8.4\)
\(3\)
\(6{,}745\)
\(984\)
\(199\)
\(2\)
\(14.0\)
\(1.5\)
\(21{,}545\)
\(3{,}532\)
\(402\)
\(4\)
\(24.9\)
\(0.75\)
\(79{,}945\)
\(14{,}180\)
\(794\)
\(8\)
\(58.0\)
\(0.375\)
\(299{,}145\)
\(55{,}108\)
\(1{,}557\)
\(16\)
\(150.0\)
\(h\)
Total
Active
Trimmed
CPUs
Wall-time [s]
\(6\)
\(2{,}145\)
\(240\)
\(96\)
\(1\)
\(27.2\)
\(3\)
\(6{,}745\)
\(984\)
\(199\)
\(2\)
\(33.1\)
\(1.5\)
\(21{,}545\)
\(3{,}532\)
\(402\)
\(4\)
\(45.5\)
\(0.75\)
\(79{,}945\)
\(14{,}180\)
\(794\)
\(8\)
\(79.3\)
\(0.375\)
\(299{,}145\)
\(55{,}108\)
\(1{,}557\)
\(16\)
\(184.1\)
Download data (CSV)
8.2
Coreform IGA for Abaqus solve
Problem size and solve-stage timing of the Abaqus job: the wall time is the Abaqus/Standard (solve) stage only — the pre and cleanup stages are excluded.
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(6\)
\(458\)
\(770\)
\(2{,}310\)
\(1\)
\(1\)
\(3\)
\(1{,}354\)
\(2{,}466\)
\(7{,}398\)
\(2\)
\(1\)
\(1.5\)
\(4{,}326\)
\(8{,}106\)
\(24{,}318\)
\(4\)
\(2\)
\(0.75\)
\(15{,}772\)
\(30{,}413\)
\(91{,}239\)
\(8\)
\(5\)
\(0.375\)
\(58{,}126\)
\(114{,}226\)
\(342{,}678\)
\(16\)
\(18\)
\(h\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(6\)
\(458\)
\(1{,}202\)
\(3{,}606\)
\(1\)
\(2\)
\(3\)
\(1{,}354\)
\(3{,}814\)
\(11{,}442\)
\(2\)
\(2\)
\(1.5\)
\(4{,}326\)
\(12{,}382\)
\(37{,}146\)
\(4\)
\(4\)
\(0.75\)
\(15{,}772\)
\(46{,}069\)
\(138{,}207\)
\(8\)
\(9\)
\(0.375\)
\(58{,}126\)
\(172{,}286\)
\(516{,}858\)
\(16\)
\(26\)
\(h\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(6\)
\(560\)
\(1{,}804\)
\(5{,}412\)
\(1\)
\(4\)
\(3\)
\(1{,}536\)
\(5{,}476\)
\(16{,}428\)
\(2\)
\(8\)
\(1.5\)
\(5{,}000\)
\(17{,}446\)
\(52{,}338\)
\(4\)
\(16\)
\(0.75\)
\(17{,}454\)
\(63{,}469\)
\(190{,}407\)
\(8\)
\(40\)
\(0.375\)
\(59{,}676\)
\(232{,}932\)
\(698{,}796\)
\(16\)
\(125\)
Download data (CSV)
8.3
Reference Abaqus solve
The body-fitted reference solve, for scale alongside the immersed cases.
\(h\)
Elements (C3D20)
Nodes
DOFs
Threads
Wall-time [s]
\(0.375\)
\(202{,}452\)
\(1{,}020{,}224\)
\(3{,}060{,}672\)
\(16\)
\(115\)
Download data (CSV)
9
Discussion
At the finest mesh Coreform IGA recovers the accepted peak stress to 0.62% of Roark’s
\(K_t\)
and the reference peak to 0.09%, while the transverse stress and displacement profiles match the reference to better than
\(10^{-3}\)
in relative
\(L^2\)
. Under uniform refinement the stress profile self-converges at rates of 1.86 (
\(p=1\)
), 2.23 (
\(p=2\)
), and 2.38 (
\(p=3\)
); the coarsest immersed meshes are pre-asymptotic, where raising the polynomial degree does not improve accuracy until the background mesh resolves the stress concentration.
Local refinement around the hole reaches the reference peak to 0.11% with 57,300 degrees of freedom, against 516,858 for the uniform mesh with the same finest cell (
Section 7.6
).
Because the reference is itself an approximate FEM solution, convergence is assessed against the finest IGA solution rather than the reference. Agreement between Coreform IGA and the independent body-fitted reference to within a fraction of a percent — with both agreeing with Roark’s analytic
\(K_t\)
— verifies the Coreform IGA immersed U-spline discretization for this stress-concentration problem.
References
[1]
W. C. Young, R. G. Budynas, and A. M. Sadegh,
Roark’s formulas for stress and strain
, 8th ed. New York: McGraw-Hill, 2011.
Appendix
9.1
Capabilities exercised
The table is generated from current passing simulation artifacts. Each row must satisfy the linked capability card’s evidence contract; declaring a capability in a test does not, by itself, create evidence.
Table 6: Capabilities exercised by the current verification artifacts.
Capability
Evidence level
Linear static analysis
Numerically validated
Static stress analysis
Numerically validated
Threaded Abaqus/Standard execution
Numerically validated
Parallel Coreform IGA mesh generation
Numerically validated
Parallel execution
Numerically validated
Pressure loading
Numerically validated
Uniform pressure loading
Numerically validated
Mechanical loads
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
Local refinement of immersed rectilinear meshes
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
Point probes
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
19 generated Abaqus/CAE journals are available.
Show 19 individual journals
immersed_h0p1–roark_17p1_7a_immersed_h0p1_model.jnl
immersed_h0p2–roark_17p1_7a_immersed_h0p2_model.jnl
immersed_h0p3–roark_17p1_7a_immersed_h0p3_model.jnl
immersed_h1p1–roark_17p1_7a_immersed_h1p1_model.jnl
immersed_h1p2–roark_17p1_7a_immersed_h1p2_model.jnl
immersed_h1p2_refined_surface_L1–roark_17p1_7a_immersed_h1p2_refined_surface_L1_model.jnl
immersed_h1p2_refined_surface_L2–roark_17p1_7a_immersed_h1p2_refined_surface_L2_model.jnl
immersed_h1p2_refined_surface_L3–roark_17p1_7a_immersed_h1p2_refined_surface_L3_model.jnl
immersed_h1p3–roark_17p1_7a_immersed_h1p3_model.jnl
immersed_h2p1–roark_17p1_7a_immersed_h2p1_model.jnl
immersed_h2p2–roark_17p1_7a_immersed_h2p2_model.jnl
immersed_h2p3–roark_17p1_7a_immersed_h2p3_model.jnl
immersed_h3p1–roark_17p1_7a_immersed_h3p1_model.jnl
immersed_h3p2–roark_17p1_7a_immersed_h3p2_model.jnl
immersed_h3p3–roark_17p1_7a_immersed_h3p3_model.jnl
immersed_h4p1–roark_17p1_7a_immersed_h4p1_model.jnl
immersed_h4p2–roark_17p1_7a_immersed_h4p2_model.jnl
immersed_h4p3–roark_17p1_7a_immersed_h4p3_model.jnl
reference–roark_17p1_7a_reference_model.jnl
9.3
Supported Abaqus keywords
The deck below is reproduced from the verification tests’ Abaqus input files with comments stripped and instance-specific names (sets, materials, files) redacted to
<placeholders>
.
Note
roark_17p1_7a.inp - Outline
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
*Elastic
1.,0.
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
9.4
Provenance
The records below identify the execution environment and verification basis for every CTest artifact reported on this page.
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 4 record-only checks (0 positive)
Test date
04-Nov-2026 19:40:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 4 record-only checks (1 positive)
Test date
04-Nov-2026 19:40:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 4 record-only checks (1 positive)
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
Test_abaqus_verification_roark_17p1_7a_immersed_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 4 record-only checks (2 positive)
Test date
04-Nov-2026 19:46:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 4 record-only checks (2 positive)
Test date
04-Nov-2026 19:46:33
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h1p2_refined_surface_L1
Parameter
Value
Evidence type
Numerical verification with record-only diagnostics
CTest result
Passed (Release)
Verification result
Passed — 5/5 enforced criteria passed; 1 record-only diagnostic
Test date
04-Nov-2026 20:14:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h1p2_refined_surface_L2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 6/6 enforced criteria passed
Test date
04-Nov-2026 20:22:22
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h1p2_refined_surface_L3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 6/6 enforced criteria passed
Test date
04-Nov-2026 20:56:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h1p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:45:29
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 4 record-only checks (3 positive)
Test date
04-Nov-2026 20:13:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h2p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:13:07
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h2p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:13:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h3p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:14:45
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h3p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
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
Test_abaqus_verification_roark_17p1_7a_immersed_h3p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:09:51
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h4p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:50:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h4p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:40:00
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_immersed_h4p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 20:01:50
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_roark_17p1_7a_reference
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:41:22
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
