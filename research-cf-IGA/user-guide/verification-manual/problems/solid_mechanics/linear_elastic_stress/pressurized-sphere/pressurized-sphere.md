# Thick-walled spherical pressure vessel – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/linear_elastic_stress/pressurized-sphere/pressurized-sphere.html

---

Verification Manual
Solid mechanics
Linear elastic stress
Thick-walled spherical pressure vessel
Thick-walled spherical pressure vessel
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
\[
\newcommand{\peakVonMises}{5.54 \times 10^{-6}}
\newcommand{\innerHoop}{4.54 \times 10^{-6}}
\]
1
Objective
Verify Coreform IGA for Abaqus on a doubly curved three-dimensional problem whose solution is known exactly: a thick-walled sphere loaded by internal pressure. The Lamé closed-form solution gives the displacement and stress at every point in the structure, so the numerical solution is checked against exact values throughout the three-dimensional domain, and there is no separate reference solve to introduce its own approximation. The quantities of interest are the volume-integrated relative global
\(L^2\)
errors of the displacement vector, radial displacement, radial stress, and tangential stress.
2
Geometry
One octant of a hollow sphere with inner radius
\(a\)
and outer radius
\(b\)
is shown in
Figure 1
, with its dimensions provided in
Table 1
. The octant is oriented so the three cut faces lie on the coordinate planes. By spherical symmetry the octant with symmetry conditions represents the full sphere exactly.
Figure 1: Schematic of the thick-walled spherical shell octant and its coordinate-plane cut faces. The wall thickness is exaggerated for clarity.
Table 1: Geometric parameters.
Quantity
Symbol
Value
Inner radius
\(a\)
0.9
Outer radius
\(b\)
1
Wall thickness
\(b - a\)
0.1
3
Material
The sphere is a homogeneous, isotropic, linear-elastic material (
Table 2
). The model is nondimensional: the modulus and pressure set the scale of the displacement solution but cancel from the relative-error metrics.
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
4
Loading & boundary conditions
A uniform pressure
\(p_i\)
acts on the inner surface; the outer surface is traction-free (
Table 3
). The three cut faces carry symmetry conditions (zero normal displacement) on the
\(x\)
,
\(y\)
, and
\(z\)
coordinate planes, which restrain rigid-body motion and complete the octant model of the full sphere.
Table 3: Applied loading.
Quantity
Symbol
Value
Inner pressure
\(p_i\)
\(1.00 \times 10^{-6}\)
Outer pressure
\(p_o\)
0
5
Mesh & discretization
Coreform IGA discretizes the octant with
immersed, trimmed U-splines
: a structured rectilinear background grid is trimmed to the curved shell rather than body-fitting a mesh to it, so the spherical surfaces cut arbitrarily through the grid cells. The background element size is
\(h = 0.1/2^{k}\)
for refinement levels
\(k = 0, 1, 2\)
, at degrees
\(p = 1\)
,
\(2\)
, and
\(3\)
with maximum continuity
\(C^{p-1}\)
.
\(h\)
=0.1
\(h\)
=0.05
\(h\)
=0.025
6
Reference & accepted solutions
The
accepted
solution is the Lamé closed form for a thick-walled sphere under internal and external pressure
[1]
. With
\[
c = \frac{p_i a^3 - p_o b^3}{b^3 - a^3},
\qquad
\kappa = \frac{(p_i - p_o)\, a^3 b^3}{b^3 - a^3},
\]
the radial and tangential stresses and the radial displacement are
\[
\sigma_{rr}(r) = c - \frac{\kappa}{r^3},
\qquad
\sigma_{\theta\theta}(r) = c + \frac{\kappa}{2 r^3},
\qquad
u_r(r) = \frac{1}{E}\left[(1 - 2\nu)\, c\, r + (1 + \nu)\, \frac{\kappa}{2 r^2}\right].
\]
The analytical reference is exact. On the represented trimmed domain, error integrals evaluated with the solver quadrature are compared with the acceptance limits: each refinement level carries its own relative
\(L^2\)
tolerance per quantity, tightening as the mesh refines (see the tolerance columns in the results table). No separate numerical reference solve is needed.
6.1
Global error evaluation
The default-off
--quadrature-output
option records physical reference positions X, reference volume weights dV, directly evaluated U, and frame-matched material S in the separate
<job>_iga_quadrature.sql
file (format version 1). Visualization remains in
<job>_iga_results.sql
(schema 5); L2 postprocessing needs only the quadrature file. Both cut and ordinary cells contribute. The postprocessor uses these solver samples to compute the L2 errors in
global_l2_metrics.json
. With
\(\mathbf{n} = \mathbf{x}_q/\lVert\mathbf{x}_q\rVert\)
, it computes
\[
u_{r,h} = \mathbf{u}_h \cdot \mathbf{n},
\qquad
\sigma_{rr,h} = \mathbf{n}^{T}\boldsymbol{\sigma}_h\mathbf{n},
\qquad
\sigma_{\theta\theta,h}
= \frac{\operatorname{tr}(\boldsymbol{\sigma}_h)-\sigma_{rr,h}}{2}.
\]
For a scalar quantity
\(q\)
, the reported relative global error is
\[
\epsilon_q
=
\frac{\lVert q_h-q\rVert_{L^2(\Omega)}}{\lVert q\rVert_{L^2(\Omega)}}
\approx
\left[
\frac{\sum_Q w_q\,\det(\mathbf{J}_q)\,[q_h(\mathbf{x}_q)-q(\mathbf{x}_q)]^2}
     {\sum_Q w_q\,\det(\mathbf{J}_q)\,q(\mathbf{x}_q)^2}
\right]^{1/2}.
\]
The vector-displacement error uses the Euclidean squared norm of
\(\mathbf{u}_h-\mathbf{u}\)
in the numerator and of
\(\mathbf{u}\)
in the denominator. The errors are integrated over the represented trimmed octant using the physical reference-volume weights. The results identify this domain as
volume_quadrature_reference
. The numerical integration uses the solver’s volume quadrature and does not independently bound quadrature error.
7
Results
7.1
Quantities of interest
For each mesh, the numerical fields for vector displacement, radial displacement, radial stress, and tangential stress are compared against the Lamé closed form in the relative global
\(L^2\)
norm. The error for each refinement level is shown beside the acceptance tolerance it is required to meet.
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h\)
Vector displacement error
Tolerance
Radial displacement error
Tolerance
Radial stress error
Tolerance
Tangential stress error
Tolerance
0.1
\(0.584\%\)
\(0.8\%\)
\(0.568\%\)
\(0.8\%\)
\(28.5\%\)
\(40\%\)
\(1.89\%\)
\(2.5\%\)
0.1
\(0.421\%\)
\(0.6\%\)
\(0.409\%\)
\(0.6\%\)
\(30.4\%\)
\(45\%\)
\(1.98\%\)
\(2.8\%\)
0.1
\(0.418\%\)
\(0.6\%\)
\(0.405\%\)
\(0.6\%\)
\(30.4\%\)
\(45\%\)
\(1.98\%\)
\(2.8\%\)
0.05
\(0.199\%\)
\(0.3\%\)
\(0.194\%\)
\(0.3\%\)
\(15.5\%\)
\(20\%\)
\(0.93\%\)
\(1.2\%\)
0.025
\(0.0649\%\)
\(0.1\%\)
\(0.0632\%\)
\(0.1\%\)
\(7.78\%\)
\(10\%\)
\(0.437\%\)
\(0.6\%\)
\(h\)
Vector displacement error
Tolerance
Radial displacement error
Tolerance
Radial stress error
Tolerance
Tangential stress error
Tolerance
0.1
\(0.094\%\)
\(0.1\%\)
\(0.093\%\)
\(0.1\%\)
\(1.85\%\)
\(2.5\%\)
\(0.16\%\)
\(0.2\%\)
0.05
\(0.0173\%\)
\(0.025\%\)
\(0.0166\%\)
\(0.025\%\)
\(0.681\%\)
\(0.8\%\)
\(0.0445\%\)
\(0.06\%\)
0.025
\(0.00365\%\)
\(0.005\%\)
\(0.00312\%\)
\(0.005\%\)
\(0.211\%\)
\(0.3\%\)
\(0.0138\%\)
\(0.02\%\)
\(h\)
Vector displacement error
Tolerance
Radial displacement error
Tolerance
Radial stress error
Tolerance
Tangential stress error
Tolerance
0.1
\(0.0317\%\)
\(0.05\%\)
\(0.0287\%\)
\(0.04\%\)
\(0.837\%\)
\(1.2\%\)
\(0.0859\%\)
\(0.12\%\)
0.05
\(0.00734\%\)
\(0.01\%\)
\(0.0064\%\)
\(0.008\%\)
\(0.207\%\)
\(0.3\%\)
\(0.0175\%\)
\(0.025\%\)
0.025
\(0.00181\%\)
\(0.0025\%\)
\(404 \times 10^{-6}\%\)
\(500 \times 10^{-6}\%\)
\(0.0716\%\)
\(0.1\%\)
\(0.00588\%\)
\(0.008\%\)
Download data (CSV)
7.2
Convergence
(a) Vector displacement
\(\mathbf{u}\)
(b) Radial displacement
\(u_r\)
(c) Radial stress
\(\sigma_{rr}\)
(d) Tangential stress
\(\sigma_{\theta\theta}\)
Figure 2: Relative global
\(L^2\)
error against the Lamé closed form under mesh refinement, per degree. Errors use the solver volume quadrature on the represented trimmed domain.
7.3
Contour plots
h0p1
h0p1
h0p1
h1p1
h2p1
h0p2
h1p2
h2p2
h0p3
h1p3
h2p3
8
Performance
Performance is reported separately for trimming (building the trimmed U-spline basis on the immersing grid) and the Abaqus solve.
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.1\)
\(1{,}728\)
\(343\)
\(343\)
\(1\)
\(11.9\)
\(0.1\)
\(1{,}728\)
\(343\)
\(343\)
\(1\)
\(40.6\)
\(0.1\)
\(1{,}728\)
\(343\)
\(343\)
\(1\)
\(12.1\)
\(0.05\)
\(10{,}648\)
\(1{,}967\)
\(1{,}963\)
\(2\)
\(47.4\)
\(0.025\)
\(74{,}088\)
\(12{,}460\)
\(7{,}994\)
\(4\)
\(74.6\)
\(h\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.1\)
\(2{,}744\)
\(343\)
\(343\)
\(1\)
\(23.4\)
\(0.05\)
\(13{,}824\)
\(1{,}967\)
\(1{,}971\)
\(2\)
\(54.5\)
\(0.025\)
\(85{,}184\)
\(12{,}460\)
\(8{,}003\)
\(4\)
\(87.4\)
\(h\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(0.1\)
\(4{,}096\)
\(343\)
\(343\)
\(1\)
\(36.1\)
\(0.05\)
\(17{,}576\)
\(1{,}967\)
\(1{,}961\)
\(2\)
\(139.5\)
\(0.025\)
\(97{,}336\)
\(12{,}460\)
\(8{,}016\)
\(4\)
\(153.9\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.1\)
\(2{,}860\)
\(2{,}195\)
\(6{,}585\)
\(1\)
\(1\)
\(0.1\)
\(6{,}818\)
\(4{,}058\)
\(12{,}174\)
\(1\)
\(4\)
\(0.1\)
\(5{,}654\)
\(3{,}278\)
\(9{,}834\)
\(1\)
\(1\)
\(0.05\)
\(11{,}960\)
\(9{,}248\)
\(27{,}744\)
\(2\)
\(5\)
\(0.025\)
\(51{,}153\)
\(40{,}266\)
\(120{,}798\)
\(4\)
\(17\)
\(h\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.1\)
\(2{,}889\)
\(2{,}557\)
\(7{,}671\)
\(1\)
\(3\)
\(0.05\)
\(13{,}130\)
\(11{,}040\)
\(33{,}120\)
\(2\)
\(9\)
\(0.025\)
\(56{,}962\)
\(47{,}604\)
\(142{,}812\)
\(4\)
\(20\)
\(h\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(0.1\)
\(4{,}678\)
\(3{,}958\)
\(11{,}874\)
\(1\)
\(7\)
\(0.05\)
\(19{,}971\)
\(15{,}790\)
\(47{,}370\)
\(2\)
\(38\)
\(0.025\)
\(84{,}539\)
\(66{,}267\)
\(198{,}801\)
\(4\)
\(92\)
Download data (CSV)
9
Discussion
The
\(p=1\)
,
\(p=2\)
, and
\(p=3\)
cases converge monotonically toward the exact solution in all four global quantities. In this problem, the tangential stress is more accurate than the radial stress because it is dominated by the smoothly varying hoop response, while the radial stress must satisfy the traction boundary conditions at both surfaces, which cut arbitrarily through the immersed grid cells. However, these results do demonstrate that a curved geometry with a pressure-loaded surface that does not conform to the background grid can reproduce expected results. We note that the convergence rates do not meet the expected rates from finite element theory, but this is not surprising given the machinery required to integrate immersed IGA into Abaqus/Standard.
References
[1]
S. P. Timoshenko and J. N. Goodier,
Theory of elasticity
, 3rd ed. New York: McGraw-Hill, 1970.
Appendix
9.1
Capabilities exercised
The table is generated from current passing simulation and study artifacts. Each row must satisfy the linked capability card’s evidence contract; declaring a capability in a test does not, by itself, create evidence.
Table 4: Capabilities exercised by the current verification artifacts.
Capability
Evidence level
Linear static analysis
Numerically validated
Static stress analysis
Numerically validated
Cartesian symmetry boundary conditions (X, Y, Z)
Numerically validated
Boundary conditions and constraints
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
Numerically validated
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
11 generated Abaqus/CAE journals are available.
Show 11 individual journals
immersed_h0p1–pressurized_sphere_h0_model.jnl
immersed_h0p2–pressurized_sphere_h0_model.jnl
immersed_h0p3–pressurized_sphere_h0_model.jnl
immersed_h1p1–pressurized_sphere_h1_model.jnl
immersed_h1p2–pressurized_sphere_h1_model.jnl
immersed_h1p3–pressurized_sphere_h1_model.jnl
immersed_h2p1–pressurized_sphere_h2_model.jnl
immersed_h2p2–pressurized_sphere_h2_model.jnl
immersed_h2p3–pressurized_sphere_h2_model.jnl
immersed_regproj_h0p1–pressurized_sphere_h0_model.jnl
immersed_wt2_h0p1–pressurized_sphere_h0_model.jnl
9.3
Supported Abaqus keywords
The deck below is reproduced from the verification test’s Abaqus input file with comments stripped and instance-specific names (sets, materials, files) redacted to
<placeholders>
.
Note
pressurized_sphere_h0.inp - outline
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
, XSYMM
*Boundary
<node-set>
, YSYMM
*Boundary
<node-set>
, ZSYMM
*Step
, name
=
<step>
, nlgeom
=NO
, inc
=1
*Static
, direct
1., 1.,
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
Test_abaqus_verification_pressurized_sphere_immersed_h0p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:37:30
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_regproj_h0p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:17:39
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_wt2_h0p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:38:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h1p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:22:50
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h2p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:21:49
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h0p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:24:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h1p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:24:24
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h2p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:20:09
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h0p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:20:23
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h1p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:24:58
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_pressurized_sphere_immersed_h2p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:18:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_study_pressurized_sphere_immersed_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:41:04 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_pressurized_sphere_immersed_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:41:04 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_pressurized_sphere_immersed_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:41:03 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
