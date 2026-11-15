# NAFEMS P18 T1: membrane with a circular hot spot – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/thermal_strain/nafems-p18-t1/nafems-p18-t1.html

---

Verification Manual
Solid mechanics
Thermal strain
NAFEMS P18 T1 — membrane with circular hot spot
NAFEMS P18 T1: membrane with a circular hot spot
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify the stress produced by the discontinuous thermal strain in Test T1 of
The Standard NAFEMS Benchmarks
[1]
. The automated study compares an immersed Coreform/Abaqus IGA discretization with the published solution.
2
Geometry
Two symmetry planes reduce the
\(20\ \text{mm} \times 20\ \text{mm}\)
membrane to the
\(10\ \text{mm} \times 10\ \text{mm} \times 0.25\ \text{mm}\)
quarter model in
Figure 1
. The quarter-circle of radius
\(1\ \text{mm}\)
is a second material region, not a hole. Point D lies at
\((x,y)=(1,0)\ \text{mm}\)
on the interface.
Figure 1: Quarter-symmetry model for NAFEMS P18 T1.
3
Material
Both regions are isotropic and linear elastic with
\(E=100{,}000\ \text{MPa}\)
and
\(\nu=0.3\)
. The circular region receives a uniform thermal strain of
\(1\times10^{-3}\)
, while the surrounding membrane has zero thermal strain.
4
Loading & boundary conditions
Symmetry conditions impose
\(u_x=0\)
on
\(x=0\)
and
\(u_y=0\)
on
\(y=0\)
. The two cells belong to one multi-cell part and use different Abaqus material sections, so the thermal expansion changes discontinuously while displacement remains continuous. The reported
\(\sigma_{yy}\)
is sampled at D on the surrounding-material side.
5
Mesh & discretization
The immersed model uses a linear,
\(C^0\)
spline space constructed from a uniform
\(44\times44\times1\)
background grid. The in-plane cell size is
\(0.227\ \text{mm}\)
, and no local or hierarchical refinement is used. The three-dimensional formulation constrains
\(u_z\)
and applies the constitutive transformation required to reproduce plane stress.
Figure 2: Exact background grid used for both immersed material regions.
6
Reference & accepted solutions
The accepted value is
\(\sigma_{yy}(D)=50.0\ \text{MPa}\)
[1]
. The Coreform IGA test enforces a 4% relative-error limit.
7
Results
7.1
Stress at point D
Table 1: NAFEMS P18 T1 verification results.
Discretization
\(\sigma_{yy}(D)\)
(MPa)
Relative error
Limit
Immersed Coreform/Abaqus IGA
48.503
2.99%
4%
Figure 3: Immersed Coreform/Abaqus IGA stress field near the circular interface.
8
Discussion
The immersed result satisfies the accepted value and verifies displacement continuity across the multi-cell material interface together with the discontinuous thermal-strain assignment.
References
[1]
The standard NAFEMS benchmarks
. NAFEMS.
Appendix
8.1
Download artifacts
Download all verification artifacts (ZIP)
8.1.1
Geometry CAD files
Download the immersed Coreform/Abaqus model
8.1.2
Abaqus/CAE journal files
1 generated Abaqus/CAE journal is available.
Show 1 individual journal
t1–nafems_p18_thermal_model.jnl
8.2
Provenance
The record below identifies the Coreform/Abaqus verification evidence used by this page.
Note
Test_abaqus_verification_nafems_p18_t1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:21:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
