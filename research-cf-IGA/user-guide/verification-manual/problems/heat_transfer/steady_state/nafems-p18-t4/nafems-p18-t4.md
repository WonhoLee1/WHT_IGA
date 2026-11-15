# NAFEMS P18 T4: two-dimensional heat transfer with convection – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/heat_transfer/steady_state/nafems-p18-t4/nafems-p18-t4.html

---

Verification Manual
Heat transfer
Steady-state heat transfer
NAFEMS P18 T4 — convection
NAFEMS P18 T4: two-dimensional heat transfer with convection
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify two-dimensional steady conduction with convection boundaries using Test T4 of
The Standard NAFEMS Benchmarks
[1]
. The study compares a body-fitted Coreform/Abaqus IGA solution with the published target.
2
Geometry
The domain is
\(0.6\ \text{m}\)
wide,
\(1.0\ \text{m}\)
high, and
\(0.05\ \text{m}\)
thick. Point E is located at
\((x,y)=(0.6,0.2)\ \text{m}\)
. Insulated front and back faces preserve the two-dimensional response.
Figure 1: Geometry and thermal boundary conditions for NAFEMS P18 T4.
3
Material
The isotropic thermal conductivity is
\(k=52\ \text{W}/(\text{m K})\)
.
4
Loading & boundary conditions
The bottom edge AB is held at
\(100\ ^\circ\text{C}\)
, and the left edge AD is insulated. The top and right edges exchange heat with a
\(0\ ^\circ\text{C}\)
environment using
\(h=750\ \text{W}/(\text{m}^2\text{K})\)
.
5
Mesh & discretization
The Coreform/Abaqus model uses a body-fitted
\(12\times20\times1\)
mesh with linear,
\(C^0\)
splines.
Figure 2: Exact body-fitted mesh used by the Coreform T4 simulation.
6
Reference & accepted solutions
The accepted temperature at E is
\(18.3\ ^\circ\text{C}\)
[1]
. The automated Coreform/Abaqus test enforces a 3% relative-error limit.
7
Results
7.1
Temperature at point E
Table 1: NAFEMS P18 T4 verification results.
Discretization
\(T_E\)
(
\(^\circ\)
C)
Relative error
Limit
Body-fitted Coreform/Abaqus IGA
18.08995
1.148%
3%
Figure 3: Coreform temperature field for NAFEMS P18 T4.
8
Discussion
The result satisfies the published target and verifies steady two-dimensional conduction and convection on a body-fitted spline discretization.
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
Download the body-fitted Coreform/Abaqus model
8.1.2
Abaqus/CAE journal files
1 generated Abaqus/CAE journal is available.
Show 1 individual journal
t4–nafems_p18_thermal_model.jnl
8.2
Provenance
The record below identifies the Coreform/Abaqus verification evidence used by this page.
Note
Test_abaqus_verification_nafems_p18_t4
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:40:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
