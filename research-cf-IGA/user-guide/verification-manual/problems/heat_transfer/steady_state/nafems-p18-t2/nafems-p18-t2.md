# NAFEMS P18 T2: one-dimensional heat transfer with radiation – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/heat_transfer/steady_state/nafems-p18-t2/nafems-p18-t2.html

---

Verification Manual
Heat transfer
Steady-state heat transfer
NAFEMS P18 T2 — radiation
NAFEMS P18 T2: one-dimensional heat transfer with radiation
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify steady one-dimensional conduction coupled to nonlinear radiation using Test T2 of
The Standard NAFEMS Benchmarks
[1]
. The study compares a body-fitted Coreform/Abaqus IGA solution with the published target.
2
Geometry
The rod has length
\(L=0.1\ \text{m}\)
. Its square cross-section is a convenient three-dimensional representation; insulated lateral surfaces make the solution one-dimensional.
Figure 1: Geometry and thermal boundary conditions for NAFEMS P18 T2.
3
Material
The thermal conductivity is
\(k=55.6\ \text{W}/(\text{m K})\)
. The specified density and specific heat do not affect this steady solution.
4
Loading & boundary conditions
End A is held at
\(T_A=1000\ \text{K}\)
. End B radiates to surroundings at
\(T_\infty=300\ \text{K}\)
according to
\[
q=\epsilon\sigma\left(T_B^4-T_\infty^4\right),
\]
where
\(\epsilon=0.98\)
and
\(\sigma=56.7\times10^{-9}\ \text{W}/(\text{m}^2\text{K}^4)\)
.
5
Mesh & discretization
The Coreform/Abaqus model uses a body-fitted
\(10\times1\times1\)
mesh with linear,
\(C^0\)
splines.
Figure 2: Exact body-fitted mesh used by the Coreform T2 simulation.
6
Reference & accepted solutions
The accepted temperature at B is
\(927\ \text{K}\)
[1]
. The automated Coreform test enforces a 1% relative-error limit.
7
Results
7.1
Temperature at point B
Table 1: NAFEMS P18 T2 verification results.
Discretization
\(T_B\)
(K)
Relative error
Limit
Body-fitted Coreform/Abaqus IGA
927.00323
0.000%
1%
Figure 3: Coreform temperature field for NAFEMS P18 T2.
8
Discussion
The result agrees with the accepted endpoint temperature and verifies conduction and nonlinear radiation in the steady heat-transfer formulation.
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
t2–nafems_p18_thermal_model.jnl
8.2
Provenance
The record below identifies the Coreform/Abaqus verification evidence used by this page.
Note
Test_abaqus_verification_nafems_p18_t2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:40:33
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
