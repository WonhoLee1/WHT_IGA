# NAFEMS P18 T3: one-dimensional transient heat transfer – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/heat_transfer/transient/nafems-p18-t3/nafems-p18-t3.html

---

Verification Manual
Heat transfer
Transient heat transfer
NAFEMS P18 T3 — transient conduction
NAFEMS P18 T3: one-dimensional transient heat transfer
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify transient thermal capacity and implicit time integration using Test T3 of
The Standard NAFEMS Benchmarks
[1]
. This verifies the Coreform IGA thermal UEL through the Abaqus heat-transfer procedure.
2
Geometry
The rod has length
\(L=0.1\ \text{m}\)
, and point C is located at
\(x=0.08\ \text{m}\)
from end A. Insulated lateral surfaces make the response one-dimensional.
Figure 1: Geometry and time-dependent boundary condition for NAFEMS P18 T3.
3
Material
The material properties are
\(k=35\ \text{W}/(\text{m K})\)
,
\(\rho=7200\ \text{kg}/\text{m}^3\)
, and
\(c=440.5\ \text{J}/(\text{kg K})\)
.
4
Loading & boundary conditions
The initial temperature is
\(0\ ^\circ\text{C}\)
. End A remains at
\(0\ ^\circ\text{C}\)
, and end B follows
\[
T_B(t)=100\sin\left(\frac{\pi t}{40}\right)\ ^\circ\text{C}.
\]
All lateral surfaces are insulated.
5
Mesh & discretization
The Coreform/Abaqus model uses the same body-fitted
\(10\times1\times1\)
linear,
\(C^0\)
mesh as T2. The Abaqus transient heat-transfer step starts with a
\(0.1\ \text{s}\)
increment, permits increments up to
\(0.5\ \text{s}\)
, and ends at
\(t=32\ \text{s}\)
.
Figure 2: Exact body-fitted mesh used by the Coreform T3 simulation.
6
Reference & accepted solutions
The accepted temperature at C and
\(t=32\ \text{s}\)
is
\(36.60\ ^\circ\text{C}\)
[1]
. The automated test enforces a 2% relative-error limit.
7
Results
7.1
Temperature at point C
Table 1: NAFEMS P18 T3 verification result at 32 s.
Discretization
\(T_C\)
(
\(^\circ\)
C)
Relative error
Limit
Body-fitted Coreform/Abaqus IGA
37.25446
1.79%
2%
Figure 3: Coreform temperature field at 32 s for NAFEMS P18 T3.
8
Discussion
The result satisfies the published benchmark and directly exercises the volumetric capacity term and Abaqus transient thermal integration.
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
t3–nafems_p18_thermal_model.jnl
8.2
Provenance
The record below identifies the Coreform/Abaqus verification evidence used by this page.
Note
Test_abaqus_verification_nafems_p18_t3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:40:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
