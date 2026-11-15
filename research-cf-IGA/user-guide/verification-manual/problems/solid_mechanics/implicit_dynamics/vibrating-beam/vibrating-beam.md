# Implicit-dynamic vibration of a three-dimensional beam – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/implicit_dynamics/vibrating-beam/vibrating-beam.html

---

Verification Manual
Solid mechanics
Implicit dynamics
Implicit-dynamic vibrating beam
Implicit-dynamic vibration of a three-dimensional beam
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify direct implicit-dynamic time integration for native Abaqus, body-fitted Coreform IGA, and immersed Coreform IGA discretizations. The six-case study exercises both the
TRANSIENT_FIDELITY
and
MODERATE_DISSIPATION
Abaqus application shortcuts
[1]
.
2
Geometry
The model is a rod of length
\(L=10\)
aligned with the
\(x\)
axis and a
\(2\times2\)
square cross-section.
Figure 1: Vibrating-beam geometry and boundary conditions.
3
Material
The rod is linear elastic with Young’s modulus
\(E=1\times10^6\)
, density
\(\rho=1\times10^{-6}\)
, and Poisson’s ratio
\(\nu=0\)
in a dimensionally consistent unit system. The zero Poisson’s ratio makes the response effectively one-dimensional.
4
Loading & boundary conditions
The
\(-x\)
face is fixed, and the lateral faces are free. A pressure on the
\(+x\)
face rises smoothly from zero to its peak during one analytical fundamental-mode period, returns smoothly to zero during the next period, and then remains zero. The free-vibration comparison starts after the load has been removed.
5
Mesh & discretization
The native reference uses structured
C3D8R
elements with nominal size
\(h=1\)
. Both Coreform IGA discretizations use quadratic,
\(C^1\)
splines with the same nominal size; one is body-fitted and the other is immersed.
Figure 2: Background mesh used by the immersed Coreform IGA model.
All cases cover
\(0\leq t\leq2\times10^{-3}\)
with an initial increment of
\(1\times10^{-6}\)
and a maximum increment of
\(2\times10^{-6}\)
. Axial displacement, velocity, and acceleration are recorded at every accepted increment and averaged over the free end.
6
Reference & accepted solutions
The fundamental axial frequency of a fixed-free rod is
\[
f_1=\frac{1}{4L}\sqrt{\frac{E}{\rho}}=25{,}000.
\]
Each case must predict
\(f_1\)
within 2%. Each Coreform IGA history is also compared with the native Abaqus history using the same application shortcut.
Quantity
Acceptance criterion
Fundamental-frequency relative error
\(<0.02\)
Normalized displacement-history RMS error
\(<0.01\)
Free-vibration phase RMS error
\(<0.2\ \text{rad}\)
Amplitude-decay-ratio difference
\(<0.03\)
Energy-decay-ratio difference
\(<0.05\)
Accepted-increment count ratio
\(\leq1.10\)
Samples per fundamental-mode period
\(\geq19.9\)
7
Results
7.1
Fundamental frequency
Table 1: Implicit-dynamic vibrating-beam verification results.
Application
Discretization
\(f_1\)
Relative error
Accepted increments
TRANSIENT_FIDELITY
Native Abaqus
24745.6
1.02%
1003
TRANSIENT_FIDELITY
Body-fitted Coreform IGA
24768.6
0.93%
1003
TRANSIENT_FIDELITY
Immersed Coreform IGA
24766.8
0.93%
1003
MODERATE_DISSIPATION
Native Abaqus
24681.0
1.28%
1002
MODERATE_DISSIPATION
Body-fitted Coreform IGA
24703.8
1.18%
1002
MODERATE_DISSIPATION
Immersed Coreform IGA
24702.1
1.19%
1002
The maximum normalized displacement-history RMS error is
\(2.21\times10^{-3}\)
, and the maximum free-vibration phase RMS error is
\(0.155\ \text{rad}\)
. The maximum amplitude- and energy-decay-ratio differences are
\(304\times10^{-6}\)
and
\(1.11\times10^{-3}\)
, respectively.
7.2
Representative displacement field
The final recorded
TRANSIENT_FIDELITY
immersed result is shown using the same CTest result database that supplies the verification evidence.
Figure 3: Axial displacement in the immersed Coreform IGA vibrating-beam model.
8
Discussion
Every spatial discretization and application shortcut predicts the analytical fundamental frequency within the enforced tolerance. The body-fitted and immersed histories track their matching native Abaqus references without using more accepted increments. The sampling criterion confirms that agreement is not obtained by taking substantially smaller time increments.
References
[1]
Abaqus analysis guide
, 2026 FD03. Dassault Systèmes, 2026.
Appendix
8.1
Download artifacts
Download all verification artifacts (ZIP)
8.1.1
Geometry CAD files
Download the moderate-dissipation body-fitted Coreform model
Download the moderate-dissipation immersed Coreform model
Download the transient-fidelity body-fitted Coreform model
Download the transient-fidelity immersed Coreform model
8.1.2
Abaqus/CAE journal files
6 generated Abaqus/CAE journals are available.
Show 6 individual journals
md_bodyfit–vibrating_beam_3d_model.jnl
md_immersed–vibrating_beam_3d_model.jnl
md_native–vibrating_beam_3d_model.jnl
tf_bodyfit–vibrating_beam_3d_model.jnl
tf_immersed–vibrating_beam_3d_model.jnl
tf_native–vibrating_beam_3d_model.jnl
8.2
Provenance
The records below identify the six simulation CTests and their cross-case study CTest.
Note
Test_abaqus_verification_vibrating_beam_3d_tf_native
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:40:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_vibrating_beam_3d_tf_bodyfit
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:15:30
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_vibrating_beam_3d_tf_immersed
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:16:49
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_vibrating_beam_3d_md_native
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:17:40
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_vibrating_beam_3d_md_bodyfit
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:17:24
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_vibrating_beam_3d_md_immersed
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 19:16:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_study_vibrating_beam_3d_core
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 26/26 enforced criteria passed
Test date
04-Nov-2026 19:41:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
