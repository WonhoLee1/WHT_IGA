# Mode-based dynamics – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/response-spectrum/index.html

---

User Guide
Mode-based dynamics
Mode-based dynamics
Estimate peak and RMS response from an Abaqus frequency extraction, outside the solver.
Mode-based dynamics
Two linear dynamics procedures reuse the eigenmodes of a frequency extraction instead of integrating a structure through time: response spectrum analysis, which estimates a peak response to base motion, and random response analysis, which gives spectral densities and RMS values for a statistically described excitation.
The
respspec
utility carries out both in Python, outside the solver. It reads the eigenmodes and participation factors that Abaqus/Standard writes during a
*FREQUENCY
step, applies your spectra or spectral densities, and reports displacement, velocity, acceleration, and stress. Because the modal data is extracted once and reused, you can evaluate many load definitions against a single eigen-solution without rerunning the solver or holding an Abaqus license.
The utility ships in the Coreform IGA for Abaqus source tree at
interop/abaqus/response_spectrum/
. It reads any Abaqus output database that contains a mass-normalized frequency extraction, whether the model uses native finite elements, Coreform IGA elements, or both.
For a Coreform IGA model the utility is not just a convenience. IGA meshes are Abaqus user elements, and Abaqus/Standard does not post-process a random response for them, so the frequency extraction is the last result the solver produces.
0.1
Topics
Start here
Running a response spectrum analysis
Prepare the frequency step, extract the modal data, combine it against a spectrum, and read the peak results.
Run an analysis →
Response spectrum case files
Define spectra, excitation directions, modal and directional summation rules, and modal damping in a case file.
Write a case file →
Random response analysis
Compute power spectral densities, RMS values, and RMS von Mises stress from a spectral density excitation.
Run a random response →
0.2
Choosing a procedure
The two procedures answer different questions, and they take different inputs.
Response spectrum
Random response
Input
A spectrum: peak response against frequency
A power spectral density of the excitation
Assumes
A transient event, such as a shock or an earthquake
A stationary, ergodic process, such as road roughness or acoustic noise
Reports
A peak estimate, unsigned
Spectral densities and RMS values
Combination
Summation rules across modes and directions
Integration of the response spectrum across the frequency range
Use response spectrum when a design code hands you a spectrum, and random response when the excitation is described statistically, as a vibration qualification test is. Neither gives a time history; if you need one, integrate through time with a modal dynamic step instead.
0.3
What the response spectrum procedure computes
The procedure follows the Abaqus response spectrum method.
First, the peak amplitude of each generalized coordinate is found for mode
\(\alpha\)
and spectrum
\(k\)
:
\[
(q_\alpha)_k = c_k \, S_k(\omega_\alpha, \xi_\alpha) \sum_j t_{kj} \Gamma_{\alpha j}
\]
Here
\(c_k\)
is the scale factor you apply to spectrum
\(k\)
,
\(S_k\)
is the spectrum interpolated at the natural frequency
\(\omega_\alpha\)
and the modal damping
\(\xi_\alpha\)
,
\(t_{kj}\)
are the direction cosines of the excitation, and
\(\Gamma_{\alpha j}\)
is the participation factor of mode
\(\alpha\)
in direction
\(j\)
.
The peak of a physical variable in that single mode follows from the mode shape
\(\Phi_\alpha\)
:
\[
(R^i_\alpha)_k = \Phi^i_\alpha \, (q_\alpha)_k
\]
These single-mode peaks are then combined twice, because peaks in different modes and different directions do not occur at the same instant. The modal summation rule combines across modes, and the directional summation rule combines across excitation directions.
0.4
Supported summation rules
Rule
Modal summation
Directional summation
Absolute values
ABS
—
Square root of the sum of the squares
SRSS
SRSS
Naval Research Laboratory
NRL
—
Ten percent
TENP
—
Complete quadratic combination
CQC
—
Grouping
GRP
—
Double sum combination
DSC
—
Algebraic
—
ALGEBRAIC
Forty percent rule
—
R40
Thirty percent rule
—
R30
ALGEBRAIC
combines the excitation directions before the modal summation, which is the correct treatment when the spectra are components of one base motion acting in an arbitrary direction. Every other directional rule combines the per-direction results after the modal summation, which suits statistically independent excitations.
0.5
What the random response procedure computes
The excitation is a cross-spectral density projected onto the modes, each mode responds through its complex transfer function, and the physical spectral density follows from the mode shapes. RMS values are the square root of that spectral density integrated over the frequency range. RMS von Mises stress needs the quadratic form of Segalman et al., because an invariant cannot be recovered from the RMS of the components. The
random response page
gives the equations and the case file format.
0.6
Limitations
For response spectrum, the rigid-response methods, Gupta and Lindley-Yow, and the missing mass method are not implemented. Results therefore contain the periodic response of the extracted modes only, with no correction for truncated high-frequency modes.
For random response, only base motion excitation is implemented. Concentrated, distributed, and connector loads, moving noise correlations, the
UPSD
and
UCORR
user subroutines, and secondary bases are not.
Reaction forces cannot be checked against a solver-computed response spectrum step, because Abaqus does not write
RF
in a SIM-based response spectrum procedure. The utility still combines
RF
when the frequency step wrote it as a mode shape.
0.7
Verification
Both procedures are verified against Abaqus/Standard itself. Three steel solids, a cube, a truncated sphere, and a cylinder, are each solved with eight response spectrum steps and six random response steps. Between them these span every supported summation rule, three spectrum types, three ways of specifying modal damping, one and three excitation directions, correlated and uncorrelated loading, superposed load cases, a displacement base spectrum, and both frequency scales.
Recombining the same modes in Python reproduces the Abaqus results to within a few parts in ten million of the peak value, which is the precision of the single-precision fields Abaqus writes. The random response frequency sweep matches the Abaqus grid point for point.
Quantity
Worst relative error
Response spectrum displacement, velocity, acceleration, stress
\(2.9 \times 10^{-6}\)
Random response
RU
,
RV
,
RA
,
RS
,
RMISES
\(320 \times 10^{-9}\)
Sampled spectral densities
\(1.3 \times 10^{-6}\)
Two details of the Abaqus response spectrum implementation differ from the published description, and the utility follows the observed behavior in both cases:
Between two tabulated damping curves, the interpolation weight is linear in the damping ratio while the magnitude is interpolated logarithmically.
In the grouping method, modal peaks inside a group are summed with their signs rather than as absolute values, so near-repeated modes of a symmetric part cancel.
