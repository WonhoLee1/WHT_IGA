# Shock and vibration – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/shock_and_vibration/index.html

---

Verification Manual
Solid mechanics
Shock and vibration
Shock and vibration
Verification of mode-based linear dynamics: frequency extraction, random response to a power spectral density, and peak response to a shock spectrum.
Coreform IGA elements are Abaqus user elements, and Abaqus/Standard does not post-process a random response or a response spectrum for them. For an immersed model the frequency extraction is the last result the solver produces, so Coreform computes the remaining procedures from that modal basis. These studies verify both halves of that arrangement: the Coreform procedures against Abaqus on an identical modal basis, and the immersed modal basis against a body-fitted Abaqus reference.
Verification studies
Production
Cantilever beam shock and vibration
Verify eigenfrequencies, RMS response to a base-acceleration spectral density, and peak response to a shock spectrum, for an immersed Coreform IGA mesh against a body-fitted Abaqus analysis.
Open study →
