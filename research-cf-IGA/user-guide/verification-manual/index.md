# Verification Manual – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/index.html

---

Verification Manual
This manual presents the engineering evidence used to assess
Coreform IGA for Abaqus
. Each study defines a problem, identifies its reference solution, and reports the numerical checks used to evaluate Coreform IGA.
The reported results are processed directly from automated verification runs. Figures, tables, and downloadable artifacts are generated from those results when the manual is built, reducing manual transcription and keeping the published evidence aligned with the tested software.
Verified capabilities
The verification hierarchy is organized by physical capability. The Coreform IGA discretization currently consists of three-dimensional solid spline elements; it does not provide one-dimensional beam or two-dimensional shell elements.
Solid mechanics
6 verification studies
Linear elastic stress
Verify stress fields, stress concentrations, load and boundary-condition types, maximally smooth spline bases, multiple materials, and curved trimmed boundaries under small-strain linear elasticity.
Browse verification studies →
1 verification study
Thermal strain
Verify stresses produced by prescribed thermal strain across a material interface.
Browse verification studies →
1 verification study
Implicit dynamics
Verify Abaqus implicit-dynamic application shortcuts with native and Coreform IGA discretizations.
Browse verification studies →
Heat transfer
3 verification studies
Steady and transient heat transfer
Verify conduction with radiation, convection, thermal capacity, and implicit time integration.
Browse verification studies →
How verification becomes documentation
Each study begins with a defined engineering problem and acceptance criteria. Automated tests run Coreform IGA, collect the relevant solution quantities, and compare them with an accepted analytic or numerical reference.
The manual is authored with
Quarto
. During a documentation build, Quarto reads the completed test results, computes the summaries used in tables, generates plots with Matplotlib, requests visualizations from ParaView where appropriate, creates download links, and combines those outputs with the engineering explanation on each page.
Figure 1: High-level workflow from product development to published verification evidence.
This workflow makes the relationship between a verification claim and its supporting evidence visible to readers while allowing the explanatory prose to remain clear and focused on engineering interpretation.
