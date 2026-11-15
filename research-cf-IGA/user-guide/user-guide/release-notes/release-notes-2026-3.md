# Coreform IGA for Abaqus 2026.3 – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/release-notes/release-notes-2026-3.html

---

User Guide
Release notes
2026.3
Coreform IGA for Abaqus 2026.3
Coreform IGA for Abaqus 2026.3 was the initial, limited-distribution release of the Coreform IGA for Abaqus product line. It introduced an integrated workflow for defining IGA discretizations in Abaqus/CAE, generating them with Coreform IGA Mesh, solving them with Abaqus/Standard, and inspecting IGA results.
Released:
17 March 2026
Note
Historical, non-public release
This release was shared with a limited audience and was not broadly announced or marketed. This page provides a concise historical summary of the capabilities available to that audience.
Version 2026.3 predates the current
Capabilities matrix
, whose support classifications begin with 2026.7. The current
Production
,
Beta
, and
Alpha
classifications are therefore not applied retroactively on this page. Consult the matrix and
Current capabilities and limitations
for the support status of a current release.
1
Initial Abaqus integration
1.1
Purpose-built Abaqus/CAE application
This release introduced the Coreform IGA for Abaqus custom application. It retained the familiar Abaqus/CAE modeling modules and added Coreform-aware mesh, probe, and job workflows.
The initial application could create and manage IGA mesh definitions on assembly instances, preview the generated discretization, define result probes, and run the Coreform IGA Mesh and Abaqus/Standard stages as one job workflow. Job controls included separate processor counts for mesh generation and the Abaqus solve.
1.2
Automated setup for Abaqus 2025
The initial setup wizard configured the Abaqus environment, installed the Python dependency used by the integration, and authenticated the Coreform IGA for Abaqus license. This release targeted Abaqus 2025; Abaqus 2024 was not formally qualified.
2
IGA mesh generation
2.1
Immersed rectilinear meshes
Coreform IGA Mesh generated trimmed, immersed IGA discretizations directly from the model geometry exported by Abaqus/CAE. Users could control the background-grid element size, choose axis-aligned, object-aligned, or user-defined orientations, and adjust the frame and padding around the geometry.
The mesh authoring interface exposed linear, quadratic, and cubic spline bases with continuous, maximally smooth, or user-defined continuity. It also included penalty enforcement for immersed boundaries and a selective-integration option for nearly incompressible analysis.
2.2
Limited body-fitted IGA workflow
An initial body-fitted IGA path converted an existing Abaqus mesh into an IGA discretization. This path was limited to fully structured, mapped, all-hexahedral meshes; bottom-up and sweep meshes were not supported by the 2026.3 interface.
2.3
Standalone and parallel preprocessing
The new
coreform_iga_mesh
executable separated IGA mesh generation and persistence from the legacy solver workflow. It processed selected model parts, generated their trimmed quadrature and boundary data, and stored the resulting IGA mesh in a SQLite database for Coreform IGA Interop and the Abaqus user subroutines. The Abaqus job workflow could assign multiple processors to this preprocessing stage independently of the solver configuration.
3
Abaqus/Standard analysis
3.1
Linear, nonlinear, and modal analysis
The initial Abaqus/Standard integration supplied the user-element and material interfaces needed to solve body-fitted and immersed IGA models. The release included workflows for linear static analysis, multistep static analysis, geometrically nonlinear static analysis, natural-frequency extraction, and early direct implicit dynamics.
Available material behavior included linear elasticity, hyperelastic response, and isotropic metal plasticity. The nonlinear implementation included large-deformation kinematics, and selective integration was available for locking-sensitive, nearly incompressible models.
3.2
Loads, boundary conditions, and constraints
Coreform IGA Interop provided initial translation paths for prescribed and pinned boundary conditions, concentrated forces, uniform pressure, surface traction, body force, gravity, and moments transferred through kinematic couplings. These definitions were authored with standard Abaqus/CAE objects and translated where they targeted an IGA discretization.
3.3
Assemblies and material interfaces
The release included assembly workflows containing multiple IGA parts, mixed native finite-element and IGA parts, tied material interfaces, and both matching and nonmatching interfaces. An initial cellular-material workflow also read cell-based section assignments from Abaqus/CAE Boolean merges created with retained intersections, allowing distinct material definitions on the retained cells.
3.4
Early surface-to-surface contact
An early Abaqus/Standard surface-to-surface contact path generated conventional boundary elements on IGA surfaces so Abaqus could apply its native contact algorithm. The release included an early two-part, finite-sliding contact workflow, but contact remained an emerging capability rather than a broadly qualified production workflow.
4
Results and probes
4.1
IGA results in Abaqus/CAE
Coreform IGA Interop created a companion
<job-name>_iga.odb
containing a visualization mesh for the IGA geometry. It mapped displacement, stress, and equivalent plastic-strain results onto that mesh for supported time-domain analyses, allowing the results to be viewed with standard Abaqus/CAE visualization tools.
4.2
Point and line probes
The custom application introduced point and line probes for sampling IGA results at physical locations independent of visualization-mesh nodes. Users could preview probe locations in Abaqus/CAE, manage probe definitions with the model, and export displacement, stress, and equivalent plastic-strain histories from the solve.
