# Coreform IGA for Abaqus 2026.11 – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/release-notes/release-notes-2026-11.html

---

User Guide
Release notes
2026.11
Coreform IGA for Abaqus 2026.11
Coreform IGA for Abaqus 2026.11 improves local refinement, job monitoring, meshing performance, and results visualization. These notes cover changes since the 2026.10 release.
Alpha
and
Beta
describe the support levels in
Current capabilities and limitations
. Check the
Capabilities matrix
for the scope and verification status of a particular workflow.
1
Abaqus/CAE application
Mesh generation and Abaqus interoperability now run in the background, keeping the interface responsive during these stages.
Revamped our
Job Monitor
window to now show progress, logs, warnings, and errors, and you can cancel jobs from the monitor.
Mesh and refinement previews are faster to redisplay and use less memory. Larger meshes can now be previewed without the previous one-million-cell limit.
2
Local refinement (Beta)
Local refinement of immersed rectilinear meshes advances to
Beta
, with verification for single-part, linear-elastic stress concentrations near named faces.
Highly local refinement near geometric features is more robust.
Leaving
Offset
blank now selects a default value based on the mesh degree and refinement level.
Refined meshes use less memory during generation and make better use of parallel processing.
3
Meshing and analysis
Meshing is more robust for complex CAD geometry, including small holes and curved surfaces.
Abaqus solves use less memory and run faster for models affected by small trimmed cells. In one large assembly benchmark, peak solver memory and solve time were both reduced by about half.
Pressure and traction loads are handled more reliably across analysis steps, including load changes and resets.
Experimental
Alpha
formulations for nearly incompressible analysis now include B-bar and F-bar options with selected Abaqus elastic, hyperelastic, and plastic materials. Material and analysis combinations remain limited.
4
Results and visualization
The newer visualization workflow is now enabled by default, providing continuous contours through the model interior. To use the previous workflow, uncheck
Use new output
box in the Job Manager and regenerate the mesh. For command-line jobs, pass
--legacy-output-mesh
to both the meshing and Abaqus interoperability stages.
Visualization meshes supported with local refinement.
An optional unified ODB mode writes IGA and native Abaqus results into the same results database. Enable it with an environment variable:
CF_ABQ_ODB_OUTPUT_MODE=native
.
The ParaView plugin reads IGA results databases directly, including results from MPI jobs. You can view deformed results, inspect mesh and CAD feature lines, and use slices and clips to explore the interior.
Strain output now follows Abaqus solid-element conventions. Mixed finite-element and IGA models gain compatible nodal output for easier comparison. The mixed-formulations’ (B-Bar, F-Bar) pressure fields are now named
PRESS_IGA
in ODB output.
Experimental (
Alpha
) support for kinetic, stored, and plastic energy histories are now available for supported material workflows . Complete energy balances are not yet available for all models.
5
Installation and upgrading
coreform_iga_mesh.bat
has been removed from Windows installations; run
coreform_iga_mesh.exe
directly. A saved Job whose
Mesh executable
points at
coreform_iga_mesh.bat
must be changed to
coreform_iga_mesh.exe
in the same
bin
directory.
New
Windows compiler setup instructions
explain how to use your own Abaqus user subroutines with Coreform IGA.
You will need to regenerate IGA input mesh databases created with earlier releases before running analyses with 2026.11.
