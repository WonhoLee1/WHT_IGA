# Coreform IGA for Abaqus 2026.10 – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/release-notes/release-notes-2026-10.html

---

User Guide
Release notes
2026.10
Coreform IGA for Abaqus 2026.10
Coreform IGA for Abaqus 2026.10 brings local refinement into the Abaqus/CAE interface, adds a beta visualization-output workflow, and improves beta MPI execution for supported jobs. These notes cover changes since the 2026.9 release.
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
1.1
Local refinement in the GUI
You can now create and manage local refinements in the Mesh module without writing a script. Attach a refinement to an immersed IGA mesh, choose a refinement level, and select a named CAD geometry set, a sphere or spherical shell, or a box as its region. The Refinement Manager supports editing, renaming, copying, suppressing, and deleting refinements.
Preview a region while editing it, or display refinement lines on the mesh preview in the viewport. The mesh preview is an approximation of the generated mesh; run the Coreform IGA Mesh stage to obtain the actual discretization. Local refinement remains an
Alpha
capability.
1.2
Mesh preview controls
The Mesh Manager now provides viewport previews of the background mesh and its boundary approximation before mesh generation. Display controls let you inspect the mesh and refinement regions while working in a part or assembly.
1.3
Improved GUI performance
Module switching and mesh and probe editing are more responsive in large projects. Viewport caching reuses unchanged geometry and mesh display data to reduce repeated calculations and drawing work.
2
Analysis and Abaqus interoperability
2.1
Improved output (Beta)
The new output workflow creates a conforming visualization mesh for
<job-name>_iga.odb
, giving continuous interior contours and field values evaluated at output-node locations. It also reduces output database size in tested models. This changes visualization output, not the underlying analysis.
This workflow is
Beta
and enabled by default for new GUI jobs, command-line meshing, and interoperability. The Job Manager’s
Use new output (Beta)
setting controls both the Mesh and Abaqus interoperability stages. To use the previous workflow, clear that setting or pass
--legacy-output-mesh
to both
coreform_iga_mesh
and
coreform_interop_abaqus
. Regenerate the mesh database after switching to the previous workflow; a database containing an immersed output mesh continues to select the new workflow. No workflow flag is needed to use the new output mesh.
2.2
MPI execution (Beta)
Multi-process Abaqus/Standard execution is available at the
Beta
level for supported Coreform IGA workflows. The runtime now discovers which IGA elements belong to each solver rank and loads their analysis data accordingly. This release also improves MPI initialization for heat transfer and contact models.
In the Job Manager’s
Parallelization
tab, set the Abaqus CPU count and select
MPI
as the Abaqus mode. The IGA mesh-generation CPU count is configured separately. Check the intended model and procedure with a single-process reference run when evaluating this beta workflow.
