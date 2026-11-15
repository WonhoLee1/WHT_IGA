# Example Problems – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/example-problems/index.html

---

Example Problems
These examples provide complete Abaqus/CAE workflows for preparing, solving, and postprocessing Coreform IGA for Abaqus simulations. Choose an example by experience level and by the capability you want to explore.
Difficulty and feature maturity are independent: an advanced example may use released functionality, while an alpha example may require an advanced activation or authoring workflow.
Getting started
Absolute beginner
Solid mechanics
Stress concentration in a thin plate with a hole
Build and solve a complete Coreform IGA model from an empty Abaqus/CAE database. This example assumes no previous Abaqus experience and introduces each module in sequence.
Start the first example →
Advanced capabilities
Advanced
Thermomechanics
Unsupported development
Thermal cycling of a hermetic glass-to-metal seal
Evaluate residual glass stress and metal plasticity through a qualification thermal cycle using temperature-dependent public material data.
Open the advanced example →
Advanced
Solid mechanics
Unreleased alpha
Nearly incompressible elasticity — Cook’s membrane
Exercise the alpha locking-free formulation in a large-deformation Cook’s membrane problem, including its environment-variable activation workflow.
Open the advanced example →
Advanced
Solid mechanics
Unreleased alpha
Surface-to-surface contact — door-lock assembly
Explore alpha surface-to-surface contact in a three-part immersed assembly with frictional sliding, large-deformation kinematics, and elastic-plastic materials.
Open the advanced example →
How examples are organized
While customers may develop their own prepared workflows, these examples will generally follow the “modules” workflow in Abaqus/CAE in a linear manner. Each example problem will begin with any relevant model database operations such as creating the model or defining physical constants. Thus the overall workflow for each example problem will generally follow:
Model definition
Create a new model database & rename the default model
Add a description to the model attributes, including the
consistent unit system
Part module
Creation of base geometry
Property module
Material and material section creation
Assignment of material sections to geometry
Assembly module
Root assembly definition
Creation of merged-geometry parts
Step module
Step definition and Abaqus outputs
[Coreform IGA]
Creation of Coreform IGA probes for accurate output
Interaction module
Load module
Mesh module
[Coreform IGA]
Creation of Coreform IGA meshes
Job module
[Coreform IGA]
Creation of Coreform IGA for Abaqus jobs
Visualization module
Visualizing field results
[Coreform IGA]
Importing Coreform IGA probe data and plotting as XY-Data
Some examples extend this workflow with additional steps — for instance,
local refinement toward a feature
adds a refinement to the mesh through the scripting interface before the job is submitted. Each such example states its extra steps explicitly.
