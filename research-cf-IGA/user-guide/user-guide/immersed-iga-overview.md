# Immersed IGA overview – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/immersed-iga-overview.html

---

User Guide
Immersed IGA overview
Immersed IGA overview
Coreform IGA for Abaqus uses immersed isogeometric analysis (IGA) to analyze models in Abaqus/Standard using spline-based approximation functions. This page introduces the immersed IGA approach and explains how geometry and discretization in immersed IGA differ from those in a conventional, body-fitted finite element mesh.
1
Traditional finite element analysis
In conventional finite element analysis (FEA), elements and nodes define the discretized geometry of the structure being analyzed. Each element represents a portion of the structure, and adjacent elements are connected through shared nodes. The nodal coordinates, element connectivity, and element types together define the finite element mesh.
A body-fitted finite element mesh places its boundary nodes and element faces on the boundary of the analysis geometry. The mesh is an approximation space for the solution and a discrete representation of the geometry. Curved surfaces, small features, and complex intersections may require smaller elements or higher-order elements to represent the geometry adequately.
The element type, element shape, element location, and total number of elements affect the computed solution. Reducing the element size generally improves the representation of the geometry and solution fields, but it also increases the number of degrees of freedom and the cost of the analysis. For a well-posed problem and an appropriate formulation, a sequence of refined meshes is expected to converge toward a unique numerical solution. A converged numerical solution still approximates the physical problem because the geometry, material behavior, boundary conditions, and loads are all model assumptions.
2
Isogeometric analysis
IGA uses smooth spline functions to represent the solution field and geometry. Conventional finite elements provide a continuous solution field across element boundaries, but the derivatives of the solution field are typically discontinuous at element boundaries. Spline basis functions, however, often extend across several elements and maintain continuous derivatives across element boundaries. As with conventional finite elements, the polynomial degree, continuity, and element size together determine the approximation space. The smoothness provided by continuous splines can provide better representation of curved geometry and yield better results when solving problems that benefit from smooth solution fields.
In body-fitted IGA, the analysis mesh conforms to the geometry boundary. The coefficients that define spline geometry are commonly called control points, while coefficients for displacement and other solution fields become analysis degrees of freedom. Unlike ordinary finite element nodes, spline control points do not generally lie on the physical geometry.
3
Immersed isogeometric analysis
Immersed IGA removes the requirement that the analysis mesh conform to the geometry. The CAD model directly defines the physical analysis domain and a simple background spline grid supplies the approximation space for the solution. For an immersed analysis, Coreform IGA for Abaqus creates a structured rectilinear background grid that encloses the selected CAD part. The background grid does not need to align with the faces, edges, holes, or fillets on the part. Each element of the background grid is classified according to its relationship to the physical CAD domain:
An element completely outside the part is excluded from the analysis.
An element completely inside the part is an untrimmed element and contributes fully to the analysis.
An element that is intersected by the boundary of the part is classified as a trimmed element. Trimmed elements are intersected with the CAD geometry, and only the interior portion contributes to the analysis.
For trimmed elements, Coreform IGA for Abaqus constructs the integration data needed to integrate the governing equations over the interior portion of the element. The resulting analysis domain follows the CAD boundary to the geometric accuracy used during trimming, even though the structured background grid does not align with the boundary.
Characteristic
Body-fitted FEA
Immersed IGA
Volumetric mesh
Conforms to the part boundary
Structured background grid encloses the part
Physical boundary
Represented by exterior element faces
Defined by trimming the background grid with the CAD geometry
Curved and small features
May require local mesh construction to conform to the feature
Intersect the background grid without requiring a conforming volumetric mesh
Solution approximation
Finite element shape functions connected through shared nodes
Smooth spline basis functions that may span several elements
Basis construction
Inferred from the element mesh and its connectivity
Defined independently of the CAD boundary, then restricted to the physical domain
Geometry preparation
Requires a valid body-fitted mesh
Requires a valid CAD domain and a background grid fine enough for the analysis
Immersed IGA can reduce the geometry cleanup and mesh-generation work required for a body-fitted analysis, especially for complex CAD models. However, difficult cut configurations, such as small slivers, can affect integration accuracy and numerical conditioning. Coreform IGA for Abaqus includes methods to address these situations, but care should be taken in areas of high interest. Small geometric features and steep solution gradients must still be resolved by the background mesh.
4
Using immersed IGA in an Abaqus model
When using Coreform IGA for Abaqus, the user defines the analysis model in Abaqus/CAE, including the CAD parts, materials, sections, assembly, steps, loads, boundary conditions, and supported interactions as usual. For each part that will use immersed IGA, the user then creates an
Immersed Rectilinear
Coreform IGA Mesh and specifies its element size, spline degree, and continuity. Coreform IGA for Abaqus then trims the background spline grid, generates the data required by its Abaqus user elements, and transfers supported Abaqus model definitions to the IGA discretization.
Some Abaqus features act on element-based boundary representations. Where required, Coreform IGA for Abaqus creates tessellated surface elements on selected CAD boundaries and couples those surface degrees of freedom to the underlying spline solution. This mechanism supports operations such as applying surface loads, boundary conditions, and supported contact interactions.
A model may contain a mix of parts discretized with Coreform IGA for Abaqus and conventionally meshed Abaqus parts. Abaqus/Standard assembles and solves the combined model, while Coreform IGA for Abaqus supplies the element residuals, tangents, and other data for the IGA user elements. The
Workflow overview
describes the preprocessing, simulation, and postprocessing stages in more detail.
5
Accuracy, refinement, and computational cost
The main discretization controls for an immersed IGA mesh are:
Element size,
\(h\)
:
Smaller background elements provide more local resolution and generally increase the number of elements and degrees of freedom.
Polynomial degree,
\(p\)
:
A higher degree increases the approximation order and the support of each basis function. The recommended degree is
\(p=2\)
.
Continuity:
Greater continuity produces smoother fields across element boundaries and changes the connectivity and support of the basis. The typical continuity is
\(p-1\)
.
Local refinement:
Refinement regions add finer levels only where they are needed, such as around a stress concentration. A refinement band around named CAD entities defaults its width to
\(p\)
times the element size at its level, the narrowest band in which that level’s basis functions touching the entities can be activated. A band much thinner than that adds elements but few or no basis functions.
Boundary tessellation and integration:
Trimmed boundaries and cut elements must be resolved accurately enough for geometry-dependent loads, constraints, contact, and integration. The optional positive
tessellationSize
mesh setting limits the physical triangle-edge length used for the part’s trimming tessellation. Leave it unset to use automatic tessellation sizing.
