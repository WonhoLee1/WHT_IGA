# coreform.abaqus.mesh – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/python-api/coreform/abaqus/mesh/index.html

---

User Guide
Scripting Reference
abaqus
mesh
coreform.abaqus.mesh
Source directory:
interop/abaqus/cf_app/coreform/abaqus/mesh
Directory tree
mesh/
1
records.py
Module path:
coreform.abaqus.mesh.records
Available API members
Name
Kind
Description
CoreformIGAMeshInfo.frame_orientation_mode
Attribute
GUI mode used to derive the Mesh frame orientation.
CoreformIGAMeshInfo.padding_mode
Attribute
GUI mode used to derive frame padding.
CoreformIGAMeshInfo.continuity_mode
Attribute
GUI mode used to derive basis continuity.
CoreformIGAMeshInfo.bc_penalty_setting
Attribute
GUI setting used to derive the boundary-condition penalty.
CoreformIGA.mesh
Attribute
The persistent Mesh record, or
None
when no Mesh is defined.
mesh.meshType
Attribute
Type of the persistent Mesh record.
mesh.visible
Attribute
Whether the Mesh is visible in Abaqus/CAE.
mesh.degree
Attribute
Polynomial degree of the Mesh basis.
mesh.continuity
Attribute
Continuity of the Mesh basis.
mesh.bcPenaltyValue
Attribute
Penalty value used to impose essential boundary conditions.
mesh.tessellationSize
Attribute
Target tessellation size, or
None
when no target is specified.
mesh.elementSize
Attribute
Target element size when size-based Mesh sizing is active.
mesh.elementIntervalCount
Attribute
Element interval counts when interval-count sizing is active.
mesh.elementIntervalCountByAxis
Attribute
Per-axis interval counts when per-axis sizing is active.
mesh.elementLayout
Attribute
Element-layout choice for each parametric axis.
mesh.smallCellVolumeRatio
Attribute
Volume-ratio threshold used to identify small immersed cells.
mesh.padding
Attribute
Additional cell padding along each frame axis.
mesh.frameOrigin
Attribute
Origin of an explicit Mesh frame, or
None
for owner resolution.
mesh.frameBasis
Attribute
Basis of an explicit Mesh frame, or
None
for owner resolution.
mesh.frameExtents
Attribute
Minimum and maximum extents in the Mesh frame.
CoreformIGAMeshInfo.to_dict
Instance method
Return the GUI sidecar values as a dictionary.
CoreformIGA.owner
Instance method
Return the Abaqus Part or PartInstance that owns this namespace.
CoreformIGA.ImmersedRectilinearMesh
Instance method
Create and persist an immersed rectilinear Mesh for this owner.
CoreformIGA.BodyFitMesh
Instance method
Create and persist a body-fit Mesh for this owner.
CoreformIGA.deleteMesh
Instance method
Delete the Mesh currently attached to this owner.
mesh.owner
Instance method
Return the Abaqus Part or PartInstance that owns this Mesh.
mesh.toSpec
Instance method
Return a portable specification representing this Mesh record.
mesh.meshInfo
Instance method
Return the Abaqus GUI sidecar associated with this Mesh.
mesh.validateForOwner
Instance method
Validate this Mesh for a specified Abaqus owner.
mesh.setValues
Instance method
Update this persistent Mesh through Abaqus-facing camelCase fields.
mesh.resolveForOwner
Instance method
Return the portable Mesh specification with its owner-dependent frame resolved.
CoreformIGAMeshInfo
Class
Store Abaqus/CAE GUI choices that derive portable Mesh values.
CoreformIGA
Class
Persistent Mesh namespace attached to an Abaqus Part or PartInstance.
mesh
Class
Persistent Abaqus-facing Mesh record.
1.1
CoreformIGAMeshInfo
Object path:
coreform.abaqus.mesh.records.CoreformIGAMeshInfo
Store Abaqus/CAE GUI choices that derive portable Mesh values.
Call syntax
result
=
CoreformIGAMeshInfo()
Signatures
CoreformIGAMeshInfo(frame_orientation_mode: Optional[
str
]
=
None
, padding_mode: Optional[
str
]
=
None
, continuity_mode: Optional[
str
]
=
None
, bc_penalty_setting: Optional[
str
]
=
None
)
Parameters
Name
Type
Required
Default
Description
frame_orientation_mode
Optional[str]
No
None
GUI mode used to derive the Mesh frame orientation.
padding_mode
Optional[str]
No
None
GUI mode used to derive frame padding.
continuity_mode
Optional[str]
No
None
GUI mode used to derive basis continuity.
bc_penalty_setting
Optional[str]
No
None
GUI setting used to derive the boundary-condition penalty.
1.1.1
Attributes
Name
Type
Value
Description
frame_orientation_mode
str or None
frame_orientation_mode
GUI mode used to derive the Mesh frame orientation.
padding_mode
str or None
padding_mode
GUI mode used to derive frame padding.
continuity_mode
str or None
continuity_mode
GUI mode used to derive basis continuity.
bc_penalty_setting
str or None
bc_penalty_setting
GUI setting used to derive the boundary-condition penalty.
1.1.2
to_dict
Defined on:
coreform.abaqus.mesh.records.CoreformIGAMeshInfo.to_dict
Description:
Return the GUI sidecar values as a dictionary.
Call syntax
result
=
coreformIGAMeshInfo.to_dict()
Returns
Type
Description
dict[str, object]
Read-only GUI-choice fields keyed by their portable names.
Signatures
to_dict()
->
dict
[
str
,
object
]
Parameters
None.
1.2
CoreformIGA
Public path:
meshOwner.customData.CoreformIGA
Available API members
Name
Kind
Description
mesh
Attribute
The persistent Mesh record, or
None
when no Mesh is defined.
owner
Instance method
Return the Abaqus Part or PartInstance that owns this namespace.
ImmersedRectilinearMesh
Instance method
Create and persist an immersed rectilinear Mesh for this owner.
BodyFitMesh
Instance method
Create and persist a body-fit Mesh for this owner.
deleteMesh
Instance method
Delete the Mesh currently attached to this owner.
1.2.1
Attributes
Name
Type
Value
Description
mesh
object
—
The persistent Mesh record, or
None
when no Mesh is defined.
1.2.2
owner
Defined on:
meshOwner.customData.CoreformIGA.owner
Description:
Return the Abaqus Part or PartInstance that owns this namespace.
Call syntax
result
=
meshes.owner()
Returns
Type
Description
object
Return the Abaqus Part or PartInstance that owns this namespace.
Signatures
owner()
->
object
Parameters
None.
1.2.3
ImmersedRectilinearMesh
Defined on:
meshOwner.customData.CoreformIGA.ImmersedRectilinearMesh
Description:
Create and persist an immersed rectilinear Mesh for this owner.
Call syntax
result
=
meshes.ImmersedRectilinearMesh(spec
=
spec)
result
=
meshes.ImmersedRectilinearMesh(elementSize
=
elementSize)
Returns
Type
Description
ImmersedRectilinearMesh record
Newly created persistent Mesh record.
Signatures
ImmersedRectilinearMesh(
*
, spec, visible
=
True
, meshInfo
=
None
, overwrite
=
False
)
->
object
ImmersedRectilinearMesh(
*
, degree
=
1
, continuity
=
0
, bcPenaltyValue
=
None
, tessellationSize
=
None
, elementSize
=
None
, elementIntervalCount
=
None
, elementIntervalCountByAxis
=
None
, elementLayout
=
None
, smallCellVolumeRatio
=
0.2
, padding
=
None
, frameOrigin
=
None
, frameBasis
=
None
, frameExtents
=
None
, visible
=
True
, meshInfo
=
None
, overwrite
=
False
)
->
object
Parameters
For
ImmersedRectilinearMesh(*, spec, visible=True, meshInfo=None, overwrite=False) -> object
:
Name
Type
Required
Default
Description
spec
ImmersedRectilinearMeshSpec or Mapping[str, object]
Yes
—
Complete portable specification or serialized specification mapping.
visible
bool
No
True
Whether the Mesh is visible in Abaqus/CAE.
meshInfo
CoreformIGAMeshInfo or None
No
None
Abaqus GUI sidecar values.
overwrite
bool
No
False
Replace an existing Mesh atomically when true.
For
ImmersedRectilinearMesh(*, degree=1, continuity=0, bcPenaltyValue=None, tessellationSize=None, elementSize=None, elementIntervalCount=None, elementIntervalCountByAxis=None, elementLayout=None, smallCellVolumeRatio=0.2, padding=None, frameOrigin=None, frameBasis=None, frameExtents=None, visible=True, meshInfo=None, overwrite=False) -> object
:
Name
Type
Required
Default
Description
degree
int
No
1
Polynomial degree of the Mesh basis.
continuity
int
No
0
Continuity of the Mesh basis.
bcPenaltyValue
float or None
No
None
Penalty value used to impose essential boundary conditions.
tessellationSize
float or None
No
None
Target tessellation size.
elementSize
Float3 or None
Exactly one of
elementSize
,
elementIntervalCount
, or
elementIntervalCountByAxis
None
Target element size.
elementIntervalCount
Tuple[int, int, int] or None
Exactly one of
elementSize
,
elementIntervalCount
, or
elementIntervalCountByAxis
None
Element interval counts.
elementIntervalCountByAxis
Tuple[int, int, int] or None
Exactly one of
elementSize
,
elementIntervalCount
, or
elementIntervalCountByAxis
None
Per-axis element interval counts.
elementLayout
Tuple[str, str, str] or None
No
None
Element-layout choice for each axis.
smallCellVolumeRatio
float
No
0.2
Volume-ratio threshold used to identify small cells.
padding
Tuple[int, int, int] or None
No
None
Additional cells along each frame axis.
frameOrigin
Vector3D or None
No
None
Origin of an explicit Mesh frame.
frameBasis
Tuple[Vector3D, Vector3D, Vector3D] or None
No
None
Basis of an explicit Mesh frame.
frameExtents
Tuple[Float3, Float3] or None
No
None
Minimum and maximum frame extents.
visible
bool
No
True
Whether the Mesh is visible in Abaqus/CAE.
meshInfo
CoreformIGAMeshInfo or None
No
None
Abaqus GUI sidecar values.
overwrite
bool
No
False
Replace an existing Mesh atomically when true.
Parameter rules
Rule
Requirement
Call forms
Use exactly one of these call forms:
spec
or
abaqus_definition
.
Required group
Supply exactly one of
elementSize
,
elementIntervalCount
, or
elementIntervalCountByAxis
.
Parameter group
Supply
frameOrigin
,
frameBasis
, and
frameExtents
together, or omit all of them.
Raises
Type
Condition
ValidationError
If the specification or owner-dependent Mesh settings are invalid.
1.2.4
BodyFitMesh
Defined on:
meshOwner.customData.CoreformIGA.BodyFitMesh
Description:
Create and persist a body-fit Mesh for this owner.
Call syntax
result
=
meshes.BodyFitMesh(spec
=
spec)
result
=
meshes.BodyFitMesh()
Returns
Type
Description
BodyFitMesh record
Newly created persistent Mesh record.
Signatures
BodyFitMesh(
*
, spec, visible
=
True
, meshInfo
=
None
, overwrite
=
False
)
->
object
BodyFitMesh(
*
, degree
=
1
, continuity
=
0
, bcPenaltyValue
=
None
, tessellationSize
=
None
, visible
=
True
, meshInfo
=
None
, overwrite
=
False
)
->
object
Parameters
For
BodyFitMesh(*, spec, visible=True, meshInfo=None, overwrite=False) -> object
:
Name
Type
Required
Default
Description
spec
BodyFitMeshSpec or Mapping[str, object]
Yes
—
Complete portable specification or serialized specification mapping.
visible
bool
No
True
Whether the Mesh is visible in Abaqus/CAE.
meshInfo
CoreformIGAMeshInfo or None
No
None
Abaqus GUI sidecar values.
overwrite
bool
No
False
Replace an existing Mesh atomically when true.
For
BodyFitMesh(*, degree=1, continuity=0, bcPenaltyValue=None, tessellationSize=None, visible=True, meshInfo=None, overwrite=False) -> object
:
Name
Type
Required
Default
Description
degree
int
No
1
Polynomial degree of the Mesh basis.
continuity
int
No
0
Continuity of the Mesh basis.
bcPenaltyValue
float or None
No
None
Penalty value used to impose essential boundary conditions.
tessellationSize
float or None
No
None
Target tessellation size.
visible
bool
No
True
Whether the Mesh is visible in Abaqus/CAE.
meshInfo
CoreformIGAMeshInfo or None
No
None
Abaqus GUI sidecar values.
overwrite
bool
No
False
Replace an existing Mesh atomically when true.
Parameter rules
Rule
Requirement
Call forms
Use exactly one of these call forms:
spec
or
abaqus_definition
.
Raises
Type
Condition
ValidationError
If the specification or owner-dependent Mesh settings are invalid.
1.2.5
deleteMesh
Defined on:
meshOwner.customData.CoreformIGA.deleteMesh
Description:
Delete the Mesh currently attached to this owner.
Call syntax
meshes.deleteMesh()
Returns
Type
Description
None
—
Signatures
deleteMesh()
->
None
Parameters
None.
1.3
mesh
Public path:
meshOwner.customData.CoreformIGA.mesh
Available API members
Name
Kind
Description
meshType
Attribute
Type of the persistent Mesh record.
visible
Attribute
Whether the Mesh is visible in Abaqus/CAE.
degree
Attribute
Polynomial degree of the Mesh basis.
continuity
Attribute
Continuity of the Mesh basis.
bcPenaltyValue
Attribute
Penalty value used to impose essential boundary conditions.
tessellationSize
Attribute
Target tessellation size, or
None
when no target is specified.
elementSize
Attribute
Target element size when size-based Mesh sizing is active.
elementIntervalCount
Attribute
Element interval counts when interval-count sizing is active.
elementIntervalCountByAxis
Attribute
Per-axis interval counts when per-axis sizing is active.
elementLayout
Attribute
Element-layout choice for each parametric axis.
smallCellVolumeRatio
Attribute
Volume-ratio threshold used to identify small immersed cells.
padding
Attribute
Additional cell padding along each frame axis.
frameOrigin
Attribute
Origin of an explicit Mesh frame, or
None
for owner resolution.
frameBasis
Attribute
Basis of an explicit Mesh frame, or
None
for owner resolution.
frameExtents
Attribute
Minimum and maximum extents in the Mesh frame.
owner
Instance method
Return the Abaqus Part or PartInstance that owns this Mesh.
toSpec
Instance method
Return a portable specification representing this Mesh record.
meshInfo
Instance method
Return the Abaqus GUI sidecar associated with this Mesh.
validateForOwner
Instance method
Validate this Mesh for a specified Abaqus owner.
setValues
Instance method
Update this persistent Mesh through Abaqus-facing camelCase fields.
resolveForOwner
Instance method
Return the portable Mesh specification with its owner-dependent frame resolved.
1.3.1
Attributes
Name
Type
Applies to
Value
Description
meshType
MeshType
All
—
Type of the persistent Mesh record.
visible
bool
All
—
Whether the Mesh is visible in Abaqus/CAE.
degree
int
All
—
Polynomial degree of the Mesh basis.
continuity
int
All
—
Continuity of the Mesh basis.
bcPenaltyValue
Optional[float]
All
—
Penalty value used to impose essential boundary conditions.
tessellationSize
Optional[float]
All
—
Target tessellation size, or
None
when no target is specified.
elementSize
Optional[Float3]
ImmersedRectilinearMesh
—
Target element size when size-based Mesh sizing is active.
elementIntervalCount
Optional[Tuple[int, int, int]]
ImmersedRectilinearMesh
—
Element interval counts when interval-count sizing is active.
elementIntervalCountByAxis
Optional[Tuple[int, int, int]]
ImmersedRectilinearMesh
—
Per-axis interval counts when per-axis sizing is active.
elementLayout
Optional[Tuple[str, str, str]]
ImmersedRectilinearMesh
—
Element-layout choice for each parametric axis.
smallCellVolumeRatio
float
ImmersedRectilinearMesh
—
Volume-ratio threshold used to identify small immersed cells.
padding
Optional[Tuple[int, int, int]]
ImmersedRectilinearMesh
—
Additional cell padding along each frame axis.
frameOrigin
Optional[Vector3D]
ImmersedRectilinearMesh
—
Origin of an explicit Mesh frame, or
None
for owner resolution.
frameBasis
Optional[Tuple[Vector3D, Vector3D, Vector3D]]
ImmersedRectilinearMesh
—
Basis of an explicit Mesh frame, or
None
for owner resolution.
frameExtents
Optional[Tuple[Float3, Float3]]
ImmersedRectilinearMesh
—
Minimum and maximum extents in the Mesh frame.
1.3.2
owner
Defined on:
meshOwner.customData.CoreformIGA.mesh.owner
Description:
Return the Abaqus Part or PartInstance that owns this Mesh.
Call syntax
result
=
mesh.owner()
Returns
Type
Description
object
Return the Abaqus Part or PartInstance that owns this Mesh.
Signatures
owner()
->
object
Parameters
None.
1.3.3
toSpec
Defined on:
meshOwner.customData.CoreformIGA.mesh.toSpec
Description:
Return a portable specification representing this Mesh record.
Call syntax
result
=
mesh.toSpec()
Returns
Type
Description
MeshSpec
Return a portable specification representing this Mesh record.
Signatures
toSpec()
->
MeshSpec
Parameters
None.
1.3.4
meshInfo
Defined on:
meshOwner.customData.CoreformIGA.mesh.meshInfo
Description:
Return the Abaqus GUI sidecar associated with this Mesh.
Call syntax
result
=
mesh.meshInfo()
Returns
Type
Description
CoreformIGAMeshInfo
Return the Abaqus GUI sidecar associated with this Mesh.
Signatures
meshInfo()
->
CoreformIGAMeshInfo
Parameters
None.
1.3.5
validateForOwner
Defined on:
meshOwner.customData.CoreformIGA.mesh.validateForOwner
Description:
Validate this Mesh for a specified Abaqus owner.
Call syntax
mesh.validateForOwner(owner
=
owner)
Returns
Type
Description
None
—
Signatures
validateForOwner(owner:
object
, allow_existing:
bool
=
False
)
->
None
Parameters
Name
Type
Required
Default
Description
owner
object
Yes
—
Abaqus Part or PartInstance against which to validate the Mesh.
allow_existing
bool
No
False
Whether validation permits an existing Mesh on the owner.
1.3.6
setValues
Defined on:
meshOwner.customData.CoreformIGA.mesh.setValues
Description:
Update this persistent Mesh through Abaqus-facing camelCase fields.
Call syntax
mesh.setValues(...)
Returns
Type
Description
None
The record is updated in place.
Signatures
ImmersedRectilinearMesh.setValues(
*
, degree
=
..., continuity
=
..., bcPenaltyValue
=
..., tessellationSize
=
..., elementSize
=
..., elementIntervalCount
=
..., elementIntervalCountByAxis
=
..., elementLayout
=
..., smallCellVolumeRatio
=
..., padding
=
..., frameOrigin
=
..., frameBasis
=
..., frameExtents
=
..., visible
=
..., meshInfo
=
...)
->
None
BodyFitMesh.setValues(
*
, degree
=
..., continuity
=
..., bcPenaltyValue
=
..., tessellationSize
=
..., visible
=
..., meshInfo
=
...)
->
None
Parameters
For
ImmersedRectilinearMesh.setValues(*, degree=..., continuity=..., bcPenaltyValue=..., tessellationSize=..., elementSize=..., elementIntervalCount=..., elementIntervalCountByAxis=..., elementLayout=..., smallCellVolumeRatio=..., padding=..., frameOrigin=..., frameBasis=..., frameExtents=..., visible=..., meshInfo=...) -> None
:
Name
Applies to
Type
Required
Default
Description
degree
All
int
No
Current value
Replace the polynomial degree.
continuity
All
int
No
Current value
Replace the basis continuity.
bcPenaltyValue
All
float or None
No
Current value
Replace the boundary-condition penalty value.
tessellationSize
All
float or None
No
Current value
Replace the target tessellation size.
elementSize
ImmersedRectilinearMesh
Float3
No
Current value
Switch to size-based rectilinear sizing.
elementIntervalCount
ImmersedRectilinearMesh
Tuple[int, int, int]
No
Current value
Switch to interval-count rectilinear sizing.
elementIntervalCountByAxis
ImmersedRectilinearMesh
Tuple[int, int, int]
No
Current value
Switch to per-axis interval-count sizing.
elementLayout
ImmersedRectilinearMesh
Tuple[str, str, str] or None
No
Current value
Replace the element-layout choices.
smallCellVolumeRatio
ImmersedRectilinearMesh
float
No
Current value
Replace the small-cell volume-ratio threshold.
padding
ImmersedRectilinearMesh
Tuple[int, int, int] or None
No
Current value
Replace the frame padding.
frameOrigin
ImmersedRectilinearMesh
Vector3D or None
No
Current value
Replace the explicit frame origin.
frameBasis
ImmersedRectilinearMesh
Tuple[Vector3D, Vector3D, Vector3D] or None
No
Current value
Replace the explicit frame basis.
frameExtents
ImmersedRectilinearMesh
Tuple[Float3, Float3] or None
No
Current value
Replace the frame extents.
visible
All
bool
No
Current value
Replace the Abaqus/CAE visibility setting.
meshInfo
All
CoreformIGAMeshInfo or None
No
Current value
Replace the Abaqus GUI sidecar.
For
BodyFitMesh.setValues(*, degree=..., continuity=..., bcPenaltyValue=..., tessellationSize=..., visible=..., meshInfo=...) -> None
:
Name
Applies to
Type
Required
Default
Description
degree
All
int
No
Current value
Replace the polynomial degree.
continuity
All
int
No
Current value
Replace the basis continuity.
bcPenaltyValue
All
float or None
No
Current value
Replace the boundary-condition penalty value.
tessellationSize
All
float or None
No
Current value
Replace the target tessellation size.
visible
All
bool
No
Current value
Replace the Abaqus/CAE visibility setting.
meshInfo
All
CoreformIGAMeshInfo or None
No
Current value
Replace the Abaqus GUI sidecar.
Parameter rules
Rule
Requirement
Optional exclusive group
If updating this group, supply exactly one of
elementSize
,
elementIntervalCount
, or
elementIntervalCountByAxis
.
Parameter group
Supply
frameOrigin
,
frameBasis
, and
frameExtents
together, or omit all of them.
Applicability
elementSize
applies only to
ImmersedRectilinearMesh
.
Applicability
elementIntervalCount
applies only to
ImmersedRectilinearMesh
.
Applicability
elementIntervalCountByAxis
applies only to
ImmersedRectilinearMesh
.
Applicability
elementLayout
applies only to
ImmersedRectilinearMesh
.
Applicability
smallCellVolumeRatio
applies only to
ImmersedRectilinearMesh
.
Applicability
padding
applies only to
ImmersedRectilinearMesh
.
Applicability
frameOrigin
applies only to
ImmersedRectilinearMesh
.
Applicability
frameBasis
applies only to
ImmersedRectilinearMesh
.
Applicability
frameExtents
applies only to
ImmersedRectilinearMesh
.
Raises
Type
Condition
ValidationError
If an update is unsupported or produces an invalid Mesh.
1.3.7
resolveForOwner
Defined on:
meshOwner.customData.CoreformIGA.mesh.resolveForOwner
Description:
Return the portable Mesh specification with its owner-dependent frame resolved.
Call syntax
result
=
mesh.resolveForOwner()
Returns
Type
Description
ImmersedRectilinearMeshSpec
Return the portable Mesh specification with its owner-dependent frame resolved.
Signatures
resolveForOwner(owner:
object
=
None
)
->
ImmersedRectilinearMeshSpec
Parameters
Name
Type
Required
Default
Description
owner
object
No
None
Abaqus owner to resolve, or
None
to use this record’s owner.
