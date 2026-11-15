# coreform.abaqus.cf_export – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/python-api/coreform/abaqus/cf_export/index.html

---

User Guide
Scripting Reference
abaqus
cf_export
coreform.abaqus.cf_export
Source directory:
interop/abaqus/cf_app/coreform/abaqus/cf_export
Directory tree
cf_export/
1
export.py
Module path:
coreform.abaqus.cf_export.export
Available API members
Name
Kind
Description
export_cf
Function
Write the Coreform
.cf
file for an Abaqus model.
1.1
export_cf
Object path:
coreform.abaqus.cf_export.export.export_cf
Description:
Write the Coreform
.cf
file for an Abaqus model.
Call syntax
result
=
coreform.abaqus.cf_export.export.export_cf(mdb_object
=
mdb_object, model_name
=
model_name, cf_path
=
cf_path)
Returns
Type
Description
str
The supplied destination path after the file has been written.
Signatures
export_cf(mdb_object:
object
, model_name:
str
, cf_path:
str
, probe_specs: Optional[Sequence[
object
]]
=
None
, part_metadata: Optional[Mapping[
str
, Mapping[
str
,
object
]]]
=
None
)
->
str
Parameters
Name
Type
Required
Default
Description
mdb_object
object
Yes
—
Abaqus model database containing the model to export.
model_name
str
Yes
—
Name of the Abaqus model to export.
cf_path
str
Yes
—
Destination path for the Coreform
.cf
file.
probe_specs
Optional[Sequence[object]]
No
None
Probe specifications to serialize. When omitted, collect them from the
model’s persistent Probe repository.
part_metadata
Optional[Mapping[str, Mapping[str, object]]]
No
None
Per-part material and locking-free overrides keyed by instance or part
name.
2
locking_free_configuration.py
Module path:
coreform.abaqus.cf_export.locking_free_configuration
Available API members
Name
Kind
Description
CutCellPressureTreatment.NONE
Attribute
Do not apply a special cut-cell pressure treatment.
CutCellPressureTreatment.PRESSURE_AGGREGATION
Attribute
Aggregate cut-cell pressure degrees of freedom.
LockingFreeStabilization.NONE
Attribute
Do not add stabilization.
LockingFreeStabilization.VMS
Attribute
Apply variational multiscale stabilization.
locking_free_parameters
Function
Build a complete composable locking-free configuration.
has_explicit_pressure_dofs
Function
Return whether a configuration retains explicit pressure degrees of freedom.
locking_free_part_metadata
Function
Build material metadata for a locking-free Part.
requires_unsymmetric_matrix_storage
Function
Return whether the locking-free tangent must retain unsymmetric blocks.
CutCellPressureTreatment
Class
Select how pressure degrees of freedom are treated in cut cells.
LockingFreeStabilization
Class
Select the stabilization used by a locking-free formulation.
2.1
CutCellPressureTreatment
Object path:
coreform.abaqus.cf_export.locking_free_configuration.CutCellPressureTreatment
Select how pressure degrees of freedom are treated in cut cells.
2.1.1
Attributes
Name
Type
Value
Description
NONE
CutCellPressureTreatment
'none'
Do not apply a special cut-cell pressure treatment.
PRESSURE_AGGREGATION
CutCellPressureTreatment
'pressure_aggregation'
Aggregate cut-cell pressure degrees of freedom.
2.2
LockingFreeStabilization
Object path:
coreform.abaqus.cf_export.locking_free_configuration.LockingFreeStabilization
Select the stabilization used by a locking-free formulation.
2.2.1
Attributes
Name
Type
Value
Description
NONE
LockingFreeStabilization
'none'
Do not add stabilization.
VMS
LockingFreeStabilization
'vms'
Apply variational multiscale stabilization.
2.3
locking_free_parameters
Object path:
coreform.abaqus.cf_export.locking_free_configuration.locking_free_parameters
Description:
Build a complete composable locking-free configuration.
Call syntax
result
=
coreform.abaqus.cf_export.locking_free_configuration.locking_free_parameters(displacement_degree
=
displacement_degree)
Returns
Type
Description
dict
Metadata for a condensed
U(p,C^(p-2))/P0
configuration by default.
Signatures
locking_free_parameters(displacement_degree:
int
, displacement_continuity: Optional[
int
]
=
None
, pressure_degree:
int
=
0
, pressure_continuity:
int
=
-
1
, macro_elements:
bool
=
False
, static_condensation:
bool
=
True
, cut_cell_pressure_treatment: Union[CutCellPressureTreatment,
str
]
=
CutCellPressureTreatment.PRESSURE_AGGREGATION, stabilization: Union[LockingFreeStabilization,
str
]
=
LockingFreeStabilization.NONE)
->
dict
[
str
,
object
]
Parameters
Name
Type
Required
Default
Description
displacement_degree
int
Yes
—
Polynomial degree of the displacement basis.
displacement_continuity
Optional[int]
No
None
Displacement continuity, or
None
to use
degree - 2
.
pressure_degree
int
No
0
Polynomial degree of the pressure basis.
pressure_continuity
int
No
-1
Continuity of the pressure basis.
macro_elements
bool
No
False
Whether to use macro-element stabilization.
static_condensation
bool
No
True
Whether to statically condense pressure degrees of freedom.
cut_cell_pressure_treatment
Union[CutCellPressureTreatment, str]
No
CutCellPressureTreatment.PRESSURE_AGGREGATION
Pressure treatment applied to cut cells.
stabilization
Union[LockingFreeStabilization, str]
No
LockingFreeStabilization.NONE
Stabilization applied to the formulation.
2.4
has_explicit_pressure_dofs
Object path:
coreform.abaqus.cf_export.locking_free_configuration.has_explicit_pressure_dofs
Description:
Return whether a configuration retains explicit pressure degrees of freedom.
Call syntax
result
=
coreform.abaqus.cf_export.locking_free_configuration.has_explicit_pressure_dofs(params
=
params)
Returns
Type
Description
bool
Return whether a configuration retains explicit pressure degrees of freedom.
Signatures
has_explicit_pressure_dofs(params: Mapping[
str
,
object
])
->
bool
Parameters
Name
Type
Required
Default
Description
params
Mapping[str, object]
Yes
—
Locking-free configuration metadata.
2.5
locking_free_part_metadata
Object path:
coreform.abaqus.cf_export.locking_free_configuration.locking_free_part_metadata
Description:
Build material metadata for a locking-free Part.
Call syntax
result
=
coreform.abaqus.cf_export.locking_free_configuration.locking_free_part_metadata(params
=
params, use_locking_free
=
use_locking_free, youngs_modulus
=
youngs_modulus, poissons_ratio
=
poissons_ratio)
Returns
Type
Description
dict
Complete exporter metadata for the Part material and formulation.
Signatures
locking_free_part_metadata(params: Mapping[
str
,
object
], use_locking_free:
bool
, youngs_modulus:
float
, poissons_ratio:
float
)
->
dict
[
str
,
object
]
Parameters
Name
Type
Required
Default
Description
params
Mapping[str, object]
Yes
—
Composable locking-free configuration and optional material settings.
use_locking_free
bool
Yes
—
Whether the Part uses the locking-free formulation.
youngs_modulus
float
Yes
—
Young’s modulus assigned to the Part material.
poissons_ratio
float
Yes
—
Poisson’s ratio assigned to the Part material.
2.6
requires_unsymmetric_matrix_storage
Object path:
coreform.abaqus.cf_export.locking_free_configuration.requires_unsymmetric_matrix_storage
Description:
Return whether the locking-free tangent must retain unsymmetric blocks.
Call syntax
result
=
coreform.abaqus.cf_export.locking_free_configuration.requires_unsymmetric_matrix_storage(params
=
params)
Returns
Type
Description
bool
Return whether the locking-free tangent must retain unsymmetric blocks.
Signatures
requires_unsymmetric_matrix_storage(params: Mapping[
str
,
object
])
->
bool
Parameters
Name
Type
Required
Default
Description
params
Mapping[str, object]
Yes
—
Locking-free configuration metadata.
