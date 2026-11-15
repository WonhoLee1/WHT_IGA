# coreform – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/python-api/coreform/index.html

---

User Guide
Scripting Reference
coreform
Source directory:
interop/abaqus/cf_app/coreform
Directory tree
coreform/
abaqus/
cf_export/
job/
mesh/
probe/
job/
mesh/
probe/
utils/
1
types.py
Module path:
coreform.types
Available API members
Name
Kind
Description
Float3
Attribute
Three floating-point components.
Point3D
Attribute
A point represented by three floating-point coordinates.
Vector3D
Attribute
A vector represented by three floating-point components.
1.1
Float3
Object path:
coreform.types.Float3
Type:
Tuple[float, float, float]
Three floating-point components.
1.2
Point3D
Object path:
coreform.types.Point3D
Type:
Tuple[float, float, float]
A point represented by three floating-point coordinates.
1.3
Vector3D
Object path:
coreform.types.Vector3D
Type:
Tuple[float, float, float]
A vector represented by three floating-point components.
2
validation.py
Module path:
coreform.validation
Available API members
Name
Kind
Description
ValidationIssue.path
Attribute
Field path locating the problem within a nested specification.
ValidationIssue.code
Attribute
Stable machine-readable issue identifier.
ValidationIssue.message
Attribute
Human-readable explanation of the problem.
ValidationError.issues
Attribute
Structured validation issues carried by this exception.
ValidationIssue
Class
Describe one validation problem with a stable path and code.
ValidationError
Class
Collect one or more structured Coreform validation issues.
2.1
ValidationIssue
Object path:
coreform.validation.ValidationIssue
Describe one validation problem with a stable path and code.
Call syntax
result
=
ValidationIssue(path
=
path, code
=
code, message
=
message)
Signatures
ValidationIssue(path: Tuple[
str
, ...], code:
str
, message:
str
)
Parameters
Name
Type
Required
Default
Description
path
Tuple[str, ...]
Yes
—
Field path locating the problem within a nested specification.
code
str
Yes
—
Stable machine-readable issue identifier.
message
str
Yes
—
Human-readable explanation of the problem.
2.1.1
Attributes
Name
Type
Value
Description
path
Tuple[str, ...]
—
Field path locating the problem within a nested specification.
code
str
—
Stable machine-readable issue identifier.
message
str
—
Human-readable explanation of the problem.
2.2
ValidationError
Object path:
coreform.validation.ValidationError
Collect one or more structured Coreform validation issues.
Call syntax
result
=
ValidationError(issues
=
issues)
Signatures
ValidationError(issues: Union[
'ValidationError'
,
str
, Iterable[ValidationIssue]])
Parameters
Name
Type
Required
Default
Description
issues
Union['ValidationError', str, Iterable[ValidationIssue]]
Yes
—
Existing error, fallback message, or nonempty structured issue collection.
2.2.1
Attributes
Name
Type
Value
Description
issues
Tuple[ValidationIssue, ...]
tuple(issues)
Structured validation issues carried by this exception.
