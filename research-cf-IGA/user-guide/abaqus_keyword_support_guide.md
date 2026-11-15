# Abaqus keyword support – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/abaqus_keyword_support_guide.html

---

User guide
Abaqus keyword support
Use this guide to see how Coreform IGA for Abaqus handles the Abaqus keyword families it actively manages in models with IGA regions. It is not a complete catalog of Abaqus keywords.
13
keyword families
49
documented cases
27
supported cases
What this guide covers
Keyword records from families not listed below remain unchanged in the translated
.inp
for Abaqus to process. Those keywords may or may not be officially supported or validated by Coreform IGA; pass-through only means that Coreform IGA did not modify them.
A listed keyword family is actively managed by Coreform IGA. Within a listed family, only the options and combinations shown have the documented behavior; an unlisted option or combination may be rejected and stop input-file translation rather than passing through. Overall Coreform IGA procedure and product limitations still apply.
How to use this page
A keyword can appear in several rows because a load type, option, output variable, model scope, or condition can change how Coreform IGA handles it. Related options with the same support and translated-input behavior are combined in one row. Search or filter the table, find the row that matches your model, then read its support status and behavior.
A supported Abaqus-native-only row means that Coreform IGA leaves that behavior to Abaqus; it does not mean that Coreform IGA implements the native Abaqus feature.
Supported
Coreform IGA supports this case under the scope and conditions shown in the table.
Not supported for IGA
Coreform IGA cannot provide this behavior for IGA regions. It may reject or skip the IGA portion, as described in the row.
Recognized
Coreform IGA reads this keyword for context or checks. This status does not mean that the keyword itself adds an IGA physics or results capability.
Output-request shortcuts
Component requests
U1
–
U3
,
V1
–
V3
,
A1
–
A3
, and
RF1
–
RF3
use the
U
,
V
,
A
, and
RF
support rows, respectively.
S11
,
S22
,
S33
,
S12
,
S13
, and
S23
use the
S
row. Coreform IGA expands
ALL
and
PRESELECT
only to the output variables listed in this guide and warns when other Abaqus variables may be implied.
Keyword support table
Search descriptions as well as keywords and output variables. Use the filters when you want to compare one keyword family or support status.
Search
Keyword
All keywords
*DLOAD (10)
*DSLOAD (7)
*SFILM (3)
*SRADIATE (3)
*SOLID SECTION (5)
*OUTPUT (1)
*FIELD OUTPUT (6)
*HISTORY OUTPUT (1)
*NODE OUTPUT (4)
*ELEMENT OUTPUT (4)
*INCLUDE (1)
*STEP (2)
*CONTACT INCLUSIONS (2)
Support
All support statuses
Supported (27)
Not supported for IGA (16)
Recognized (6)
Clear filters
Showing
49
of
49
cases in
13
keyword families
.
No matching cases.
Try a broader search or clear one of the filters.
*DLOAD
10 cases
Support cases for *DLOAD
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Gravity
GRAV
IGA only
When the direction vector is nonzero
Supported
Coreform IGA replaces the source
GRAV
data line in the translated
.inp
with internal
U994
,
U995
, and
U996
component-load records. The IGA solver applies those components to the targeted IGA elements.
GRAV
model-wide
When the direction vector is nonzero
Supported
Coreform IGA keeps the model-wide
GRAV
record in the translated
.inp
for Abaqus-native elements and adds internal
U994
,
U995
, and
U996
component-load records for every IGA volume element. Abaqus applies the original gravity load to native elements, and the IGA solver applies the generated components to IGA elements.
GRAV
mixed Abaqus-native and IGA
When the direction vector is nonzero
When Abaqus-native and IGA membership can be separated
Supported
Coreform IGA replaces the mixed source record with separate records on generated native-only and IGA-only element sets. The native record retains the original
GRAV
data; internal
U994
,
U995
, and
U996
records apply the gravity components to the IGA elements.
Body forces
BX
,
BY
, or
BZ
IGA only
Supported
Coreform IGA replaces the source
BX
,
BY
, or
BZ
data line in the translated
.inp
with an internal
U997
,
U998
, or
U999
load record, respectively. The IGA solver applies the corresponding x-, y-, or z-direction body force to the targeted IGA elements.
BX
,
BY
, or
BZ
mixed Abaqus-native and IGA
When Abaqus-native and IGA membership can be separated
Supported
Coreform IGA replaces the mixed source record with separate records on generated native-only and IGA-only element sets. The native record retains the original
BX
,
BY
, or
BZ
data; an internal
U997
,
U998
, or
U999
record applies the corresponding body force to the IGA elements.
Reserved load types
U994
–
U996
IGA only
Not supported for IGA
Source input files must use the standard Abaqus
GRAV
load type. Coreform IGA reserves
U994
,
U995
, and
U996
for the generated x-, y-, and z-components of IGA gravity loads.
Result:
Coreform IGA stops input-file translation for this case.
U997
–
U999
IGA only
Not supported for IGA
Source input files must use the standard Abaqus
BX
,
BY
, or
BZ
load type. Coreform IGA reserves
U997
,
U998
, and
U999
for the generated x-, y-, and z-components of IGA body-force loads.
Result:
Coreform IGA stops input-file translation for this case.
Other targets
Any load type
Abaqus-native only
Supported
Coreform IGA leaves native-only
*DLOAD
records unchanged in the translated
.inp
for Abaqus to apply.
Any load type
mixed Abaqus-native and IGA
When Abaqus-native and IGA membership cannot be separated
Not supported for IGA
Coreform IGA cannot safely divide this mixed
*DLOAD
target into native-only and IGA-only element sets.
Result:
Coreform IGA stops input-file translation for this case.
Any load type
unresolved
Not supported for IGA
Coreform IGA cannot determine which elements an unresolved
*DLOAD
record would load, so it cannot apply the load to IGA elements.
Result:
Coreform IGA stops input-file translation for this case.
*DSLOAD
7 cases
Support cases for *DSLOAD
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Pressure
P
IGA only
Supported
Coreform IGA removes the IGA pressure data line from the translated
.inp
, stores the pressure and any named amplitude internally, and applies the load to the targeted IGA surface at runtime.
Result:
OP=NEW resets previous surface loads once per step; retained loads must be respecified. Removed loads follow the step amplitude variation to zero. Supported single-definition replacements retain ramp continuity. Modifying additive load groups remains unsupported.
P
mixed Abaqus-native and IGA
When Abaqus-native and IGA membership can be separated
Supported
Coreform IGA splits the target into generated native-only and IGA-only surfaces. The translated
.inp
retains an Abaqus pressure record for the native surface; Coreform IGA removes the IGA pressure data from the deck, stores it and any named amplitude internally, and applies it to the IGA surface at runtime.
Result:
OP=NEW resets previous surface loads once per step; retained loads must be respecified. Removed loads follow the step amplitude variation to zero. Supported single-definition replacements retain ramp continuity. Modifying additive load groups remains unsupported.
Vector tractions
TRSHR
or
TRVEC
IGA only
When the traction direction vector is nonzero
Supported
Coreform IGA removes the IGA
TRSHR
or
TRVEC
traction data line from the translated
.inp
, stores the traction vector and any named amplitude internally, and applies it to the targeted IGA surface at runtime.
Result:
OP=NEW resets previous surface loads once per step; retained loads must be respecified. Removed loads follow the step amplitude variation to zero. Supported single-definition replacements retain ramp continuity. Modifying additive load groups remains unsupported.
TRSHR
or
TRVEC
mixed Abaqus-native and IGA
When the traction direction vector is nonzero
When Abaqus-native and IGA membership can be separated
Supported
Coreform IGA splits the target into generated native-only and IGA-only surfaces. The translated
.inp
retains an Abaqus
TRSHR
or
TRVEC
record for the native surface; Coreform IGA removes the IGA traction data from the deck, stores it and any named amplitude internally, and applies it to the IGA surface at runtime.
Result:
OP=NEW resets previous surface loads once per step; retained loads must be respecified. Removed loads follow the step amplitude variation to zero. Supported single-definition replacements retain ramp continuity. Modifying additive load groups remains unsupported.
Other targets
Any load type
Abaqus-native only
Supported
Coreform IGA leaves native-only
*DSLOAD
records unchanged in the translated
.inp
for Abaqus to apply.
Any load type
mixed Abaqus-native and IGA
When Abaqus-native and IGA membership cannot be separated
Not supported for IGA
Coreform IGA cannot safely divide this mixed
*DSLOAD
target into native-only and IGA-only surfaces.
Result:
Coreform IGA stops input-file translation for this case.
Any load type
model-wide or unresolved
Not supported for IGA
Coreform IGA cannot extract an IGA surface load when
*DSLOAD
is model-wide or its target cannot be resolved to a supported surface.
Result:
Coreform IGA stops input-file translation for this case.
*SFILM
3 cases
Support cases for *SFILM
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Surface film condition
F
IGA only
During a heat-transfer step
Without a named amplitude
Supported
Coreform IGA removes the IGA film condition from the translated deck, stores its sink temperature and film coefficient, and integrates it on the tessellated physical boundary at runtime.
Any film condition
Abaqus-native only
Supported
Coreform IGA leaves native-only
*SFILM
records unchanged for Abaqus to apply.
Unsupported target or option
IGA, mixed, or unresolved
Not supported for IGA
Coreform IGA stops translation when it cannot safely own the film condition.
Result:
Input-file translation stops for these cases.
*SRADIATE
3 cases
Support cases for *SRADIATE
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Radiation to ambient
R
IGA only
During a heat-transfer step
Without a named amplitude
Supported
Coreform IGA removes the IGA radiation condition from the translated deck, stores its ambient temperature and emissivity, and integrates it on the tessellated physical boundary at runtime.
Any radiation condition
Abaqus-native only
Supported
Coreform IGA leaves native-only
*SRADIATE
records unchanged for Abaqus to apply.
Unsupported target or option
IGA, mixed, or unresolved
Not supported for IGA
Coreform IGA stops translation when it cannot safely own the radiation condition.
Result:
Input-file translation stops for these cases.
*SOLID SECTION
5 cases
Support cases for *SOLID SECTION
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
IGA sections
With
MATERIAL
IGA only
When
MATERIAL
is specified
Supported
Coreform IGA removes the IGA
*SOLID SECTION
record from the translated setup
.inp
and transfers its
MATERIAL
assignment to generated
*UEL PROPERTY
data.
Without
MATERIAL
IGA only
When
MATERIAL
is missing
Not supported for IGA
Coreform IGA cannot create the IGA
*UEL PROPERTY
data without a
MATERIAL
assignment on the source
*SOLID SECTION
record.
Result:
Coreform IGA stops input-file translation for this case.
Other targets
—
Abaqus-native only
Supported
Coreform IGA leaves native-only
*SOLID SECTION
records unchanged in the translated setup
.inp
for Abaqus to use.
—
mixed Abaqus-native and IGA
When the section spans Abaqus-native and IGA elements
Not supported for IGA
Coreform IGA cannot split one
*SOLID SECTION
assignment between Abaqus-native and IGA elements.
Result:
Coreform IGA stops input-file translation for this case.
—
model-wide or unresolved
Supported
Coreform IGA leaves an unresolved
*SOLID SECTION
record unchanged in the translated setup
.inp
; it is not used to create IGA properties.
*OUTPUT
1 case
Support cases for *OUTPUT
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Field or history request context
Not target-specific
Recognized
Coreform IGA reads
*OUTPUT
to apply its scope and schedule to the requests that follow;
*OUTPUT
alone does not add an IGA result variable. Variables specified directly on
*OUTPUT
follow the
*FIELD OUTPUT
rows.
*FIELD OUTPUT
6 cases
Support cases for *FIELD OUTPUT
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Request context
Field request context
Not target-specific
Recognized
Coreform IGA uses
*FIELD OUTPUT
to identify field-output requests; the variable rows below determine what Coreform IGA writes for IGA regions.
Displacement, stress, and strain
U
,
S
,
PEEQ
,
E
,
NE
,
LE
, or
NT11
IGA only or model-wide
Supported
Coreform IGA leaves the source output request in the translated
.inp
, records the IGA displacement, stress, equivalent-plastic-strain, total-strain, nominal-strain, logarithmic-strain, or nodal-temperature request internally, and adds any helper output needed. Abaqus produces any native-region output, and Coreform IGA writes the requested IGA values to the IGA results database.
Result:
Requesting PEEQ automatically includes the binary AC YIELD (actively yielding) companion in the IGA ODB, with the same regions and output schedule. AC YIELD is not a standalone input variable.
Result:
For supported general solid analyses, E/NE/LE requests resolve to total infinitesimal E in small strain; finite-strain E requests resolve to LE while NE remains nominal strain. Linear perturbation strain output is not implemented for IGA.
U
,
S
,
PEEQ
,
E
,
NE
, or
LE
mixed Abaqus-native and IGA
Supported
Coreform IGA preserves the Abaqus-native output request in the translated
.inp
and records the IGA portion internally. Abaqus produces the native-region output; Coreform IGA adds any needed helper output and writes the requested IGA values to the IGA results database.
Result:
Requesting PEEQ automatically includes the binary AC YIELD (actively yielding) companion in the IGA ODB, with the same regions and output schedule. AC YIELD is not a standalone input variable.
Result:
For supported general solid analyses, E/NE/LE requests resolve to total infinitesimal E in small strain; finite-strain E requests resolve to LE while NE remains nominal strain. Linear perturbation strain output is not implemented for IGA.
Other field variables
RMISES
,
MISES
,
RE
,
PE
,
V
,
A
, or
RF
Includes IGA entities
Not supported for IGA
Coreform IGA leaves any Abaqus-native portion of the request in the translated
.inp
for Abaqus to produce, but does not write these variables for IGA regions.
Result:
Coreform IGA reports a warning and omits the requested variable from IGA results; translation of the remaining model continues.
Abaqus-native results
Any other Abaqus output variable
All targets
Not supported for IGA
Coreform IGA leaves native-only field-output requests unchanged. For IGA regions, every other Abaqus output variable is recognized by this fallback row and skipped with a warning when no IGA extractor is implemented.
Result:
No IGA result is written for variables matched by this fallback.
Any output variable
Abaqus-native only
Supported
Coreform IGA leaves native-only field-output requests unchanged in the translated
.inp
for Abaqus to produce.
*HISTORY OUTPUT
1 case
Support cases for *HISTORY OUTPUT
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
History request context
Not target-specific
Recognized
Coreform IGA reads
*HISTORY OUTPUT
to interpret the requests that follow and leaves them in the translated
.inp
for Abaqus to process, but Coreform IGA does not write history results for IGA regions.
*NODE OUTPUT
4 cases
Support cases for *NODE OUTPUT
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Nodal displacement
U
or
NT11
IGA only or model-wide
Supported
Coreform IGA leaves the source nodal-output request in the translated
.inp
, records the IGA nodal-displacement or nodal-temperature request internally, and adds any helper output needed. Abaqus produces any native-region output, and Coreform IGA writes the requested nodal field to the IGA results database.
Nodal motion and reaction variables
V
,
A
, or
RF
Includes IGA entities
Not supported for IGA
Coreform IGA leaves any Abaqus-native portion of the request in the translated
.inp
for Abaqus to produce, but does not write nodal velocity, acceleration, or reaction-force values for IGA regions.
Result:
Coreform IGA reports a warning and omits the requested variable from IGA results; translation of the remaining model continues.
Abaqus-native results
Any other Abaqus nodal output variable
All targets
Not supported for IGA
Coreform IGA leaves native-only
*NODE OUTPUT
requests unchanged. For IGA regions, every other Abaqus nodal variable—including thermal output such as
NT
—is recognized and skipped with a warning when no IGA extractor is implemented.
Result:
No IGA nodal result is written for variables matched by this fallback.
Any output variable
Abaqus-native only
Supported
Coreform IGA leaves native-only
*NODE OUTPUT
requests unchanged in the translated
.inp
for Abaqus to produce.
*ELEMENT OUTPUT
4 cases
Support cases for *ELEMENT OUTPUT
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
Stress and strain
S
,
PEEQ
,
E
,
NE
, or
LE
IGA only or model-wide
Supported
Coreform IGA leaves the source element-output request in the translated
.inp
, records the IGA element-stress, equivalent-plastic-strain, total-strain, nominal-strain, or logarithmic-strain request internally, and adds any helper output needed. Abaqus produces any native-region output, and Coreform IGA writes the requested IGA element values to the IGA results database.
Result:
Requesting PEEQ automatically includes the binary AC YIELD (actively yielding) companion in the IGA ODB, with the same regions and output schedule. AC YIELD is not a standalone input variable.
Result:
For supported general solid analyses, E/NE/LE requests resolve to total infinitesimal E in small strain; finite-strain E requests resolve to LE while NE remains nominal strain. Linear perturbation strain output is not implemented for IGA.
Other stress and strain variables
RMISES
,
MISES
,
RE
, or
PE
Includes IGA entities
Not supported for IGA
Coreform IGA leaves any Abaqus-native portion of the request in the translated
.inp
for Abaqus to produce, but does not write these element variables for IGA regions.
Result:
Coreform IGA reports a warning and omits the requested variable from IGA results; translation of the remaining model continues.
Abaqus-native results
Any other Abaqus element output variable
All targets
Not supported for IGA
Coreform IGA leaves native-only
*ELEMENT OUTPUT
requests unchanged. For IGA regions, every other Abaqus element variable is recognized and skipped with a warning when no IGA extractor is implemented.
Result:
No IGA element result is written for variables matched by this fallback.
Any output variable
Abaqus-native only
Supported
Coreform IGA leaves native-only
*ELEMENT OUTPUT
requests unchanged in the translated
.inp
for Abaqus to produce.
*INCLUDE
1 case
Support cases for *INCLUDE
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
—
Not target-specific
Recognized
Coreform IGA follows input files referenced by
*INCLUDE
. In copy mode, it rewrites include paths when needed so translated setup copies and the generated IGA geometry file are referenced from their new locations.
*STEP
2 cases
Support cases for *STEP
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
UNSYMM
Not target-specific
When
UNSYMM
is enabled
Recognized
Coreform IGA reads
UNSYMM
when validating solver settings and otherwise leaves the Abaqus step option unchanged.
NLGEOM
without
UNSYMM
Not target-specific
When
NLGEOM
is enabled without
UNSYMM
Recognized
Coreform IGA leaves the step unchanged and reports a warning because some nonlinear IGA workflows require unsymmetric matrix storage.
Result:
Coreform IGA reports a warning and leaves the step unchanged; input-file translation continues.
*CONTACT INCLUSIONS
2 cases
Support cases for *CONTACT INCLUSIONS
Options, load types, or output variables
Applies to
Support
What Coreform IGA does
ALL EXTERIOR
ALL EXTERIOR
model-wide
When IGA exterior surfaces are available
Supported
After generating IGA exterior surfaces, Coreform IGA inserts exclusion and initialization records for them into the translated
.inp
. This covers this setup only, not every Abaqus general-contact workflow.
ALL EXTERIOR
model-wide
When no IGA exterior surfaces were generated
Supported
Coreform IGA inserts no IGA exterior-surface initialization records because no IGA exterior surfaces were generated. It may still add general-contact exclusion records to the translated
.inp
.
Result:
Coreform IGA reports a warning; translation continues without IGA exterior-surface initialization records.
Need help deciding whether a keyword not listed in this guide is compatible with your model? Contact Coreform support with the keyword, options, and target region from your input file.
Coreform IGA for Abaqus documentation
