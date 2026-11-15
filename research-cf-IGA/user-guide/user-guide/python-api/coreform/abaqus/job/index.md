# coreform.abaqus.job – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/user-guide/python-api/coreform/abaqus/job/index.html

---

User Guide
Scripting Reference
abaqus
job
coreform.abaqus.job
Source directory:
interop/abaqus/cf_app/coreform/abaqus/job
Directory tree
job/
1
errors.py
Module path:
coreform.abaqus.job.errors
Available API members
Name
Kind
Description
JobFileCollisionError.files
Attribute
Existing output paths in deterministic reporting order.
JobRuntimeError
Class
A Coreform IGA job build/submit failure with a user-actionable message.
JobFileCollisionError
Class
Report recognized Abaqus output files that block Job submission.
1.1
JobRuntimeError
Object path:
coreform.abaqus.job.errors.JobRuntimeError
A Coreform IGA job build/submit failure with a user-actionable message.
1.2
JobFileCollisionError
Object path:
coreform.abaqus.job.errors.JobFileCollisionError
Report recognized Abaqus output files that block Job submission.
Call syntax
result
=
JobFileCollisionError(files
=
files)
Signatures
JobFileCollisionError(files: Union[
'JobFileCollisionError'
, Iterable[
str
]])
Parameters
Name
Type
Required
Default
Description
files
Union['JobFileCollisionError', Iterable[str]]
Yes
—
Existing output paths that would be overwritten by submission.
1.2.1
Attributes
Name
Type
Value
Description
files
Tuple[str, ...]
tuple(files)
Existing output paths in deterministic reporting order.
2
records.py
Module path:
coreform.abaqus.job.records
Available API members
Name
Kind
Description
CoreformIGA.jobs
Attribute
Persistent Coreform IGA Job records indexed by name.
jobs[jobName].name
Attribute
Repository key and Job name.
jobs[jobName].model
Attribute
Abaqus model used by the Job, when applicable.
jobs[jobName].description
Attribute
User-facing Job description.
jobs[jobName].meshExecutable
Attribute
Coreform IGA Mesh executable override.
jobs[jobName].mpiexec
Attribute
MPI launcher override for the Mesh stage.
jobs[jobName].numIgaCpus
Attribute
Number of CPUs assigned to the Mesh stage.
jobs[jobName].interopExecutable
Attribute
Abaqus interop executable override, when applicable.
jobs[jobName].numCpus
Attribute
Number of CPUs assigned to the solve stage.
jobs[jobName].numGPUs
Attribute
Number of GPUs assigned to an Abaqus solve stage.
jobs[jobName].multiprocessingMode
Attribute
Parallelization mode for an Abaqus solve stage.
jobs[jobName].threadsPerProcess
Attribute
Threads assigned to each Abaqus solver process.
jobs[jobName].numDomains
Attribute
Number of domains used by the Abaqus solver.
jobs[jobName].memory
Attribute
Memory limit for an Abaqus solve stage.
jobs[jobName].memoryUnits
Attribute
Units used to interpret the Abaqus solver memory limit.
CoreformIGA.deleteJob
Instance method
Delete a named persistent Job from this namespace.
CoreformIGA.Job
Instance method
Create and persist a Coreform IGA Job in the MDB namespace.
CoreformIGA.renameJob
Instance method
Rename a Job that has not been prepared or submitted.
jobs[jobName].toSpec
Instance method
Return the lossless portable specification for this Job record.
jobs[jobName].setValues
Instance method
Update this persistent Job through Abaqus-facing camelCase fields.
jobs[jobName].buildMesh
Instance method
Run the Mesh and interop stages for this Job.
jobs[jobName].writeInput
Instance method
Write the raw Abaqus input deck consumed by the interop stage.
jobs[jobName].submit
Instance method
Prepare and submit this Job through its configured pipeline.
jobs[jobName].waitForCompletion
Instance method
Block until a submitted Job finishes.
jobs[jobName].kill
Instance method
Terminate the running process associated with this Job.
jobs[jobName].getRuntimeStatus
Instance method
Return a snapshot of this Job’s current runtime status.
jobs[jobName].writeRunScript
Instance method
Write this Job pipeline as a runnable script file.
jobs[jobName].openOdb
Instance method
Open and return this Job’s Abaqus output database.
CoreformIGA
Class
The MDB-level Coreform IGA namespace (
mdb.customData.CoreformIGA
) that owns
the
jobs
repository and the public
Job(...)
constructor.
jobs[jobName]
Class
Persistent job record: the nested spec is stored losslessly as JSON
(
_spec_json
, the source of truth for
toSpec
), and the common scalar
fields are loaded as read-only camelCase inspection attrs. Edits go through
setValues(...)
so validation, the stored spec, and journaling stay in step.
2.1
CoreformIGA
Public path:
mdb.customData.CoreformIGA
Available API members
Name
Kind
Description
jobs
Attribute
Persistent Coreform IGA Job records indexed by name.
deleteJob
Instance method
Delete a named persistent Job from this namespace.
Job
Instance method
Create and persist a Coreform IGA Job in the MDB namespace.
renameJob
Instance method
Rename a Job that has not been prepared or submitted.
2.1.1
Attributes
Name
Type
Value
Description
jobs
object
—
Persistent Coreform IGA Job records indexed by name.
2.1.2
deleteJob
Defined on:
mdb.customData.CoreformIGA.deleteJob
Description:
Delete a named persistent Job from this namespace.
Call syntax
jobs.deleteJob(name
=
name)
Returns
Type
Description
None
—
Signatures
deleteJob(name:
str
)
->
None
Parameters
Name
Type
Required
Default
Description
name
str
Yes
—
Repository key of the Job to delete.
2.1.3
Job
Defined on:
mdb.customData.CoreformIGA.Job
Description:
Create and persist a Coreform IGA Job in the MDB namespace.
Call syntax
result
=
jobs.Job(spec
=
spec)
result
=
jobs.Job(name
=
name, model
=
model)
Returns
Type
Description
Job record
Newly created persistent Job record.
Signatures
Job(
*
, spec)
->
object
Job(
*
, name, model, description
=
''
, meshExecutable
=
''
, mpiexec
=
''
, numIgaCpus
=
1
, interopExecutable
=
''
, numCpus
=
1
, numGPUs
=
0
, multiprocessingMode
=
ParallelMode.THREADS, threadsPerProcess
=
None
, numDomains
=
None
, memory
=
90
, memoryUnits
=
MemoryUnits.PERCENTAGE, queue
=
None
, meshQueue
=
None
, interopQueue
=
None
, solveQueue
=
None
)
->
object
Parameters
For
Job(*, spec) -> object
:
Name
Type
Required
Default
Description
spec
CoreformIGAJobSpec or Mapping[str, object]
Yes
—
Complete portable specification or serialized specification mapping.
For
Job(*, name, model, description='', meshExecutable='', mpiexec='', numIgaCpus=1, interopExecutable='', numCpus=1, numGPUs=0, multiprocessingMode=ParallelMode.THREADS, threadsPerProcess=None, numDomains=None, memory=90, memoryUnits=MemoryUnits.PERCENTAGE, queue=None, meshQueue=None, interopQueue=None, solveQueue=None) -> object
:
Name
Type
Required
Default
Description
name
str
Yes
—
Unique, nonempty Job name.
model
str
Yes
—
Abaqus model used by the solve pipeline.
description
str
No
''
Optional Job description.
meshExecutable
str
No
''
Coreform IGA Mesh executable override.
mpiexec
str
No
''
MPI launcher override for the Mesh stage.
numIgaCpus
int
No
1
Number of CPUs for the Mesh stage.
interopExecutable
str
No
''
Abaqus interop executable override.
numCpus
int
No
1
Number of CPUs for the solve stage.
numGPUs
int
No
0
Number of GPUs for the solve stage.
multiprocessingMode
ParallelMode or str
No
ParallelMode.THREADS
Abaqus solver parallelization mode.
threadsPerProcess
int or None
No
None
Threads allocated to each solver process.
numDomains
int or None
No
None
Number of Abaqus solver domains.
memory
int
No
90
Abaqus solver memory limit.
memoryUnits
MemoryUnits or str
No
MemoryUnits.PERCENTAGE
Units used to interpret
memory
.
queue
Local, RemoteQueue, Mapping[str, object], or None
No
None
Default submission settings for every stage.
meshQueue
Local, RemoteQueue, Mapping[str, object], or None
No
None
Mesh-stage submission override.
interopQueue
Local, RemoteQueue, Mapping[str, object], or None
No
None
Interop-stage submission override.
solveQueue
Local, RemoteQueue, Mapping[str, object], or None
No
None
Solve-stage submission override.
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
If the specification, model, name, or execution settings are invalid.
2.1.4
renameJob
Defined on:
mdb.customData.CoreformIGA.renameJob
Description:
Rename a Job that has not been prepared or submitted.
Call syntax
jobs.renameJob(fromName
=
fromName, toName
=
toName)
Returns
Type
Description
None
—
Signatures
renameJob(fromName:
str
, toName:
str
)
->
None
Parameters
Name
Type
Required
Default
Description
fromName
str
Yes
—
Existing Job repository key.
toName
str
Yes
—
New unique Job name and repository key.
2.2
jobs[jobName]
Public path:
mdb.customData.CoreformIGA.jobs[jobName]
Available API members
Name
Kind
Description
name
Attribute
Repository key and Job name.
model
Attribute
Abaqus model used by the Job, when applicable.
description
Attribute
User-facing Job description.
meshExecutable
Attribute
Coreform IGA Mesh executable override.
mpiexec
Attribute
MPI launcher override for the Mesh stage.
numIgaCpus
Attribute
Number of CPUs assigned to the Mesh stage.
interopExecutable
Attribute
Abaqus interop executable override, when applicable.
numCpus
Attribute
Number of CPUs assigned to the solve stage.
numGPUs
Attribute
Number of GPUs assigned to an Abaqus solve stage.
multiprocessingMode
Attribute
Parallelization mode for an Abaqus solve stage.
threadsPerProcess
Attribute
Threads assigned to each Abaqus solver process.
numDomains
Attribute
Number of domains used by the Abaqus solver.
memory
Attribute
Memory limit for an Abaqus solve stage.
memoryUnits
Attribute
Units used to interpret the Abaqus solver memory limit.
toSpec
Instance method
Return the lossless portable specification for this Job record.
setValues
Instance method
Update this persistent Job through Abaqus-facing camelCase fields.
buildMesh
Instance method
Run the Mesh and interop stages for this Job.
writeInput
Instance method
Write the raw Abaqus input deck consumed by the interop stage.
submit
Instance method
Prepare and submit this Job through its configured pipeline.
waitForCompletion
Instance method
Block until a submitted Job finishes.
kill
Instance method
Terminate the running process associated with this Job.
getRuntimeStatus
Instance method
Return a snapshot of this Job’s current runtime status.
writeRunScript
Instance method
Write this Job pipeline as a runnable script file.
openOdb
Instance method
Open and return this Job’s Abaqus output database.
2.2.1
Attributes
Name
Type
Value
Description
name
str
—
Repository key and Job name.
model
Optional[str]
—
Abaqus model used by the Job, when applicable.
description
str
—
User-facing Job description.
meshExecutable
str
—
Coreform IGA Mesh executable override.
mpiexec
str
—
MPI launcher override for the Mesh stage.
numIgaCpus
int
—
Number of CPUs assigned to the Mesh stage.
interopExecutable
Optional[str]
—
Abaqus interop executable override, when applicable.
numCpus
int
—
Number of CPUs assigned to the solve stage.
numGPUs
Optional[int]
—
Number of GPUs assigned to an Abaqus solve stage.
multiprocessingMode
Optional[Union[ParallelMode, str]]
—
Parallelization mode for an Abaqus solve stage.
threadsPerProcess
Optional[int]
—
Threads assigned to each Abaqus solver process.
numDomains
Optional[int]
—
Number of domains used by the Abaqus solver.
memory
Optional[int]
—
Memory limit for an Abaqus solve stage.
memoryUnits
Optional[Union[MemoryUnits, str]]
—
Units used to interpret the Abaqus solver memory limit.
2.2.2
toSpec
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].toSpec
Description:
Return the lossless portable specification for this Job record.
Call syntax
result
=
job.toSpec()
Returns
Type
Description
CoreformIGAJobSpec
Return the lossless portable specification for this Job record.
Signatures
toSpec()
->
CoreformIGAJobSpec
Parameters
None.
2.2.3
setValues
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].setValues
Description:
Update this persistent Job through Abaqus-facing camelCase fields.
Call syntax
job.setValues(...)
Returns
Type
Description
None
The record is updated in place.
Signatures
setValues(
*
, model
=
..., description
=
..., meshExecutable
=
..., mpiexec
=
..., numIgaCpus
=
..., interopExecutable
=
..., numCpus
=
..., numGPUs
=
..., multiprocessingMode
=
..., threadsPerProcess
=
..., numDomains
=
..., memory
=
..., memoryUnits
=
..., meshQueue
=
..., interopQueue
=
..., solveQueue
=
...)
->
None
Parameters
Name
Type
Required
Default
Description
model
str
No
Current value
Replace the Abaqus model name.
description
str
No
Current value
Replace the Job description.
meshExecutable
str
No
Current value
Replace the Mesh executable override.
mpiexec
str
No
Current value
Replace the MPI launcher override.
numIgaCpus
int
No
Current value
Replace the Mesh-stage CPU count.
interopExecutable
str
No
Current value
Replace the interop executable override.
numCpus
int
No
Current value
Replace the solve-stage CPU count.
numGPUs
int
No
Current value
Replace the solve-stage GPU count.
multiprocessingMode
ParallelMode or str
No
Current value
Replace the solver parallelization mode.
threadsPerProcess
int or None
No
Current value
Replace the threads-per-process setting.
numDomains
int or None
No
Current value
Replace the Abaqus solver domain count.
memory
int
No
Current value
Replace the solver memory limit.
memoryUnits
MemoryUnits or str
No
Current value
Replace the solver memory units.
meshQueue
Local, RemoteQueue, Mapping[str, object], or None
No
Current value
Replace the Mesh-stage submission settings.
interopQueue
Local, RemoteQueue, Mapping[str, object], or None
No
Current value
Replace the interop-stage submission settings.
solveQueue
Local, RemoteQueue, Mapping[str, object], or None
No
Current value
Replace the solve-stage submission settings.
Raises
Type
Condition
ValidationError
If an update is unsupported or produces an invalid Job.
2.2.4
buildMesh
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].buildMesh
Description:
Run the Mesh and interop stages for this Job.
Call syntax
result
=
job.buildMesh()
Returns
Type
Description
JobExecutionResult
Run the Mesh and interop stages for this Job.
Signatures
buildMesh(workDirectory:
object
=
None
)
->
JobExecutionResult
Parameters
Name
Type
Required
Default
Description
workDirectory
object
No
None
Existing operation directory, or
None
to use the Job default.
2.2.5
writeInput
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].writeInput
Description:
Write the raw Abaqus input deck consumed by the interop stage.
Call syntax
job.writeInput()
Returns
Type
Description
None
—
Signatures
writeInput(workDirectory:
object
=
None
)
->
None
Parameters
Name
Type
Required
Default
Description
workDirectory
object
No
None
Existing operation directory, or
None
to use the Job default.
2.2.6
submit
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].submit
Description:
Prepare and submit this Job through its configured pipeline.
Call syntax
result
=
job.submit()
Returns
Type
Description
JobExecutionResult
Prepare and submit this Job through its configured pipeline.
Signatures
submit(workDirectory:
object
=
None
, overwriteExistingFiles:
bool
=
False
, consistencyChecking:
str
=
'ON'
)
->
JobExecutionResult
Parameters
Name
Type
Required
Default
Description
workDirectory
object
No
None
Existing operation directory, or
None
to use the Job default.
overwriteExistingFiles
bool
No
False
Whether recognized prior Job outputs may be replaced.
consistencyChecking
str
No
'ON'
Abaqus input-generation consistency-checking value.
2.2.7
waitForCompletion
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].waitForCompletion
Description:
Block until a submitted Job finishes.
Call syntax
job.waitForCompletion()
Returns
Type
Description
None
—
Signatures
waitForCompletion()
->
None
Parameters
None.
2.2.8
kill
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].kill
Description:
Terminate the running process associated with this Job.
Call syntax
job.kill()
Returns
Type
Description
None
—
Signatures
kill()
->
None
Parameters
None.
2.2.9
getRuntimeStatus
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].getRuntimeStatus
Description:
Return a snapshot of this Job’s current runtime status.
Call syntax
result
=
job.getRuntimeStatus()
Returns
Type
Description
dict[str, object]
Return a snapshot of this Job’s current runtime status.
Signatures
getRuntimeStatus()
->
dict
[
str
,
object
]
Parameters
None.
2.2.10
writeRunScript
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].writeRunScript
Description:
Write this Job pipeline as a runnable script file.
Call syntax
job.writeRunScript(fileName
=
fileName)
Returns
Type
Description
None
—
Signatures
writeRunScript(fileName:
object
, dialect:
str
=
'POWERSHELL'
, prepareInputs:
bool
=
True
)
->
None
Parameters
Name
Type
Required
Default
Description
fileName
object
Yes
—
Destination script path accepted by the Abaqus file API.
dialect
str
No
'POWERSHELL'
Uppercase script dialect accepted by
render_script
.
prepareInputs
bool
No
True
Whether to write required Abaqus inputs before rendering the script.
2.2.11
openOdb
Defined on:
mdb.customData.CoreformIGA.jobs[jobName].openOdb
Description:
Open and return this Job’s Abaqus output database.
Call syntax
result
=
job.openOdb()
Returns
Type
Description
object
Open and return this Job’s Abaqus output database.
Signatures
openOdb()
->
object
Parameters
None.
