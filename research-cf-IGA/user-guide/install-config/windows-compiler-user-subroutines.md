# Windows Compiler for User Subroutines – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/install-config/windows-compiler-user-subroutines.html

---

Installation & Configuration
Windows Compiler for User Subroutines
Windows Compiler for User Subroutines
To compile Fortran user subroutines for Coreform IGA for Abaqus on Windows, you need Intel Fortran Compiler (
ifx
) and Microsoft Visual Studio Build Tools. This compiler setup is optional when you use only the libraries included with Coreform IGA.
1
Software and Downloads
The following combination passed our Windows user-subroutine test with Abaqus 2026 Hot Fix 5 and Coreform IGA 2026.11:
Software
Tested version
Where to get it
Intel Fortran Compiler (
ifx
)
2026.2.0
Windows offline installer
Microsoft Visual Studio Build Tools
2026, version 18.11.1
Visual Studio Build Tools download
In the Visual Studio Installer, select
Desktop development with C++
, including the x64 MSVC tools and a Windows SDK. Our test used MSVC
14.51.36231
and Windows SDK
10.0.26100.0
. Build Tools supplies the linker and Windows libraries;
ifx
compiles the Fortran source.
Note
These are tested versions, not a minimum-version range or a Dassault Systèmes qualification statement.
2
Environment Prerequisites
Submit the job from
Windows Command Prompt (
cmd.exe
)
with the compiler environment initialized in that same terminal. For the tested Build Tools and Intel oneAPI installation paths:
call
"C:\Program
Files (x86)\Microsoft Visual Studio\18\BuildTools\Common7\Tools\VsDevCmd.bat
" -arch=amd64 -host_arch=amd64
call
"C:\Program
Files (x86)\Intel\oneAPI\setvars.bat
"
These scripts configure
PATH
, include directories, and library directories for the current session. Adjust the paths to your installation; other Intel layouts may use a different
environment initialization script
.
The terminal must be able to find
ifx
,
cl
,
link
, and
abaqus
. Use these commands to verify:
where.exe ifx
ifx -
-version
where.exe cl
where.exe link
where.exe abaqus
Also make the installed Coreform
bin
directory available on
PATH
, for example:
set
"
PATH
=
C:\Program Files\Coreform IGA 2026.11\bin;
%
PATH
%
"
Any directory containing customer DLLs loaded through
uniquelibs
must also be discoverable by the solver. Run from a writable job directory containing the prepared Coreform input files, with access to the Fortran source and the required
licenses
.
Abaqus must be configured to select
ifx
in its
compile_fortran
setting in
abaqus_v6.env
. Our tested configuration replaced
ifort
with
ifx
, removed
/Qprec-sqrt
and
/Qfp-stack-check
, and retained the other compiler options. Keep any job-specific library staging separate from ordinary verification commands.
3
User Subroutine Submit Command
Once compiler selection and library loading are configured, the command has this form:
abaqus job=MyJob input=MyJob.inp user=my_subroutine.for cpus=1 uniquelibs interactive
user=
supplies the Fortran source for Abaqus to compile and link during submission.
uniquelibs
enables the individual user-subroutine libraries used by Coreform IGA. Use the basename of the prepared Coreform job, including its matching
.sql
database;
cpus=1
matches the tested configuration.
