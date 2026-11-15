# Consistent unit systems – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/theory/consistent-unit-systems.html

---

Theory
Consistent unit systems
Consistent unit systems
Mechanical quantities and dimensions in SI, CGS, MMTS, IPS, and FPS.
1
Consistent unit systems
Neither Abaqus nor Coreform IGA for Abaqus prescribes units. They operate on the numerical values provided by the analyst, so every input and every interpreted result must belong to one self-consistent unit system. A length in millimeters, for example, cannot be combined with a modulus in pascals unless the modulus is first converted to the corresponding force-per-square-millimeter unit.
This chapter uses
\(M\)
for mass,
\(L\)
for length,
\(T\)
for time, and
\(\Theta\)
for thermodynamic temperature. A dimensionless quantity is written
\(1\)
. The
defining units for mechanics
tables identify the independent units used to construct the remaining quantities; this is deliberately narrower than the formal list of base units in a standards document. The
derived units
tables show a common coherent expression for each derived quantity.
Energy and moment have the same dimensionality,
\(M L^2 T^{-2}\)
, but use different unit notation to preserve their distinct physical meanings: for example,
\(\mathrm{J}\)
for energy and
\(\mathrm{N\,m}\)
for moment.
The
BIPM SI Brochure
and the
NIST Guide to the SI
provide the standards basis for the SI definitions and conversions used here.
1.1
SI
For mechanics, the International System of Units provides a coherent meter-kilogram-second system. It is usually the clearest choice when geometry is modeled in meters.
1.1.1
Defining units for mechanics
Quantity
Unit
Dimensionality
Mass
\(\mathrm{kg}\)
\(M\)
Length
\(\mathrm{m}\)
\(L\)
Time
\(\mathrm{s}\)
\(T\)
Temperature
\(\mathrm{K}\)
\(\Theta\)
1.1.2
Derived units
Quantity
Unit
Dimensionality
Area
\(\mathrm{m^2}\)
\(L^2\)
Volume
\(\mathrm{m^3}\)
\(L^3\)
Velocity
\(\mathrm{m/s}\)
\(L T^{-1}\)
Acceleration
\(\mathrm{m/s^2}\)
\(L T^{-2}\)
Frequency
\(\mathrm{Hz} = \mathrm{s^{-1}}\)
\(T^{-1}\)
Force
\(\mathrm{N}\)
\(M L T^{-2}\)
Energy or work
\(\mathrm{J}\)
\(M L^2 T^{-2}\)
Moment or torque
\(\mathrm{N\,m}\)
\(M L^2 T^{-2}\)
Power
\(\mathrm{W}\)
\(M L^2 T^{-3}\)
Pressure or stress
\(\mathrm{Pa} = \mathrm{N/m^2}\)
\(M L^{-1} T^{-2}\)
Mass density
\(\mathrm{kg/m^3}\)
\(M L^{-3}\)
Body force per volume
\(\mathrm{N/m^3}\)
\(M L^{-2} T^{-2}\)
1.2
CGS
The centimeter-gram-second system is common in some scientific literature. Its named mechanical units include the dyne for force, the erg for energy, and the barye for pressure.
1.2.1
Defining units for mechanics
Quantity
Unit
Dimensionality
Mass
\(\mathrm{g}\)
\(M\)
Length
\(\mathrm{cm}\)
\(L\)
Time
\(\mathrm{s}\)
\(T\)
Temperature
\(\mathrm{K}\)
\(\Theta\)
1.2.2
Derived units
Quantity
Unit
Dimensionality
Area
\(\mathrm{cm^2}\)
\(L^2\)
Volume
\(\mathrm{cm^3}\)
\(L^3\)
Velocity
\(\mathrm{cm/s}\)
\(L T^{-1}\)
Acceleration
\(\mathrm{cm/s^2}\)
\(L T^{-2}\)
Frequency
\(\mathrm{s^{-1}}\)
\(T^{-1}\)
Force
\(\mathrm{dyn}\)
\(M L T^{-2}\)
Energy or work
\(\mathrm{erg}\)
\(M L^2 T^{-2}\)
Moment or torque
\(\mathrm{dyn\,cm}\)
\(M L^2 T^{-2}\)
Power
\(\mathrm{erg/s}\)
\(M L^2 T^{-3}\)
Pressure or stress
\(\mathrm{Ba} = \mathrm{dyn/cm^2}\)
\(M L^{-1} T^{-2}\)
Mass density
\(\mathrm{g/cm^3}\)
\(M L^{-3}\)
Body force per volume
\(\mathrm{dyn/cm^3}\)
\(M L^{-2} T^{-2}\)
1.3
MMTS
The millimeter-tonne-second system is especially convenient for mechanical finite-element models whose geometry is dimensioned in millimeters:
\[
1\ \mathrm{tonne}\,\frac{\mathrm{mm}}{\mathrm{s^2}}
=
1\ \mathrm{N}.
\]
Consequently, stress is naturally expressed in
\(\mathrm{N/mm^2}=\mathrm{MPa}\)
. Gravitational acceleration near Earth’s surface is approximately
\(9810\ \mathrm{mm/s^2}\)
.
1.3.1
Defining units for mechanics
Quantity
Unit
Dimensionality
Mass
\(\mathrm{tonne}\)
\(M\)
Length
\(\mathrm{mm}\)
\(L\)
Time
\(\mathrm{s}\)
\(T\)
Temperature
\(\mathrm{K}\)
\(\Theta\)
1.3.2
Derived units
Quantity
Unit
Dimensionality
Area
\(\mathrm{mm^2}\)
\(L^2\)
Volume
\(\mathrm{mm^3}\)
\(L^3\)
Velocity
\(\mathrm{mm/s}\)
\(L T^{-1}\)
Acceleration
\(\mathrm{mm/s^2}\)
\(L T^{-2}\)
Frequency
\(\mathrm{s^{-1}}\)
\(T^{-1}\)
Force
\(\mathrm{N}\)
\(M L T^{-2}\)
Energy or work
\(\mathrm{mJ} = \mathrm{N\,mm}\)
\(M L^2 T^{-2}\)
Moment or torque
\(\mathrm{N\,mm}\)
\(M L^2 T^{-2}\)
Power
\(\mathrm{mW} = \mathrm{N\,mm/s}\)
\(M L^2 T^{-3}\)
Pressure or stress
\(\mathrm{MPa} = \mathrm{N/mm^2}\)
\(M L^{-1} T^{-2}\)
Mass density
\(\mathrm{tonne/mm^3}\)
\(M L^{-3}\)
Body force per volume
\(\mathrm{N/mm^3}\)
\(M L^{-2} T^{-2}\)
1.4
IPS
The coherent inch-pound-second system used here takes the pound-force as the familiar force unit. Its corresponding mass unit is the
slinch
:
\[
1\ \mathrm{slinch}
\equiv
1\ \frac{\mathrm{lbf\,s^2}}{\mathrm{in}}.
\]
1.4.1
Defining units for mechanics
Quantity
Unit
Dimensionality
Mass
\(\mathrm{slinch}\)
\(M\)
Length
\(\mathrm{in}\)
\(L\)
Time
\(\mathrm{s}\)
\(T\)
Temperature
\({}^\circ\mathrm{R}\)
\(\Theta\)
1.4.2
Derived units
Quantity
Unit
Dimensionality
Area
\(\mathrm{in^2}\)
\(L^2\)
Volume
\(\mathrm{in^3}\)
\(L^3\)
Velocity
\(\mathrm{in/s}\)
\(L T^{-1}\)
Acceleration
\(\mathrm{in/s^2}\)
\(L T^{-2}\)
Frequency
\(\mathrm{s^{-1}}\)
\(T^{-1}\)
Force
\(\mathrm{lbf}\)
\(M L T^{-2}\)
Energy or work
\(\mathrm{lbf\,in}\)
\(M L^2 T^{-2}\)
Moment or torque
\(\mathrm{lbf\,in}\)
\(M L^2 T^{-2}\)
Power
\(\mathrm{lbf\,in/s}\)
\(M L^2 T^{-3}\)
Pressure or stress
\(\mathrm{psi} = \mathrm{lbf/in^2}\)
\(M L^{-1} T^{-2}\)
Mass density
\(\mathrm{slinch/in^3}\)
\(M L^{-3}\)
Body force per volume
\(\mathrm{lbf/in^3}\)
\(M L^{-2} T^{-2}\)
1.5
FPS
The coherent foot-pound-second system likewise uses pound-force for force, but its mass unit is the
slug
:
\[
1\ \mathrm{slug}
\equiv
1\ \frac{\mathrm{lbf\,s^2}}{\mathrm{ft}},
\qquad
1\ \mathrm{slinch}=12\ \mathrm{slug}.
\]
1.5.1
Defining units for mechanics
Quantity
Unit
Dimensionality
Mass
\(\mathrm{slug}\)
\(M\)
Length
\(\mathrm{ft}\)
\(L\)
Time
\(\mathrm{s}\)
\(T\)
Temperature
\({}^\circ\mathrm{R}\)
\(\Theta\)
1.5.2
Derived units
Quantity
Unit
Dimensionality
Area
\(\mathrm{ft^2}\)
\(L^2\)
Volume
\(\mathrm{ft^3}\)
\(L^3\)
Velocity
\(\mathrm{ft/s}\)
\(L T^{-1}\)
Acceleration
\(\mathrm{ft/s^2}\)
\(L T^{-2}\)
Frequency
\(\mathrm{s^{-1}}\)
\(T^{-1}\)
Force
\(\mathrm{lbf}\)
\(M L T^{-2}\)
Energy or work
\(\mathrm{ft\,lbf}\)
\(M L^2 T^{-2}\)
Moment or torque
\(\mathrm{lbf\,ft}\)
\(M L^2 T^{-2}\)
Power
\(\mathrm{ft\,lbf/s}\)
\(M L^2 T^{-3}\)
Pressure or stress
\(\mathrm{psf} = \mathrm{lbf/ft^2}\)
\(M L^{-1} T^{-2}\)
Mass density
\(\mathrm{slug/ft^3}\)
\(M L^{-3}\)
Body force per volume
\(\mathrm{lbf/ft^3}\)
\(M L^{-2} T^{-2}\)
1.6
Common mistakes
1.6.1
On the pound
A common error—perhaps the most common error made by new analysts—is entering density incorrectly in the IPS unit system. This likely stems from confusion about what “lb” means: is it weight, mass, or force? The situation is not helped by the many names inherited by the Imperial and United States customary systems. In contrast to metric systems, their unit names were adopted from many sources. In the pre-industrial era, several conveyed useful pragmatic information.
1
Common names for weights, forces, and masses include pound (lb), poundal (pdl), pound-force (lbf), pound-mass (lbm), short and long ton, slug, and slinch. Whenever one of these terms appears, take extra care to determine its dimensionality and convert it into the selected coherent unit system. The
NIST definition of pound-force
, for example, explicitly relates it to the force produced by standard gravity on one avoirdupois pound.
The first example starts with a density expressed using the FPS mass unit and converts the mass unit from slugs to slinches:
\[
\begin{aligned}
\rho
  &= 0.00898\ \frac{\mathrm{slug}}{\mathrm{in}^3},\\
1\ \mathrm{slug}
  &\equiv \frac{1\ \mathrm{lbf}}{1\ \mathrm{ft/s^2}},\\
1\ \mathrm{slinch}
  &\equiv \frac{1\ \mathrm{lbf}}{1\ \mathrm{in/s^2}},\\
\rho
  &= 0.00898\ \frac{\mathrm{slug}}{\mathrm{in}^3}
     \left(\frac{1\ \mathrm{ft}}{12\ \mathrm{in}}\right)\\
  &= 0.000748\ \frac{\mathrm{lbf\,s^2}}{\mathrm{in}^4}\\
  &= 0.000748\ \frac{\mathrm{slinch}}{\mathrm{in}^3}.
\end{aligned}
\]
The second example begins with the familiar tabulated steel density
\(0.289\ \mathrm{lbm/in^3}\)
. Using standard gravitational acceleration rounded to
\(32.2\ \mathrm{ft/s^2}\)
gives the coherent mass equivalence
\[
1\ \mathrm{lbm}
\equiv
\frac{1\ \mathrm{lbf}}{32.2\ \mathrm{ft/s^2}}
=
\frac{1\ \mathrm{lbf\,s^2}}{32.2\ \mathrm{ft}}.
\]
The complete conversion is therefore
\[
\begin{aligned}
\rho
  &= 0.289\ \frac{\mathrm{lbm}}{\mathrm{in}^3}
     \left(
       \frac{1\ \mathrm{lbf\,s^2}}
            {32.2\ \mathrm{lbm\,ft}}
     \right)
     \left(\frac{1\ \mathrm{ft}}{12\ \mathrm{in}}\right)\\
  &= 0.000748\ \frac{\mathrm{lbf\,s^2}}{\mathrm{in}^4}\\
  &= 0.000748\ \frac{\mathrm{slinch}}{\mathrm{in}^3}.
\end{aligned}
\]
The number
\(0.289\)
is not itself wrong; the mistake is entering it as though
\(\mathrm{lbm/in^3}\)
were already the coherent IPS mass-density unit. If a load is specified as weight density instead, write and convert it explicitly as force per volume,
\(\mathrm{lbf/in^3}\)
, rather than silently treating it as mass density.
Footnotes
For instance, an
acre
was approximately the area a team of oxen could plow in a day, and a
ton of refrigeration
was based on the heat-transfer rate required to melt a short ton (2000 lb) of ice at its freezing point in 24 hours.
↩︎
