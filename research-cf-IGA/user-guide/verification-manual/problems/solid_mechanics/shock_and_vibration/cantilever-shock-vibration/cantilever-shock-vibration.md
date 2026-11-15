# Cantilever beam shock and vibration – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/shock_and_vibration/cantilever-shock-vibration/cantilever-shock-vibration.html

---

Verification Manual
Solid mechanics
Shock and vibration
Cantilever beam shock and vibration
Cantilever beam shock and vibration
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify mode-based linear dynamics on a cantilever beam: a frequency extraction, a random response to a base-acceleration power spectral density, and a peak response to a base-acceleration shock spectrum.
Coreform IGA elements are Abaqus user elements. Abaqus/Standard will extract eigenmodes for them, but it does not post-process a random response or a response spectrum from those elements, so for an immersed model the frequency extraction is the last result the solver produces. Coreform computes the two remaining procedures from that modal basis instead.
That arrangement has two things to verify, and this study checks both.
The implementation.
On the body-fitted mesh Abaqus computes all three steps itself, so the Coreform random response and response spectrum can be compared against the solver’s own answers on an identical modal basis. Any difference here is implementation error alone.
The discretization.
The immersed Coreform IGA modal basis, post-processed by Coreform, is compared against the body-fitted Abaqus analysis end to end. Any difference here is the difference between two meshes.
2
Geometry
Table 1: Cantilever beam geometry.
Parameter
Value
Width
\(b\)
\(0.1\ m\)
Height
\(h\)
\(0.1\ m\)
Length
\(L\)
\(1\ m\)
The beam is the same one used by the static
cantilever beam bending
study, so the two chapters describe one physical beam under different analyses. The section is square, which makes every bending mode a degenerate pair and puts a torsional mode inside the analysis band.
3
Material
Table 2: Linear elastic material.
Property
Value
Young’s modulus
\(E\)
\(200 \times 10^{9}\ Pa\)
Poisson’s ratio
\(\nu\)
0
Mass density
\(\rho\)
\(8.00 \times 10^{3}\ kg/m^3^\)
Poisson’s ratio is zero so that the three-dimensional solid reproduces the one-dimensional beam theory the accepted frequencies come from.
4
Loading & boundary conditions
The root face at
\(z = 0\)
is fully clamped. Both dynamic steps excite the structure through a primary base motion in the global
\(y\)
direction, the weak bending direction of the section.
Table 3: Excitation and analysis settings, shared by both discretizations.
Quantity
Value
Analysis band
\(10\ Hz\)
to
\(2.00 \times 10^{3}\ Hz\)
Modes extracted
\(20\)
Modal damping
\(\xi\)
20.0 ^{-3}
Base acceleration PSD
10.0 ^{-3} g²/Hz, flat
Frequency points per interval
\(20\)
Bias parameter
3
Modal summation
SRSS
The random response is driven by a flat base-acceleration spectral density across the band, the shape a vibration qualification test specifies. The response spectrum is driven by a base-acceleration shock spectrum that rises to a plateau and rolls off.
Table 4: Base-acceleration shock spectrum, tabulated at the modal damping ratio.
Magnitude (m/s²)
Frequency (Hz)
20
10
200
100
200
600
60
2.00 ^{3}
5
Mesh & discretization
Two discretizations of the same beam. The body-fitted case is a structured C3D8 mesh that Abaqus owns end to end. The immersed case is a Coreform IGA mesh on a rectilinear background grid, trimmed to the beam.
Table 5: Discretizations compared.
Case
Discretization
Element size
Dynamics computed by
abaqus_h0p1
Body-fitted Abaqus (C3D8)
\(50.0 \times 10^{-3}\ m\)
Abaqus/Standard
immersed_h0p1
Immersed Coreform IGA
\(25.0 \times 10^{-3}\ m\)
Coreform, from the modal basis
The immersed case uses a finer background cell than the body-fitted element size. A linear immersed basis needs four cells across the section to resolve bending; at two cells the fundamental frequency comes out several percent stiff.
6
Reference & accepted solutions
The fundamental bending frequency of a clamped-free Euler-Bernoulli beam is
\[
f_1 = \frac{\beta_1^2}{2 \pi L^2} \sqrt{\frac{E I}{\rho A}}, \qquad \beta_1 = 1.8751
\]
with
\(I = b h^3 / 12\)
the second moment of area and
\(A = b h\)
the cross-sectional area. A three-dimensional solid sits slightly above this value because beam theory omits shear flexibility.
The random response and response spectrum quantities have no closed form, so the body-fitted Abaqus analysis is the accepted reference for them. Because Abaqus also computes both procedures on the body-fitted modal basis, the Coreform implementations are additionally checked against the solver directly, which is a stronger statement than the mesh comparison.
7
Results
7.1
Implementation against Abaqus
On a given modal basis the Coreform procedures and Abaqus solve the same problem, so they must agree to solver precision. Abaqus produces displacement-based dynamics output even for user elements, so this comparison is available on both discretizations.
The reported quantity is the component the base motion excites, which is what this chapter’s results are about and the only component whose magnitude is set by the excitation rather than by cross-coupling.
Table 6: Coreform against Abaqus/Standard on an identical modal basis, worst relative difference in the excited component. Accepted below 5.00 ^{-6}, and the immersed error is additionally required to stay within 4 times the body-fitted one.
Procedure
Body-fitted
Immersed
Ratio
Random response
212 ^{-9}
213 ^{-9}
1
Response spectrum
97.1 ^{-9}
121 ^{-9}
1.25
Download data (CSV)
Both discretizations agree with the solver at the precision of the single-precision fields Abaqus writes, and they agree with each other: the immersed modal basis costs the implementation nothing. That parity is asserted by the study rather than left to the reader, because it is the claim that lets the immersed workflow stand in for a solver-computed one.
Note
Why the comparison is per component
Judging the whole displacement field instead would let a small secondary component decide the number. The axial component here peaks near eight percent of the field peak, so points that clear a field-wide floor still sit in that component’s own round-off, and the worst of them grows with however many output points a mesh happens to carry. The immersed mesh has roughly eight times as many output points as the body-fitted one, which on a field-wide measure inflates its figure to around
\(1.4 \times 10^{-6}\)
while the excited component stays at
\(210 \times 10^{-9}\)
on both. The field-wide numbers are recorded in the downloadable data for transparency, but nothing is gated on them.
7.2
Eigenfrequencies
Table 7: Eigenfrequencies of the two discretizations, with the participation factor in the excitation direction. A mode participates when its factor reaches 1% of the largest.
Mode
Body-fitted (Hz)
Immersed (Hz)
Relative difference
\(|\Gamma_y|\)
Participating
\(1\)
81.6
81.8
2.10 ^{-3}
0.116
yes
\(2\)
81.6
81.8
2.25 ^{-3}
6.99
yes
\(3\)
492
494
5.20 ^{-3}
0.2
yes
\(4\)
492
494
5.50 ^{-3}
3.89
yes
\(5\)
722
822
0.139
816 ^{-15}
no
\(6\)
1.25 ^{3}
1.25 ^{3}
292 ^{-6}
31.4 ^{-15}
no
\(7\)
1.30 ^{3}
1.32 ^{3}
8.46 ^{-3}
9.96 ^{-3}
no
\(8\)
1.30 ^{3}
1.32 ^{3}
8.75 ^{-3}
2.29
yes
Download data (CSV)
\[
\newcommand{\analyticalFOne}{80.8}
\newcommand{\nativeFOne}{81.6}
\newcommand{\immersedFOne}{81.8}
\]
Against the closed form, the body-fitted mesh gives
\(81\.6\\\\\\\\ Hz\)
and the immersed mesh
\(81\.8\\\\\\\\ Hz\)
, where Euler-Bernoulli theory predicts
\(80\.8\\\\\\\\ Hz\)
. Both sit above the beam-theory value, as a shear-flexible solid should.
7.3
Response
Table 8: Peak response in the excitation direction, body-fitted Abaqus against immersed Coreform. Accepted below 50.0 ^{-3}.
Quantity
Body-fitted Abaqus
Immersed Coreform
Relative difference
RMS displacement, random response
\(340 \times 10^{-6}\ m\)
\(346 \times 10^{-6}\ m\)
17.9 ^{-3}
Peak displacement, shock spectrum
\(965 \times 10^{-6}\ m\)
\(983 \times 10^{-6}\ m\)
18.3 ^{-3}
Download data (CSV)
8
Discussion
The two halves of the arrangement behave differently, and separating them is the point of this study.
The implementation is exact to solver precision.
Given the same modal basis, the Coreform random response and response spectrum reproduce Abaqus to a few parts in ten million, which is the precision at which Abaqus stores its fields. Nothing in the combination introduces error.
The discretizations differ by the amount two meshes differ.
The participating eigenfrequencies agree to under one percent, and the RMS and peak responses to under two percent, on a comparison between a body-fitted linear hexahedral mesh and a trimmed immersed spline basis.
Not every mode needs to agree for the response to agree.
The square section puts a torsional mode inside the analysis band, and a linear immersed basis represents torsional warping far less accurately than bending: that mode is stiff by an order of ten percent. It carries no participation in a transverse base motion, so it never reaches the response, and the study gates on the participating modes for that reason. The full spectrum is tabulated above so the departure is visible rather than hidden, and it is a useful reminder that agreement mode by mode is neither necessary nor sufficient for agreement in the response.
RMS values depend on the frequency sweep.
The random response integrates a spectral density that is sharply peaked at lightly damped resonances, so the reported RMS is only as converged as the sweep beneath it. At the settings used here both codes integrate the same points and agree with each other, but neither is the converged integral. A production run should raise the number of points per interval until the RMS stops moving.
9
References
10
Appendix
10.1
Capabilities exercised
No qualifying capability evidence was generated for these run artifacts.
10.2
Supported Abaqus keywords
Note
Analysis deck keywords
*Heading
*Preprint, echo=NO, model=NO, history=NO, contact=NO
*Part, name=<part>
*Node
1,  0.100000001,  0.100000001,           1.
2,  0.100000001, 0.0500000007,           1.
3,  0.100000001,           0.,           1.
4,  0.100000001,  0.100000001,  0.949999988
5,  0.100000001, 0.0500000007,  0.949999988
6,  0.100000001,           0.,  0.949999988
7,  0.100000001,  0.100000001,  0.899999976
8,  0.100000001, 0.0500000007,  0.899999976
9,  0.100000001,           0.,  0.899999976
10,  0.100000001,  0.100000001,  0.850000024
11,  0.100000001, 0.0500000007,  0.850000024
12,  0.100000001,           0.,  0.850000024
13,  0.100000001,  0.100000001,  0.800000012
14,  0.100000001, 0.0500000007,  0.800000012
15,  0.100000001,           0.,  0.800000012
16,  0.100000001,  0.100000001,         0.75
17,  0.100000001, 0.0500000007,         0.75
18,  0.100000001,           0.,         0.75
19,  0.100000001,  0.100000001,  0.699999988
20,  0.100000001, 0.0500000007,  0.699999988
21,  0.100000001,           0.,  0.699999988
22,  0.100000001,  0.100000001,  0.649999976
23,  0.100000001, 0.0500000007,  0.649999976
24,  0.100000001,           0.,  0.649999976
25,  0.100000001,  0.100000001,  0.600000024
26,  0.100000001, 0.0500000007,  0.600000024
27,  0.100000001,           0.,  0.600000024
28,  0.100000001,  0.100000001,  0.550000012
29,  0.100000001, 0.0500000007,  0.550000012
30,  0.100000001,           0.,  0.550000012
31,  0.100000001,  0.100000001,          0.5
32,  0.100000001, 0.0500000007,          0.5
33,  0.100000001,           0.,          0.5
34,  0.100000001,  0.100000001,  0.449999988
35,  0.100000001, 0.0500000007,  0.449999988
36,  0.100000001,           0.,  0.449999988
37,  0.100000001,  0.100000001,  0.400000006
38,  0.100000001, 0.0500000007,  0.400000006
39,  0.100000001,           0.,  0.400000006
40,  0.100000001,  0.100000001,  0.349999994
41,  0.100000001, 0.0500000007,  0.349999994
42,  0.100000001,           0.,  0.349999994
43,  0.100000001,  0.100000001,  0.300000012
44,  0.100000001, 0.0500000007,  0.300000012
45,  0.100000001,           0.,  0.300000012
46,  0.100000001,  0.100000001,         0.25
47,  0.100000001, 0.0500000007,         0.25
48,  0.100000001,           0.,         0.25
49,  0.100000001,  0.100000001,  0.200000003
50,  0.100000001, 0.0500000007,  0.200000003
51,  0.100000001,           0.,  0.200000003
52,  0.100000001,  0.100000001,  0.150000006
53,  0.100000001, 0.0500000007,  0.150000006
54,  0.100000001,           0.,  0.150000006
55,  0.100000001,  0.100000001,  0.100000001
56,  0.100000001, 0.0500000007,  0.100000001
57,  0.100000001,           0.,  0.100000001
58,  0.100000001,  0.100000001, 0.0500000007
59,  0.100000001, 0.0500000007, 0.0500000007
60,  0.100000001,           0., 0.0500000007
61,  0.100000001,  0.100000001,           0.
62,  0.100000001, 0.0500000007,           0.
63,  0.100000001,           0.,           0.
64, 0.0500000007,  0.100000001,           1.
65, 0.0500000007, 0.0500000007,           1.
66, 0.0500000007,           0.,           1.
67, 0.0500000007,  0.100000001,  0.949999988
68, 0.0500000007, 0.0500000007,  0.949999988
69, 0.0500000007,           0.,  0.949999988
70, 0.0500000007,  0.100000001,  0.899999976
71, 0.0500000007, 0.0500000007,  0.899999976
72, 0.0500000007,           0.,  0.899999976
73, 0.0500000007,  0.100000001,  0.850000024
74, 0.0500000007, 0.0500000007,  0.850000024
75, 0.0500000007,           0.,  0.850000024
76, 0.0500000007,  0.100000001,  0.800000012
77, 0.0500000007, 0.0500000007,  0.800000012
78, 0.0500000007,           0.,  0.800000012
79, 0.0500000007,  0.100000001,         0.75
80, 0.0500000007, 0.0500000007,         0.75
81, 0.0500000007,           0.,         0.75
82, 0.0500000007,  0.100000001,  0.699999988
83, 0.0500000007, 0.0500000007,  0.699999988
84, 0.0500000007,           0.,  0.699999988
85, 0.0500000007,  0.100000001,  0.649999976
86, 0.0500000007, 0.0500000007,  0.649999976
87, 0.0500000007,           0.,  0.649999976
88, 0.0500000007,  0.100000001,  0.600000024
89, 0.0500000007, 0.0500000007,  0.600000024
90, 0.0500000007,           0.,  0.600000024
91, 0.0500000007,  0.100000001,  0.550000012
92, 0.0500000007, 0.0500000007,  0.550000012
93, 0.0500000007,           0.,  0.550000012
94, 0.0500000007,  0.100000001,          0.5
95, 0.0500000007, 0.0500000007,          0.5
96, 0.0500000007,           0.,          0.5
97, 0.0500000007,  0.100000001,  0.449999988
98, 0.0500000007, 0.0500000007,  0.449999988
99, 0.0500000007,           0.,  0.449999988
100, 0.0500000007,  0.100000001,  0.400000006
101, 0.0500000007, 0.0500000007,  0.400000006
102, 0.0500000007,           0.,  0.400000006
103, 0.0500000007,  0.100000001,  0.349999994
104, 0.0500000007, 0.0500000007,  0.349999994
105, 0.0500000007,           0.,  0.349999994
106, 0.0500000007,  0.100000001,  0.300000012
107, 0.0500000007, 0.0500000007,  0.300000012
108, 0.0500000007,           0.,  0.300000012
109, 0.0500000007,  0.100000001,         0.25
110, 0.0500000007, 0.0500000007,         0.25
111, 0.0500000007,           0.,         0.25
112, 0.0500000007,  0.100000001,  0.200000003
113, 0.0500000007, 0.0500000007,  0.200000003
114, 0.0500000007,           0.,  0.200000003
115, 0.0500000007,  0.100000001,  0.150000006
116, 0.0500000007, 0.0500000007,  0.150000006
117, 0.0500000007,           0.,  0.150000006
118, 0.0500000007,  0.100000001,  0.100000001
119, 0.0500000007, 0.0500000007,  0.100000001
120, 0.0500000007,           0.,  0.100000001
121, 0.0500000007,  0.100000001, 0.0500000007
122, 0.0500000007, 0.0500000007, 0.0500000007
123, 0.0500000007,           0., 0.0500000007
124, 0.0500000007,  0.100000001,           0.
125, 0.0500000007, 0.0500000007,           0.
126, 0.0500000007,           0.,           0.
127,           0.,  0.100000001,           1.
128,           0., 0.0500000007,           1.
129,           0.,           0.,           1.
130,           0.,  0.100000001,  0.949999988
131,           0., 0.0500000007,  0.949999988
132,           0.,           0.,  0.949999988
133,           0.,  0.100000001,  0.899999976
134,           0., 0.0500000007,  0.899999976
135,           0.,           0.,  0.899999976
136,           0.,  0.100000001,  0.850000024
137,           0., 0.0500000007,  0.850000024
138,           0.,           0.,  0.850000024
139,           0.,  0.100000001,  0.800000012
140,           0., 0.0500000007,  0.800000012
141,           0.,           0.,  0.800000012
142,           0.,  0.100000001,         0.75
143,           0., 0.0500000007,         0.75
144,           0.,           0.,         0.75
145,           0.,  0.100000001,  0.699999988
146,           0., 0.0500000007,  0.699999988
147,           0.,           0.,  0.699999988
148,           0.,  0.100000001,  0.649999976
149,           0., 0.0500000007,  0.649999976
150,           0.,           0.,  0.649999976
151,           0.,  0.100000001,  0.600000024
152,           0., 0.0500000007,  0.600000024
153,           0.,           0.,  0.600000024
154,           0.,  0.100000001,  0.550000012
155,           0., 0.0500000007,  0.550000012
156,           0.,           0.,  0.550000012
157,           0.,  0.100000001,          0.5
158,           0., 0.0500000007,          0.5
159,           0.,           0.,          0.5
160,           0.,  0.100000001,  0.449999988
161,           0., 0.0500000007,  0.449999988
162,           0.,           0.,  0.449999988
163,           0.,  0.100000001,  0.400000006
164,           0., 0.0500000007,  0.400000006
165,           0.,           0.,  0.400000006
166,           0.,  0.100000001,  0.349999994
167,           0., 0.0500000007,  0.349999994
168,           0.,           0.,  0.349999994
169,           0.,  0.100000001,  0.300000012
170,           0., 0.0500000007,  0.300000012
171,           0.,           0.,  0.300000012
172,           0.,  0.100000001,         0.25
173,           0., 0.0500000007,         0.25
174,           0.,           0.,         0.25
175,           0.,  0.100000001,  0.200000003
176,           0., 0.0500000007,  0.200000003
177,           0.,           0.,  0.200000003
178,           0.,  0.100000001,  0.150000006
179,           0., 0.0500000007,  0.150000006
180,           0.,           0.,  0.150000006
181,           0.,  0.100000001,  0.100000001
182,           0., 0.0500000007,  0.100000001
183,           0.,           0.,  0.100000001
184,           0.,  0.100000001, 0.0500000007
185,           0., 0.0500000007, 0.0500000007
186,           0.,           0., 0.0500000007
187,           0.,  0.100000001,           0.
188,           0., 0.0500000007,           0.
189,           0.,           0.,           0.
*Element, type=C3D8
1,  64,  65,  68,  67,   1,   2,   5,   4
2,  65,  66,  69,  68,   2,   3,   6,   5
3,  67,  68,  71,  70,   4,   5,   8,   7
4,  68,  69,  72,  71,   5,   6,   9,   8
5,  70,  71,  74,  73,   7,   8,  11,  10
6,  71,  72,  75,  74,   8,   9,  12,  11
7,  73,  74,  77,  76,  10,  11,  14,  13
8,  74,  75,  78,  77,  11,  12,  15,  14
9,  76,  77,  80,  79,  13,  14,  17,  16
10,  77,  78,  81,  80,  14,  15,  18,  17
11,  79,  80,  83,  82,  16,  17,  20,  19
12,  80,  81,  84,  83,  17,  18,  21,  20
13,  82,  83,  86,  85,  19,  20,  23,  22
14,  83,  84,  87,  86,  20,  21,  24,  23
15,  85,  86,  89,  88,  22,  23,  26,  25
16,  86,  87,  90,  89,  23,  24,  27,  26
17,  88,  89,  92,  91,  25,  26,  29,  28
18,  89,  90,  93,  92,  26,  27,  30,  29
19,  91,  92,  95,  94,  28,  29,  32,  31
20,  92,  93,  96,  95,  29,  30,  33,  32
21,  94,  95,  98,  97,  31,  32,  35,  34
22,  95,  96,  99,  98,  32,  33,  36,  35
23,  97,  98, 101, 100,  34,  35,  38,  37
24,  98,  99, 102, 101,  35,  36,  39,  38
25, 100, 101, 104, 103,  37,  38,  41,  40
26, 101, 102, 105, 104,  38,  39,  42,  41
27, 103, 104, 107, 106,  40,  41,  44,  43
28, 104, 105, 108, 107,  41,  42,  45,  44
29, 106, 107, 110, 109,  43,  44,  47,  46
30, 107, 108, 111, 110,  44,  45,  48,  47
31, 109, 110, 113, 112,  46,  47,  50,  49
32, 110, 111, 114, 113,  47,  48,  51,  50
33, 112, 113, 116, 115,  49,  50,  53,  52
34, 113, 114, 117, 116,  50,  51,  54,  53
35, 115, 116, 119, 118,  52,  53,  56,  55
36, 116, 117, 120, 119,  53,  54,  57,  56
37, 118, 119, 122, 121,  55,  56,  59,  58
38, 119, 120, 123, 122,  56,  57,  60,  59
39, 121, 122, 125, 124,  58,  59,  62,  61
40, 122, 123, 126, 125,  59,  60,  63,  62
41, 127, 128, 131, 130,  64,  65,  68,  67
42, 128, 129, 132, 131,  65,  66,  69,  68
43, 130, 131, 134, 133,  67,  68,  71,  70
44, 131, 132, 135, 134,  68,  69,  72,  71
45, 133, 134, 137, 136,  70,  71,  74,  73
46, 134, 135, 138, 137,  71,  72,  75,  74
47, 136, 137, 140, 139,  73,  74,  77,  76
48, 137, 138, 141, 140,  74,  75,  78,  77
49, 139, 140, 143, 142,  76,  77,  80,  79
50, 140, 141, 144, 143,  77,  78,  81,  80
51, 142, 143, 146, 145,  79,  80,  83,  82
52, 143, 144, 147, 146,  80,  81,  84,  83
53, 145, 146, 149, 148,  82,  83,  86,  85
54, 146, 147, 150, 149,  83,  84,  87,  86
55, 148, 149, 152, 151,  85,  86,  89,  88
56, 149, 150, 153, 152,  86,  87,  90,  89
57, 151, 152, 155, 154,  88,  89,  92,  91
58, 152, 153, 156, 155,  89,  90,  93,  92
59, 154, 155, 158, 157,  91,  92,  95,  94
60, 155, 156, 159, 158,  92,  93,  96,  95
61, 157, 158, 161, 160,  94,  95,  98,  97
62, 158, 159, 162, 161,  95,  96,  99,  98
63, 160, 161, 164, 163,  97,  98, 101, 100
64, 161, 162, 165, 164,  98,  99, 102, 101
65, 163, 164, 167, 166, 100, 101, 104, 103
66, 164, 165, 168, 167, 101, 102, 105, 104
67, 166, 167, 170, 169, 103, 104, 107, 106
68, 167, 168, 171, 170, 104, 105, 108, 107
69, 169, 170, 173, 172, 106, 107, 110, 109
70, 170, 171, 174, 173, 107, 108, 111, 110
71, 172, 173, 176, 175, 109, 110, 113, 112
72, 173, 174, 177, 176, 110, 111, 114, 113
73, 175, 176, 179, 178, 112, 113, 116, 115
74, 176, 177, 180, 179, 113, 114, 117, 116
75, 178, 179, 182, 181, 115, 116, 119, 118
76, 179, 180, 183, 182, 116, 117, 120, 119
77, 181, 182, 185, 184, 118, 119, 122, 121
78, 182, 183, 186, 185, 119, 120, 123, 122
79, 184, 185, 188, 187, 121, 122, 125, 124
80, 185, 186, 189, 188, 122, 123, 126, 125
*Nset, nset=<node-set>, generate
1,  189,    1
*Elset, elset=<element-set>, generate
1,  80,   1
*Solid Section, elset=<element-set>, material=<material>
,
*End Part
*Assembly, name=<name>
*Instance, name=<instance>, part=<part>
*End Instance
*Nset, nset=<node-set>, instance=<instance>
61,  62,  63, 124, 125, 126, 187, 188, 189
*Elset, elset=<element-set>, instance=<instance>
39, 40, 79, 80
*Nset, nset=<node-set>, instance=<instance>
127,
*End Assembly
*Psd-definition, name=<name>, g=9.81, type=BASE
0.01, 0., 10.
0.01, 0., 2000.
*Spectrum, name=<name>, type=ACCELERATION
20., 10., 0.02
200., 100., 0.02
200., 600., 0.02
60., 2000., 0.02
*Material, name=<material>
*Density
8000.,
*Elastic
2e+11,0.
*Boundary
<node-set>, ENCASTRE
*Step, name=<step>, nlgeom=NO, perturbation
*Frequency, eigensolver=Lanczos, sim=NO, acoustic coupling=on, normalization=mass
20, , , , ,
*Restart, write, frequency=0
*Output, field
*Node Output, variable=PRESELECT
*Element Output, directions=YES
S,
*End Step
*Step, name=<step>, nlgeom=NO, perturbation
*Random Response
10., 2000., 20, 3.,
*Modal Damping, definition=FREQUENCY RANGE
10., 0.02
2000., 0.02
*Base Motion, type=ACCELERATION, dof=2, scale=1., load case=1
*Correlation, type=CORRELATED, psd=base-psd, complex=NO
1, 1.
*Output, field
*Node Output
RU, U
*Element Output, directions=YES
RMISES, RS
*Output, history, frequency=0
*End Step
*Step, name=<step>, nlgeom=NO, perturbation
*Response Spectrum, sum=SRSS
shock-spectrum,0., 1., 0., 1.
*Modal Damping, definition=FREQUENCY RANGE
10., 0.02
2000., 0.02
*Output, field
*Node Output
U,
*Element Output, directions=YES
S,
*Output, history, frequency=0
*End Step
10.3
Download artifacts
Download all verification artifacts (ZIP)
10.3.1
Geometry CAD files
No separate CAD file is required for this problem; the geometry is fully defined in the Abaqus/CAE journal files.
10.3.2
Abaqus/CAE journal files
2 generated Abaqus/CAE journals are available.
Show 2 individual journals
abaqus_h0p1–cantilever_shock_vibration_model.jnl
immersed_h0p1–cantilever_shock_vibration_model.jnl
10.4
Provenance
Note
Test_abaqus_verification_cantilever_shock_vibration_abaqus_h0p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:40:38
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus Unofficial Packaging Version
Note
Test_abaqus_verification_cantilever_shock_vibration_immersed_h0p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 4/4 enforced criteria passed
Test date
04-Nov-2026 19:40:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus Unofficial Packaging Version
