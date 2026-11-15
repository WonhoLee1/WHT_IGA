# Multi-planar material interfaces – Coreform IGA for Abaqus

Source: https://docs.coreform.com/cifa/latest/verification-manual/problems/solid_mechanics/linear_elastic_stress/multiplanar-material-interfaces/multiplanar-material-interfaces.html

---

Verification Manual
Solid mechanics
Linear elastic stress
Multi-planar material interfaces
Multi-planar material interfaces
\[
\newcommand{\ltwonorm}{L^2}
\newcommand{\degree}[1]{p = #1}
\newcommand{\Degree}[1]{\mathcal{P}^{#1}}
\]
Current evidence (exact required source revision e48e5098b2f869612b491032a1ca2b17aa1b6354). Execution status is not itself numerical verification.
1
Objective
Verify that internal material-cell boundaries do not materially change the response of a homogeneous trimmed solid, including where multiple planar interface faces meet along longitudinal edges or at cross-section corners. Six partitioned square cross-sections are compared with an unpartitioned monolithic beam using identical elastic properties in every material definition. The quantities of interest are the relative differences in tip displacement and in displacement and stress profiles sampled along the beam and across its section.
2
Geometry
The specimen is a prismatic cantilever with square cross-section
\(b \times b\)
and length
\(L\)
. Every internal cross-section boundary is extruded through the full beam length, producing planar material-interface faces parallel to the beam axis.
Figure 1
shows the monolithic reference and the six partition patterns.
Figure 1: Normalized cross-sections used by the homogeneous material-interface study. Each color denotes a distinct retained cell and material definition; all material properties are identical.
Table 1: Geometric parameters.
Quantity
Symbol
Value
Cross-section width
\(b\)
10
Cross-section height
\(b\)
10
Beam length
\(L\)
100
The first three partitions exercise rectangular inclusions that respectively touch one exterior face, touch two exterior faces, or remain completely enclosed. The zig-zag and nested-L partitions divide the square into equal-area cells and introduce non-collinear interface planes. The Tetris partition uses a
\(3 \times 3\)
construction: two equal four-square L cells surround a one-square central cell, so three materials meet at the central square’s corners.
3
Material
The partitioned models use distinct Abaqus material and solid-section definitions for each colored cell: two definitions in the first five cases and three in the Tetris case. Every definition has the same isotropic, linear-elastic properties, so the internal boundaries are physically invisible in the continuum problem.
Table 2: Properties assigned independently to every material definition.
Property
Symbol
Value
Young’s modulus
\(E\)
\(100 \times 10^{9}\)
Poisson’s ratio
\(\nu\)
0.3
4
Loading & boundary conditions
The
\(z=0\)
end face is pinned in all three translational directions. A uniform surface traction of magnitude
\(t_y\)
acts in the negative
\(y\)
direction over every cell face at
\(z=L\)
. Applying the same traction to the full partitioned end face produces resultant shear force
\(P=t_y b^2\)
and exercises bending and transverse shear without introducing density as another material parameter.
Table 3: Applied loading.
Quantity
Symbol
Value
End shear traction
\(t_y\)
\(10.0 \times 10^{-3}\)
Resultant force
\(P\)
\(1\)
Direction
—
Negative
\(y\)
5
Mesh & discretization
Each model uses the same immersed rectilinear background grid, outer geometry, grid origin, and
\(h/p\)
coordinates. The background element sizes are
\((5,5,10)/2^k\)
in
\((x,y,z)\)
for
\(k=0,\ldots,3\)
, at degrees
\(p=1\)
,
\(2\)
, and
\(3\)
with maximum continuity
\(C^{p-1}\)
. Only the retained internal cell topology changes between the monolithic and partitioned models. The complete study contains 84 independently schedulable simulations: seven cross-sections, four refinement levels, and three degrees. The following views look along the beam axis and use the
\(p=2\)
results databases. Retained cells are colored categorically by
part_id
using the Coreform primary palette in index order; each 2-by-2 figure progresses from the coarsest to the finest cross-sectional background mesh.
5.1
Monolithic reference
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 2: Monolithic reference cross-section meshes.
5.2
Bottom-centered one-third square
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 3: Bottom-centered one-third-square partition meshes.
5.3
Bottom-left one-half square
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 4: Bottom-left one-half-square partition meshes.
5.4
Centered one-third square
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 5: Centered one-third-square partition meshes.
5.5
Zig-zag interface
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 6: Zig-zag partition meshes.
5.6
Nested L interface
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 7: Nested-L partition meshes.
5.7
Tetris partition
(a)
\(h_{xy}=5\)
(b)
\(h_{xy}=2.5\)
(c)
\(h_{xy}=1.25\)
(d)
\(h_{xy}=0.625\)
Figure 8: Three-cell Tetris partition meshes.
6
Reference & accepted solutions
The primary numerical reference for every partitioned result is the monolithic beam at the identical degree, refinement level, background-grid placement, loading, and material properties. This same-
\(h\)
, same-
\(p\)
pairing removes ordinary discretization error from the comparison and isolates the effect of creating and coupling retained material cells.
Six response differences are recorded for every pairing: tip
\(u_y\)
; relative
\(L^2\)
differences in
\(u_y\)
along an axial line; and relative
\(L^2\)
differences in
\(\sigma_{zz}\)
and
\(\tau_{yz}\)
along vertical and horizontal midspan probes. The study-level CTest takes the maximum of those six quantities at the finest level and enforces the calibrated acceptance tolerance. The homogeneous Euler–Bernoulli cantilever solution
[1]
is also recorded by each simulation as a secondary, record-only diagnostic; it is not used to judge interface invariance because it carries independent beam-theory modeling error.
7
Results
Each partition geometry below is presented independently. The case-tabbed table reports the measured partitioned response, the same-
\(h\)
, same-
\(p\)
monolithic reference response, and their relative difference for every mesh size and degree. For the tip displacement these are scalar values. For each line probe the measured and reference values are the profile
\(L^2\)
norms, while the relative difference is the
\(L^2\)
norm of the pointwise profile difference normalized by the reference norm. The three probe figures use one subfigure per degree and overlay every partitioned mesh size against the finest matching-degree monolithic result; the final figure tracks the maximum same-mesh relative difference across all six verification quantities under refinement.
7.1
Bottom-centered one-third square
h=10, p=1
h=5, p=1
h=2.5, p=1
h=1.25, p=1
h=10, p=2
h=5, p=2
h=2.5, p=2
h=1.25, p=2
h=10, p=3
h=5, p=3
h=2.5, p=3
h=1.25, p=3
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-1.89 \times 10^{-9}\)
\(-2.09 \times 10^{-9}\)
\(9.44\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(9.28 \times 10^{-9}\)
\(10.2 \times 10^{-9}\)
\(9.26\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.25\)
\(0.279\)
\(10.4\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.3 \times 10^{-3}\)
\(11.9\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(48.4 \times 10^{-3}\)
\(53.6 \times 10^{-3}\)
\(9.76\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(38.9 \times 10^{-3}\)
\(36.1 \times 10^{-3}\)
\(8.09\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.32 \times 10^{-9}\)
\(-3.41 \times 10^{-9}\)
\(2.67\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(16.2 \times 10^{-9}\)
\(16.6 \times 10^{-9}\)
\(2.63\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.457\)
\(0.471\)
\(2.96\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.9 \times 10^{-3}\)
\(30.6 \times 10^{-3}\)
\(6.39\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(91.0 \times 10^{-3}\)
\(94.4 \times 10^{-3}\)
\(3.64\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(43.9 \times 10^{-3}\)
\(43.9 \times 10^{-3}\)
\(0.967\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.87 \times 10^{-9}\)
\(-3.88 \times 10^{-9}\)
\(0.238\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(18.8 \times 10^{-9}\)
\(18.8 \times 10^{-9}\)
\(0.236\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.531\)
\(0.533\)
\(0.46\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(1.08\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.108\)
\(0.108\)
\(0.292\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.8 \times 10^{-3}\)
\(0.375\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.97 \times 10^{-9}\)
\(-3.97 \times 10^{-9}\)
\(0.036\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.3 \times 10^{-9}\)
\(19.3 \times 10^{-9}\)
\(0.0361\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.544\)
\(0.544\)
\(0.219\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.936\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(90.1 \times 10^{-3}\)
\(90.1 \times 10^{-3}\)
\(0.0398\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.7 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.0783\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.99 \times 10^{-9}\)
\(-3.99 \times 10^{-9}\)
\(0.0107\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.0128\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0948\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(26.4 \times 10^{-3}\)
\(26.7 \times 10^{-3}\)
\(2.37\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.2 \times 10^{-3}\)
\(93.8 \times 10^{-3}\)
\(0.561\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(36.2 \times 10^{-3}\)
\(36.5 \times 10^{-3}\)
\(1.11\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(88.9 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(200 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.017\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.8 \times 10^{-3}\)
\(30.8 \times 10^{-3}\)
\(1.32\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0823\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(45.3 \times 10^{-3}\)
\(45.3 \times 10^{-3}\)
\(0.0602\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(311 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(398 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00184\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.192\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0235\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.0625\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(73.1 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(85.1 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00646\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.0412\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(540 \times 10^{-6}\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.00315\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00695\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00819\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.073\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.0 \times 10^{-3}\)
\(0.808\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.7 \times 10^{-3}\)
\(94.7 \times 10^{-3}\)
\(0.244\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.547\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00255\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00308\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0238\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.527\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0236\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.0 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.233\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(510 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(608 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00305\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0572\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0127\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.00754\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(87.6 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(103 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00209\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.00742\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.00157\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.00207\%\)
Download data (CSV)
7.1.1
Centerline displacement profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 9: Bottom-centered one-third square: centerline displacement profiles under mesh refinement.
7.1.2
Transverse normal-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 10: Bottom-centered one-third square: transverse normal-stress profiles under mesh refinement.
7.1.3
Transverse shear-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 11: Bottom-centered one-third square: transverse shear-stress profiles under mesh refinement.
7.1.4
Relative differences under refinement
Figure 12: Bottom-centered one-third square: maximum homogeneous-response difference under refinement.
7.2
Bottom-left one-half square
h=10, p=1
h=5, p=1
h=2.5, p=1
h=1.25, p=1
h=10, p=2
h=5, p=2
h=2.5, p=2
h=1.25, p=2
h=10, p=3
h=5, p=3
h=2.5, p=3
h=1.25, p=3
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-1.95 \times 10^{-9}\)
\(-2.09 \times 10^{-9}\)
\(6.52\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(9.59 \times 10^{-9}\)
\(10.2 \times 10^{-9}\)
\(6.26\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.256\)
\(0.279\)
\(8.35\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.0 \times 10^{-3}\)
\(31.3 \times 10^{-3}\)
\(7.75\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(45.4 \times 10^{-3}\)
\(53.6 \times 10^{-3}\)
\(18.6\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(37.9 \times 10^{-3}\)
\(36.1 \times 10^{-3}\)
\(5.92\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.34 \times 10^{-9}\)
\(-3.41 \times 10^{-9}\)
\(2.13\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(16.2 \times 10^{-9}\)
\(16.6 \times 10^{-9}\)
\(2.1\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.459\)
\(0.471\)
\(2.74\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.1 \times 10^{-3}\)
\(30.6 \times 10^{-3}\)
\(8.78\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(87.2 \times 10^{-3}\)
\(94.4 \times 10^{-3}\)
\(10.2\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(43.7 \times 10^{-3}\)
\(43.9 \times 10^{-3}\)
\(2.8\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.87 \times 10^{-9}\)
\(-3.88 \times 10^{-9}\)
\(0.186\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(18.8 \times 10^{-9}\)
\(18.8 \times 10^{-9}\)
\(0.184\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.531\)
\(0.533\)
\(0.429\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.828\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.104\)
\(0.108\)
\(6.01\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.7 \times 10^{-3}\)
\(46.8 \times 10^{-3}\)
\(0.698\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.97 \times 10^{-9}\)
\(-3.97 \times 10^{-9}\)
\(0.0198\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.3 \times 10^{-9}\)
\(19.3 \times 10^{-9}\)
\(0.0197\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.544\)
\(0.544\)
\(0.119\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.48\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(88.6 \times 10^{-3}\)
\(90.1 \times 10^{-3}\)
\(2.29\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.43\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.99 \times 10^{-9}\)
\(-3.99 \times 10^{-9}\)
\(0.00292\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00333\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.184\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(27.0 \times 10^{-3}\)
\(26.7 \times 10^{-3}\)
\(4.54\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(95.0 \times 10^{-3}\)
\(93.8 \times 10^{-3}\)
\(1.65\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(37.6 \times 10^{-3}\)
\(36.5 \times 10^{-3}\)
\(3.66\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00203\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00245\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.022\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.8 \times 10^{-3}\)
\(30.8 \times 10^{-3}\)
\(0.375\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(95.0 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.241\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(45.3 \times 10^{-3}\)
\(45.3 \times 10^{-3}\)
\(0.322\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(491 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(574 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0205\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.398\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0478\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.357\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(326 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(393 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00425\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.121\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0258\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.199\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00656\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00807\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.149\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.0 \times 10^{-3}\)
\(0.661\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(95.3 \times 10^{-3}\)
\(94.7 \times 10^{-3}\)
\(0.752\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.2 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(1.57\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00192\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00232\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0527\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.591\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.102\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.0 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.222\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(219 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(258 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00366\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0783\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0313\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.106\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(158 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(188 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00174\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0158\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0205\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0104\%\)
Download data (CSV)
7.2.1
Centerline displacement profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 13: Bottom-left one-half square: centerline displacement profiles under mesh refinement.
7.2.2
Transverse normal-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 14: Bottom-left one-half square: transverse normal-stress profiles under mesh refinement.
7.2.3
Transverse shear-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 15: Bottom-left one-half square: transverse shear-stress profiles under mesh refinement.
7.2.4
Relative differences under refinement
Figure 16: Bottom-left one-half square: maximum homogeneous-response difference under refinement.
7.3
Centered one-third square
h=10, p=1
h=5, p=1
h=2.5, p=1
h=1.25, p=1
h=10, p=2
h=5, p=2
h=2.5, p=2
h=1.25, p=2
h=10, p=3
h=5, p=3
h=2.5, p=3
h=1.25, p=3
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-1.95 \times 10^{-9}\)
\(-2.09 \times 10^{-9}\)
\(6.8\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(9.57 \times 10^{-9}\)
\(10.2 \times 10^{-9}\)
\(6.46\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.252\)
\(0.279\)
\(9.72\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.3 \times 10^{-3}\)
\(7.82\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(48.3 \times 10^{-3}\)
\(53.6 \times 10^{-3}\)
\(10\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(38.1 \times 10^{-3}\)
\(36.1 \times 10^{-3}\)
\(5.9\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.17 \times 10^{-9}\)
\(-3.41 \times 10^{-9}\)
\(7.06\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(15.4 \times 10^{-9}\)
\(16.6 \times 10^{-9}\)
\(6.99\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.434\)
\(0.471\)
\(7.85\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.0 \times 10^{-3}\)
\(30.6 \times 10^{-3}\)
\(3.11\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(87.2 \times 10^{-3}\)
\(94.4 \times 10^{-3}\)
\(7.72\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(44.1 \times 10^{-3}\)
\(43.9 \times 10^{-3}\)
\(1.8\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.87 \times 10^{-9}\)
\(-3.88 \times 10^{-9}\)
\(0.293\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(18.8 \times 10^{-9}\)
\(18.8 \times 10^{-9}\)
\(0.291\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.531\)
\(0.533\)
\(0.424\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(1.59\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.108\)
\(0.108\)
\(0.561\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.8 \times 10^{-3}\)
\(0.401\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.97 \times 10^{-9}\)
\(-3.97 \times 10^{-9}\)
\(0.0566\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.3 \times 10^{-9}\)
\(19.3 \times 10^{-9}\)
\(0.0567\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.544\)
\(0.544\)
\(0.306\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(1.13\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(90.1 \times 10^{-3}\)
\(90.1 \times 10^{-3}\)
\(0.161\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.7 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.101\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.99 \times 10^{-9}\)
\(-3.99 \times 10^{-9}\)
\(0.00503\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00667\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.142\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(26.7 \times 10^{-3}\)
\(26.7 \times 10^{-3}\)
\(3.36\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(93.1 \times 10^{-3}\)
\(93.8 \times 10^{-3}\)
\(0.801\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(36.6 \times 10^{-3}\)
\(36.5 \times 10^{-3}\)
\(0.526\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00156\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.0017\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0216\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.8 \times 10^{-3}\)
\(30.8 \times 10^{-3}\)
\(0.647\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.7 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.488\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(45.3 \times 10^{-3}\)
\(45.3 \times 10^{-3}\)
\(0.503\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(175 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(199 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00963\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.157\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0264\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.275\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(92.1 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(114 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00756\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.108\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.022\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.0529\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00316\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00339\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.204\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.0 \times 10^{-3}\)
\(32.0 \times 10^{-3}\)
\(3.01\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.1 \times 10^{-3}\)
\(94.7 \times 10^{-3}\)
\(0.918\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.0 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(1.92\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(677 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(822 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0231\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0878\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(95.2 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.369\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.293\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(117 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(134 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0063\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.167\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0333\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0263\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(86.9 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(101 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00253\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0382\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0204\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0259\%\)
Download data (CSV)
7.3.1
Centerline displacement profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 17: Centered one-third square: centerline displacement profiles under mesh refinement.
7.3.2
Transverse normal-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 18: Centered one-third square: transverse normal-stress profiles under mesh refinement.
7.3.3
Transverse shear-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 19: Centered one-third square: transverse shear-stress profiles under mesh refinement.
7.3.4
Relative differences under refinement
Figure 20: Centered one-third square: maximum homogeneous-response difference under refinement.
7.4
Zig-zag interface
h=10, p=1
h=5, p=1
h=2.5, p=1
h=1.25, p=1
h=10, p=2
h=5, p=2
h=2.5, p=2
h=1.25, p=2
h=10, p=3
h=5, p=3
h=2.5, p=3
h=1.25, p=3
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-1.99 \times 10^{-9}\)
\(-2.09 \times 10^{-9}\)
\(4.56\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(9.79 \times 10^{-9}\)
\(10.2 \times 10^{-9}\)
\(4.33\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.264\)
\(0.279\)
\(5.48\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.6 \times 10^{-3}\)
\(31.3 \times 10^{-3}\)
\(4.29\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(50.7 \times 10^{-3}\)
\(53.6 \times 10^{-3}\)
\(6.81\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(37.2 \times 10^{-3}\)
\(36.1 \times 10^{-3}\)
\(3.22\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.34 \times 10^{-9}\)
\(-3.41 \times 10^{-9}\)
\(2.18\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(16.2 \times 10^{-9}\)
\(16.6 \times 10^{-9}\)
\(2.14\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.459\)
\(0.471\)
\(2.75\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.9 \times 10^{-3}\)
\(30.6 \times 10^{-3}\)
\(6.14\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(92.8 \times 10^{-3}\)
\(94.4 \times 10^{-3}\)
\(3.48\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(44.6 \times 10^{-3}\)
\(43.9 \times 10^{-3}\)
\(4.02\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.87 \times 10^{-9}\)
\(-3.88 \times 10^{-9}\)
\(0.192\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(18.8 \times 10^{-9}\)
\(18.8 \times 10^{-9}\)
\(0.19\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.532\)
\(0.533\)
\(0.295\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(1.13\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.107\)
\(0.108\)
\(2.02\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.8 \times 10^{-3}\)
\(0.627\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.97 \times 10^{-9}\)
\(-3.97 \times 10^{-9}\)
\(0.0271\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.3 \times 10^{-9}\)
\(19.3 \times 10^{-9}\)
\(0.0269\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.544\)
\(0.544\)
\(0.133\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(1.16\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(90.1 \times 10^{-3}\)
\(90.1 \times 10^{-3}\)
\(0.838\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.499\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.99 \times 10^{-9}\)
\(-3.99 \times 10^{-9}\)
\(0.0533\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.0643\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0948\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(26.9 \times 10^{-3}\)
\(26.7 \times 10^{-3}\)
\(6.88\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.2 \times 10^{-3}\)
\(93.8 \times 10^{-3}\)
\(1.06\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(39.0 \times 10^{-3}\)
\(36.5 \times 10^{-3}\)
\(7.19\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00719\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00852\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0179\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.8 \times 10^{-3}\)
\(30.8 \times 10^{-3}\)
\(0.611\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.8 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.463\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(45.2 \times 10^{-3}\)
\(45.3 \times 10^{-3}\)
\(0.793\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00153\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00182\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00351\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.153\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0678\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.255\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(45.2 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(55.8 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00252\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.141\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0508\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.166\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.0157\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.0181\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0638\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.0 \times 10^{-3}\)
\(0.788\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.8 \times 10^{-3}\)
\(94.7 \times 10^{-3}\)
\(0.56\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(1.32\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(427 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(494 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0309\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.438\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.7 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.208\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.282\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(148 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(172 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00354\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.174\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.8 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0881\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.0 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.263\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(48.9 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(52.8 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00146\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.464\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0186\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.263\%\)
Download data (CSV)
7.4.1
Centerline displacement profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 21: Zig-zag interface: centerline displacement profiles under mesh refinement.
7.4.2
Transverse normal-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 22: Zig-zag interface: transverse normal-stress profiles under mesh refinement.
7.4.3
Transverse shear-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 23: Zig-zag interface: transverse shear-stress profiles under mesh refinement.
7.4.4
Relative differences under refinement
Figure 24: Zig-zag interface: maximum homogeneous-response difference under refinement.
7.5
Nested L interface
h=10, p=1
h=5, p=1
h=2.5, p=1
h=1.25, p=1
h=10, p=2
h=5, p=2
h=2.5, p=2
h=1.25, p=2
h=10, p=3
h=5, p=3
h=2.5, p=3
h=1.25, p=3
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-1.97 \times 10^{-9}\)
\(-2.09 \times 10^{-9}\)
\(5.7\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(9.67 \times 10^{-9}\)
\(10.2 \times 10^{-9}\)
\(5.46\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.259\)
\(0.279\)
\(7.17\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.7 \times 10^{-3}\)
\(31.3 \times 10^{-3}\)
\(6.61\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(49.8 \times 10^{-3}\)
\(53.6 \times 10^{-3}\)
\(8.16\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(37.9 \times 10^{-3}\)
\(36.1 \times 10^{-3}\)
\(5.33\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.27 \times 10^{-9}\)
\(-3.41 \times 10^{-9}\)
\(4.25\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(15.9 \times 10^{-9}\)
\(16.6 \times 10^{-9}\)
\(4.2\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.451\)
\(0.471\)
\(4.5\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.1 \times 10^{-3}\)
\(30.6 \times 10^{-3}\)
\(2.72\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(91.7 \times 10^{-3}\)
\(94.4 \times 10^{-3}\)
\(5.81\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(44.2 \times 10^{-3}\)
\(43.9 \times 10^{-3}\)
\(3.01\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.87 \times 10^{-9}\)
\(-3.88 \times 10^{-9}\)
\(0.236\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(18.8 \times 10^{-9}\)
\(18.8 \times 10^{-9}\)
\(0.235\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.531\)
\(0.533\)
\(0.268\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(2.08\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.108\)
\(0.108\)
\(1.21\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.8 \times 10^{-3}\)
\(0.705\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.97 \times 10^{-9}\)
\(-3.97 \times 10^{-9}\)
\(0.042\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.3 \times 10^{-9}\)
\(19.3 \times 10^{-9}\)
\(0.0421\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.544\)
\(0.544\)
\(0.153\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.394\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(90.1 \times 10^{-3}\)
\(90.1 \times 10^{-3}\)
\(0.22\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.7 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.145\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.99 \times 10^{-9}\)
\(-3.99 \times 10^{-9}\)
\(0.0274\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.0334\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.106\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(26.7 \times 10^{-3}\)
\(26.7 \times 10^{-3}\)
\(3.33\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(95.3 \times 10^{-3}\)
\(93.8 \times 10^{-3}\)
\(1.76\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(37.8 \times 10^{-3}\)
\(36.5 \times 10^{-3}\)
\(4.26\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(728 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(873 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.07\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.8 \times 10^{-3}\)
\(30.8 \times 10^{-3}\)
\(1.16\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.7 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.509\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(45.3 \times 10^{-3}\)
\(45.3 \times 10^{-3}\)
\(0.513\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(261 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(272 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0109\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.328\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0756\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.368\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(14.7 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(18.8 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00653\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.13\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0391\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.128\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00807\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00971\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0789\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.0 \times 10^{-3}\)
\(32.0 \times 10^{-3}\)
\(0.931\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(96.4 \times 10^{-3}\)
\(94.7 \times 10^{-3}\)
\(2.32\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.1 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(2.69\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00115\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00131\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0313\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.31\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.7 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.246\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.0 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.279\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(179 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(209 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00402\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0637\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.071\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0807\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(41.2 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(45.4 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00142\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0313\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0136\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0193\%\)
Download data (CSV)
7.5.1
Centerline displacement profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 25: Nested L interface: centerline displacement profiles under mesh refinement.
7.5.2
Transverse normal-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 26: Nested L interface: transverse normal-stress profiles under mesh refinement.
7.5.3
Transverse shear-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 27: Nested L interface: transverse shear-stress profiles under mesh refinement.
7.5.4
Relative differences under refinement
Figure 28: Nested L interface: maximum homogeneous-response difference under refinement.
7.6
Tetris partition
h=10, p=1
h=5, p=1
h=2.5, p=1
h=1.25, p=1
h=10, p=2
h=5, p=2
h=2.5, p=2
h=1.25, p=2
h=10, p=3
h=5, p=3
h=2.5, p=3
h=1.25, p=3
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-1.93 \times 10^{-9}\)
\(-2.09 \times 10^{-9}\)
\(7.42\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(9.50 \times 10^{-9}\)
\(10.2 \times 10^{-9}\)
\(7.21\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.255\)
\(0.279\)
\(8.83\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.7 \times 10^{-3}\)
\(31.3 \times 10^{-3}\)
\(12.1\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(48.9 \times 10^{-3}\)
\(53.6 \times 10^{-3}\)
\(8.97\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(38.9 \times 10^{-3}\)
\(36.1 \times 10^{-3}\)
\(8.2\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.17 \times 10^{-9}\)
\(-3.41 \times 10^{-9}\)
\(7.15\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(15.4 \times 10^{-9}\)
\(16.6 \times 10^{-9}\)
\(7.08\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.435\)
\(0.471\)
\(7.83\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.2 \times 10^{-3}\)
\(30.6 \times 10^{-3}\)
\(4.9\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(87.0 \times 10^{-3}\)
\(94.4 \times 10^{-3}\)
\(8.27\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(43.8 \times 10^{-3}\)
\(43.9 \times 10^{-3}\)
\(1.5\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.86 \times 10^{-9}\)
\(-3.88 \times 10^{-9}\)
\(0.394\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(18.8 \times 10^{-9}\)
\(18.8 \times 10^{-9}\)
\(0.392\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.531\)
\(0.533\)
\(0.616\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(1.35\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.107\)
\(0.108\)
\(1.1\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.8 \times 10^{-3}\)
\(0.0943\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.97 \times 10^{-9}\)
\(-3.97 \times 10^{-9}\)
\(0.0768\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.3 \times 10^{-9}\)
\(19.3 \times 10^{-9}\)
\(0.077\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.544\)
\(0.544\)
\(0.311\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(1.13\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(90.1 \times 10^{-3}\)
\(90.1 \times 10^{-3}\)
\(0.223\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.7 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.107\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-3.99 \times 10^{-9}\)
\(-3.99 \times 10^{-9}\)
\(0.00813\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00971\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0858\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(26.7 \times 10^{-3}\)
\(26.7 \times 10^{-3}\)
\(7.39\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.4 \times 10^{-3}\)
\(93.8 \times 10^{-3}\)
\(1.26\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(37.8 \times 10^{-3}\)
\(36.5 \times 10^{-3}\)
\(5.03\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00345\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00435\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0439\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(30.8 \times 10^{-3}\)
\(30.8 \times 10^{-3}\)
\(1.7\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(95.1 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.643\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(45.7 \times 10^{-3}\)
\(45.3 \times 10^{-3}\)
\(1.08\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00226\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00269\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0093\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(31.8 \times 10^{-3}\)
\(31.8 \times 10^{-3}\)
\(0.503\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0642\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.8 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.48\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(80.6 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(91.5 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00646\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.1 \times 10^{-3}\)
\(32.1 \times 10^{-3}\)
\(0.0936\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.8 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.0575\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.9 \times 10^{-3}\)
\(46.9 \times 10^{-3}\)
\(0.0918\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00278\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00316\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.079\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.0 \times 10^{-3}\)
\(32.0 \times 10^{-3}\)
\(0.814\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.8 \times 10^{-3}\)
\(94.7 \times 10^{-3}\)
\(0.354\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(46.6 \times 10^{-3}\)
\(46.7 \times 10^{-3}\)
\(0.71\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(0.00121\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(0.00145\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.0167\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.22\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.7 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.213\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.112\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(65.9 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(84.2 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00902\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.112\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.171\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0339\%\)
Quantity
Measured value
Reference value
Relative difference
Tip displacement
\(u_y\)
\(-4.00 \times 10^{-9}\)
\(-4.00 \times 10^{-9}\)
\(267 \times 10^{-6}\%\)
Centerline displacement
\(\lVert u_y\rVert_{L^2}\)
\(19.4 \times 10^{-9}\)
\(19.4 \times 10^{-9}\)
\(319 \times 10^{-6}\%\)
Vertical normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(0.548\)
\(0.548\)
\(0.00243\%\)
Vertical shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(32.2 \times 10^{-3}\)
\(32.2 \times 10^{-3}\)
\(0.0387\%\)
Horizontal normal stress
\(\lVert\sigma_{zz}\rVert_{L^2}\)
\(94.9 \times 10^{-3}\)
\(94.9 \times 10^{-3}\)
\(0.00883\%\)
Horizontal shear stress
\(\lVert\tau_{yz}\rVert_{L^2}\)
\(47.1 \times 10^{-3}\)
\(47.1 \times 10^{-3}\)
\(0.0189\%\)
Download data (CSV)
7.6.1
Centerline displacement profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 29: Tetris partition: centerline displacement profiles under mesh refinement.
7.6.2
Transverse normal-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 30: Tetris partition: transverse normal-stress profiles under mesh refinement.
7.6.3
Transverse shear-stress profile
(a)
\(p=1\)
(b)
\(p=2\)
(c)
\(p=3\)
Figure 31: Tetris partition: transverse shear-stress profiles under mesh refinement.
7.6.4
Relative differences under refinement
Figure 32: Tetris partition: maximum homogeneous-response difference under refinement.
8
Performance
Trimming and Abaqus solve statistics are reported for every mesh level and degree, grouped by geometry.
8.1
Monolithic reference
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(99\)
\(90\)
\(1\)
\(5.3\)
\(5\)
\(1{,}127\)
\(525\)
\(354\)
\(2\)
\(11.2\)
\(2.5\)
\(5{,}203\)
\(3{,}321\)
\(1{,}410\)
\(4\)
\(33.8\)
\(1.25\)
\(29{,}963\)
\(23{,}409\)
\(5{,}634\)
\(8\)
\(155.2\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(99\)
\(90\)
\(1\)
\(6.7\)
\(5\)
\(2{,}025\)
\(525\)
\(354\)
\(2\)
\(13.9\)
\(2.5\)
\(7{,}605\)
\(3{,}321\)
\(1{,}410\)
\(4\)
\(45.3\)
\(1.25\)
\(37{,}485\)
\(23{,}409\)
\(5{,}634\)
\(8\)
\(156.3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(99\)
\(90\)
\(1\)
\(26.5\)
\(5\)
\(3{,}267\)
\(525\)
\(354\)
\(2\)
\(32.0\)
\(2.5\)
\(10{,}575\)
\(3{,}321\)
\(1{,}410\)
\(4\)
\(65.3\)
\(1.25\)
\(46{,}023\)
\(23{,}409\)
\(5{,}634\)
\(8\)
\(180.4\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(579\)
\(482\)
\(1{,}446\)
\(1\)
\(1.0\)
\(5\)
\(1{,}137\)
\(1{,}152\)
\(3{,}456\)
\(2\)
\(0.0\)
\(2.5\)
\(4{,}345\)
\(4{,}778\)
\(14{,}334\)
\(4\)
\(2.0\)
\(1.25\)
\(26{,}605\)
\(28{,}296\)
\(84{,}888\)
\(8\)
\(4.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(571\)
\(607\)
\(1{,}821\)
\(1\)
\(1.0\)
\(5\)
\(1{,}137\)
\(1{,}487\)
\(4{,}461\)
\(2\)
\(2.0\)
\(2.5\)
\(4{,}345\)
\(5{,}781\)
\(17{,}343\)
\(4\)
\(2.0\)
\(1.25\)
\(27{,}055\)
\(31{,}916\)
\(95{,}748\)
\(8\)
\(11.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(571\)
\(786\)
\(2{,}358\)
\(1\)
\(2.0\)
\(5\)
\(1{,}137\)
\(1{,}896\)
\(5{,}688\)
\(2\)
\(5.0\)
\(2.5\)
\(4{,}993\)
\(7{,}274\)
\(21{,}822\)
\(4\)
\(12.0\)
\(1.25\)
\(29{,}817\)
\(37{,}002\)
\(111{,}006\)
\(8\)
\(55.0\)
Download data (CSV)
8.2
Bottom-centered one-third square
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(132\)
\(110\)
\(1\)
\(7.1\)
\(5\)
\(1{,}127\)
\(630\)
\(390\)
\(2\)
\(17.1\)
\(2.5\)
\(5{,}203\)
\(3{,}690\)
\(1{,}638\)
\(4\)
\(56.2\)
\(1.25\)
\(29{,}963\)
\(24{,}786\)
\(6{,}374\)
\(8\)
\(206.6\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(110\)
\(1\)
\(8.3\)
\(5\)
\(2{,}025\)
\(630\)
\(390\)
\(2\)
\(18.6\)
\(2.5\)
\(7{,}605\)
\(3{,}690\)
\(1{,}638\)
\(4\)
\(63.0\)
\(1.25\)
\(37{,}485\)
\(24{,}786\)
\(6{,}374\)
\(8\)
\(202.7\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(110\)
\(1\)
\(33.0\)
\(5\)
\(3{,}267\)
\(630\)
\(390\)
\(2\)
\(39.6\)
\(2.5\)
\(10{,}575\)
\(3{,}690\)
\(1{,}638\)
\(4\)
\(85.2\)
\(1.25\)
\(46{,}023\)
\(24{,}786\)
\(6{,}374\)
\(8\)
\(255.7\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}556\)
\(3{,}696\)
\(11{,}088\)
\(1\)
\(1.0\)
\(5\)
\(7{,}062\)
\(4{,}974\)
\(14{,}922\)
\(2\)
\(2.0\)
\(2.5\)
\(12{,}406\)
\(10{,}172\)
\(30{,}516\)
\(4\)
\(3.0\)
\(1.25\)
\(48{,}538\)
\(42{,}586\)
\(127{,}758\)
\(8\)
\(6.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}578\)
\(3{,}938\)
\(11{,}814\)
\(1\)
\(2.0\)
\(5\)
\(7{,}054\)
\(5{,}518\)
\(16{,}554\)
\(2\)
\(2.0\)
\(2.5\)
\(12{,}406\)
\(11{,}750\)
\(35{,}250\)
\(4\)
\(5.0\)
\(1.25\)
\(48{,}580\)
\(47{,}759\)
\(143{,}277\)
\(8\)
\(16.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}644\)
\(4{,}299\)
\(12{,}897\)
\(1\)
\(4.0\)
\(5\)
\(7{,}470\)
\(6{,}438\)
\(19{,}314\)
\(2\)
\(6.0\)
\(2.5\)
\(15{,}174\)
\(15{,}028\)
\(45{,}084\)
\(4\)
\(17.0\)
\(1.25\)
\(64{,}942\)
\(62{,}052\)
\(186{,}156\)
\(8\)
\(67.0\)
Download data (CSV)
8.3
Bottom-left one-half square
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(132\)
\(88\)
\(1\)
\(7.2\)
\(5\)
\(1{,}127\)
\(630\)
\(346\)
\(2\)
\(17.8\)
\(2.5\)
\(5{,}203\)
\(3{,}690\)
\(1{,}378\)
\(4\)
\(56.0\)
\(1.25\)
\(29{,}963\)
\(24{,}786\)
\(5{,}506\)
\(8\)
\(207.5\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(88\)
\(1\)
\(8.8\)
\(5\)
\(2{,}025\)
\(630\)
\(346\)
\(2\)
\(19.4\)
\(2.5\)
\(7{,}605\)
\(3{,}690\)
\(1{,}378\)
\(4\)
\(60.4\)
\(1.25\)
\(37{,}485\)
\(24{,}786\)
\(5{,}506\)
\(8\)
\(213.8\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(88\)
\(1\)
\(47.7\)
\(5\)
\(3{,}267\)
\(630\)
\(346\)
\(2\)
\(37.8\)
\(2.5\)
\(10{,}575\)
\(3{,}690\)
\(1{,}378\)
\(4\)
\(84.0\)
\(1.25\)
\(46{,}023\)
\(24{,}786\)
\(5{,}506\)
\(8\)
\(246.8\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}534\)
\(3{,}505\)
\(10{,}515\)
\(1\)
\(1.0\)
\(5\)
\(6{,}918\)
\(4{,}738\)
\(14{,}214\)
\(2\)
\(1.0\)
\(2.5\)
\(11{,}522\)
\(9{,}576\)
\(28{,}728\)
\(4\)
\(1.0\)
\(1.25\)
\(47{,}298\)
\(41{,}780\)
\(125{,}340\)
\(8\)
\(6.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}548\)
\(3{,}746\)
\(11{,}238\)
\(1\)
\(2.0\)
\(5\)
\(6{,}918\)
\(5{,}292\)
\(15{,}876\)
\(2\)
\(3.0\)
\(2.5\)
\(11{,}522\)
\(11{,}158\)
\(33{,}474\)
\(4\)
\(5.0\)
\(1.25\)
\(47{,}334\)
\(46{,}960\)
\(140{,}880\)
\(8\)
\(17.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(5{,}584\)
\(4{,}092\)
\(12{,}276\)
\(1\)
\(7.0\)
\(5\)
\(7{,}238\)
\(6{,}162\)
\(18{,}486\)
\(2\)
\(8.0\)
\(2.5\)
\(14{,}542\)
\(14{,}572\)
\(43{,}716\)
\(4\)
\(16.0\)
\(1.25\)
\(63{,}590\)
\(61{,}092\)
\(183{,}276\)
\(8\)
\(67.0\)
Download data (CSV)
8.4
Centered one-third square
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(110\)
\(99\)
\(1\)
\(6.9\)
\(5\)
\(1{,}127\)
\(693\)
\(504\)
\(2\)
\(20.7\)
\(2.5\)
\(5{,}203\)
\(3{,}649\)
\(1{,}720\)
\(4\)
\(55.8\)
\(1.25\)
\(29{,}963\)
\(25{,}353\)
\(7{,}480\)
\(8\)
\(231.0\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(110\)
\(99\)
\(1\)
\(8.4\)
\(5\)
\(2{,}025\)
\(693\)
\(504\)
\(2\)
\(25.5\)
\(2.5\)
\(7{,}605\)
\(3{,}649\)
\(1{,}720\)
\(4\)
\(62.3\)
\(1.25\)
\(37{,}485\)
\(25{,}353\)
\(7{,}480\)
\(8\)
\(232.6\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(110\)
\(99\)
\(1\)
\(40.2\)
\(5\)
\(3{,}267\)
\(693\)
\(504\)
\(2\)
\(73.5\)
\(2.5\)
\(10{,}575\)
\(3{,}649\)
\(1{,}720\)
\(4\)
\(89.5\)
\(1.25\)
\(46{,}023\)
\(25{,}353\)
\(7{,}480\)
\(8\)
\(487.0\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}366\)
\(5{,}010\)
\(15{,}030\)
\(1\)
\(4.0\)
\(5\)
\(9{,}993\)
\(6{,}672\)
\(20{,}016\)
\(2\)
\(1.0\)
\(2.5\)
\(14{,}513\)
\(11{,}442\)
\(34{,}326\)
\(4\)
\(2.0\)
\(1.25\)
\(57{,}941\)
\(48{,}232\)
\(144{,}696\)
\(8\)
\(8.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}196\)
\(5{,}135\)
\(15{,}405\)
\(1\)
\(5.0\)
\(5\)
\(9{,}993\)
\(7{,}218\)
\(21{,}654\)
\(2\)
\(2.0\)
\(2.5\)
\(14{,}513\)
\(12{,}848\)
\(38{,}544\)
\(4\)
\(5.0\)
\(1.25\)
\(57{,}815\)
\(53{,}568\)
\(160{,}704\)
\(8\)
\(19.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}196\)
\(5{,}437\)
\(16{,}311\)
\(1\)
\(15.0\)
\(5\)
\(10{,}505\)
\(8{,}192\)
\(24{,}576\)
\(2\)
\(8.0\)
\(2.5\)
\(17{,}625\)
\(16{,}130\)
\(48{,}390\)
\(4\)
\(18.0\)
\(1.25\)
\(78{,}465\)
\(70{,}226\)
\(210{,}678\)
\(8\)
\(76.0\)
Download data (CSV)
8.5
Zig-zag interface
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(132\)
\(66\)
\(1\)
\(7.2\)
\(5\)
\(1{,}127\)
\(630\)
\(258\)
\(2\)
\(18.3\)
\(2.5\)
\(5{,}203\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(57.1\)
\(1.25\)
\(29{,}963\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(203.6\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(66\)
\(1\)
\(8.6\)
\(5\)
\(2{,}025\)
\(630\)
\(258\)
\(2\)
\(20.1\)
\(2.5\)
\(7{,}605\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(55.2\)
\(1.25\)
\(37{,}485\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(213.5\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(66\)
\(1\)
\(31.6\)
\(5\)
\(3{,}267\)
\(630\)
\(258\)
\(2\)
\(39.5\)
\(2.5\)
\(10{,}575\)
\(3{,}690\)
\(1{,}026\)
\(4\)
\(79.4\)
\(1.25\)
\(46{,}023\)
\(24{,}786\)
\(4{,}098\)
\(8\)
\(270.2\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(6{,}200\)
\(4{,}026\)
\(12{,}078\)
\(1\)
\(1.0\)
\(5\)
\(7{,}442\)
\(5{,}212\)
\(15{,}636\)
\(2\)
\(1.0\)
\(2.5\)
\(13{,}514\)
\(10{,}830\)
\(32{,}490\)
\(4\)
\(3.0\)
\(1.25\)
\(51{,}154\)
\(44{,}414\)
\(133{,}242\)
\(8\)
\(6.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(6{,}220\)
\(4{,}270\)
\(12{,}810\)
\(1\)
\(3.0\)
\(5\)
\(7{,}406\)
\(5{,}790\)
\(17{,}370\)
\(2\)
\(2.0\)
\(2.5\)
\(13{,}514\)
\(12{,}498\)
\(37{,}494\)
\(4\)
\(4.0\)
\(1.25\)
\(51{,}026\)
\(50{,}004\)
\(150{,}012\)
\(8\)
\(17.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(6{,}278\)
\(4{,}629\)
\(13{,}887\)
\(1\)
\(5.0\)
\(5\)
\(8{,}150\)
\(6{,}930\)
\(20{,}790\)
\(2\)
\(8.0\)
\(2.5\)
\(16{,}002\)
\(15{,}734\)
\(47{,}202\)
\(4\)
\(18.0\)
\(1.25\)
\(69{,}778\)
\(66{,}026\)
\(198{,}078\)
\(8\)
\(72.0\)
Download data (CSV)
8.6
Nested L interface
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(132\)
\(66\)
\(1\)
\(7.7\)
\(5\)
\(1{,}127\)
\(672\)
\(298\)
\(2\)
\(19.6\)
\(2.5\)
\(5{,}203\)
\(3{,}772\)
\(1{,}106\)
\(4\)
\(59.2\)
\(1.25\)
\(29{,}963\)
\(25{,}272\)
\(4{,}578\)
\(8\)
\(223.9\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(132\)
\(66\)
\(1\)
\(9.4\)
\(5\)
\(2{,}025\)
\(672\)
\(298\)
\(2\)
\(21.3\)
\(2.5\)
\(7{,}605\)
\(3{,}772\)
\(1{,}106\)
\(4\)
\(56.2\)
\(1.25\)
\(37{,}485\)
\(25{,}272\)
\(4{,}578\)
\(8\)
\(229.1\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(132\)
\(66\)
\(1\)
\(31.1\)
\(5\)
\(3{,}267\)
\(672\)
\(298\)
\(2\)
\(41.1\)
\(2.5\)
\(10{,}575\)
\(3{,}772\)
\(1{,}106\)
\(4\)
\(87.6\)
\(1.25\)
\(46{,}023\)
\(25{,}272\)
\(4{,}578\)
\(8\)
\(363.1\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}202\)
\(4{,}539\)
\(13{,}617\)
\(1\)
\(3.0\)
\(5\)
\(9{,}148\)
\(6{,}056\)
\(18{,}168\)
\(2\)
\(2.0\)
\(2.5\)
\(14{,}108\)
\(11{,}170\)
\(33{,}510\)
\(4\)
\(3.0\)
\(1.25\)
\(55{,}176\)
\(46{,}522\)
\(139{,}566\)
\(8\)
\(7.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}160\)
\(4{,}752\)
\(14{,}256\)
\(1\)
\(2.0\)
\(5\)
\(9{,}148\)
\(6{,}652\)
\(19{,}956\)
\(2\)
\(2.0\)
\(2.5\)
\(14{,}108\)
\(12{,}840\)
\(38{,}520\)
\(4\)
\(4.0\)
\(1.25\)
\(55{,}088\)
\(52{,}136\)
\(156{,}408\)
\(8\)
\(18.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(7{,}206\)
\(5{,}105\)
\(15{,}315\)
\(1\)
\(5.0\)
\(5\)
\(9{,}596\)
\(7{,}640\)
\(22{,}920\)
\(2\)
\(8.0\)
\(2.5\)
\(17{,}956\)
\(16{,}770\)
\(50{,}310\)
\(4\)
\(19.0\)
\(1.25\)
\(75{,}792\)
\(69{,}142\)
\(207{,}426\)
\(8\)
\(72.0\)
Download data (CSV)
8.7
Tetris partition
Coreform IGA Mesh
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(325\)
\(143\)
\(66\)
\(1\)
\(15.0\)
\(5\)
\(1{,}127\)
\(777\)
\(294\)
\(2\)
\(24.1\)
\(2.5\)
\(5{,}203\)
\(3{,}977\)
\(1{,}102\)
\(4\)
\(68.5\)
\(1.25\)
\(29{,}963\)
\(26{,}325\)
\(4{,}542\)
\(8\)
\(274.8\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(735\)
\(143\)
\(66\)
\(1\)
\(16.6\)
\(5\)
\(2{,}025\)
\(777\)
\(294\)
\(2\)
\(26.6\)
\(2.5\)
\(7{,}605\)
\(3{,}977\)
\(1{,}102\)
\(4\)
\(68.2\)
\(1.25\)
\(37{,}485\)
\(26{,}325\)
\(4{,}542\)
\(8\)
\(268.4\)
\(h_z\)
Background
Active
Trimmed
CPUs
Wall-time [s]
\(10\)
\(1{,}377\)
\(143\)
\(66\)
\(1\)
\(32.3\)
\(5\)
\(3{,}267\)
\(777\)
\(294\)
\(2\)
\(115.0\)
\(2.5\)
\(10{,}575\)
\(3{,}977\)
\(1{,}102\)
\(4\)
\(97.3\)
\(1.25\)
\(46{,}023\)
\(26{,}325\)
\(4{,}542\)
\(8\)
\(390.8\)
Download data (CSV)
Coreform IGA for Abaqus solve
\(p = 1\)
\(p = 2\)
\(p = 3\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(10{,}513\)
\(6{,}895\)
\(20{,}685\)
\(1\)
\(6.0\)
\(5\)
\(13{,}721\)
\(9{,}176\)
\(27{,}528\)
\(2\)
\(3.0\)
\(2.5\)
\(20{,}069\)
\(15{,}330\)
\(45{,}990\)
\(4\)
\(3.0\)
\(1.25\)
\(72{,}117\)
\(57{,}780\)
\(173{,}340\)
\(8\)
\(8.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(10{,}419\)
\(7{,}157\)
\(21{,}471\)
\(1\)
\(8.0\)
\(5\)
\(13{,}709\)
\(9{,}975\)
\(29{,}925\)
\(2\)
\(3.0\)
\(2.5\)
\(20{,}069\)
\(17{,}399\)
\(52{,}197\)
\(4\)
\(6.0\)
\(1.25\)
\(71{,}727\)
\(64{,}654\)
\(193{,}962\)
\(8\)
\(22.0\)
\(h_z\)
Elements
Nodes
DOFs
Threads
Wall-time [s]
\(10\)
\(10{,}483\)
\(7{,}638\)
\(22{,}914\)
\(1\)
\(9.0\)
\(5\)
\(14{,}485\)
\(11{,}426\)
\(34{,}278\)
\(2\)
\(12.0\)
\(2.5\)
\(24{,}801\)
\(22{,}296\)
\(66{,}888\)
\(4\)
\(21.0\)
\(1.25\)
\(101{,}505\)
\(88{,}218\)
\(264{,}654\)
\(8\)
\(158.0\)
Download data (CSV)
9
Discussion
Because all material definitions are identical, the exact continuum response is independent of the internal cell decomposition. Differences from the monolithic result therefore measure numerical sensitivity to retained interfaces, including closed inclusions, re-entrant partition corners, diagonal interfaces, and three-material junctions. The full refinement history is retained so that a passing finest-level result can be distinguished from a coarse-grid coincidence and so regressions at individual interface topologies remain visible.
References
[1]
J. M. Gere and S. P. Timoshenko,
Mechanics of materials
, 4th ed. Boston: PWS Publishing Company, 1997.
Appendix
9.1
Capabilities exercised
The table is generated from current passing simulation and study artifacts. Each row must satisfy the linked capability card’s evidence contract; declaring a capability in a test does not, by itself, create evidence.
Table 4: Capabilities exercised by the current verification artifacts.
Capability
Evidence level
Linear static analysis
Exercised
Static stress analysis
Exercised
Boundary conditions and constraints
Exercised
Pinned boundary condition
Exercised
Threaded Abaqus/Standard execution
Exercised
Parallel Coreform IGA mesh generation
Exercised
Parallel execution
Exercised
Surface traction loading
Exercised
Mechanical loads
Exercised
Cell-based material section assignment
Exercised
Multi-cell material interfaces
Numerically validated
Linear elastic material behavior
Exercised
Materials and section assignments
Numerically validated
Meshing
Exercised
Linear C0 spline basis
Exercised
Quadratic C1 spline basis
Exercised
Cubic C2 spline basis
Exercised
Immersed rectilinear meshes
Exercised
Spline basis degree and continuity
Exercised
Explicit field-output selection
Exercised
Field output requests
Exercised
Stress field output (S)
Exercised
Displacement field output (U)
Exercised
Individual field-output variables
Exercised
Line probes
Exercised
Coreform IGA result probes
Exercised
Boolean merge with retained boundaries
Exercised
Dependent part instances
Exercised
9.2
Download artifacts
Download all verification artifacts (ZIP)
9.2.1
Geometry CAD files
No separate CAD file is required for this problem; the geometry is fully defined in the Abaqus/CAE journal files.
9.2.2
Abaqus/CAE journal files
84 generated Abaqus/CAE journals are available.
Show 84 individual journals
bottom_center_third_h0p1–multiplanar_material_interfaces_bottom_center_third_h0_model.jnl
bottom_center_third_h0p2–multiplanar_material_interfaces_bottom_center_third_h0_model.jnl
bottom_center_third_h0p3–multiplanar_material_interfaces_bottom_center_third_h0_model.jnl
bottom_center_third_h1p1–multiplanar_material_interfaces_bottom_center_third_h1_model.jnl
bottom_center_third_h1p2–multiplanar_material_interfaces_bottom_center_third_h1_model.jnl
bottom_center_third_h1p3–multiplanar_material_interfaces_bottom_center_third_h1_model.jnl
bottom_center_third_h2p1–multiplanar_material_interfaces_bottom_center_third_h2_model.jnl
bottom_center_third_h2p2–multiplanar_material_interfaces_bottom_center_third_h2_model.jnl
bottom_center_third_h2p3–multiplanar_material_interfaces_bottom_center_third_h2_model.jnl
bottom_center_third_h3p1–multiplanar_material_interfaces_bottom_center_third_h3_model.jnl
bottom_center_third_h3p2–multiplanar_material_interfaces_bottom_center_third_h3_model.jnl
bottom_center_third_h3p3–multiplanar_material_interfaces_bottom_center_third_h3_model.jnl
bottom_left_half_h0p1–multiplanar_material_interfaces_bottom_left_half_h0_model.jnl
bottom_left_half_h0p2–multiplanar_material_interfaces_bottom_left_half_h0_model.jnl
bottom_left_half_h0p3–multiplanar_material_interfaces_bottom_left_half_h0_model.jnl
bottom_left_half_h1p1–multiplanar_material_interfaces_bottom_left_half_h1_model.jnl
bottom_left_half_h1p2–multiplanar_material_interfaces_bottom_left_half_h1_model.jnl
bottom_left_half_h1p3–multiplanar_material_interfaces_bottom_left_half_h1_model.jnl
bottom_left_half_h2p1–multiplanar_material_interfaces_bottom_left_half_h2_model.jnl
bottom_left_half_h2p2–multiplanar_material_interfaces_bottom_left_half_h2_model.jnl
bottom_left_half_h2p3–multiplanar_material_interfaces_bottom_left_half_h2_model.jnl
bottom_left_half_h3p1–multiplanar_material_interfaces_bottom_left_half_h3_model.jnl
bottom_left_half_h3p2–multiplanar_material_interfaces_bottom_left_half_h3_model.jnl
bottom_left_half_h3p3–multiplanar_material_interfaces_bottom_left_half_h3_model.jnl
centered_third_h0p1–multiplanar_material_interfaces_centered_third_h0_model.jnl
centered_third_h0p2–multiplanar_material_interfaces_centered_third_h0_model.jnl
centered_third_h0p3–multiplanar_material_interfaces_centered_third_h0_model.jnl
centered_third_h1p1–multiplanar_material_interfaces_centered_third_h1_model.jnl
centered_third_h1p2–multiplanar_material_interfaces_centered_third_h1_model.jnl
centered_third_h1p3–multiplanar_material_interfaces_centered_third_h1_model.jnl
centered_third_h2p1–multiplanar_material_interfaces_centered_third_h2_model.jnl
centered_third_h2p2–multiplanar_material_interfaces_centered_third_h2_model.jnl
centered_third_h2p3–multiplanar_material_interfaces_centered_third_h2_model.jnl
centered_third_h3p1–multiplanar_material_interfaces_centered_third_h3_model.jnl
centered_third_h3p2–multiplanar_material_interfaces_centered_third_h3_model.jnl
centered_third_h3p3–multiplanar_material_interfaces_centered_third_h3_model.jnl
monolithic_h0p1–multiplanar_material_interfaces_monolithic_h0_model.jnl
monolithic_h0p2–multiplanar_material_interfaces_monolithic_h0_model.jnl
monolithic_h0p3–multiplanar_material_interfaces_monolithic_h0_model.jnl
monolithic_h1p1–multiplanar_material_interfaces_monolithic_h1_model.jnl
monolithic_h1p2–multiplanar_material_interfaces_monolithic_h1_model.jnl
monolithic_h1p3–multiplanar_material_interfaces_monolithic_h1_model.jnl
monolithic_h2p1–multiplanar_material_interfaces_monolithic_h2_model.jnl
monolithic_h2p2–multiplanar_material_interfaces_monolithic_h2_model.jnl
monolithic_h2p3–multiplanar_material_interfaces_monolithic_h2_model.jnl
monolithic_h3p1–multiplanar_material_interfaces_monolithic_h3_model.jnl
monolithic_h3p2–multiplanar_material_interfaces_monolithic_h3_model.jnl
monolithic_h3p3–multiplanar_material_interfaces_monolithic_h3_model.jnl
nested_l_h0p1–multiplanar_material_interfaces_nested_l_h0_model.jnl
nested_l_h0p2–multiplanar_material_interfaces_nested_l_h0_model.jnl
nested_l_h0p3–multiplanar_material_interfaces_nested_l_h0_model.jnl
nested_l_h1p1–multiplanar_material_interfaces_nested_l_h1_model.jnl
nested_l_h1p2–multiplanar_material_interfaces_nested_l_h1_model.jnl
nested_l_h1p3–multiplanar_material_interfaces_nested_l_h1_model.jnl
nested_l_h2p1–multiplanar_material_interfaces_nested_l_h2_model.jnl
nested_l_h2p2–multiplanar_material_interfaces_nested_l_h2_model.jnl
nested_l_h2p3–multiplanar_material_interfaces_nested_l_h2_model.jnl
nested_l_h3p1–multiplanar_material_interfaces_nested_l_h3_model.jnl
nested_l_h3p2–multiplanar_material_interfaces_nested_l_h3_model.jnl
nested_l_h3p3–multiplanar_material_interfaces_nested_l_h3_model.jnl
tetris_h0p1–multiplanar_material_interfaces_tetris_h0_model.jnl
tetris_h0p2–multiplanar_material_interfaces_tetris_h0_model.jnl
tetris_h0p3–multiplanar_material_interfaces_tetris_h0_model.jnl
tetris_h1p1–multiplanar_material_interfaces_tetris_h1_model.jnl
tetris_h1p2–multiplanar_material_interfaces_tetris_h1_model.jnl
tetris_h1p3–multiplanar_material_interfaces_tetris_h1_model.jnl
tetris_h2p1–multiplanar_material_interfaces_tetris_h2_model.jnl
tetris_h2p2–multiplanar_material_interfaces_tetris_h2_model.jnl
tetris_h2p3–multiplanar_material_interfaces_tetris_h2_model.jnl
tetris_h3p1–multiplanar_material_interfaces_tetris_h3_model.jnl
tetris_h3p2–multiplanar_material_interfaces_tetris_h3_model.jnl
tetris_h3p3–multiplanar_material_interfaces_tetris_h3_model.jnl
zig_zag_h0p1–multiplanar_material_interfaces_zig_zag_h0_model.jnl
zig_zag_h0p2–multiplanar_material_interfaces_zig_zag_h0_model.jnl
zig_zag_h0p3–multiplanar_material_interfaces_zig_zag_h0_model.jnl
zig_zag_h1p1–multiplanar_material_interfaces_zig_zag_h1_model.jnl
zig_zag_h1p2–multiplanar_material_interfaces_zig_zag_h1_model.jnl
zig_zag_h1p3–multiplanar_material_interfaces_zig_zag_h1_model.jnl
zig_zag_h2p1–multiplanar_material_interfaces_zig_zag_h2_model.jnl
zig_zag_h2p2–multiplanar_material_interfaces_zig_zag_h2_model.jnl
zig_zag_h2p3–multiplanar_material_interfaces_zig_zag_h2_model.jnl
zig_zag_h3p1–multiplanar_material_interfaces_zig_zag_h3_model.jnl
zig_zag_h3p2–multiplanar_material_interfaces_zig_zag_h3_model.jnl
zig_zag_h3p3–multiplanar_material_interfaces_zig_zag_h3_model.jnl
9.3
Supported Abaqus keywords
The deck below is reproduced from the monolithic
\(h_0p_2\)
test with comments stripped and instance-specific names redacted.
Note
multiplanar_material_interfaces_monolithic_h0.inp - outline
*Heading
*Preprint
, echo
=NO
, model
=NO
, history
=NO
, contact
=NO
*Include
, input
=
<iga-mesh-inp>
*System
*Material
, name
=
<material>
*Elastic
1e+11, 0.3
*User Output Variables
1,
*Boundary
<node-set>
, PINNED
*Step
, name
=
<step>
, nlgeom
=NO
*Static
1., 1., 1e-05, 1.
*Restart
, write, frequency
=0
*Output
, field
*Node Output
U,
*Element Output
, position
=AVERAGED AT NODES
, directions
=YES
S,
*Node Output
U
*Output
, history, variable
=PRESELECT
*End Step
9.4
Provenance
The records below identify the execution environment and verification basis for every simulation and study CTest reported on this page.
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:29
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:45:01
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:12:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:10:06
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:38:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:44:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:03:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:00:29
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:21:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:40:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:52:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_center_third_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:51:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:45:26
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:12:20
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:03:15
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:38:23
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:44:19
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:05:05
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:03:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:23:23
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:43:00
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:55:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_bottom_left_half_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:41:59
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:00
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:42:50
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:11:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:56:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:36:11
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:41:10
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:11:05
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:55:41
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:22:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:35:30
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:50:16
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_centered_third_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:36:31
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:40:08
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:46:28
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:13:16
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:09:53
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:39:58
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:46:25
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:13:40
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:09:55
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:38:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:46:11
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:11:16
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_monolithic_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:58:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:37:56
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:43:17
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:11:42
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:55:40
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:35:48
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:43:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:03:40
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:55:12
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:20:22
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:38:26
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:51:53
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_nested_l_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:37:57
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:25:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:40:04
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:01:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:54:13
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:24:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:39:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:01:28
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:52:35
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:20:21
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:35:09
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:50:00
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_tetris_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:23:25
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h0p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:39:05
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h1p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (0 positive)
Test date
04-Nov-2026 19:45:01
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h2p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:13:18
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h3p1
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:09:32
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h0p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:37:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h1p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:43:24
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h2p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 20:02:34
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h3p2
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:59:11
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h0p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:21:16
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h1p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:40:36
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h2p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:57:26
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_multiplanar_material_interfaces_zig_zag_h3p3
Parameter
Value
Evidence type
Record-only diagnostic
CTest result
Passed (Release)
Verification result
Not asserted — 1 record-only check (1 positive)
Test date
04-Nov-2026 19:45:14
Operating system
Ubuntu 24.04.5 LTS
Coreform IGA version
2026.11+72575
Source revision
e48e5098b2f8
Abaqus version
Abaqus 6.26-4
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_bottom_center_third_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:13 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_bottom_center_third_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_bottom_center_third_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:13 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_bottom_left_half_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_bottom_left_half_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_bottom_left_half_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_centered_third_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_centered_third_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_centered_third_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_nested_l_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_nested_l_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_nested_l_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_tetris_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_tetris_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_tetris_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:15 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_zig_zag_p1
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_zig_zag_p2
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:14 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
Note
Test_abaqus_verification_study_multiplanar_material_interfaces_zig_zag_p3
Parameter
Value
Evidence type
Numerical verification
CTest result
Passed (Release)
Verification result
Passed — 1/1 enforced criterion passed
Test date
04-Nov-2026 20:14:13 UTC
Operating system
—
Coreform IGA version
Not applicable — postprocessing-only CTest
Source revision
e48e5098b2f8
Abaqus version
Not applicable — postprocessing-only CTest
