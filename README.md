<center><h2>Exterior Poisson Reconstruction (Version 1.00)</h2></center>
<center>
<a href="#LINKS">links</a>
<!--
<a href="#COMPILATION">compilation</a>
-->
<a href="#EXECUTABLES">executables</a>
<a href="#USAGE">usage</a>
<a href="#CHANGES">changes</a>
</center>
<hr>
This software supports reconstruction of co-dimension two manifolds.
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
<b>Papers:</b>
<a href="https://www.cs.jhu.edu/~misha/MyPapers/SGP06.pdf">[Kazhdan, Bolitho, and Hoppe, 2006]</a>,
<br>
<b>Executables: </b>
<a href="https://www.cs.jhu.edu/~misha/Code/ExteriorPoissonRecon/Version1.00/ExteriorPoissonRecon.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="https://www.cs.jhu.edu/~misha/Code/ExteriorPoissonRecon/Version1.00/ExteriorPoissonRecon.zip">ZIP</a> <a href="https://github.com/mkazhdan/ExteriorPoissonRecon">GitHub</a><br>
<b>Older Versions:</b>
<!--
<a href="https://www.cs.jhu.edu/~misha/Code/ExteriorPoissonRecon/Version1.00/">V1.00</a>,
-->
</ul>

<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Sample</b></font>:
Generates a sampling of points and frames (encoded as the entries of a skew-symmetric matrix) from co-dimension two manifolds in 3D and 4D.
</SUMMARY>

<dt><b>--type</b> &lt;<i>input geometry type</i>&gt;</dt>
<dd>
This string specifies the type of geometry the points should be sampled from. Supported types include:
<UL>
<LI><code>line_segment</code>: Points lie on a (straight) line segment
<LI><code>circle</code>: Points lie on a circle
<LI><code>link</code>: Points lie on two interlocking circles
<LI><code>spiral:&lt;r&gt;</code>: Points lie on a spiral with <code>r</code> rotations
<LI><code>torus_knot:&lt;p&gt;:&lt;q&gt;</code>: Points lie on a (<code>p</code>,<code>q</code>) torus-knot 
<LI><code>borromean_rings</code>: Points lie on interlocking Borromean rings
<LI><code>clifford_torus</code>: Points lie on the Clifford torus
<LI><code>hopf_torus:&lt;n&gt;:&lt;a&gt;</code>: Points lie on the Hopf torus with <code>n</code> nodes and amplitude <code>a</code>.<BR>
Reasonable values for amplitude are in the range [0.1,0.5].
</UL>
</dd>

<dt>[<b>--out</b> &lt;<i>output file name</i>&gt;]</dt>
<dd>
This optional string value specifies the name of the file to which the samples will be written.<br>
The file will be written out in <a href="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format,
with x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i> and the orientation of the sample given by the (linearized) coefficients of a skew-symmetric matrix,
encoded by the properties <i>skew0</i>,...,<i>skew&lt;n&gt;</i> with <i>n=2</i> for curves in 3D and <i>n=5</i> for surfaces in 4D.<br>
If this argument is not provided, no output is generated.
</dd>

<dt>[<b>--res</b> &lt;<i>sample resolution</i>&gt;]</dt>
<dd> This optional integer value specifies the resolution of the sampling.<BR>
The default value for this parameter is 1024.
</dd>

<dt>[<b>--aNoise</b> &lt;<i>angular noise</i>&gt;]</dt>
<dd> This optional floating point value specifies the maximum amount of noise in the samples' orientations (in units of radians).<BR>
The default value for this parameter is 0.
</dd>

<dt>[<b>--pNoise</b> &lt;<i>positional noise</i>&gt;]</dt>
<dd> This optional floating point value specifies the maximum amount of noise in the samples' positions (in units of voxels).<BR>
The default value for this parameter is 0.
</dd>

<dt>[<b>--regular</b>]</dt>
<dd>If enabled, samples will be obtained by regularly sampling in parameter space.</dd>

</DETAILS>
</dl>
</ul>

<!--------------------->

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>ExteriorPoissonRecon</b></font>:
Reconstructs an implicit function, sampled over a regular grid with two values per grid cell, from a sampling of points and frames (encoded as the entries of a skew-symmetric matrix).
</SUMMARY>

<dt><b>--in</b> &lt;<i>dimension, input points and frames</i>&gt;</dt>
<dd>
This integer/string pair value specifies the dimension in which the points are embedded and the name of the file containing the points.
</dd>

<dt>[<b>--out</b> &lt;<i>grid header</i>&gt;]</dt>
<dd>
This optional string value specifies the header for the grid files describing the estimated density distribution and the reconstructed implicit function.<br>
The density will be output to the file <code>&lt;grid header&gt;.density.grid</code> and the reconstructed implicit function will be output to the file <code>&lt;grid header&gt;.grid</code>.<br>
If this argument is not provided, no output is generated.
</dd>

<dt>[<b>--depth</b> &lt;<i>reconstruction depth</i>&gt;]</dt>
<dd>
This optional integer value is the depth of the grid that will be used for reconstruction.
Running at depth <i>d</i> corresponds to solving on a grid whose resolution is than <i>2^D x 2^d x ... </i>.<br>
The default value for this parameter is 5.
</dd>

<dt>[<b>--sWeight</b> &lt;<i>screening weight</i>&gt;]</dt>
<dd>
This optional floating point value is the screening weight used for reconstruction.<br>
The default value for this parameter is 50.
</dd>

<dt>[<b>--dWeight</b> &lt;<i>Dirichlet weight</i>&gt;]</dt>
<dd>
This optional floating point value is the Dirichlet weight used for reconstruction.<br>
The default value for this parameter is 0.003125.
</dd>

<dt>[<b>--scale</b> &lt;<i>scale factor</i>&gt;]</dt>
<dd>
This optional floating point value specifies the ratio between the diameter of the cube used for reconstruction and the diameter of the samples' bounding cube.<br>
The default value is 1.1.
</dd>

<dt>[<b>--verbosity</b> &lt;<i>verbosity</i>&gt;</b>]
<dd>
This optional integer value specifies the level of verbosity of the executable's output to the command prompt.
<UL>
<LI>0: No ooutput
<LI>1: Global residual error
<LI>2: Residual error after each level of the multigrid hierarchy 
</UL>
</dd>

</DETAILS>
</dl>
</ul>

<!--------------------->

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Visualize3D</b></font>:
Extracts a curve from the reconstructed implicit function.
</SUMMARY>

<dt><b>--in</b> &lt;<i>input implicit grid</i>&gt;</dt>
<dd>
This string value is the file-name of the grid sampling the reconstructed implicit function.
</dd>

<dt>[<b>--density</b> &lt;<i>input density grid</i>&gt;]</dt>
<dd>
This optional string value is the file-name of the grid sampling the density values.<br>
If this argument is not provided, no density-based trimming is performed.
</dd>

<dt>[<b>--out</b> &lt;<i>output curve</i>&gt;]</dt>
<dd>
This optional string value specifies the file to which the extracted level-set will be written.<br>
The file will be written out in <a href="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, 
with x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i>.<br>
If this argument is not provided, no output is generated.
</dd>

<dt>[<b>--trimDensity</b> &lt;<i>trimming density</i>&gt;]</dt>
<dd>
This optional non-negative floating point value specifies the density that must be met by some point on a connected component of the reconstruction for the connected component to be kept.<br>
The default value for this argument is 0.0.
</dd>

</DETAILS>
</dl>
</ul>

<!--------------------->

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Visualize4D</b></font>:
Extracts a surface from the reconstructed implicit function.
</SUMMARY>

<dt><b>--in</b> &lt;<i>input implicit grid</i>&gt;</dt>
<dd>
This string value is the file-name of the grid sampling the reconstructed implicit function.
</dd>

<dt>[<b>--density</b> &lt;<i>input density grid</i>&gt;]</dt>
<dd>
This optional string value is the file-name of the grid sampling the density values.<br>
If this argument is not provided, no density-based trimming is performed.
</dd>

<dt>[<b>--out</b> &lt;<i>output mesh</i>&gt;]</dt>
<dd>
This optional string value specifies the file to which the extracted level-set will be written.<br>
The file will be written out in <a href="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, 
with x-, y-, z-, and w-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, <i>z</i>, and <i>w</i>.<br>
If this argument is not provided, no output is generated.
</dd>

<dt>[<b>--trimDensity</b> &lt;<i>trimming density</i>&gt;]</dt>
<dd>
This optional non-negative floating point value specifies the density that must be met by some point on a connected component of the reconstruction for the connected component to be kept.<br>
The default value for this argument is 0.0.
</dd>

</DETAILS>
</dl>
</ul>

<!--------------------->

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Stereo</b></font>:
Replaces the vertices of a mesh in 4D with their positions in 3D after stereographic projection.
</SUMMARY>

<dt><b>--in</b> &lt;<i>input 4D mesh</i>&gt;</dt>
<dd>
This string value is the file-name of the input (4D) mesh.<br>
The file is assumed to be in <a href="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, 
with x-, y-, z-, and w-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, <i>z</i>, and <i>w</i>.
</dd>

<dt>[<b>--out</b> &lt;<i>output 3D mesh</i>&gt;]</dt>
<dd>
This string value is the file-name of the output (3D) mesh.<br>
The will be written in <a href="https://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, 
with x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i>.<br>
If this argument is not provided, no output is generated.
</dd>

<dt>[<b>--stereo</b> &lt;<i>x, y, z, and w</i>&gt;]</dt>
<dd>
This optional quadruple of floating point values specifies the 4D axis of stereographic projection.<br>
The default value for this argument is 0, 0, 0, 1.
</dd>

</DETAILS>
</dl>
</ul>


<hr>
<a name="USAGE"><b>USAGE EXAMPLES</b></a><br>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Sample / ExteriorPoissonRecon / Visualize3D</b></font>
</SUMMARY>
For testing purposes, four point sets are provided:
<ol>

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/horse.npts"><b>Horse</b></a>:
A set of 100,000 oriented point samples (represented in ASCII format) was obtained by sampling a virtual horse model with a sampling density proportional to curvature, giving a set of non-uniformly distributed points.<br>
The surface of the model can be reconstructed by calling the either Poisson surface reconstructor:
<blockquote><code>% PoissonRecon --in horse.npts --out horse.ply --depth 10</code></blockquote>
or the SSD surface reconstructor
<blockquote><code>% SSDRecon --in horse.npts --out horse.ply --depth 10</code></blockquote>
</li>

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/bunny.points.ply"><b>Bunny</b></a>:
A set of 362,271 oriented point samples (represented in PLY format) was obtained by merging the data from the original Stanford Bunny
<a href="ftp://graphics.stanford.edu/pub/3Dscanrep/bunny.tar.gz">range scans</a>. The orientation of the sample points was estimated
using the connectivity information within individual range scans.<br>
The original (unscreened) Poisson reconstruction can be obtained by setting the point interpolation weight to zero:
<blockquote><code>% PoissonRecon --in bunny.points.ply --out bunny.ply --depth 10 --pointWeight 0</code></blockquote>
By default, the Poisson surface reconstructor uses degree-2 B-splines. A more efficient reconstruction can be obtained using degree-1 B-splines:
<blockquote><code>% PoissonRecon --in bunny.points.ply --out bunny.ply --depth 10 --pointWeight 0 --degree 1</code></blockquote>
(The SSD reconstructor requires B-splines of degree at least 2 since second derivatives are required to formulate the bi-Laplacian energy.)
</li>

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/eagle.points.ply"><b>Eagle</b></a>:
A set of 796,825 oriented point samples with color (represented in PLY format) was obtained in the EPFL <a href="https://lgg.epfl.ch/statues.php">Scanning 3D Statues from Photos</a> course.<br>
A reconstruction of the eagle can be obtained by calling:
<blockquote><code>% PoissonRecon --in eagle.points.ply --out eagle.pr.ply --depth 10</code></blockquote>
(with the RGBA color properties automatically detected from the .ply header).<BR>
A reconstruction of the eagle that does not close up the holes can be obtained by first calling:
<blockquote><code>% SSDRecon --in eagle.points.ply --out eagle.ssd.ply --depth 10 --density</code></blockquote>
using the <b>--density</b> flag to indicate that density estimates should be output with the vertices of the mesh, and then calling:
<blockquote><code>% SurfaceTrimmer --in eagle.ssd.ply --out eagle.ssd.trimmed.ply --trim 7</code></blockquote>
to remove all subsets of the surface where the sampling density corresponds to a depth smaller than 7.<BR>
This reconstruction can be chunked into cubes of size 4&times;4&times;4 by calling:
<blockquote><code>% ChunkPly --in 1 eagle.ssd.trimmed.ply --out eagle.ssd.trimmed.chnks --width 4</code></blockquote>
which partitions the reconstruction into 11 pieces.

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/torso.zip"><b>Torso</b></a>:
A set of 3,488,432 (torso.points.ply) and an envelope (torso.envelope.ply).<br>
A reconstruction of the torso that constrains the reconstruction to be contained within the envelope can be obtained by calling:
<blockquote><code>% PoissonRecon --in torso.points.ply --envelope torso.envelope.ply --out torso.pr.ply --depth 10</code></blockquote>
using the <b>--envelope</b> flag to specify the water-tight mesh constraining the reconstruction.<BR>
</li>

</ol>

</DETAILS>
</dl>
</ul>

<!--------------------->

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Sample / ExteriorPoissonRecon / Visualize4D / Stereo</b></font>
</SUMMARY>
For testing purposes, four point sets are provided:
<ol>

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/horse.npts"><b>Horse</b></a>:
A set of 100,000 oriented point samples (represented in ASCII format) was obtained by sampling a virtual horse model with a sampling density proportional to curvature, giving a set of non-uniformly distributed points.<br>
The surface of the model can be reconstructed by calling the either Poisson surface reconstructor:
<blockquote><code>% PoissonRecon --in horse.npts --out horse.ply --depth 10</code></blockquote>
or the SSD surface reconstructor
<blockquote><code>% SSDRecon --in horse.npts --out horse.ply --depth 10</code></blockquote>
</li>

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/bunny.points.ply"><b>Bunny</b></a>:
A set of 362,271 oriented point samples (represented in PLY format) was obtained by merging the data from the original Stanford Bunny
<a href="ftp://graphics.stanford.edu/pub/3Dscanrep/bunny.tar.gz">range scans</a>. The orientation of the sample points was estimated
using the connectivity information within individual range scans.<br>
The original (unscreened) Poisson reconstruction can be obtained by setting the point interpolation weight to zero:
<blockquote><code>% PoissonRecon --in bunny.points.ply --out bunny.ply --depth 10 --pointWeight 0</code></blockquote>
By default, the Poisson surface reconstructor uses degree-2 B-splines. A more efficient reconstruction can be obtained using degree-1 B-splines:
<blockquote><code>% PoissonRecon --in bunny.points.ply --out bunny.ply --depth 10 --pointWeight 0 --degree 1</code></blockquote>
(The SSD reconstructor requires B-splines of degree at least 2 since second derivatives are required to formulate the bi-Laplacian energy.)
</li>

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/eagle.points.ply"><b>Eagle</b></a>:
A set of 796,825 oriented point samples with color (represented in PLY format) was obtained in the EPFL <a href="https://lgg.epfl.ch/statues.php">Scanning 3D Statues from Photos</a> course.<br>
A reconstruction of the eagle can be obtained by calling:
<blockquote><code>% PoissonRecon --in eagle.points.ply --out eagle.pr.ply --depth 10</code></blockquote>
(with the RGBA color properties automatically detected from the .ply header).<BR>
A reconstruction of the eagle that does not close up the holes can be obtained by first calling:
<blockquote><code>% SSDRecon --in eagle.points.ply --out eagle.ssd.ply --depth 10 --density</code></blockquote>
using the <b>--density</b> flag to indicate that density estimates should be output with the vertices of the mesh, and then calling:
<blockquote><code>% SurfaceTrimmer --in eagle.ssd.ply --out eagle.ssd.trimmed.ply --trim 7</code></blockquote>
to remove all subsets of the surface where the sampling density corresponds to a depth smaller than 7.<BR>
This reconstruction can be chunked into cubes of size 4&times;4&times;4 by calling:
<blockquote><code>% ChunkPly --in 1 eagle.ssd.trimmed.ply --out eagle.ssd.trimmed.chnks --width 4</code></blockquote>
which partitions the reconstruction into 11 pieces.

<li> <a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/torso.zip"><b>Torso</b></a>:
A set of 3,488,432 (torso.points.ply) and an envelope (torso.envelope.ply).<br>
A reconstruction of the torso that constrains the reconstruction to be contained within the envelope can be obtained by calling:
<blockquote><code>% PoissonRecon --in torso.points.ply --envelope torso.envelope.ply --out torso.pr.ply --depth 10</code></blockquote>
using the <b>--envelope</b> flag to specify the water-tight mesh constraining the reconstruction.<BR>
</li>

</ol>

</DETAILS>
</dl>
</ul>




<hr>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><font size="+1"><b><B>HISTORY OF CHANGES</B></b></font></A>
</SUMMARY>
<a href="https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version1.00/">Version 1.00</a>:
<ol>
<li>Initial release
</ol>

</DETAILS>

<!--
<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
<UL>
<LI>This work was genersouly supported by the National Science Foundation (NSF) grant numbers <A HREF="https://www.nsf.gov/awardsearch/showAward?AWD_ID=0746039">0746039</A> and <A HREF="https://www.nsf.gov/awardsearch/showAward?AWD_ID=1422325">1422325</A>.
<LI>We are extremely grateful to the EPFL <a href="https://lgg.epfl.ch/statues.php">Scanning 3D Statues from Photos</a> course, the <A HREF="http://graphics.stanford.edu/data/3Dscanrep/">Stanford 3D Scanning Repository</A>, and <A HREF="https://www.resonai.com/">Resonai</A> for sharing their data.
<LI>This work was carried out at the <A HREF="https://www.arch.jhu.edu/">Advanced Research Computing at Hopkins (ARCH)</A> core facility, which is supported by the National Science Foundation (NSF) grant number <A HREF="https://www.nsf.gov/awardsearch/showAward?AWD_ID=1920103">1920103</A>.
</UL>
--<
