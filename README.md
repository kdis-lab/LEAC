# LEAC - Library Evolutionary Algorithms for Clustering <img align="right" width="100" height="100" src="leac_logo.png">

Library Evolutionary Algorithms for Clustering (**LEAC**) is a library for the implementation
of evolutionary and genetic algorithms to solve the problem of *partition clustering*.
The tools is free and open source, under the GNU General Public
[GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) 
[\[License file\]](../../LICENSE),
for Windows&reg;, GNU/Linux&reg; and Mac OS X&reg;.
For more details see <a href="https://github.com/kdis-lab/leac/tree/master/leac-userManual.pdf">user
manual</a> and <a href="https://hbrobles.github.io/APILeac/index.html">API</a> or locally, the API can be found in the document directory

Clustering is useful in several exploratory *pattern-analysis*,
*grouping*, *decision-making*, and *machine-learning situations*,
including *data mining*, *document retrieval*, *image segmentation*,
and *pattern classification* <a href="http://doi.acm.org/10.1145/331499.331504">Jain et al</a>.

![](../master/leac_cluster.svg)

LEAC is built in C ++ with the current standards of C++11, C++14.
Take advantage of the STL library, through the use of 
data structures and algorithms, for example,
Numerics library of Pseudo-random number generation which is 
fundamental in the genetic and evolutionary algorithms.
An additional feature that LEAC has is a low-level software layer
based on Streaming SIMD Extensions (SSE) and in the 
[OpenBLAS](http://www.openblas.net) library ,in order to
increase performance. For now only for Linux x86_64, for future
versions it is intended to port to other architectures and
increase the number of genetic operators implemented with
these software layers

LEAC also includes several implementations of evolutionary algorithms for
partitional clustering which are based on state-of-art of evolutionary
proposals in this area:

<table>
	<tr>
	  <th>Encode</th>
	  <th>Fixed-K</th>
	  <th>Variable-K</th>
	</tr>
	<tr>
	  <td>Label</td>
	  <td> <a href="http://dx.doi.org/10.1016/0167-8655(96)00043-8">gaclustering_fklabel</a>,
	    <a href="http://dx.doi.org/10.1109/3477.764879">gka_fklabel</a>,
	    <a href="http://dx.doi.org/10.1186/1471-2105-5-172">igka_fklabel</a>,
	    <a href="http://doi.acm.org/10.1145/967900.968029">fgka_fklabel</a>
	  </td>
	  <td> <a href="http://dx.doi.org/10.1016/j.eswa.2012.02.149">gga_vklabeldbindex and gga_vklabelsilhouette</a>,
	    <a href="http://dl.acm.org/citation.cfm?id=1293920.1293922">cga_vklabel</a>,
	    <a href="http://dx.doi.org/10.1016/j.ins.2005.07.015">eac_vklabel</a>,
	    <a href="http://dx.doi.org/10.1109/CEC.2006.1688522">eaci_vklabel</a>,
	    <a href="http://dx.doi.org/10.1109/CEC.2006.1688522">eacii_vklabel</a>,
	    <a href="http://dx.doi.org/10.1109/CEC.2006.1688522">eaciii_vklabel</a>,
	    <a href="http://dx.doi.org/10.1109/CEC.2006.1688522">feac_vklabel</a>
	  </td>
	</tr>
	<tr>
	  <td>Crisp-matrix</td>
	  <td> <a href="http://dx.doi.org/10.1109/ICEC.1994.350046">gaclustering_fkcrispmatrix</a>
	  </td>
	  <td></td>
	</tr>
	<tr>
	  <td>Centroid</td>
	  <td> <a href="http://dx.doi.org/10.1016/S0031-3203(99)00137-5">gas_fkcentroid</a>,
	    <a href="http://dx.doi.org/10.1016/S0020-0255(02)00208-6">kga_fkcentroid</a>,
	    <a href="http://dx.doi.org/10.1016/j.patcog.2008.11.006">gagr_fkcentroid</a>,
	    <a href="http://dx.doi.org/10.1093/comjnl/40.9.547">cbga_fkcentroid</a>
	  </td>
	  <td> <a href="http://dx.doi.org/10.1016/S0031-3203(01)00108-X">gcuk_vkcentroid</a>,
	    <a href="http://dx.doi.org/10.1016/j.neucom.2011.11.001">tgca_vkcentroid</a>
	  </td>
	</tr>
	<tr>
	  <td>Medoid</td>
	  <td> gaprototypes_fkmedoid,
	    <a href="http://dx.doi.org/10.1109/CEC.2004.1330840">hka_fkmedoid</a>,
	    <a href="https://doi.org/10.1016/0003-2670(93)80130-D">gca_fkmedoid</a>
	  </td>
	  <td></td>
	</tr>
	<tr>
	  <td>Tree</td>
	  <td></td>
	  <td>  <a href="http://dx.doi.org/10.1007/978-3-540-39398-6_7">gaclustering_vktreebinary</a> 
	</tr>
	<tr>
	  <td>Sub-cluster</td>
	  <td></td>
	  <td> <a href="http://dx.doi.org/10.1016/S0031-3203(00)00005-4">clustering_vksubclusterbinary</a> 
	</tr>
    </table>

