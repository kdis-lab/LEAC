# LEAC - Library Evolutionary Algorithms for Clustering <img align="right" width="100" height="100" src="leac_logo.png">

* [What is LEAC library?](https://github.com/kdis-lab/LEAC/blob/master/README.md#what-is-leac-library)
* [Repository description](https://github.com/kdis-lab/LEAC/blob/master/README.md#repository-description)
* [Getting the library](https://github.com/kdis-lab/LEAC/blob/master/README.md#getting-the-library)
* [Tutorials and documentation](https://github.com/kdis-lab/LEAC/blob/master/README.md#tutorials-and-documentation)
* [Methods included](https://github.com/kdis-lab/LEAC/blob/master/README.md#methods-included)
* [LEAC library data format](https://github.com/kdis-lab/LEAC/blob/master/README.md#leac-library-data-format)
* [Citation](https://github.com/kdis-lab/LEAC/blob/master/README.md#citation)
* [License](https://github.com/kdis-lab/LEAC/blob/master/README.md#license)
* [Reporting bugs](https://github.com/kdis-lab/LEAC/blob/master/README.md#reporting-bugs)


## What is LEAC library?
Library Evolutionary Algorithms for Clustering (**LEAC**) is a library of genetic algorithms to solve the problem of *partition clustering*. It includes 22 classification genetic algorithms for solving partitional clustering. These are considered the state-of-the-art of mono-objective genetic algorithms in the area.

LEAC is a modular library which make easier to develop new evolutionary algorithm proposals for solving the partitional clustering. Besides including the most representative proposals of Evolutionary Algorithms for partitional clustering. LEAC allows you to implement easily new evolutionary algorithm proposals for clustering using the classes developed in the library. Thus, new algorithms can use the diversity of proposed strategies, procedures and genetic operators included in the library to evolve the population according to the flowchart of these algorithms.

LEAC is built in C ++ with the current standards of C++11, C++14. Take advantage of the STL library, through the use of  data structures and algorithms, for example, Numerics library of Pseudo-random number generation which is  fundamental in the genetic and evolutionary algorithms.
An additional feature that LEAC has is a low-level software layer based on Streaming SIMD Extensions (SSE) and in the  [OpenBLAS](http://www.openblas.net) library ,in order to increase performance. For now only for Linux x86_64, for future versions it is intended to port to other architectures and increase the number of genetic operators implemented with these software layers.

## LEAC description

LEAC library is based on a layered software architecture composed of four layers: algorithms, EA, Clustering and Performance. Each layer consists of a set of related packets as shown in the image.
<img width="686" alt="image" src="https://github.com/kdis-lab/LEAC/assets/37608799/8171595c-595e-4199-9765-03f7d0e019f7">

* Algorithm: It contains final implementations of several evolutionary algorithms for clustering.
* EA: It contains several packets with operators and strategies to configure EAs, such as: encoding criteria, initialization methods, selection methods, crossover and mutation operators and updating and replacement strategies.
* Clustering: It contains several packets with specific clustering operators, such as, supervised and unsupervised performance measures and clustering operators based on centroids, crispmatrix or medoids.
* Performance: It consists of low-level programmed functions under the current CPU architectures.


## Getting the library


## Tutorials and documentation
For more details see <a href="https://github.com/kdis-lab/leac/tree/master/leac-userManual.pdf">user
manual</a> and <a href="https://hbrobles.github.io/APILeac/index.html">API</a> or locally, the API can be found in the 'docs' directory


## Methods included
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
	    <a href="http://dx.doi.org/10.1109/CEC.2006.1688522">feac_vklabelssilhouette</a>,
	    <a href="http://dx.doi.org/10.1109/CEC.2006.1688522">feac_vklabelrandindex</a>
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

## LEAC library data format

The format of data is based on the Weka's format.


## Citation

Robles-Berumen, H., Zafra, A., Fardoun, H. M., & Ventura, S. (2019). LEAC: An efficient library for clustering with evolutionary algorithms. *Knowledge-Based Systems*, *179*, 117-119. https://doi.org/10.1016/j.knosys.2019.05.008

## License

The tools is free and open source, under the GNU General Public [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) 
[\[License file\]](../../LICENSE), for Windows&reg;, GNU/Linux&reg; and Mac OS X&reg;.

## Reporting bugs
Feel free to open an [issue](https://github.com/kdis-lab/LEAC/issues) at Github if anything is not working as expected. Merge request are also encouraged, it will be carefully reviewed and merged if everything is all right.


