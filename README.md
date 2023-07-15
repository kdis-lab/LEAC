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

LEAC is a modular library that facilitates the development of new proposals for evolutionary algorithms to solve partitional clustering.  Besides including the most representative proposals of Evolutionary Algorithms for partitional clustering. LEAC allows you to implement easily new evolutionary algorithm proposals for clustering using the classes developed in the library. Thus, new algorithms can use the diversity of proposed strategies, procedures, and genetic operators included in the library to evolve the population according to the flowchart of these algorithms.

LEAC is built in C ++ with the current standards of C++11, C++14. Take advantage of the STL library, through the use of data structures and algorithms, for example, the Numerics library of Pseudo-random number generation which is fundamental in the genetic and evolutionary algorithms. An additional feature that LEAC has is a low-level software layer based on Streaming SIMD Extensions (SSE) and in the  [OpenBLAS](http://www.openblas.net) library, in order to increase performance. For now, only for Linux x86_64, for future versions, it is intended to port to other architectures and increase the number of genetic operators implemented with these software layers.

## LEAC description

LEAC library is based on a layered software architecture composed of four layers: algorithms, EA, Clustering, and Performance. Each layer consists of a set of related packets as shown in the image.
<img width="686" alt="image" src="https://github.com/kdis-lab/LEAC/assets/37608799/8171595c-595e-4199-9765-03f7d0e019f7">

* Algorithm: It contains the final implementations of several evolutionary algorithms for clustering.
* EA: It contains several packets with operators and strategies to configure EAs, such as encoding criteria, initialization methods, selection methods, crossover, and mutation operators, and updating and replacement strategies.
* Clustering: It contains several packets with specific clustering operators, such as supervised and unsupervised performance measures and clustering operators based on centroids, crisp matrix or medoids.
* Performance: It consists of low-level programmed functions under the current CPU architectures.


## Getting the library

The project can be download in this repository. Detailed information about getting and running the library can be consulted in [user's manual](https://github.com/kdis-lab/LEAC/blob/master/leac-userManual.pdf) where all dependencies and libraries necessaries are commented.

## Tutorials and documentation

[LEAC user's manual](https://github.com/kdis-lab/LEAC/blob/master/leac-userManual.pdf) can be found in the main directory of this repository and includes:
* Detailed steps for getting and running the library.
* A description of LEAC library architecture.
* Examples for developing a new Evolutionary Algorithm in the library.
* Examples for running an Evolutionary Algorithm included in the library
* Examples for carrying out an experimental study with Evolutionary algorithms included in the library.


[LEAC API](https://github.com/kdis-lab/LEAC/tree/master/docs) can be found in the [*docs*](https://github.com/kdis-lab/LEAC/tree/master/docs) folder.


## Methods included
LEAC includes several implementations of evolutionary algorithms for partitional clustering which are based on state-of-art of evolutionary
proposals in this area. Following the taxonomy given by Hruschka et al., a first classification of these algorithms is carried out
according to the use of fixed or variable number initial of clusters. Concretely, the algorithms included are:

* Proposals that need to set the number of clusters: [GA](http://dx.doi.org/10.1016/0167-8655(96)00043-8), [GKA](http://dx.doi.org/10.1109/3477.764879), [IGKA](http://dx.doi.org/10.1186/1471-2105-5-172), [FGKA](http://doi.acm.org/10.1145/967900.968029), [GA_B](http://dx.doi.org/10.1007/978-3-540-39398-6_7), [GCA](https://doi.org/10.1016/0003-2670(93)80130-D), [HKA](http://dx.doi.org/10.1109/CEC.2004.1330840), [CBGA](http://dx.doi.org/10.1093/comjnl/40.9.547), [GAGR](http://dx.doi.org/10.1016/j.patcog.2008.11.006), [GAs](http://dx.doi.org/10.1016/S0031-3203(99)00137-5) and [KGA](http://dx.doi.org/10.1016/S0020-0255(02)00208-6).
* Proposals that discover the number of clusters: [CGA](http://dl.acm.org/citation.cfm?id=1293920.1293922), [GGA_S](http://dx.doi.org/10.1016/j.eswa.2012.02.149), [GGA_DB](http://dx.doi.org/10.1016/j.eswa.2012.02.149), [EAC](http://dx.doi.org/10.1016/j.ins.2005.07.015), [FEAC_SS](http://dx.doi.org/10.1109/CEC.2006.1688522), [FEAC_RI](http://dx.doi.org/10.1109/CEC.2006.1688522), [GCUK](http://dx.doi.org/10.1016/S0031-3203(01)00108-X), [VGA](http://dx.doi.org/10.1109/5326.923275), [TGCA](http://dx.doi.org/10.1016/j.neucom.2011.11.001), [GA_C](http://dx.doi.org/10.1109/ICEC.1994.350046) and [Clustering](http://dx.doi.org/10.1016/S0031-3203(00)00005-4). 


## LEAC library data format

The format of arff data based on the [Weka's format](https://www.cs.waikato.ac.nz/ml/weka/arff.html)


## Citation

Robles-Berumen, H., Zafra, A., Fardoun, H. M., & Ventura, S. (2019). LEAC: An efficient library for clustering with evolutionary algorithms. *Knowledge-Based Systems*, *179*, 117-119. https://doi.org/10.1016/j.knosys.2019.05.008

## License

The tools is free and open source, under the GNU General Public [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) 
[\[License file\]](../../LICENSE), for Windows&reg;, GNU/Linux&reg; and Mac OS X&reg;.

## Reporting bugs
Feel free to open an [issue](https://github.com/kdis-lab/LEAC/issues) at Github if anything is not working as expected. Merge request are also encouraged, it will be carefully reviewed and merged if everything is all right.


