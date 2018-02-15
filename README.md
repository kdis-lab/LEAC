# LEAC - Library Evolutionary Algorithms for Clustering

Library Evolutionary Algorithms for Clustering (**LEAC**) is a library for the implementation
of evolutionary and genetic algorithms to solve the problem of *partition clustering*.
The tools is free and open source, under the GNU General Public
[GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) 
[\[License file\]](../../LICENSE),
for Windows&reg;, GNU/Linux&reg; and OS X&reg;.

Clustering is useful in several exploratory *pattern-analysis*,
*grouping*, *decision-making*, and *machine-learning situations*,
including *data mining*, *document retrieval*, *image segmentation*,
and *pattern classification*
\cite Jain:Murty:Flynn:ClusteringSurvey:1999

![](../master/doc/leac_cluster.svg)

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
    <td> <a href="http://dx.doi.org/10.1016/0167-8655(96)00043-8">GA gaclustering_fklabel</a>
         GKA	gka_fklabel	\cite Krishna:Murty:GAClustering:GKA:1999
	 IGKA 	igka_fklabel 	\cite Lu:etal:GAclusteringLabel:IGKA:2004
	 FGKA 	fgka_fklabebl 	\cite Lu:etal:GAclusteringLabel:FGKA:2004
    </td>
    <td> <a href="http://dx.doi.org/10.1016/j.eswa.2012.02.149">GGA gga_vklabeldbindex gga_vklabelsilhouette</a>,
         CGA [cga_vklabel](http://iospress.metapress.com/content/adhnkma5h48f1l0q/)
	 EAC eac_vklabel	\cite Hruschka:etal:GAClusteringLabelKVar:EAC:2006
	 EAC I			eaci_vklabel	\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
	 EAC II			eacii_vklabel	\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
	 EAC III		eaciii_vklabel	\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
	 FEAC			feac_vklabel	\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006
    </td>
   </tr>
   <tr>
    <td>Crisp-matrix</td>
     <td> GA	gaclustering_fkcrispmatrix	\cite Bezdek:etal:GAclustering:GA:1994
     </td>
     <td></td>
   </tr>
  <tr>
    <td>Centroid</td>
    <td> GAS	gas_fkcentroid	\cite Maulik:Bandyopadhyay:GAclustering:GAS:2000
         KGA	kga_fkcentroid	\cite Bandyopadhyay:Maulik:GAclustering:KGA:2002
	 GAGR	gagr_fkcentroid	\cite Chang:etal:GAclustering:GAGR:2009
	 CBGA	cbga_fkcentroid_int and cbga_fkcentroid	\cite Franti:etal:GAclustering:gafranti:1997
    </td>
    <td> GCUK	gcuk_vkcentroid	\cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002,
         TGCA	tgca_vkcentroid	\cite He:Tan:GAclusteringVarK:TGCA:2012
    </td>
    </tr>
    <tr>
    <td>Medoid</td>
     <td> GA-Prototypes	gaprototypes_fkmedoid	\cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997
          HKA	hka_fkmedoid	\cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004
          GCA	gca_fkmedoid	\cite Lucasius:etal:GAclusteringMedoid:GCA:1993
     </td>
     <td></td>
     </tr>
     <tr>
     <td>Tree</td>
     <td></td>
     <td>GA	gaclustering_vktreebinary	\cite Casillas:etal:GAclusteringVarK:GA:2003
     </tr>
     <tr>
     <td>Sub-cluster</td>
     <td></td>
     <td> GA	gaclustering_vktreebinary	\cite Casillas:etal:GAclusteringVarK:GA:2003</td>
     </tr>
     <tr>
</table>

