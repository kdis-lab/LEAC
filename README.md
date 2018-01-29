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

<img src="../leac_cluster.svg" alt="Clustering">


LEAC is built in C ++ with the current standards of C++11, C++14.
Take advantage of the STL library, through the use of 
data structures and algorithms, for example,
Numerics library of Pseudo-random number generation which is 
fundamental in the genetic and evolutionary algorithms.
An additional feature that LEAC has is a low-level software layer
based on Streaming SIMD Extensions (SSE)\cite progguide:intel10
and in the 
[OpenBLAS](http://www.openblas.net) library.
\href{http://www.openblas.net}{OpenBLAS} 
library,in order to increase performance.

LEAC demonstrates its ease of implementation of algorithms, 
by including a collection of genetic and evolutionary programs (EAC), 
which can be used to perform cluster analysis of a data set.

<table>
  <tr>
    <th>Fixed-K</th>
  </tr>
   <tr>
    <th>Algorithms</th>
    <th>Program</th>
    <th>Based on</th>
  </tr>
  <tr>
    <td>GA</td>
    <td>gaclustering_fklabel</td>
    <td>\cite Murthy:Chowdhury:GAclustering:GA:1996</td>
  </tr>
  <tr>
    <td>GKA</td>
    <td>gka_fklabel</td>
    <td>\cite Krishna:Murty:GAClustering:GKA:1999</td>
  </tr>
  <tr>
    <td>IGKA</td>
    <td>igka_fklabel</td>
    <td>\cite Lu:etal:GAclusteringLabel:IGKA:2004</td>
  </tr>
  <tr>
    <td>FGKA</td>
    <td>fgka_fklabebl</td>
    <td>\cite Lu:etal:GAclusteringLabel:FGKA:2004</td>
  </tr>
  <tr>
    <td>GA</td>
    <td>gaclustering_fkcrispmatrix</td>
    <td>\cite Bezdek:etal:GAclustering:GA:1994</td>
  </tr>
  <tr>
    <td>GAS</td>
    <td>gas_fkcentroid</td>
    <td>\cite Maulik:Bandyopadhyay:GAclustering:GAS:2000</td>
  </tr>
  <tr>
    <td>KGA</td>
    <td>kga_fkcentroid</td>
    <td>\cite Bandyopadhyay:Maulik:GAclustering:KGA:2002</td>
  </tr>
  <tr>
    <td>GAGR</td>
    <td>gagr_fkcentroid</td>
    <td>\cite Chang:etal:GAclustering:GAGR:2009</td>
  </tr>
  <tr>
    <td>CBGA</td>
    <td>cbga_fkcentroid_int and cbga_fkcentroid</td>
    <td>\cite Franti:etal:GAclustering:gafranti:1997</td> 
  </tr>
  <tr>
    <td>GA-Prototypes</td>
    <td>gaprototypes_fkmedoid</td>
    <td>\cite Kuncheva:Bezdek:GAMedoid:GAPrototypes:1997</td>
  </tr>
  <tr>
    <td>HKA</td>
    <td>hka_fkmedoid</td>
    <td>\cite Sheng:Xiaohui:GAclusteringMedoid:HKA:2004</td>
  </tr>
  <tr>
    <td>GCA</td>
    <td>gca_fkmedoid</td>
    <td>\cite Lucasius:etal:GAclusteringMedoid:GCA:1993</td>
  </tr>
  <tr>
    <th>Variable-K</th>
  </tr>
  <tr>
    <th>Algorithms</th>
    <th>Program</th>
    <th>Based on</th>
  </tr>
  <tr>
    <td>GGA</td>
    <td>gga_vklabeldbindex and gga_vklabelsilhouette</td>
    <td>\cite Agustin:etal:GAclusteringVarK:GGA:2012</td>
  </tr>
   <tr>
    <td>CGA</td>
    <td>cga_vklabel</td>
    <td>\cite Hruschka:etal:GAclusteringLabelKVar:CGAII:2004 \cite Hruschka:etal:GAClusteringLabelKVar:EAC:2006</td>
  </tr>
  <tr>
    <td>EAC</td>
    <td>eac_vklabel</td>
    <td>\cite Hruschka:etal:GAClusteringLabelKVar:EAC:2006</td>
  </tr>
  <tr>
    <td>EAC I</td>
    <td>eaci_vklabel</td>
    <td>\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006</td>
  </tr>
  <tr>
    <td>EAC II</td>
    <td>eacii_vklabel</td>
    <td>\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006</td>
  </tr>
  <tr>
    <td>EAC III</td>
    <td>eaciii_vklabel</td>
    <td>\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006</td>
  </tr>
  <tr>
    <td>FEAC</td>
    <td>feac_vklabel</td>
    <td>\cite Alves:etal:GAclusteringLabelKVar:FEAC:2006</td>
  </tr>
  <tr>
    <td>GCUK</td>
    <td>gcuk_vkcentroid</td>
    <td>\cite Bandyopadhyay:Maulik:GACVarK:GCUK:2002</td>
  </tr>
  <tr>
    <td>TGCA</td>
    <td>tgca_vkcentroid</td>
    <td>\cite He:Tan:GAclusteringVarK:TGCA:2012</td>
  </tr>
  <tr>
    <td>CLUSTERING</td>
    <td>clustering_vksubclusterbinary</td>
    <td>\cite Tseng:Yang:GAclusteringVarK:CLUSTERING:2001</td>
  </tr>
  <tr>
    <td>GA</td>
    <td>gaclustering_vktreebinary</td>
    <td>\cite Casillas:etal:GAclusteringVarK:GA:2003</td>
  </tr>
</table>

