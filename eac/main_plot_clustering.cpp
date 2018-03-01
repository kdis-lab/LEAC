/*! \file main_plot_clustering.cpp
 *
 * \brief Main program for plot clustering
 *
 * \details This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */
#include <iostream>
#include <istream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include "matrix.hpp"
#include "matrix_read.hpp"
#include "inparamclustering_getparameter.hpp"
#include "instances_read.hpp"
#include "gnuplot_i.h"
#include "transformation.hpp"
#include "datatype_instance_real.hpp"
#include "stats_instances.hpp"

#include "verbose.hpp"

std::pair<bool,std::string>
plotclustering_getDataFile(std::ifstream& lfp_file)
{
  std::string lostr_lineData;
  bool        lob_noteof = true;
  std::string lstr_lineWhithFormat;
  if ( std::getline(lfp_file, lstr_lineWhithFormat) ) {
    std::size_t lst_found =  lstr_lineWhithFormat.find_last_of('>');
    if (lst_found != std::string::npos) {
      lostr_lineData = lstr_lineWhithFormat.substr(lst_found+1);
    }
  }
  else {
    lob_noteof = false;
  }
  return std::make_pair(lob_noteof,lostr_lineData);
}


/*---< main() >-------------------------------------------------------------
 */
int main(int argc, char **argv)
{

#ifdef __VERBOSE_YES
  const char* lpc_labeMain = "main_plot_clustering";
  geverbosepc_labelstep = lpc_labeMain;
#endif /*__VERBOSE_YES*/

  inout::InParamPlotClustering
    <DATATYPE_FEATURE,
     DATATYPE_INSTANCES_CLUSTER_K,
     DATATYPE_CLUSTERIDX
     >
    lipplotclustering_inParam
    ("GPC","GNU Plot Clustering",inout::OTHER);

  /*READ PARAMETER
   */
  inparamclustering_getParameter
    (lipplotclustering_inParam, argc, argv);

#ifdef __VERBOSE_YES

  ++geiinparam_verbose;
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "main:  IN"
	      << '(' << geiinparam_verbose << ')'
	      << "\nargv: ";
    std::cout << argv[0] << ' ';
    for (int li_i = 1; li_i < argc;  li_i++ ) {
      std::cout << ' '  << argv[li_i];
    }
    std::cout << std::endl;
  }
#endif /*__VERBOSE_YES*/

  using namespace std;
  istringstream   liss_stringstream;

  int li_numPtCluster      = 16;
  int liarray_ptCluster[]  = {1,2,3,7,9,11,13,15,17,18,20,24,28,49,50,52};
  int li_idxPtClusterCurrent     = 0;
  int li_idxColorClusterCurrent  = 0;
  int li_numPtClass  = 16;
  int liarray_ptClass[] = {4,6,8,10,12,14,33,34,36,40,64,65,66,67,68,69};
  int li_idxPtClassCurrent       = 0;
  int li_idxColorClassCurrent    = 1;
  const int li_numColor          = 6;
  const char *lpt_pointColor[] = {"red", "pink", "green", "purple", "orange", "blue"}; 
  const char*  lptc_pointColorCentroid = "black";
  const char*  lptc_pointColorCentroid2 = "black";
  const int    li_pointLabelSize = 12;

  if ( lipplotclustering_inParam.getOutFileGraphics() == NULL ) {
    li_numPtCluster      = 5;
    li_numPtClass        = 5;
    liarray_ptCluster[0]  = 1;
    liarray_ptCluster[1]  = 2;
    liarray_ptCluster[2]  = 3;
    liarray_ptCluster[3]  = 5;
    liarray_ptCluster[4]  = 7;
    liarray_ptClass[0] = 4;
    liarray_ptClass[1] = 6;
    liarray_ptClass[2] = 8;
    liarray_ptClass[3] = 10;
    liarray_ptClass[4] = 12;
    if ( lipplotclustering_inParam.getPointTypeCentroids() == 65 )
      lipplotclustering_inParam.setPointTypeCentroids(6);
  }

  
  /*BEGIN READ: INSTANCES
   */
  std::vector<data::Instance<DATATYPE_FEATURE>* >  lvectorptinst_instances;
  std::vector<data::Instance<DATATYPE_FEATURE>* >  lvectorptinst_instancesTest;
  std::vector<std::pair<uintidx,DATATYPE_CLUSTERIDX> > lvectorpair_idxInstClass;
  data::Instance<DATATYPE_FEATURE>::setHomogeneousCoord(true);

  if ( lipplotclustering_inParam.getClassInstanceColumn() ) {

    data::InstanceIterfazClass
      <DATATYPE_INSTANCES_CLUSTER_K,
       DATATYPE_CLUSTERIDX
       >
      ::initialize();

    lvectorptinst_instances =
      inout::instancesReadWithClass
      <DATATYPE_FEATURE,
       DATATYPE_INSTANCES_CLUSTER_K,
       DATATYPE_CLUSTERIDX>
      (lipplotclustering_inParam,
       false
       );

    if ( lipplotclustering_inParam.getNumFilesInstanceTest() > 0) {

      lvectorptinst_instancesTest =
	inout::instancesReadWithClass
	<DATATYPE_FEATURE,
	 DATATYPE_INSTANCES_CLUSTER_K,
	 DATATYPE_CLUSTERIDX>
	(lipplotclustering_inParam,
	 true
	 );

    }

    lvectorpair_idxInstClass.reserve(lvectorptinst_instances.size());

    uintidx luintidx_idxInst = 0;
    std::for_each
      (lvectorptinst_instances.begin(),
       lvectorptinst_instances.end(),
       [&](const data::Instance<DATATYPE_FEATURE> *linst_iter)
       {

	 data::InstanceClass<DATATYPE_FEATURE,DATATYPE_INSTANCES_CLUSTER_K,DATATYPE_CLUSTERIDX>
	   *linstclass_iter =
	   (data::InstanceClass
	    <DATATYPE_FEATURE,DATATYPE_INSTANCES_CLUSTER_K,DATATYPE_CLUSTERIDX>*)
	   linst_iter;

	 lvectorpair_idxInstClass.emplace_back
	   (luintidx_idxInst++,
	    linstclass_iter->getClassIdx()
	    );
       }
       );

    std::sort
      (lvectorpair_idxInstClass.begin(), lvectorpair_idxInstClass.end(),
       [](const std::pair<uintidx,DATATYPE_CLUSTERIDX> &left,
	  const std::pair<uintidx,DATATYPE_CLUSTERIDX> &right)
       {
	 return left.second < right.second; /*descend*/
       }
       );
  }
  else { /*Instances class*/
    lvectorptinst_instances =
      inout::instancesRead<DATATYPE_FEATURE>
      (lipplotclustering_inParam,
       false
       );

    if ( lipplotclustering_inParam.getNumFilesInstanceTest() > 0) {
      lvectorptinst_instancesTest =
	inout::instancesRead<DATATYPE_FEATURE>
	(lipplotclustering_inParam,
	 true
	 );
    }
  } /*END READ INSTANCES
     */

  mat::MatrixRow<DATATYPE_FEATURE> lmatrixt_trans;

  if ( lipplotclustering_inParam.getOpProjection() == INPARAMCLUSTERING_PLOTPLOJECTION_PCA ) {

    lmatrixt_trans =
      mat::getMatrixTransformPCA
      (lvectorptinst_instances.begin(),
       lvectorptinst_instances.end()
       );

  }
  else { //INPARAMCLUSTERING_PLOTPLOJECTION_IDENTITY

    lmatrixt_trans =
      mat::getIdentity<DATATYPE_FEATURE>
      (data::Instance<DATATYPE_FEATURE>::getNumDimensions());

  }

  if ( data::Instance<DATATYPE_FEATURE>::getNumDimensions() != lmatrixt_trans.getNumColumns() )
    throw std::invalid_argument
      ("main_plot_clustering.cpp: number of dimensions is equal to zero");

  if ( data::Instance<DATATYPE_FEATURE>::getNumDimensions() == 0 )
    throw std::invalid_argument
      ("main_plot_clustering.cpp: dimensiones de la instancias es diferente a la matriz de transformacion");
  // }

  long int  lli_idxCoordX =
    (lipplotclustering_inParam.getCoordX() != 0)?
    (long int) lipplotclustering_inParam.getCoordX() -1:-1;

  long int  lli_idxCoordY =
    (lipplotclustering_inParam.getCoordY() != 0)?
    (long int) lipplotclustering_inParam.getCoordY() -1:-1;

  long int  lli_idxCoordZ  =
    (lipplotclustering_inParam.getCoordZ() != 0)?
    (long int) lipplotclustering_inParam.getCoordZ() -1:-1;

  gnuplot_ctrl  *h1;
  const char    *lpc_tmpFileName;
  FILE          *lfile_tmpCoord;

  std::ifstream lfp_fileCentroids(lipplotclustering_inParam.getInFileCentroids());
  std::ifstream lfp_fileCentroids2(lipplotclustering_inParam.getInFileCentroids2());
  std::ifstream lfp_fileMemberCluster(lipplotclustering_inParam.getInFileMemberCluster());
  std::ifstream lfp_fileGraph(lipplotclustering_inParam.getInFileGraph());

  /*OPEN CENTROIDS
   */
  if ( (lipplotclustering_inParam.getInFileCentroids() != NULL) && !lfp_fileCentroids ) {
    std::string lstr_error("main_plot_clustering.cpp: unable to open file centroids: ");
    lstr_error += lipplotclustering_inParam.getInFileCentroids();
    throw  std::invalid_argument(lstr_error);
  }

  if ( (lipplotclustering_inParam.getInFileCentroids2() != NULL) && !lfp_fileCentroids2 ) {
    std::string lstr_error("main_plot_clustering.cpp: unable to open file centroids2: ");
    lstr_error += lipplotclustering_inParam.getInFileCentroids2();
    throw  std::invalid_argument(lstr_error);
  }

  /*OPEN MEMBER CLUSTER
   */
  if ( (lipplotclustering_inParam.getInFileMemberCluster() != NULL) && !lfp_fileMemberCluster ) {
    std::string lstr_error("main_plot_clustering.cpp: unable to open file member cluster: ");
    lstr_error += lipplotclustering_inParam.getInFileMemberCluster();
    throw  std::invalid_argument(lstr_error);
  }

  if ( (lipplotclustering_inParam.getInFileGraph() != NULL) && !lfp_fileGraph ) {
    std::string lstr_error("main_plot_clustering.cpp: unable to open file graph");
    lstr_error += lipplotclustering_inParam.getInFileGraph();
    throw  std::invalid_argument(lstr_error);
  }

  bool  lb_plotMemberCluster = false;
  bool  lb_plotMemberClusterFile =
    lipplotclustering_inParam.getInFileMemberCluster() != NULL;
  bool  lb_plotCentroidsFile =
    lipplotclustering_inParam.getInFileCentroids() != NULL;
  bool  lb_plotCentroidsFile2 =
    lipplotclustering_inParam.getInFileCentroids2() != NULL;
  bool  lb_plotGraphFile =
    lipplotclustering_inParam.getInFileGraph() != NULL;

  bool  lb_plotInstances  = true;
  int   li_numPlots = 0;

  while ( lb_plotMemberClusterFile || lb_plotCentroidsFile || lb_plotCentroidsFile2 || lb_plotGraphFile || lb_plotInstances  ) {

    ++li_numPlots;
    h1 = gnuplot_init();
    std::stringstream lss_cmdPlot//("",ios_base::app | ios_base::out);
      (( -1 < lli_idxCoordZ) ? "splot" : "plot",ios_base::app | ios_base::out);

    li_idxPtClusterCurrent = 0;
    li_idxPtClassCurrent   = 0;

    /*MEMBER CLUSTER
     */
    std::string lstr_memberCluster;
    if ( lb_plotMemberClusterFile ) {
      std::tie(lb_plotMemberClusterFile,lstr_memberCluster )
	=  plotclustering_getDataFile(lfp_fileMemberCluster);
    }
    /*GRAPH COMPONENT
     */
    std::string lstr_lineGraph;
    if ( lb_plotGraphFile ) {
      std::tie(lb_plotGraphFile,lstr_lineGraph )
	=  plotclustering_getDataFile(lfp_fileGraph);
    }
    /*WITH CENTROIDS
     */
    std::string lstr_centroids;
    if ( lb_plotCentroidsFile  ) {
      std::tie(lb_plotCentroidsFile,lstr_centroids )
	=  plotclustering_getDataFile(lfp_fileCentroids);
    }
    std::string lstr_centroids2;
    if ( lb_plotCentroidsFile2  ) {
      std::tie(lb_plotCentroidsFile2,lstr_centroids2 )
	=  plotclustering_getDataFile(lfp_fileCentroids2);
    }

    if ( !(lb_plotMemberClusterFile || lb_plotCentroidsFile || lb_plotCentroidsFile2 || lb_plotGraphFile || lb_plotInstances) )
      break;

    /*MEMBER CLUSTER
     */
    if ( lb_plotMemberClusterFile ) { /* IF plotMemberCluster */

      inout::LineSplit lls_memberCluster(",","");
      uintidx luintidx_numInstMemberCluster =
	lls_memberCluster.split(lstr_memberCluster);

      if ( luintidx_numInstMemberCluster > 0 ) {
	std::vector<std::string>& lvector_memberCluster =
	  lls_memberCluster.getVectorRawString();

	std::vector<std::pair<uintidx,DATATYPE_CLUSTERIDX> > lvectorpair_idxInstCluster;
	lvectorpair_idxInstCluster.reserve(luintidx_numInstMemberCluster);

	uintidx luintidx_idxInst = 0;

	std::for_each
	  (lvector_memberCluster.begin(),
	   lvector_memberCluster.end(),
	   [&](std::string& lstr_iterMemberCluster)
	   {
	     DATATYPE_CLUSTERIDX lmcidxt_clusterK;

	     liss_stringstream.clear();
	     liss_stringstream.str(lstr_iterMemberCluster);
	     liss_stringstream >> lmcidxt_clusterK;

	     lvectorpair_idxInstCluster.emplace_back
	       (luintidx_idxInst++,
		lmcidxt_clusterK
		);
	   }
	   );

	std::sort
	  (lvectorpair_idxInstCluster.begin(), lvectorpair_idxInstCluster.end(),
	   [](const std::pair<uintidx,DATATYPE_CLUSTERIDX> &left,
	      const std::pair<uintidx,DATATYPE_CLUSTERIDX> &right)
	   {
	     return left.second < right.second; /*descend*/
	   }
	   );

	DATATYPE_CLUSTERIDX lmcidxt_clusterKLast = lvectorpair_idxInstCluster.at(0).second;
	lpc_tmpFileName = gnuplot_tmpfile(h1);
	lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");

	if (lfile_tmpCoord == NULL) {
	  fprintf
	    (stderr,
	     "main_plot_clustering.cpp:: cannot create temporary file for plot: %s",
	     lpc_tmpFileName
	     );
	  return  -1;
	}

	if ( lipplotclustering_inParam.getIDInstanceColumn() ) {

	  mat::point
	    (lmatrixt_trans,
	     lvectorptinst_instances.at(lvectorpair_idxInstCluster.at(0).first)->getFeatures(),
	     lvectorptinst_instances.at(lvectorpair_idxInstCluster.at(0).first)->getId(),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );
	}
	else {
	  mat::point
	    (lmatrixt_trans,
	     lvectorptinst_instances.at(lvectorpair_idxInstCluster.at(0).first)->getFeatures(),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );
	}

	uintidx luintidx_i = 1;

	while ( luintidx_i < lvectorpair_idxInstCluster.size() ) {
	  if ( lvectorpair_idxInstCluster.at(luintidx_i).second != lmcidxt_clusterKLast )  {
	    fclose(lfile_tmpCoord) ;

	    if (h1->nplots > 0)
	      lss_cmdPlot << ",";

	    lss_cmdPlot << ' ' << '\"' << lpc_tmpFileName << '\"';
	    if ( lipplotclustering_inParam.getIDInstanceColumn()  ) {

	      if (-1 < lli_idxCoordZ)
		lss_cmdPlot << " using 2:3:4 ";
	      else
		lss_cmdPlot << " using 2:3 ";
	    }
	    else {
	      if (-1 < lli_idxCoordZ)
		lss_cmdPlot << " using 1:2:3 with points ";
	      else
		lss_cmdPlot << " using 1:2 with points ";
	    }
	    lss_cmdPlot
	      << " pt " << liarray_ptCluster[li_idxPtClusterCurrent];
	    if ( lipplotclustering_inParam.getOutFileGraphics() == NULL ) {
	      lss_cmdPlot
		<< " lc " << "'" << lpt_pointColor[li_idxColorClusterCurrent] << "'";
	    }
	    lss_cmdPlot << "  ps " << lipplotclustering_inParam.getPointSizeMemberCluster();
	    lss_cmdPlot
	      << " title \'" << lipplotclustering_inParam.getLabelMemberCluster()
	      << (lmcidxt_clusterKLast + lipplotclustering_inParam.getIDincrement())
	      << '\'';

	    if ( lipplotclustering_inParam.getIDInstanceColumn() &&
		 (lipplotclustering_inParam.getClassInstanceColumn() == 0 ) ) {
	      lss_cmdPlot
		<< ", " << '\"' << lpc_tmpFileName << '\"'
		<< " using 2:3:1 "
		<< " with labels offset char "
		<< lipplotclustering_inParam.getLabelOffset()
		<< " notitle";
	    }

	    h1->nplots++;
	    li_idxPtClusterCurrent = (li_idxPtClusterCurrent+1) % li_numPtCluster;
	    li_idxColorClusterCurrent = (li_idxColorClusterCurrent+1) % li_numColor;
	    lmcidxt_clusterKLast = lvectorpair_idxInstCluster.at(luintidx_i).second;
	    lpc_tmpFileName = gnuplot_tmpfile(h1);
	    lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");

	    if (lfile_tmpCoord == NULL) {
	      fprintf
		(stderr,
		 "main_plot_clustering.cpp:: cannot create temporary file for plot: %s",
		 lpc_tmpFileName
		 );
	      return  -1;
	    }
	  }


	  if ( lipplotclustering_inParam.getIDInstanceColumn() ) {

	    mat::point
	      (lmatrixt_trans,
	       lvectorptinst_instances.at
	       (lvectorpair_idxInstCluster.at(luintidx_i).first)->getFeatures(),
	       lvectorptinst_instances.at
	       (lvectorpair_idxInstCluster.at(luintidx_i).first)->getId(),
	       lli_idxCoordX,
	       lli_idxCoordY,
	       lli_idxCoordZ,
	       lfile_tmpCoord
	       );

	  }
	  else {

	    mat::point
	      (lmatrixt_trans,
	       lvectorptinst_instances.at
	       (lvectorpair_idxInstCluster.at(luintidx_i).first)->getFeatures(),
	       lli_idxCoordX,
	       lli_idxCoordY,
	       lli_idxCoordZ,
	       lfile_tmpCoord
	       );
	  }
	  ++luintidx_i;
	}
	fclose(lfile_tmpCoord);

	if (h1->nplots > 0)
	  lss_cmdPlot << ",";
	lss_cmdPlot << ' ' << '\"' << lpc_tmpFileName << '\"';
	if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	  if (-1 < lli_idxCoordZ)
	    lss_cmdPlot << " using 2:3:4 ";
	  else
	    lss_cmdPlot << " using 2:3 ";
	}
	else {
	  if (-1 < lli_idxCoordZ)
	    lss_cmdPlot << " using 1:2:3 with points ";
	  else
	    lss_cmdPlot << " using 1:2 with points ";
	}
	if ( lipplotclustering_inParam.getOutFileGraphics() == NULL ) {
	      lss_cmdPlot
		<< " lc " << "'" << lpt_pointColor[li_idxColorClusterCurrent] << "'";
        }
	lss_cmdPlot
	  << " pt " << liarray_ptCluster[li_idxPtClusterCurrent]
	  << " ps " << lipplotclustering_inParam.getPointSizeMemberCluster();
	lss_cmdPlot
	  << " title \'" << lipplotclustering_inParam.getLabelMemberCluster()
	  << (lmcidxt_clusterKLast + lipplotclustering_inParam.getIDincrement() )
	  << '\'';

	if ( lipplotclustering_inParam.getIDInstanceColumn() &&  (lipplotclustering_inParam.getClassInstanceColumn() == 0 ) ) {
	  lss_cmdPlot
	    << ", " << '\"' << lpc_tmpFileName << '\"'
	    << " using 2:3:1 "
	    << " with labels offset char "
	    << lipplotclustering_inParam.getLabelOffset()
	    << " notitle";
	}

	h1->nplots++;
	li_idxPtClusterCurrent = (li_idxPtClusterCurrent +1) % li_numPtCluster;
	lb_plotMemberCluster = true;
      }
      else {
	lb_plotMemberCluster = false;
      }
    } //END MEMBER CLUSTER

    /*WITH CLASS-----------------------------------------------------------------
     */
    if ( (lipplotclustering_inParam.getClassInstanceColumn() > 0 ) ) {

      const std::vector<data::CountLabel<DATATYPE_CLUSTERIDX,DATATYPE_INSTANCES_CLUSTER_K>* >&
	lstcvector_classLabel = data::InstanceIterfazClass
	<DATATYPE_INSTANCES_CLUSTER_K,DATATYPE_CLUSTERIDX>
	::getVectorClassLabel();
      DATATYPE_CLUSTERIDX lmcidx_clusterBegin(0);
      DATATYPE_CLUSTERIDX lmcidx_clusterEnd(lstcvector_classLabel.size()-1);
      uintidx luintidx_i = 0;

#ifdef __VERBOSE_YES
      ++geiinparam_verbose;
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "main_plot_clustering.cpp:WITH CLASS  IN"
		  << '(' << geiinparam_verbose << ")\n"
		  << "lipplotclustering_inParam.getClassInstanceColumn() = "
		  << lipplotclustering_inParam.getClassInstanceColumn()
		  << "\nlmcidx_clusterBegin = " << lmcidx_clusterBegin
		  << "\nlmcidx_clusterEnd = " << lmcidx_clusterEnd
		  << "\n)"
		  << std::endl;
      }
#endif //__VERBOSE_YES

      for (DATATYPE_CLUSTERIDX lmcidx_k = lmcidx_clusterBegin;
	   lmcidx_k <= lmcidx_clusterEnd;
	   ++lmcidx_k)
	{

	  lpc_tmpFileName = gnuplot_tmpfile(h1);
	  lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");

	  if (lfile_tmpCoord == NULL) {
	    fprintf
	      (stderr,
	       "main_plot_clustering.cpp:: cannot create temporary file for plot: %s",
	       lpc_tmpFileName
	       );
	    return  -1;
	  }

	  DATATYPE_INSTANCES_CLUSTER_K   linstk_numeroInst =
	    lstcvector_classLabel.at(lmcidx_k)->getNumLabel();
	  DATATYPE_INSTANCES_CLUSTER_K   linstk_countInst = 0;

	  while ( linstk_countInst < linstk_numeroInst) {

	    data::InstanceClass<DATATYPE_FEATURE,DATATYPE_INSTANCES_CLUSTER_K,DATATYPE_CLUSTERIDX>
	      *linstclass_iter =
	      (data::InstanceClass
	       <DATATYPE_FEATURE,DATATYPE_INSTANCES_CLUSTER_K,DATATYPE_CLUSTERIDX>*)
	      lvectorptinst_instances.at(lvectorpair_idxInstClass.at(luintidx_i++).first);

	    if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	      mat::point
		(lmatrixt_trans,
		 linstclass_iter->getFeatures(),
		 linstclass_iter->getId(),
		 lli_idxCoordX,
		 lli_idxCoordY,
		 lli_idxCoordZ,
		 lfile_tmpCoord
		 );
	    }
	    else {
	      mat::point
		(lmatrixt_trans,
		 linstclass_iter->getFeatures(),
		 lli_idxCoordX,
		 lli_idxCoordY,
		 lli_idxCoordZ,
		 lfile_tmpCoord
		 );
	    }

	    ++linstk_countInst;

	  }

	  fclose(lfile_tmpCoord) ;

	  if (h1->nplots > 0)
	    lss_cmdPlot << ",";

	  lss_cmdPlot
	    << ' ' << '\"' << lpc_tmpFileName << '\"';
	  if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	    if (-1 < lli_idxCoordZ)
	      lss_cmdPlot << " using 2:3:4 ";
	    else
	      lss_cmdPlot << " using 2:3 ";
	  }
	  else {
	    if (-1 < lli_idxCoordZ)
	      lss_cmdPlot << " using 1:2:3 with points ";
	    else
	      lss_cmdPlot << " using 1:2 with points ";
	  }
	  lss_cmdPlot << " pt " <<  liarray_ptClass[li_idxPtClassCurrent];
	  if ( lipplotclustering_inParam.getOutFileGraphics() == NULL ) {
	      lss_cmdPlot
		<< " lc " << "'" << lpt_pointColor[li_idxColorClassCurrent] << "'";
	  }
	  lss_cmdPlot << " ps " <<  lipplotclustering_inParam.getSizeInstance();
	  lss_cmdPlot << " title "
		      << '\"' << (lstcvector_classLabel.at(lmcidx_k)->getLabel()).c_str() << '\"';
	  if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	    lss_cmdPlot
	      << ", " << '\"' << lpc_tmpFileName << '\"'
	      << " using 2:3:1 "
	      << " with labels offset char "
	      << lipplotclustering_inParam.getLabelOffset()
	      << " notitle";
	  }

	  h1->nplots++;
	  li_idxPtClassCurrent = (li_idxPtClassCurrent+1) % li_numPtClass;
          li_idxColorClassCurrent = (li_idxColorClassCurrent+1) % li_numColor;
	}

#ifdef __VERBOSE_YES
      if ( geiinparam_verbose <= geiinparam_verboseMax ) {
	std::cout << "main_plot_clustering.cpp:WITH CLASS  IN OUT"
		  << '(' << geiinparam_verbose << ")\n"
		  << lss_cmdPlot.str()
		  << std::endl;
      }
      --geiinparam_verbose;
#endif //__VERBOSE_YES


    }
    /*WITHOUT CLASS---------------------------------------------------------------
     */
    else if ( !lb_plotMemberCluster ) {

      lpc_tmpFileName = gnuplot_tmpfile(h1);
      lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");
      if (lfile_tmpCoord == NULL) {
	fprintf
	  (stderr,
	   "main_plot_clustering.cpp:: cannot create temporary file for plot: %s",
	   lpc_tmpFileName
	   );
	return  -1;
      }

      for (uintidx luintidx_i = 0; luintidx_i < lvectorptinst_instances.size(); ++luintidx_i) {

	data::Instance<DATATYPE_FEATURE>* liter_iInstance =
	  lvectorptinst_instances.at(luintidx_i);

	if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	  mat::point
	    (lmatrixt_trans,
	     liter_iInstance->getFeatures(),
	     liter_iInstance->getId(),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );
	}
	else {
	  mat::point
	    (lmatrixt_trans,
	     liter_iInstance->getFeatures(),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );
	}

      }
      fclose(lfile_tmpCoord) ;

      if (h1->nplots > 0)
	lss_cmdPlot << ",";
      lss_cmdPlot << ' ' << '\"' << lpc_tmpFileName << '\"';

      if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	if (-1 < lli_idxCoordZ)
	  lss_cmdPlot << " using 2:3:4 ";
	else
	  lss_cmdPlot << " using 2:3 ";
      }
      else {
	if (-1 < lli_idxCoordZ)
	  lss_cmdPlot << " using 1:2:3 with points ";
	else
	  lss_cmdPlot << " using 1:2 with points ";
      }
      lss_cmdPlot << " ps " <<  lipplotclustering_inParam.getSizeInstance();
      if ( lipplotclustering_inParam.getIDInstanceColumn() ) {
	lss_cmdPlot
	  << ", " << '\"' << lpc_tmpFileName << '\"'
	  << " using 2:3:1 "
	  << " with labels offset char "
	  << lipplotclustering_inParam.getLabelOffset()
	  << " notitle";
      }


      h1->nplots++;
    }
    lb_plotInstances   = false;

    /*GRAPH COMPONENT
     */
    if ( lb_plotGraphFile ) {
      using namespace std;
      istringstream liss_stringstream;

      inout::LineSplit lls_graphListAdj(";","");
      inout::LineSplit lls_graphVertexAdj(",","");

      lpc_tmpFileName = gnuplot_tmpfile(h1);
      lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");
      if (lfile_tmpCoord == NULL) {
	fprintf
	  (stderr,
	   "main_plot_clustering.cpp: cannot create temporary file for plot graph"
	   );
	return  -1;
      }

      uintidx luintidx_numListAdjGraph =
	lls_graphListAdj.split(lstr_lineGraph);

      for ( uintidx luintidx_i = 1; luintidx_i <= luintidx_numListAdjGraph; ++luintidx_i ) {

	uintidx luintidx_numVertexAdj =
	  lls_graphVertexAdj.split(lls_graphListAdj.getItem(luintidx_i));
	uintidx luintidx_vertexFrom;

	liss_stringstream.clear();
	liss_stringstream.str(lls_graphVertexAdj.getItem(1));
	liss_stringstream >>  luintidx_vertexFrom;

	for ( uintidx luintidx_j = 2; luintidx_j <= luintidx_numVertexAdj; ++luintidx_j ) {
	  uintidx luintidx_vertexTo;
	  liss_stringstream.clear();
	  liss_stringstream.str(lls_graphVertexAdj.getItem(luintidx_j));
	  liss_stringstream >>  luintidx_vertexTo;

	  mat::point
	    (lmatrixt_trans,
	     lvectorptinst_instances.at(luintidx_vertexFrom)->getFeatures(),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );

	  mat::point
	    (lmatrixt_trans,
	     lvectorptinst_instances.at(luintidx_vertexTo)->getFeatures(),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );

	  if (-1 < lli_idxCoordZ)
	    fprintf(lfile_tmpCoord, "\n\n");
	  else
	    fprintf(lfile_tmpCoord, "\n");
	}
      }

      fclose(lfile_tmpCoord);

      if (h1->nplots > 0)
	lss_cmdPlot << ",";
      lss_cmdPlot
	<< ' ' << '\"' << lpc_tmpFileName << '\"'
	<< " notitle  with lines lt -1 lc rgb \"black\" lw 0.5";
      h1->nplots++;
    }

    /*WITH CENTROIDS
     */
    if ( lb_plotCentroidsFile  ) {
      mat::MatrixRow<DATATYPE_FEATURE> lmatrixt_centroids =
	mat::matrixread_get<DATATYPE_FEATURE>(lstr_centroids,true);

      lpc_tmpFileName = gnuplot_tmpfile(h1);
      lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");
      if (lfile_tmpCoord == NULL) {
	fprintf
	  (stderr,
	   "main_plot_clustering.cpp: cannot create temporary file for plot centroids"
	   );
	return  -1;
      }

      if (lmatrixt_centroids.getNumRows() > 0 ) {

	for ( uintidx lst_i = 0; lst_i < lmatrixt_centroids.getNumRows(); ++lst_i) {

	  mat::point
	    (lmatrixt_trans,
	     lmatrixt_centroids.getRow(lst_i),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );
	}

	fclose(lfile_tmpCoord);

	char lpc_using[20];
	if(-1 < lli_idxCoordZ) {
	  sprintf(lpc_using,"1:2:3:($0+%d)",lipplotclustering_inParam.getIDincrement());
	}
	else {
	  sprintf(lpc_using,"1:2:($0+%d)",lipplotclustering_inParam.getIDincrement());
	}

	if (h1->nplots > 0)
	  lss_cmdPlot << ",";
	lss_cmdPlot
	  << ' ' << '\"' << lpc_tmpFileName << '\"'
	  << " title \'" << lipplotclustering_inParam.getTitleCentroids() << "\'"
	  << " with points "
	  << " pt " << lipplotclustering_inParam.getPointTypeCentroids()
	  << " ps " << lipplotclustering_inParam.getPointSizeCentroids()
	  << " lc rgb " << '\"' << lptc_pointColorCentroid << '\"';
	h1->nplots++;

	if ( lipplotclustering_inParam.getWithIDCentroids() ) {
	  if (h1->nplots > 0)
	    lss_cmdPlot << ",";
	  lss_cmdPlot
	    << ' ' << '\"' << lpc_tmpFileName << '\"'
	    << " using "  << lpc_using
	    << " with labels font 'Arial Bold," <<  li_pointLabelSize
	    << " offset char " << lipplotclustering_inParam.getLabelOffset()
	    << "' notitle ";
	  h1->nplots++;
	}

      } /*IS NOT NULL CENTROIDS*/

    } /*END if */

    if ( lb_plotCentroidsFile2  ) {
      mat::MatrixRow<DATATYPE_FEATURE> lmatrixt_centroids2 =
	mat::matrixread_get<DATATYPE_FEATURE>(lstr_centroids2,true);

      lpc_tmpFileName = gnuplot_tmpfile(h1);
      lfile_tmpCoord  = fopen(lpc_tmpFileName, "w");
      if (lfile_tmpCoord == NULL) {
	fprintf
	  (stderr,
	   "main_plot_clustering.cpp: cannot create temporary file for plot centroids"
	   );
	return  -1;
      }

      if (lmatrixt_centroids2.getNumRows() > 0 ) {  /*IS NOT NULL CENTROIDS*/

	for ( uintidx lst_i = 0; lst_i < lmatrixt_centroids2.getNumRows(); ++lst_i) {

	  mat::point
	    (lmatrixt_trans,
	     lmatrixt_centroids2.getRow(lst_i),
	     lli_idxCoordX,
	     lli_idxCoordY,
	     lli_idxCoordZ,
	     lfile_tmpCoord
	     );
	}

	fclose(lfile_tmpCoord) ;

	char lpc_usinglabel[20];
	if(-1 < lli_idxCoordZ) {
	  sprintf(lpc_usinglabel,"1:2:3:($0+%d)",lipplotclustering_inParam.getIDincrement());
	}
	else {
	  sprintf(lpc_usinglabel,"1:2:($0+%d)",lipplotclustering_inParam.getIDincrement());
	}

	if (h1->nplots > 0)
	  lss_cmdPlot << ",";
	lss_cmdPlot
	  << ' ' << '\"' << lpc_tmpFileName << '\"'
	  << " title \'" << lipplotclustering_inParam.getTitleCentroids2() << "\'"
	  << " with points "
	  << " pt " << lipplotclustering_inParam.getPointTypeCentroids2()
	  << " ps " << lipplotclustering_inParam.getPointSizeCentroids2()
	  << " lc rgb " << '\"' << lptc_pointColorCentroid2 << '\"';
	h1->nplots++;

	if ( lipplotclustering_inParam.getWithIDCentroids2() ) {
	  if (h1->nplots > 0)
	    lss_cmdPlot << ",";
	  lss_cmdPlot
	    << ' ' << '\"' << lpc_tmpFileName << '\"'
	    << " using "  << lpc_usinglabel
	    << " with labels font 'Arial Bold," <<  li_pointLabelSize
	    << " offset char " << lipplotclustering_inParam.getLabelOffset()
	    << "' notitle ";
	  h1->nplots++;
	}

      } /*IS NOT NULL CENTROIDS*/

    } /*END if   */


    if ( lipplotclustering_inParam.getTitle() != NULL ) {
      gnuplot_cmd(h1, "set title %s",lipplotclustering_inParam.getTitle());
    }
    /*set key
     */
    if ( lipplotclustering_inParam.getLegendOption() != NULL ) {
      gnuplot_cmd
	(h1,
	 "set key %s",
	 lipplotclustering_inParam.getLegendOption()
	 );
    }

    if ( lipplotclustering_inParam.getOutFileGraphics() != NULL ) {
      if ( lipplotclustering_inParam.getColorGraphics() == false )
	gnuplot_cmd
	  (h1,
	   "set terminal postscript eps enhanced monochrome solid font 'Helvetica,10'"
	   );
      else
	gnuplot_cmd
	  (h1,
	   "set terminal postscript eps enhanced color font 'Helvetica,10'"
	   );
      gnuplot_cmd
	(h1,
	 "set output '%s%d.eps'",
	 lipplotclustering_inParam.getOutFileGraphics(),
	 li_numPlots
	 );

    }
    else {
#if defined(_WIN32) || defined(_WIN64)

#else
      gnuplot_cmd(h1, "set terminal x11"); //for mac an linux
#endif // #ifdef WIN32
    }

    std::string lstr_cmdPlot = lss_cmdPlot.str();
    gnuplot_cmd(h1,"set datafile sep ','");


    if ( lipplotclustering_inParam.getXRange() != NULL ) {
      std::stringstream lss_cmdPlotXRange;
      lss_cmdPlotXRange <<  "set xrange [" << lipplotclustering_inParam.getXRange() << "]";
      std::string lstr_cmdPlotXRange = lss_cmdPlotXRange.str();
      gnuplot_cmd(h1,lstr_cmdPlotXRange.c_str());
    }
    if ( lipplotclustering_inParam.getYRange() !=  NULL ) {
      std::stringstream lss_cmdPlotYRange;
      lss_cmdPlotYRange << "set yrange [" << lipplotclustering_inParam.getYRange() << "]";
      std::string lstr_cmdPlotYRange = lss_cmdPlotYRange.str();
      gnuplot_cmd(h1,lstr_cmdPlotYRange.c_str());

    }

#ifdef __VERBOSE_YES
    ++geiinparam_verbose;
    if ( geiinparam_verbose <= geiinparam_verboseMax ) {
      std::cout
	<< "#-------------------------------------------------------------------------------------\n"
	<< lstr_cmdPlot << std::endl
	<< "#--------------------------------------------------------------------------------------"
	<< std::endl;
    }
    --geiinparam_verbose;
#endif /*__VERBOSE_YES*/
    /*VERBOSE THE COMMAND
     */
    gnuplot_cmd(h1,lstr_cmdPlot.c_str());
    if ( lipplotclustering_inParam.getOutFileGraphics() == NULL ) {
      sleep(lipplotclustering_inParam.getSleep());
    }

    gnuplot_close(h1);

  };

  /*DELETE INSTANCES
   */
  for ( auto  liter_instance: lvectorptinst_instances )
    delete liter_instance;

  if ( lipplotclustering_inParam.getInFileCentroids() != NULL ) {
    lfp_fileCentroids.close();
  }

  if ( lipplotclustering_inParam.getInFileMemberCluster()  != NULL ) {
    lfp_fileMemberCluster.close();
  }

  if ( lipplotclustering_inParam.getInFileGraph()  != NULL ) {
    lfp_fileGraph.close();
  }

#ifdef __VERBOSE_YES
  if ( geiinparam_verbose <= geiinparam_verboseMax ) {
    std::cout << "main: OUT"
	      << '(' << geiinparam_verbose << ")\n";
  }
  --geiinparam_verbose;
#endif /*__VERBOSE_YES*/

  return 0;
}

