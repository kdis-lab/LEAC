/*! \file inparam_plotclustering.hpp
 *
 * \brief Definition of plot clustering program parameters
 *
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef IN_PARAM_PLOT_CLUSTERING_HPP
#define IN_PARAM_PLOT_CLUSTERING_HPP

#include "inparam.hpp"
#include "inparam_readinst.hpp"

#define INPARAMCLUSTERING_PLOTPLOJECTION_PCA       0
#define INPARAMCLUSTERING_PLOTPLOJECTION_IDENTITY  1
  
#define INPARAMCLUSTERING_PLOT_PLOJECTION {"pca", "identity", (char *) NULL }

/*! \namespace inout
  \brief Module for input and output parameters
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace  inout {
  
#define SLEEP_LGTH  1800

/*! \class InParamPlotClustering
  \brief Input parameter for plot clustering 
*/
template < typename T_FEATURE,
	   typename T_INSTANCES_CLUSTER_K,
	   typename T_CLUSTERIDX
	  > 
class InParamPlotClustering 
  : public InParamAlgorithmo
  , public InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>
{
public: 
  InParamPlotClustering
  (std::string        ais_algorithmoName,
   std::string        ais_algorithmoAuthor,
   InParam_algTypeOut aiato_algTypeOut)
    : InParamAlgorithmo(ais_algorithmoName,ais_algorithmoAuthor,aiato_algTypeOut)
    , InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>()
    , _i_opProjection(0)
    , _ps_inFileCentroids(NULL)
    , _str_titleCentroids("Centroids")
    , _b_withIDCentroids(true)
    , _i_pointTypeCentroid(65) /*Point type Cicle=65,Diamond=69, */
    , _d_pointSizeCentroid(2.0)  /*Point size */
    , _ps_inFileCentroids2(NULL)
    , _b_withIDCentroids2(true)
    , _i_pointTypeCentroid2(69) /*Point type Cicle=65,Diamond=69, */
    , _d_pointSizeCentroid2(2.0)  /*Point size */

    , _d_labelOffset(1.8)   
       
    , _ps_inFileMemberCluster(NULL)
    , _str_labelMemberCluster("Cluster-")
    , _d_pointSizeMemberCluster(1.0)   
    , _ps_inFileGraph(NULL)
    , _ps_outFileGraphics(NULL)
    , _ps_title(NULL)
    , _ps_legendOption(NULL)
    , _i_IDincrement(0)
    , _uintidx_coordX(1) 
    , _uintidx_coordY(2) 
    , _uintidx_coordZ(0)
    , _d_sizeInstance(1.0)
    , _b_colorGraphics(true)
    , _ps_xRange(NULL)
    , _ps_yRange(NULL)
    , _uint_sleep(SLEEP_LGTH)
  {}

  ~InParamPlotClustering() {}

  inline void setOpProjection(int ai_opProjection) 
  {
    this->_i_opProjection = ai_opProjection;
  }

  inline int getOpProjection() 
  {
    return this->_i_opProjection;
  }
  
  inline void setInFileCentroids(char* aips_fileName) 
  {
    this->_ps_inFileCentroids = aips_fileName;
  }

  inline char* getInFileCentroids() 
  {
    return this->_ps_inFileCentroids;
  }

  inline void setWithIDCentroids(bool aib_withIDCentroids) 
  {
    this->_b_withIDCentroids = aib_withIDCentroids;
  }

  inline bool getWithIDCentroids() 
  {
    return this->_b_withIDCentroids;
  }

  inline void setIDincrement(int aii_IDincrement) 
  {
    this->_i_IDincrement = aii_IDincrement;
  }

  inline int getIDincrement() 
  {
    return this->_i_IDincrement;
  }

  inline void setPointTypeCentroids(int aii_pointTypeCentroid) 
  {
    this->_i_pointTypeCentroid = aii_pointTypeCentroid;
  }

  inline int getPointTypeCentroids() 
  {
    return this->_i_pointTypeCentroid;
  }

  inline void setPointSizeCentroids(double aid_pointSizeCentroid) 
  {
    this->_d_pointSizeCentroid = aid_pointSizeCentroid;
  }

  inline double getPointSizeCentroids() 
  {
    return this->_d_pointSizeCentroid;
  }

  inline void  setTitleCentroids(char* aips_titleCentroids)
  {
    _str_titleCentroids = aips_titleCentroids;
  }
  

  inline void setPointTypeCentroids2(int aii_pointTypeCentroid) 
  {
    this->_i_pointTypeCentroid2 = aii_pointTypeCentroid;
  }

  inline int getPointTypeCentroids2() 
  {
    return this->_i_pointTypeCentroid2;
  }

  inline void setPointSizeCentroids2(double aid_pointSizeCentroid) 
  {
    this->_d_pointSizeCentroid2 = aid_pointSizeCentroid;
  }

  inline double getPointSizeCentroids2() 
  {
    return this->_d_pointSizeCentroid2;
  }

  inline void setLabelOffset(double aid_labelOffset) 
  {
    this->_d_labelOffset = aid_labelOffset;
  }

  inline double getLabelOffset() 
  {
    return this->_d_labelOffset;
  }

  inline void setWithIDCentroids2(bool aib_withIDCentroids) 
  {
    this->_b_withIDCentroids2 = aib_withIDCentroids;
  }

  inline bool getWithIDCentroids2() 
  {
    return this->_b_withIDCentroids2;
  }

  inline std::string& getTitleCentroids()
  {
    return this->_str_titleCentroids;
  }

  //--
  inline void setInFileCentroids2(char* aips_fileName2) 
  {
    this->_ps_inFileCentroids2 = aips_fileName2;
  }

  inline char* getInFileCentroids2() 
  {
    return this->_ps_inFileCentroids2;
  }

  inline void  setTitleCentroids2(char* aips_titleCentroids2)
  {
    _str_titleCentroids2 = aips_titleCentroids2;
  }
  
  inline std::string& getTitleCentroids2()
  {
    return this->_str_titleCentroids2;
  }
  
  inline void setInFileMemberCluster(char* aips_fileName) 
  {
    this->_ps_inFileMemberCluster = aips_fileName;
  }

  inline char* getInFileMemberCluster() 
  {
    return this->_ps_inFileMemberCluster;
  }

  inline void  setLabelMemberCluster(char* aips_labelMemberCluster)
  {
    _str_labelMemberCluster = aips_labelMemberCluster;
  }
  
  inline std::string& getLabelMemberCluster()
  {
    return this->_str_labelMemberCluster;
  }

  inline void setPointSizeMemberCluster(double aid_pointSizeMemberCluster) 
  {
    this->_d_pointSizeMemberCluster = aid_pointSizeMemberCluster;
  }

  inline double getPointSizeMemberCluster() 
  {
    return this->_d_pointSizeMemberCluster;
  }

  inline void setInFileGraph(char* aips_fileName) 
  {
    this->_ps_inFileGraph = aips_fileName;
  }

  inline char* getInFileGraph() 
  {
    return this->_ps_inFileGraph;
  }

  inline void setOutFileGraphics(char* aips_fileName) 
  {
    this->_ps_outFileGraphics = aips_fileName;
  }

  inline char* getOutFileGraphics() 
  {
    return this->_ps_outFileGraphics;
  }

 inline void  setTitle(char* aips_title)
  {
    _ps_title = aips_title;
  }
  
  inline char*  getTitle()
  {
    return this->_ps_title;
  }

  inline void  setLegendOption(char* aips_legendOption)
  {
    _ps_legendOption = aips_legendOption;
  }
  
  inline char*  getLegendOption()
  {
    return this->_ps_legendOption;
  }

  inline void setSleep(unsigned int aiuint_sleep)
  {
    this->_uint_sleep = aiuint_sleep;
  }
  
  inline unsigned int getSleep()
  {
    return this->_uint_sleep;
  }
  
  inline void setCoordX(uintidx aiuintidx_coordX) 
  {
    this->_uintidx_coordX = aiuintidx_coordX;
  }

  inline uintidx getCoordX() 
  {
    return this->_uintidx_coordX;
  }

  inline void setCoordY(uintidx aiuintidx_coordY) 
  {
    this->_uintidx_coordY = aiuintidx_coordY;
  }

  inline uintidx getCoordY() 
  {
    return this->_uintidx_coordY;
  }

  inline void setCoordZ(uintidx aiuintidx_coordZ) 
  {
    this->_uintidx_coordZ = aiuintidx_coordZ;
  }

  inline uintidx getCoordZ() 
  {
    return this->_uintidx_coordZ;
  }


  inline void setSizeInstance(double aid_sizeInstance) 
  {
    this->_d_sizeInstance = aid_sizeInstance;
  }

  inline double getSizeInstance() 
  {
    return this->_d_sizeInstance;
  }


  inline void setColorGraphics(bool aib_colorGraphics) 
  {
    this->_b_colorGraphics = aib_colorGraphics;
  }

  inline bool getColorGraphics() 
  {
    return this->_b_colorGraphics;
  }

  inline void setXRange(char* aips_xRange) 
  {
    this->_ps_xRange = aips_xRange;
  }

  inline char* getXRange() 
  {
    return this->_ps_xRange;
  }

  inline void setYRange(char* aips_yRange) 
  {
    this->_ps_yRange = aips_yRange;
  }

  inline char* getYRange() 
  {
    return this->_ps_yRange;
  }
    
  virtual void print(std::ostream& aipf_outFile=std::cout, const char aic_separator=',') const
  {   
    InParamAlgorithmo::print(aipf_outFile,aic_separator);
    InParamReadInst<T_FEATURE,T_INSTANCES_CLUSTER_K,T_CLUSTERIDX>::print(aipf_outFile,aic_separator);
  }

protected:

  int         _i_opProjection;
  char        *_ps_inFileCentroids;
  std::string _str_titleCentroids;
  bool        _b_withIDCentroids;
  int         _i_pointTypeCentroid;
  double      _d_pointSizeCentroid;
  char        *_ps_inFileCentroids2;
  std::string _str_titleCentroids2;
  bool        _b_withIDCentroids2;
  int         _i_pointTypeCentroid2;
  double      _d_pointSizeCentroid2;
  
  double      _d_labelOffset;
  
  char        *_ps_inFileMemberCluster;
  std::string _str_labelMemberCluster;
  double      _d_pointSizeMemberCluster;
  
  char        *_ps_inFileGraph;
  char        *_ps_outFileGraphics;

  char        *_ps_title;
  char        *_ps_legendOption;

  int         _i_IDincrement;
 
  uintidx _uintidx_coordX; //0...N
  uintidx _uintidx_coordY; //0...N
  uintidx _uintidx_coordZ; //0...N

  double  _d_sizeInstance;
  bool    _b_colorGraphics;
  
  char   *_ps_xRange;
  char   *_ps_yRange;

  unsigned int _uint_sleep;
  
}; /*InParamPlotClustering*/

} /* END namespace inout*/
  
#endif /*IN_PARAM_PLOT_CLUSTERING_HPP*/
