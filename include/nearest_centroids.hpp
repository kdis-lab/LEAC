/*! \file nearest_centroids.hpp
 *
 * \brief  Function for search nearest centroids
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */

#ifndef  __NEAREST_CENTROIDS_HPP
#define  __NEAREST_CENTROIDS_HPP


/*! \namespace nearest
  \brief Function to find nearest instances
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/
namespace nearest {

/*! \class NearestCentroids
  \brief Stores the relationship between centroids, the closest and its distance
*/  
template < class T_CLUSTERIDX,
           class T_DIST
	   >
class NearestCentroids
{
public:
  NearestCentroids() : 
    _cidx_nearestClusterK(-1), 
    _rt_distCentroidCentroid(0), 
    _b_distRecalculate(true) 
  { }
  NearestCentroids
  (const T_CLUSTERIDX  aimidx_nearestClusterK,
   T_DIST              airt_distCentroidCentroid,
   bool                aib_distRecalculate
   ):
    _cidx_nearestClusterK(aimidx_nearestClusterK),
    _rt_distCentroidCentroid(airt_distCentroidCentroid),
    _b_distRecalculate(aib_distRecalculate)
  { }

  //copy constructor
  NearestCentroids
  (const NearestCentroids<T_CLUSTERIDX,T_DIST> &ainearcent_b):
    _cidx_nearestClusterK(ainearcent_b._cidx_nearestClusterK), 
    _rt_distCentroidCentroid(ainearcent_b._rt_distCentroidCentroid), 
    _b_distRecalculate(ainearcent_b._b_distRecalculate)
  { }

  //move constructor
  NearestCentroids
  (NearestCentroids<T_CLUSTERIDX,T_DIST> &&ainearcent_b):
    _cidx_nearestClusterK(ainearcent_b._cidx_nearestClusterK), 
    _rt_distCentroidCentroid(ainearcent_b._rt_distCentroidCentroid), 
    _b_distRecalculate(ainearcent_b._b_distRecalculate)
  {
    ainearcent_b._cidx_nearestClusterK = -1;
    ainearcent_b._rt_distCentroidCentroid = 0;
    ainearcent_b._b_distRecalculate = true;
  }
  
  ~NearestCentroids() { }

  NearestCentroids<T_CLUSTERIDX,T_DIST>&
  operator=(const NearestCentroids<T_CLUSTERIDX,T_DIST> &ainearcent_b)
  {
    if( this != &ainearcent_b ) { 
      _cidx_nearestClusterK = ainearcent_b._cidx_nearestClusterK; 
      _rt_distCentroidCentroid = ainearcent_b._rt_distCentroidCentroid; 
      _b_distRecalculate = ainearcent_b._b_distRecalculate;
    }
    
    return *this;
  }

  NearestCentroids<T_CLUSTERIDX,T_DIST>&
  operator=(NearestCentroids<T_CLUSTERIDX,T_DIST> &&ainearcent_b)
  {
    if( this != &ainearcent_b ) { 
      _cidx_nearestClusterK = ainearcent_b._cidx_nearestClusterK; 
      _rt_distCentroidCentroid = ainearcent_b._rt_distCentroidCentroid; 
      _b_distRecalculate = ainearcent_b._b_distRecalculate;

      ainearcent_b._cidx_nearestClusterK = -1;
      ainearcent_b._rt_distCentroidCentroid = 0;
      ainearcent_b._b_distRecalculate = true;
    }

    return *this;
  }
  
  inline 
  void 
  setDistCentroidCentroid(const T_DIST airt_distCentroidCentroid)
  {
    this->_rt_distCentroidCentroid = airt_distCentroidCentroid;
  }

  inline 
  void 
  setNearestClusterK(const T_CLUSTERIDX aimidx_midx_nearestClusterK)
  {
    this->_cidx_nearestClusterK = aimidx_midx_nearestClusterK;
  }

  inline 
  void 
  setDistRecalculate(const bool aib_distRecalculate)
  {
    this->_b_distRecalculate = aib_distRecalculate;
  }

  inline 
  T_DIST
  getDistCentroidCentroid() const
  {
    return this->_rt_distCentroidCentroid;
  }

  inline 
  T_CLUSTERIDX
  getNearestClusterK() const 
  {
    return this->_cidx_nearestClusterK;
  }

  inline 
  bool
  getDistRecalculate()
  {
    return this->_b_distRecalculate;
  }

  bool menor(const NearestCentroids<T_DIST, T_CLUSTERIDX> &ainearcent_b) 
  {
    return  ainearcent_b._rt_distCentroidCentroid < this->_rt_distCentroidCentroid;
  }
  
  void  print(std::ostream &os=std::cout, char delim=',')
  {
    os << '('
       << _cidx_nearestClusterK
       << delim
       << _rt_distCentroidCentroid
       << delim << _b_distRecalculate << ')';
  }

protected:
  T_CLUSTERIDX    _cidx_nearestClusterK;
  T_DIST          _rt_distCentroidCentroid;
  bool            _b_distRecalculate; 
}; /*NearestCentroids*/

} /*END namespace nearest*/

#endif  /*__NEAREST_CENTROIDS_HPP*/
