/*! \file random_ext.hpp
 *
 * \brief Operation for the generation of pseudo random numbers
 *
 * \details  This file is part of the LEAC.\n\n
 * \version 1.0
 * \date 2015-2017
 * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
 * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
 */


#ifndef __RANDOM_EXT_HPP
#define __RANDOM_EXT_HPP

#include <algorithm>  //std::generate_n
#include <functional> //std::ref
#include <iterator>   //ostream_iterator
#include <random>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>  //std::generate_n
#include <functional> //std::ref
#include <iterator>   //ostream_iterator

//BEGIN WINDOWS
#if defined(_WIN32) || defined(_WIN64)

std::random_device grd_randomdevice;

//END WINDOWS
#elif defined(__unix__) || defined(__APPLE__) && defined(__MACH__)

std::random_device grd_randomdevice;

#endif

/*! A global variable of type mersenne_twister_engine is a random number engine based on Mersenne Twister algorithm. It produces high quality unsigned integer random numbers of type UIntType on the interval \f$[0, 2w-1]\f$.*/
//std::mersenne_twister_engine gmt19937_eng;

#if defined(__LP64__) || defined(_WIN64) || (defined(__x86_64__) && !defined(__ILP32__) ) || defined(_M_X64) || defined(__ia64) || defined (_M_IA64) || defined(__aarch64__) || defined(__powerpc64__)
#define __STD_MT19937_64
typedef std::mt19937_64 StdMT19937;
#else
typedef std::mt19937 StdMT19937;
#endif

StdMT19937 gmt19937_eng;

/*! \namespace randomext
  \brief Functions to generate random numbers
  \details
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace randomext {

std::string
setSeed(const unsigned int aiu_numSeed = 8)
{
  std::vector<unsigned int> lvectorui_seeddata(aiu_numSeed);
  std::generate_n
    (lvectorui_seeddata.data(),
     lvectorui_seeddata.size(),
     std::ref(grd_randomdevice)
     );
  std::seed_seq lseedseq_s(std::begin(lvectorui_seeddata), std::end(lvectorui_seeddata));
  
  gmt19937_eng.seed(lseedseq_s);

  std::ostringstream lostrstream_secuencia;
  lseedseq_s.param(std::ostream_iterator<unsigned int>(lostrstream_secuencia, " "));
  std::string lstr_seed_seq = lostrstream_secuencia.str().c_str();
  if (lstr_seed_seq.size () > 0)  lstr_seed_seq.resize (lstr_seed_seq.size () - 1);
      
  return lstr_seed_seq;
      
}

void
setSeed(std::string aistr_seed_seq)
{
  std::istringstream lstrstream_buffer(aistr_seed_seq);
  std::istream_iterator<unsigned int> beg(lstrstream_buffer), end;
  std::seed_seq lseedseq_s(beg, end);

  gmt19937_eng.seed(lseedseq_s);
}

}

#endif /*__RANDOM_EXT_HPP*/
