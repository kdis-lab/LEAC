//gcc saasxpa_kh.c -O3 -funroll-loops -msse3  -S
//Flags: -msse, -msse2, -msse3, -march=, -mfpmath=sse
#include <stdint.h> //int64_t

void
saasxpa_kh
(const float    aif_alpha,
 float          *aimatrixrowfloat_a,
 const int64_t  aiint64_numRows,
 const int64_t  aiint64_numColumns,
 const float    *aiarrayfloat_x
 )
{
  const float *last1 = aiarrayfloat_x + aiint64_numColumns;
  int64_t liu_i;
  
  for ( liu_i = 0; liu_i < aiint64_numRows; liu_i++) {
    const float *larrayt_xIter = aiarrayfloat_x;
    while (larrayt_xIter != last1) {
      *aimatrixrowfloat_a += aif_alpha * ((*aimatrixrowfloat_a) - (*larrayt_xIter));
      ++aimatrixrowfloat_a;
      ++larrayt_xIter;
    }
    
  }
}


