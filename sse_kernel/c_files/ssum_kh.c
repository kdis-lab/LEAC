//gcc ssum_kh.c -O3 -funroll-loops -msse3  -S
//Flags: -msse, -msse2, -msse3, -march=, -mfpmath=sse
#include <stdint.h> //int64_t

float
ssum_kh
(const float   *aiarrayfloat_y,
 const int64_t  aiint64_t_n
)
{
  int64_t lint64_i;
  float lf_sum = 0.0f;
  
  for ( lint64_i = 0; lint64_i < aiint64_t_n; lint64_i++) {
    lf_sum += aiarrayfloat_y[lint64_i];
  }
  return lf_sum;
  
}
