//gcc iriscal_kh.c -O3 -funroll-loops  -S
#include <stdint.h> //int64_t
int   iriscal_kh(const int64_t n, int64_t a, int64_t b, const double alpha, int  *x, const int64_t incx, int* c, int64_t d, int* e, int64_t f)
{
  int64_t i;
  int64_t j = 0;
  if ( incx ==  1 ) {
    for ( i = 0; i < n; i++) { 
      x[i] = (int) (alpha * (double)x[i]);
    } 
  }
  else {
    for ( i = 0; i < n; i++) { 
      x[j] = (int) (alpha * (double)x[j]);
      j += incx;
    }
  }

    return 1;
}
  
