//gcc daysxpy_kh.c -O3 -funroll-loops -msse3  -S
//Flags: -msse, -msse2, -msse3, -mavx -march=, -mfpmath=sse
#include <stdint.h> //int64_t

int   daysxpy_kh(const int64_t n, int64_t A, int64_t B, const double alpha, const double *x, const int64_t incx, double *y, const int64_t incy, double *C, int64_t D)
{
  int64_t j = 0;
  int64_t k = 0;
  int64_t i;

  if ( incx  == 1 && incy == 1 ) {
    for ( i = 0; i < n; i++) {
      /*y = y + a(y-x)*/
      y[i] += alpha * (y[i]-x[i]);
    }
  }
  else {
    for ( i = 0; i < n; i++) { 
      y[k] += alpha * (y[k]-x[j]);
      j += incx;
      k += incy;
    }
  }
  return 1;
}
