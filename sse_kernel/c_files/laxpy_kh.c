//gcc laxpy_kh.c -O3 -funroll-loops  -S
#include <stdint.h> //int64_t

int   laxpy_kh(const int64_t n, int64_t a, int64_t b, const int alpha, const int *x, const int64_t incx, long int *y, const int64_t incy, long int* c , int64_t d)
{
  int64_t j = 0;
  int64_t k = 0;
  int64_t i;

  if ( incx  == 1 && incy == 1 ) {
    for ( i = 0; i < n; i++) { 
      y[i] += alpha * x[i];
    }
  }
  else {
    for ( i = 0; i < n; i++) { 
      y[k] += alpha * x[j];
      j += incx;
      k += incy;
    }
  }
  return 1;
}
