//gcc inrm2_k.c -O3 -funroll-loops  -S
#include <stdint.h> //int64_t

unsigned long int inrm2_k(const int64_t n, const int *x, const int64_t incx, const int *y, const int64_t incy)
{
  unsigned long int  lul_sumDiff =  (unsigned long int) 0;
  int li_diff;
  int64_t i;
  int64_t j = 0;
  int64_t k = 0;

  if ( (incx ==  1) && (incy ==  1) ) {   
    for ( i = 0; i < n; i++ ) {
      li_diff  = x[i] - y[i];
      lul_sumDiff += (unsigned long int) li_diff*li_diff;
    }
  }
  else {
    for ( i = 0; i < n; i++ ) {
      li_diff  = x[j] - y[k];
      lul_sumDiff += (unsigned long int) li_diff*li_diff;
      j += incx;
      k += incy;
    }
  }
  
  return lul_sumDiff;
}
