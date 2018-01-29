//gcc shnrm2_k.c -Wall -O3 -funroll-loops  -S
#include <stdint.h> //int64_t


unsigned int shnrm2_k(const int64_t n, const short *x, const int64_t incx, const short *y, const int64_t incy)
{
  unsigned int  lul_sumDiff =  (unsigned int) 0;
  short li_diff;
  int64_t i;
  int64_t j=0;
  int64_t k=0;

  if ( (incx ==  1) && (incy ==  1) ) {   
    for ( i = 0; i < n; i++ ) {
      li_diff  = x[i] - y[i];
      lul_sumDiff += (unsigned int) li_diff*li_diff;
    }
  }
  else {
    for ( i = 0; i < n; i++ ) {
      li_diff  = x[j] - y[k];
      lul_sumDiff += (unsigned int) li_diff*li_diff;
      j += incx;
      k += incy;
    }
  }
  
  return lul_sumDiff;
}
