/*! \file jacobi_eigenvalue.hpp
 *
 * \brief Function jacobi_eigenvalue are part of
 *
 * \details  This file is part of John Burkardt 2013 
 * <a href="https://people.sc.fsu.edu/~jburkardt/cpp_src/jacobi_eigenvalue/">Jacobi Eigenvalue</a>
 *
 * \copyright John Burkardt 2013, JACOBI_EIGENVALUE library is distributed under the GNU LGPL license.
 */

#ifndef JACOBI_EIGENVALUE_HPP
#define JACOBI_EIGENVALUE_HPP


/*! \namespace mat
  \brief Matrix module and associated operations
  \details Implementation of the data type Matrix and operations, uses OpenBLAS when compiling with this option, otherwise functions that are not based in the Interface to Streaming SIMD Extensions (SSE).
  
  \version 1.0
  \date   2015-2017
  \copyright GPLv3 license
*/

namespace mat {
  
template <class T>
void r8mat_identity ( int n, T a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IDENTITY sets the square matrix A to the identity.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, T A[N*N], the N by N identity matrix.
//
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}
//****************************************************************************80

  

template <class T>
void r8mat_diag_get_vector ( int n, T a[], T v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, T A[N*N], the N by N matrix.
//
//    Output, T V[N], the diagonal entries
//    of the matrix.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return;
}
//****************************************************************************80


template <class T>
T r8mat_is_eigen_right ( int n, int k, T a[], T x[],
  T lambda[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.
//
//  Discussion:
//
//    An R8MAT is a matrix of Ts.
//
//    This routine computes the Frobenius norm of
//
//      A * X - X * LAMBDA
//
//    where
//
//      A is an N by N matrix,
//      X is an N by K matrix (each of K columns is an eigenvector)
//      LAMBDA is a K by K diagonal matrix of eigenvalues.
//
//    This routine assumes that A, X and LAMBDA are all real.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int K, the number of eigenvectors.
//    K is usually 1 or N.
//
//    Input, T A[N*N], the matrix.
//
//    Input, T X[N*K], the K eigenvectors.
//
//    Input, T LAMBDA[K], the K eigenvalues.
//
//    Output, T R8MAT_IS_EIGEN_RIGHT, the Frobenius norm
//    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
//    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
//
{
  T *c;
  T error_frobenius;
  int i;
  int j;
  int l;

  c = new T[n*k];

  for ( j = 0; j < k; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
      for ( l = 0; l < n; l++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+l*n] * x[l+j*n];
      }
    }
  }

  for ( j = 0; j < k; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] - lambda[j] * x[i+j*n];
    }
  }

  error_frobenius = r8mat_norm_fro ( n, k, c );

  delete [] c;

  return error_frobenius;
}
  
//****************************************************************************80

template <class T>
T r8mat_norm_fro ( int m, int n, T a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DIST_FRO returns the Frobenius norm of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_DIST_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, T A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, T R8MAT_DIST_FRO, the Frobenius norm of A.
//
{
  int i;
  int j;
  T value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

  



//****************************************************************************80
template <class T>
void jacobi_eigenvalue ( int n, T a[], int it_max, T v[], 
  T d[], int &it_num, int &rot_num )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
//
//  Discussion:
//
//    This function computes the eigenvalues and eigenvectors of a
//    real symmetric matrix, using Rutishauser's modfications of the classical
//    Jacobi rotation method with threshold pivoting. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2013
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, T A[N*N], the matrix, which must be square, real,
//    and symmetric.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, T V[N*N], the matrix of eigenvectors.
//
//    Output, T D[N], the eigenvalues, in descending order.
//
//    Output, int &IT_NUM, the total number of iterations.
//
//    Output, int &ROT_NUM, the total number of rotations.
//
{
  T *bw;
  T c;
  T g;
  T gapq;
  T h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int p;
  int q;
  T s;
  T t;
  T tau;
  T term;
  T termp;
  T termq;
  T theta;
  T thresh;
  T w;
  T *zw;

  r8mat_identity ( n, v );

  r8mat_diag_get_vector ( n, a, d );

  bw = new T[n];
  zw = new T[n];

  for ( i = 0; i < n; i++ )
  {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0;
  rot_num = 0;

  while ( it_num < it_max )
  {
    it_num = it_num + 1;
//
//  The convergence threshold is based on the size of the elements in
//  the strict upper triangle of the matrix.
//
    thresh = 0.0;
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < j; i++ )
      {
        thresh = thresh + a[i+j*n] * a[i+j*n];
      }
    }

    thresh = sqrt ( thresh ) / ( T ) ( 4 * n );

    if ( thresh == 0.0 )
    {
      break;
    }

    for ( p = 0; p < n; p++ )
    {
      for ( q = p + 1; q < n; q++ )
      {
        gapq = 10.0 * fabs ( a[p+q*n] );
        termp = gapq + fabs ( d[p] );
        termq = gapq + fabs ( d[q] );
//
//  Annihilate tiny offdiagonal elements.
//
        if ( 4 < it_num &&
             termp == fabs ( d[p] ) &&
             termq == fabs ( d[q] ) )
        {
          a[p+q*n] = 0.0;
        }
//
//  Otherwise, apply a rotation.
//
        else if ( thresh <= fabs ( a[p+q*n] ) )
        {
          h = d[q] - d[p];
          term = fabs ( h ) + gapq;

          if ( term == fabs ( h ) )
          {
            t = a[p+q*n] / h;
          }
          else
          {
            theta = 0.5 * h / a[p+q*n];
            t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
            if ( theta < 0.0 )
            {
              t = - t;
            }
          }
          c = 1.0 / sqrt ( 1.0 + t * t );
          s = t * c;
          tau = s / ( 1.0 + c );
          h = t * a[p+q*n];
//
//  Accumulate corrections to diagonal elements.
//
          zw[p] = zw[p] - h;                 
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          a[p+q*n] = 0.0;
//
//  Rotate, using information from the upper triangle of A only.
//
          for ( j = 0; j < p; j++ )
          {
            g = a[j+p*n];
            h = a[j+q*n];
            a[j+p*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = p + 1; j < q; j++ )
          {
            g = a[p+j*n];
            h = a[j+q*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = q + 1; j < n; j++ )
          {
            g = a[p+j*n];
            h = a[q+j*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[q+j*n] = h + s * ( g - h * tau );
          }
//
//  Accumulate information in the eigenvector matrix.
//
          for ( j = 0; j < n; j++ )
          {
            g = v[j+p*n];
            h = v[j+q*n];
            v[j+p*n] = g - s * ( h + g * tau );
            v[j+q*n] = h + s * ( g - h * tau );
          }
          rot_num = rot_num + 1;
        }
      }
    }

    for ( i = 0; i < n; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
//
//  Restore upper triangle of input matrix.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[i+j*n] = a[j+i*n];
    }
  }
//
//  Ascending sort the eigenvalues and eigenvectors.
//
  for ( k = 0; k < n - 1; k++ )
  {
    m = k;
    for ( l = k + 1; l < n; l++ )
    {
      if ( d[l] < d[m] )
      {
        m = l;
      }
    }

    if ( m != k )
    {
      t    = d[m];
      d[m] = d[k];
      d[k] = t;
      for ( i = 0; i < n; i++ )
      {
        w        = v[i+m*n];
        v[i+m*n] = v[i+k*n];
        v[i+k*n] = w;
      }
    }
  }

  delete [] bw;
  delete [] zw;

  return;
}
//****************************************************************************80

} /*END namespace mat*/

#endif /* JACOBI_EIGENVALUE_HPP */
