#include "MatMult.h"

template <>
void wrap_gemm(
   const CBLAS_ORDER Order, 
   const CBLAS_TRANSPOSE TransA, 
   const CBLAS_TRANSPOSE TransB, 
   const int M, 
   const int N,
   const int K, 
   const double alpha, 
   const double* A, 
   const int lda, 
   const double* B, 
   const int ldb, 
   const double beta, 
   double* C, 
   const int ldc)
{
   cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <>
void wrap_gemm(
   const CBLAS_ORDER Order, 
   const CBLAS_TRANSPOSE TransA, 
   const CBLAS_TRANSPOSE TransB, 
   const int M, 
   const int N,
   const int K, 
   const double alpha, 
   const complex* A, 
   const int lda, 
   const complex* B, 
   const int ldb, 
   const double beta, 
   complex* C, 
   const int ldc)
{
   complex al(alpha,0.0);
   complex be(beta,0.0);
   cblas_zgemm(Order, TransA, TransB, M, N, K, &al, A, lda, B, ldb, &be, C, ldc);
}


template <>
void wrap_gbmv(
   const CBLAS_ORDER Order, 
   const CBLAS_TRANSPOSE TransA, 
   const int M, 
   const int N, 
   const int KL, 
   const int KU, 
   const double alpha, 
   const double* A, 
   const int lda, 
   const double* X, 
   const int incX, 
   const double beta, 
   double* Y, 
   const int incY)
{
   cblas_dgbmv(Order, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}

template <>
void wrap_gbmv(
   const CBLAS_ORDER Order, 
   const CBLAS_TRANSPOSE TransA, 
   const int M, 
   const int N, 
   const int KL, 
   const int KU, 
   const double alpha, 
   const complex* A, 
   const int lda, 
   const complex* X, 
   const int incX, 
   const double beta, 
   complex* Y, 
   const int incY)
{
   complex al(alpha,0.0);
   complex be(beta,0.0);
   cblas_zgbmv(Order, TransA, M, N, KL, KU, &al, A, lda, X, incX, &be, Y, incY);
}


template <>
void wrap_gemv(
   const CBLAS_ORDER Order, 
   const CBLAS_TRANSPOSE TransA, 
   const int M, 
   const int N, 
   const double alpha, 
   const double* A, 
   const int lda, 
   const double* X, 
   const int incX, 
   const double beta, 
   double* Y, 
   const int incY)
{
   cblas_dgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


template <>
void wrap_gemv(
   const CBLAS_ORDER Order, 
   const CBLAS_TRANSPOSE TransA, 
   const int M, 
   const int N, 
   const double alpha, 
   const complex* A, 
   const int lda, 
   const complex* X, 
   const int incX, 
   const double beta, 
   complex* Y, 
   const int incY)
{
   complex al(alpha,0.0);
   complex be(beta,0.0);
   cblas_zgemv(Order, TransA, M, N, &al, A, lda, X, incX, &be, Y, incY);
}
