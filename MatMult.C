#include "MatMult.h"


template<>
void gemm(VMatrix<double,RowMajor> const& A,  VMatrix<double,RowMajor> const& B,  
   VMatrix<double,RowMajor>& C)
{
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nCols(), 
      B.data(), B.nCols(), 1.0, C.data(), C.nCols());
}


template<>
void gemm(VMatrix<double,ColumnMajor> const& A,  VMatrix<double,ColumnMajor> const& B,  
   VMatrix<double,ColumnMajor>& C)
{
   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nCols(), 
      B.data(), B.nRows(), 1.0, C.data(), C.nRows());
}


template<>
void gemm(VMatrix<complex,RowMajor> const& A,  VMatrix<complex,RowMajor> const& B,  
   VMatrix<complex,RowMajor>& C)
{
   complex one(1.0,0.0);
   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), &one, A.data(), A.nCols(), 
      B.data(), B.nCols(), &one, C.data(), C.nCols());
}


template<>
void gemm(VMatrix<complex,ColumnMajor> const& A,  VMatrix<complex,ColumnMajor> const& B,  
   VMatrix<complex,ColumnMajor>& C)
{
   complex one(1.0,0.0);
   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), &one, A.data(), A.nCols(), 
      B.data(), B.nRows(), &one, C.data(), C.nRows());
}



template<>
void gbmv(VMatrix<double,RowMajor> const& A,  VMatrix<double,RowMajor> const& B,  
   VMatrix<double,RowMajor>& C, int const kl, int const ku, int const band, int const vec)
{
   cblas_dgbmv(CblasRowMajor, CblasNoTrans,
      A.nRows(), A.nCols(), kl, ku, 1.0, A.data(), band,
      B.data()+vec, B.nCols(), 1.0, C.data()+vec, C.nCols());
}


template<>
void gbmv(VMatrix<double, ColumnMajor> const& A,  VMatrix<double, ColumnMajor> const& B,  
   VMatrix<double, ColumnMajor>& C, int const kl, int const ku, int const band, int const vec)
{
   cblas_dgbmv(CblasColMajor, CblasNoTrans,
      A.nRows(), A.nCols(), kl, ku, 1.0, A.data(), band,
      B.data()+vec*B.nRows(), 1, 1.0, C.data()+vec*B.nRows(), 1);
}


template<>
void gbmv(VMatrix<complex, RowMajor> const& A,  VMatrix<complex, RowMajor> const& B,  
   VMatrix<complex,RowMajor>& C, int const kl, int const ku, int const band, int const vec)
{
   complex one(1.0,0.0);
   cblas_zgbmv(CblasRowMajor, CblasNoTrans,
      A.nRows(), A.nCols(), kl, ku, &one, A.data(), band,
      B.data()+vec, B.nCols(), &one, C.data()+vec, C.nCols());
}


template<>
void gbmv(VMatrix<complex, ColumnMajor> const& A,  VMatrix<complex, ColumnMajor> const& B,  
   VMatrix<complex,ColumnMajor>& C, int const kl, int const ku, int const band, int const vec)
{
   complex one(1.0,0.0);
   cblas_zgbmv(CblasColMajor, CblasNoTrans,
      A.nRows(), A.nCols(), kl, ku, &one, A.data(), band,
      B.data()+vec*B.nRows(), 1, &one, C.data()+vec*B.nRows(), 1);
}
