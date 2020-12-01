#include "CMTile.h"
#include "DiagonalTile.h"


template <>
void tile_product(CMTile<double> const& A, CMTile<double> const& B, 
   double const beta, CMTile<double>& C)
{
   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), beta, C.data(), C.leadingDim());
}


template <>
void tile_product(CMTile<complex> const& A, CMTile<complex> const& B, 
   complex const beta, CMTile<complex>& C)
{
   complex one(1.0,0);
   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), &one, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), &beta, C.data(), C.leadingDim());
}


template <>
void tile_product(StripedTile<double> const& A, CMTile<double> const& B,
   double const beta, CMTile<double>& C)
{
   unsigned rowsA(A.nRows());
   unsigned colsA(A.nCols());
   unsigned colsC(C.nCols());

   unsigned ldb(B.nRows());
   unsigned ldc(C.nRows());

   unsigned m(std::min(rowsA,colsA));

   double const* a(A.data());
   double const* b(B.data());
   double*       c(C.data());

   std::vector<int> const& stripes(A.stripes());

   for (unsigned s = 0; s < stripes.size(); ++s) {
       int offset(stripes[s]);
       int offC(std::max(0,-offset));
       int offB(std::max(0, offset));
       // Contraction length
       int len = (offset < 0) ? std::min(rowsA+offset, colsA)
                              : std::min(rowsA, colsA-offset);
       double const* a0(&a[s*m]);
#pragma omp parallel for
       for (unsigned j = 0; j < colsC; ++j) {
           double const* b0(&b[offB+j*ldb]);
           double*       c0(&c[offC+j*ldc]);
           for (unsigned i = 0; i < len; ++i) {
               c0[i] += a0[i]*b0[i];
           }   
       }   
   }   
}

