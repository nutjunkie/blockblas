#include "RMTile.h"

void tile_product(RMTile<double> const& A, RMTile<double> const& B, 
   double const beta, RMTile<double>& C)
{
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), beta, C.data(), C.leadingDim());
}


void tile_product(RMTile<complex> const& A, RMTile<complex> const& B, 
   complex const beta, RMTile<complex>& C)
{
   complex one(1.0,0);
   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
      A.nRows(), B.nCols(), A.nCols(), &one, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), &beta, C.data(), C.leadingDim());
}
