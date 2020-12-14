#include "TileProduct.h"


template <>
void tile_product(CMTile<double> const& A, CMTile<double> const& B, 
   double const beta, CMTile<double>& C, CBLAS_TRANSPOSE const Atrans)
{
   unsigned aRows(A.nRows());
   unsigned aCols(A.nCols());

   if (Atrans == CblasTrans) {
      aRows = A.nCols();
      aCols = A.nRows();
   }

   cblas_dgemm(CblasColMajor, Atrans, CblasNoTrans, 
      aRows, B.nCols(), aCols, 1.0, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), beta, C.data(), C.leadingDim());
}



template <>
void tile_product(CMTile<complex> const& A, CMTile<complex> const& B, 
   complex const beta, CMTile<complex>& C, CBLAS_TRANSPOSE const Atrans)
{
   complex one(1.0,0);
   unsigned aRows(A.nRows());
   unsigned aCols(A.nCols());

   if (Atrans == CblasTrans) {
      aRows = A.nCols();
      aCols = A.nRows();
   }

   cblas_zgemm(CblasColMajor, Atrans, CblasNoTrans, 
      aRows, B.nCols(), aCols, &one, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), &beta, C.data(), C.leadingDim());
}

