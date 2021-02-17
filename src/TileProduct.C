#include "mkl.h"
#include "TileProduct.h"


// double * double = double
template <>
void tile_product(CMTile<double> const& A, CMTile<double> const& B, 
   double const beta, CMTile<double>& C, CBLAS_TRANSPOSE const Atrans)
{
   unsigned aRows(A.nRows());
   unsigned aCols(A.nCols());

   if (Atrans == CblasTrans) {
      std::swap(aRows,aCols);
      //std::cout << "Trans trigger" << std::endl;
   }

   //std::cout << "A rows/cols " << aRows << "/" << aCols << std::endl;
   //std::cout << "B rows/cols " << B.nRows() << "/" << B.nCols() << std::endl;
   //std::cout << "C rows/cols " << C.nRows() << "/" << C.nCols() << std::endl;

   //if (A.nRows() == 1 && A.nCols() == 1 &&
   //    B.nRows() == 1 && B.nCols() == 1 &&
   //    C.nRows() == 1 && C.nCols() == 1) {
   //    C.set(0,0, A(0,0)*B(0,0) + beta*C(0,0));
   //}else {
      cblas_dgemm(CblasColMajor, Atrans, CblasNoTrans, 
         aRows, B.nCols(), aCols, 1.0, A.data(), A.leadingDim(),
         B.data(), B.leadingDim(), beta, C.data(), C.leadingDim());
   //}
}



// complex * complex = complex 
template <>
void tile_product(CMTile<complex> const& A, CMTile<complex> const& B, 
   complex const beta, CMTile<complex>& C, CBLAS_TRANSPOSE const Atrans)
{
   complex one(1.0,0);
   unsigned aRows(A.nRows());
   unsigned aCols(A.nCols());

   if (Atrans == CblasTrans) {
      std::swap(aRows,aCols);
   }

   cblas_zgemm(CblasColMajor, Atrans, CblasNoTrans, 
      aRows, B.nCols(), aCols, &one, A.data(), A.leadingDim(),
      B.data(), B.leadingDim(), &beta, C.data(), C.leadingDim());
}


// double * complex = complex 
template <>
void tile_product(CMTile<double> const& A, CMTile<complex> const& B, 
   complex const beta, CMTile<complex>& C, CBLAS_TRANSPOSE const Atrans)
{
   double zero(0.0);
   double one(1.0);
   unsigned aRows(A.nRows());
   unsigned aCols(A.nCols());

   if (Atrans == CblasTrans) {
      std::swap(aRows,aCols);
   }

   CMTile<double> Bt, Ct(C.nRows(), C.nCols());
   Ct.alloc();
   C.scale(beta);

   B.getReal(Bt);
   cblas_dgemm(CblasColMajor, Atrans, CblasNoTrans, 
      aRows, Bt.nCols(), aCols, one, A.data(), A.leadingDim(),
      Bt.data(), Bt.leadingDim(), zero, Ct.data(), Ct.leadingDim());
   C.addReal(Ct);

   B.getImag(Bt);
   cblas_dgemm(CblasColMajor, Atrans, CblasNoTrans, 
      aRows, Bt.nCols(), aCols, one, A.data(), A.leadingDim(),
      Bt.data(), Bt.leadingDim(), zero, Ct.data(), Ct.leadingDim());
   C.addImag(Ct);
}
