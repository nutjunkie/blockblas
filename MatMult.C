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
      A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nRows(), 
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
      A.nRows(), B.nCols(), A.nCols(), &one, A.data(), A.nRows(), 
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



template<>
void gsmm(VMatrix<double,RowMajor> const& A,  VMatrix<double,RowMajor> const& B,  
   VMatrix<double,RowMajor>& C)
{
   unsigned rowsC(C.nRows());
   unsigned colsC(C.nCols());

   std::vector<int> const& stripes(A.stripes());
   unsigned rowsA(A.nRows());
   unsigned colsA(A.nCols());
   unsigned m(std::min(rowsA,colsA));

   double const* a(A.data());
   double const* b(B.data());
   double*       c(C.data());

   for (unsigned s = 0; s < stripes.size(); ++s) {
       int offset(stripes[s]);
       int offC(std::max(0,-offset));
       int offB(std::max(0, offset));
       // Contraction length
       int len = (offset < 0) ? std::min(rowsA+offset, colsA)
                                 : std::min(rowsA, colsA-offset);

//     std::cout << "offset: " << offset << " offC: " << offC << " offB: " << offB 
//               << " contraction: " <<  len  << std::endl;

#pragma omp parallel for
       for (unsigned i = 0; i < len; ++i) {
           double        a0(a[i + s*m]);
           double const* b0(&b[(offB+i)*colsC]);  // Is this colsC or colsB?
           double*       c0(&c[(offC+i)*colsC]);
//         cblas_daxpy(colsC, a0, b0, 1, c0, 1);
           for (unsigned j = 0; j < colsC; ++j) {
               c0[j] += a0*b0[j];
           }   
       }   
   } 
}


template<>
void gsmm(VMatrix<double, ColumnMajor> const& A,  VMatrix<double, ColumnMajor> const& B,  
   VMatrix<double, ColumnMajor>& C)
{
   unsigned rowsB(B.nRows());
   unsigned rowsC(C.nRows());
   unsigned colsC(C.nCols());

   std::vector<int> const& stripes(A.stripes());
   unsigned rowsA(A.nRows());
   unsigned colsA(A.nCols());
   unsigned m(std::min(rowsA,colsA));

   double const* a(A.data());
   double const* b(B.data());
   double*       c(C.data());

   for (unsigned s = 0; s < stripes.size(); ++s) {
       int offset(stripes[s]);
       int offC(std::max(0,-offset));
       int offB(std::max(0, offset));
       // Contraction length
       //
       int len = (offset < 0) ? std::min(rowsA+offset, colsA)
                              : std::min(rowsA, colsA-offset);

//     std::cout << "offset: " << offset << " offC: " << offC << " offB: " << offB 
//               << " contraction: " <<  len  << std::endl;

       double const* a0(&a[s*m]);
#pragma omp parallel for
       for (unsigned j = 0; j < colsC; ++j) {
           double const* b0(&b[offB+j*rowsB]);
           double*       c0(&c[offC+j*rowsC]);
           for (unsigned i = 0; i < len; ++i) {
               c0[i] += a0[i]*b0[i];
           }   
       }   
   }
}



template<>
void gsmm(VMatrix<complex, RowMajor> const& A,  VMatrix<complex, RowMajor> const& B,  
   VMatrix<complex, RowMajor>& C)
{
    std::cout << "Complex gsmm NYI" << std::endl;
}

template<>
void gsmm(VMatrix<complex, ColumnMajor> const& A,  VMatrix<complex, ColumnMajor> const& B,  
   VMatrix<complex, ColumnMajor>& C)
{
    std::cout << "Complex gsmm NYI" << std::endl;
}
