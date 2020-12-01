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

 

template <>
void CMTile<double>::invert()
{
   size_t nRows(this->m_nRows);
   size_t nCols(this->m_nCols);
   int    lda(this->m_leadingDim);

    if (nRows != nCols || !this->isBound()) {
       std::cerr << "ERROR: invert() called on invalid Tile (" 
                 << nRows << "," << nCols << ") -> " << this->storage() << std::endl;
       return;
    }   

    int n(nRows);
    int *ipiv = new int[n+1];
    int lwork = n*n;
    int info;
    double* work = new double[lwork];

#ifdef __INTEL_COMPILER
    dgetrf(&n,&n,this->m_data,&lda,ipiv,&info);
    dgetri(&n,   this->m_data,&lda,ipiv,work,&lwork,&info);
#else
    dgetrf_(&n,&n,this->m_data,&lda,ipiv,&info);
    dgetri_(&n,   this->m_data,&lda,ipiv,work,&lwork,&info);
#endif

    delete ipiv;
    delete work;
}



template <>
void CMTile<complex>::invert()
{
   size_t nRows(this->m_nRows);
   size_t nCols(this->m_nCols);
   int    lda(this->m_leadingDim);

   if (nRows != nCols || !this->isBound()) {
      std::cerr << "invert() called on invalid matrix (" 
                << nRows << "," << nCols << ") -> " << this->storage() << std::endl;
      return;
   }   

   int n(nRows);
   int *ipiv = new int[n+1];
   int lwork = n*n;
   int info;
   complex* work = new complex[lwork];

#ifdef __INTEL_COMPILER
   zgetrf(&n,&n,this->m_data,&lda,ipiv,&info);
   zgetri(&n,   this->m_data,&lda,ipiv,work,&lwork,&info);
#else
   zgetrf_(&n,&n,this->m_data,&lda,ipiv,&info);
   zgetri_(&n,   this->m_data,&lda,ipiv,work,&lwork,&info);
#endif

   delete ipiv;
   delete work;
}


template <>
double CMTile<double>::norm2() const
{
   double n2(0.0);
   unsigned lda(m_leadingDim);

   for (unsigned j = 0; j < this->m_nCols; ++j) {
       double* a0(this->m_data+j*lda);
       for (unsigned i = 0; i < this->m_nRows; ++i) {
           n2 += a0[i] * a0[i];
       }
   }

   return n2;
}


template <>
double CMTile<complex>::norm2() const
{
   double n(0.0);
   double n2(0.0);
   unsigned lda(m_leadingDim);

   for (unsigned j = 0; j < this->m_nCols; ++j) {
       complex* a0(this->m_data+j*lda);
       for (unsigned i = 0; i < this->m_nRows; ++i) {
           n   = std::norm(a0[i]);
           n2 += n*n;
       }
   }

   return n2;
}
