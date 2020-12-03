#include "Log.h"
#include "CMTile.h"


template <>
void CMTile<double>::invert()
{
   size_t nRows(this->m_nRows);
   size_t nCols(this->m_nCols);
   int    lda(this->m_leadingDim);

    if (nRows != nCols || !this->isBound()) {
       std::stringstream ss("ERROR: invert() called on invalid Tile (");
       ss << nRows << "," << nCols << ") -> " << this->storage();
       Log::error(ss.str());
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

    delete [] ipiv;
    delete [] work;
}



template <>
void CMTile<complex>::invert()
{
   size_t nRows(this->m_nRows);
   size_t nCols(this->m_nCols);
   int    lda(this->m_leadingDim);

   if (nRows != nCols || !this->isBound()) {
      std::stringstream ss("ERROR: invert() called on invalid Tile (");
      ss << nRows << "," << nCols << ") -> " << this->storage();
      Log::error(ss.str());
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

   delete [] ipiv;
   delete [] work;
}



template <>
void CMTile<double>::factorLU(int* ipiv)
{
   int nRows(this->m_nRows);
   int nCols(this->m_nCols);
   int lda(this->m_leadingDim);

    if (nRows != nCols || !this->isBound()) {
       std::stringstream ss("factorLU() called on invalid Tile (");
       ss << nRows << "," << nCols << ") -> " << this->storage();
       Log::error(ss.str());
       return;
    }   

    int info = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'U', nRows, this->m_data, lda, ipiv);
    if (info != 0) {
       std::string s("Bad return from dsytrf: ");
       s += std::to_string(info);
       Log::error(s);
    }

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
