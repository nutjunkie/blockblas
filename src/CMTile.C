#include "Log.h"
#include "CMTile.h"
#include "TileArray.h"


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
void CMTile<complex>::factorLU(int* ipiv)
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

    int info = LAPACKE_zsytrf(LAPACK_COL_MAJOR, 'U', nRows, this->m_data, lda, ipiv);
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



template <>
CMTile<double>::CMTile(TileArray<double> const& TA) : Tile<double>(TA.nRows(), TA.nCols())
{
   m_leadingDim = TA.nRows();

   this->alloc();

   size_t nRowTiles(TA.nRowTiles());
   size_t nColTiles(TA.nColTiles());

   for (size_t bj = 0; bj < nColTiles; ++bj) {
       size_t joff(TA.colOffset(bj));
       
       for (size_t bi = 0; bi < nRowTiles; ++bi) {
           size_t ioff(TA.rowOffset(bi));

           Tile<double>const& t(TA.tile(bi,bj));

           for (size_t j = 0; j < t.nCols(); ++j) {
               for (size_t i = 0; i < t.nRows(); ++i) {
                   this->set(ioff+i, joff+j, t(i,j));
               }
           }
       }
   }
}



template <>
CMTile<complex>::CMTile(TileArray<complex> const& TA) : Tile<complex>(TA.nRows(), TA.nCols())
{
   m_leadingDim = TA.nRows();

   this->alloc();

   size_t nRowTiles(TA.nRowTiles());
   size_t nColTiles(TA.nColTiles());

   for (size_t bj = 0; bj < nColTiles; ++bj) {
       size_t joff(TA.colOffset(bj));
       
       for (size_t bi = 0; bi < nRowTiles; ++bi) {
           size_t ioff(TA.rowOffset(bi));
           Tile<complex>const& t(TA.tile(bi,bj));

           for (size_t j = 0; j < t.nCols(); ++j) {
               for (size_t i = 0; i < t.nRows(); ++i) {
                   this->set(ioff+i, joff+j, t(i,j));
               }
           }
       }
   }
}



template <>
void CMTile<complex>::getReal(CMTile<double>& that) const
{
   size_t nRows(this->nRows());
   size_t nCols(this->nCols());

   that.resize(nRows, nCols);
   that.alloc();

   double const* a0 = reinterpret_cast<double const*>(this->data());
   double*       b0(that.data());

   size_t lda(m_leadingDim);
   size_t ldb(that.leadingDim());

   for (unsigned j = 0; j < nCols; ++j) {
       double const* a(a0+2*j*lda);            
       double*       b(b0+  j*ldb);            
       for (unsigned i = 0; i < nRows; ++i) {
           b[i] = a[2*i];
       }
   }
}


template <>
void CMTile<complex>::getImag(CMTile<double>& that) const
{
   size_t nRows(this->nRows());
   size_t nCols(this->nCols());

   that.resize(nRows, nCols);
   that.alloc();

   double const* a0 = reinterpret_cast<double const*>(this->data());
   double*       b0(that.data());

   size_t lda(m_leadingDim);
   size_t ldb(that.leadingDim());

   for (unsigned j = 0; j < nCols; ++j) {
       double const* a(a0+1+2*j*lda);            
       double*       b(b0+  j*ldb);            
       for (unsigned i = 0; i < nRows; ++i) {
           b[i] = a[2*i];
       }
   }
}


template <>
void CMTile<complex>::addReal(CMTile<double> const& that)
{
   //asserts

   size_t nRows(this->nRows());
   size_t nCols(this->nCols());
   size_t lda(m_leadingDim);
   size_t ldb(that.leadingDim());

   double*       a0 = reinterpret_cast<double*>(this->data());
   double const* b0(that.data());

   for (unsigned j = 0; j < nCols; ++j) {
       double*       a(a0 +2*j*lda);            
       double const* b(b0+  j*ldb);            
       for (unsigned i = 0; i < nRows; ++i) {
           a[2*i] += b[i];
       }
   }
}

template <>
void CMTile<complex>::addImag(CMTile<double> const& that)
{
   //asserts

   size_t nRows(this->nRows());
   size_t nCols(this->nCols());
   size_t lda(m_leadingDim);
   size_t ldb(that.leadingDim());

   double*       a0 = reinterpret_cast<double*>(this->data());
   double const* b0(that.data());

   for (unsigned j = 0; j < nCols; ++j) {
       double*       a(a0+1+2*j*lda);            
       double const* b(b0+  j*ldb);            
       for (unsigned i = 0; i < nRows; ++i) {
           a[2*i] += b[i];
       }
   }
}
