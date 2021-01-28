#include "Log.h"
#include "DiagonalTile.h"
#include "CMTile.h"


template <>
void DiagonalTile<double>::takeDiagonal(CMTile<double> const& that)
{
   resize(that.nRows(), that.nCols());

   if (that.isBound()) {
      this->alloc();
      double* a(this->data());
      double const* b(that.data());
      size_t ldb(that.leadingDim());

      for (size_t i = 0; i < this->m_nData; ++i) {
          a[2*i] = b[i+i*ldb];
      } 
   }
}


template <>
void DiagonalTile<complex>::takeDiagonal(CMTile<double> const& that)
{
   resize(that.nRows(), that.nCols());

   if (that.isBound()) {
      this->alloc();
      double* a(reinterpret_cast<double*>(this->data()));
      double const* b(that.data());
      size_t ldb(that.leadingDim());

      for (size_t i = 0; i < this->m_nData; ++i) {
          a[2*i] = b[i+i*ldb];
      } 
   }
}


template <>
void DiagonalTile<complex>::takeDiagonal(CMTile<complex> const& that)
{
   resize(that.nRows(), that.nCols());

   if (that.isBound()) {
      this->alloc();
      complex* a(this->data());
      complex const* b(that.data());
      size_t ldb(that.leadingDim());

      for (size_t i = 0; i < this->m_nData; ++i) {
          a[i] = b[i+i*ldb];
      } 
   }
}
