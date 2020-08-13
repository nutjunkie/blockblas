#include "VMatrix.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <veclib/veclib.h>
/******************************************************************************
 * 
 *  Virtual matrix class
 * 
 *  
 * 
 *****************************************************************************/


void VMatrix::bind()
{
   release();

   switch (m_storage) {
      case Zero:
         fillZero();
         break;
      case Dense:
         fillDense();
         break;
      case Diagonal:
         fillDiagonal();
         break;
      case Tridiagonal:
         fillTridiagonal();
         break;
      case Pentadiagonal:
         fillPentadiagonal();
         break;
      default:
         fillDense();
         break;
   }
}


void VMatrix::fillZero()
{
   // This represents a zero block matrix where the entries are not
   // explicitly stored.  To initialize a non-zero block matrix use
   // the appropriate storage type and the ZeroFunctor.
   m_data = new double[1];
   m_data[0] = 0.0;
}


void VMatrix::fillDiagonal()
{
   unsigned m(std::min(m_nRows,m_nCols));
   m_data = new double[m];
   if (!m_functor) return;

   for (unsigned i = 0; i < m; ++i) {
       m_data[i] = (*m_functor)(i,i);
   }
}


// Packed band structure using row-major layout
void VMatrix::fillTridiagonal()
{
   unsigned  m(std::min(m_nRows,m_nCols));
   m_data = new double[3*m+1];  // +1 if m_nNows > m
   if (!m_functor) return;

   // First row
   unsigned k(0);
   m_data[k++] = 0.0;
   m_data[k++] = (*m_functor)(0,0);
   m_data[k++] = (*m_functor)(0,1);

   for (unsigned i = 1; i < m-1; ++i) {
       m_data[k++] = (*m_functor)(i,i-1);
       m_data[k++] = (*m_functor)(i,i  );
       m_data[k++] = (*m_functor)(i,i+1);
   }

   // Clean up the tassels
   m_data[k++] = (*m_functor)(m-1,m-2);
   m_data[k++] = (*m_functor)(m-1,m-1);

   if (m_nCols > m) {
      m_data[k++] = (*m_functor)(m-1,m);
   }else if (m_nRows > m) {
      m_data[k++] = 0.0;
      m_data[k++] = (*m_functor)(m,m-1);
   }

   //for (unsigned i=0; i < k; ++i) std::cout<< i << " " << m_data[i] << std::endl;

   if (k > 3*m+1)  std::cout << "Range check: " << 3*m+1 << " < " << k << std::endl;
}


// Packed band structure using row-major layout
void VMatrix::fillPentadiagonal()
{
   unsigned m(std::min(m_nRows,m_nCols)), k(0);
   m_data = new double[5*m+10];  // +2 if m_nNows > m
   if (!m_functor) return;

   // First row
   m_data[k++] = 0.0;
   m_data[k++] = 0.0;
   m_data[k++] = (*m_functor)(0,0);
   m_data[k++] = (*m_functor)(0,1);
   m_data[k++] = (*m_functor)(0,2);

   // Second row
   m_data[k++] = 0.0;
   m_data[k++] = (*m_functor)(1,0);
   m_data[k++] = (*m_functor)(1,1);
   m_data[k++] = (*m_functor)(1,2);
   m_data[k++] = (*m_functor)(1,3);

   for (unsigned i = 2; i < m-2; ++i) {
       m_data[k++] = (*m_functor)(i,i-2);
       m_data[k++] = (*m_functor)(i,i-1);
       m_data[k++] = (*m_functor)(i,i  );
       m_data[k++] = (*m_functor)(i,i+1);
       m_data[k++] = (*m_functor)(i,i+2);
   }

   //seond to last row
   m_data[k++] = (*m_functor)(m-2,m-4);
   m_data[k++] = (*m_functor)(m-2,m-3);
   m_data[k++] = (*m_functor)(m-2,m-2);
   m_data[k++] = (*m_functor)(m-2,m-1);

   if (m_nCols <= m_nRows) m_data[k++] = 0.0;
   if (m_nCols > m_nRows ) m_data[k++] = (*m_functor)(m-2,m);
      
   // last row
   m_data[k++] = (*m_functor)(m-1,m-3);
   m_data[k++] = (*m_functor)(m-1,m-2);
   m_data[k++] = (*m_functor)(m-1,m-1);

   if (m_nCols > m_nRows  ) m_data[k++] = (*m_functor)(m-1,m);
   if (m_nCols > m_nRows+1) m_data[k++] = (*m_functor)(m-1,m+1);

   if (m_nRows > m_nCols) {
      m_data[k++] = 0.0;
      m_data[k++] = 0.0;
      m_data[k++] = (*m_functor)(m,m-2);
      m_data[k++] = (*m_functor)(m,m-1);
   }
   if (m_nRows > m_nCols+1) {
      m_data[k++] = 0.0;
      m_data[k++] = 0.0;
      m_data[k++] = 0.0;
      m_data[k++] = (*m_functor)(m+1,m-1);
   }



   for (unsigned i=0; i < k; ++i) std::cout<< i << " " << m_data[i] << std::endl;

   if (k > 5*m+2)  std::cout << "Range check: " << 5*m+2 << " < " << k << std::endl;
}



void VMatrix::fillDense()
{
   m_data = new double[m_nRows*m_nCols];
   if (!m_functor) return;

   unsigned k(0);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j, ++k) {
           m_data[k] = (*m_functor)(i,j); 
       }
   }
}


double VMatrix::operator()(int const i, int const j) const
{
   double value(0.0);

   switch (m_storage) {
      case Zero:
         break;
      case Diagonal:
         if (i==j) value = m_data[i];
         break;
      case Tridiagonal:
         if (std::abs(j-i) <= 1) {
            value = m_data[i*3+(j-i)+1];
         }
         break;
      case Pentadiagonal:
         if (std::abs(j-i) <= 2) {
            value = m_data[i*5+(j-i)+2];
         }
         break;
      case Dense:
         value = m_data[i*m_nCols + j];
         break;
      default:
         break;
   }

   return value;
}


void VMatrix::set(int const i, int const j, double value)
{
   switch (m_storage) {
      case Zero:
         break;
      case Diagonal:
         if (i==j) m_data[i] = value;
         break;
      case Tridiagonal:
         if (std::abs(j-i) <= 1) {
            m_data[i*3+(j-i)+1] = value;
         }
         break;
      case Pentadiagonal:
         if (std::abs(j-i) <= 2) {
            m_data[i*5+(j-i)+2] = value;
         }
         break;
      case Dense:
         m_data[i*m_nCols + j] = value;
         break;
      default:
         break;
   }
}


void VMatrix::print() const
{
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j) {
           std::cout << std::setw(5) << (*this)(i,j) << " ";
       }
       std::cout << std::endl;
   }
   std::cout << std::endl;
}


// Accumulates the A.B product into C:
//    C += A.B 
// Also handles the case when B and C are vectors.
void VMatrix::matrix_product(VMatrix& C, VMatrix& A, VMatrix& B)
{
   if (A.nCols() != B.nRows() ||
       B.nCols() != C.nCols() ||
       A.nRows() != C.nRows()) {

       std::cerr << "Barf on the dimensions:" << std::endl;
       std::cout << A.nCols() << " != " << B.nRows() << " || " << std::endl;
       std::cout << B.nCols() << " != " << C.nCols() << " || " << std::endl;
       std::cout << A.nRows() << " != " << C.nRows() << std::endl;
   }

   switch (A.storage()) {
      case VMatrix::Zero:
         // Nothing to do
         break;
      case VMatrix::Diagonal: 
         std::cerr << "matrix_product for VMatrix::Diagonal NYI" << std::endl;
         break;
      case VMatrix::Tridiagonal:
         std::cerr << "matrix_product for VMatrix::Tridiagonal NYI" << std::endl;
         break;
      case VMatrix::Pentadiagonal:
         std::cerr << "matrix_product for VMatrix::Pentadiagonal NYI" << std::endl;
         break;
      case VMatrix::Dense:
         if (C.nCols() == 1) {
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
              A.nRows(), A.nCols(), 1.0, A.data(), A.nCols(),
              B.data(), 1, 0.0, C.data(), 1);
         }else {
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nCols(),
              B.data(), B.nCols(), 1.0, C.data(), C.nCols());
         }
         break;
   }
}


/*
VMatrix& VMatrix::operator*=(VMatrix const& rhs)
{

. . .  . .  . . . . .     .
. . .  . .  . . . . .     .
. . .  . .  . . . . .     .

. . .  . .  . . . . .     .
. . .  . .  . . . . .     .

. . .  . .  . . . . .     .
. . .  . .  . . . . .     .
. . .  . .  . . . . .     .
. . .  . .  . . . . .     .
. . .  . .  . . . . .     .

   }
}
*/
