#include "Matrix.h"
#include <algorithm>
#include <iostream>
/******************************************************************************
 * 
 *  Class for managing block matrices.
 * 
 *****************************************************************************/

ZeroFunctor Matrix::s_zero = ZeroFunctor();

void Matrix::bind()
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


void Matrix::fillZero()
{
   m_data = new double[1];
   m_data[0] = 0.0;
}


void Matrix::fillDiagonal()
{
   unsigned m(std::min(m_nRows,m_nCols));
   m_data = new double[m];
   for (unsigned i = 0; i < m; ++i) {
       m_data[i] = (*m_functor)(i,i);
   }
}


/*
x  x  .  .  .  .
x  x  x  .  .  .
.  x  x  x  .  .
.  .  x  x  x  .
.  .  .  x  x  x
.  .  .  .  x  x

.  .  .  .  .  x

.  .  .  .  .  .
.  .  .  .  .  .
.  .  .  .  .  .


0  1  2  3  4  5   6  7  8
x  x  .  .  .  .   .  .  .
x  x  x  .  .  .   .  .  .
.  x  x  x  .  .   .  .  .
.  .  x  x  x  .   .  .  .
.  .  .  x  x  x   .  .  .
.  .  .  .  x  x   x  .  .

* x x
x x x

*/


// Packed band structure using row-major layout
void Matrix::fillTridiagonal()
{
   unsigned  m(std::min(m_nRows,m_nCols));
   m_data = new double[3*m+1];  // +1 if m_nNows > m

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
void Matrix::fillPentadiagonal()
{
   unsigned m(std::min(m_nRows,m_nCols)), k(0);
   m_data = new double[5*m+10];  // +2 if m_nNows > m

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



void Matrix::fillDense()
{
   m_data = new double[m_nRows*m_nCols];

   unsigned k(0);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j, ++k) {
           m_data[k] = (*m_functor)(i,j); 
       }
   }
}


double Matrix::operator()(int const i, int const j) const
{
   // Tmp, this should be replaced by array lookup methods
   //return (*m_functor)(i,j);
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


