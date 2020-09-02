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


std::string VMatrix::toString(StorageT storage)
{
   std::string s;
   switch (storage) {
      case Zero:           s = "Zero";           break;
      case Diagonal:       s = "Diagonal";       break;
      case Tridiagonal:    s = "Tridiagonal";    break;
      case Pentadiagonal:  s = "Pentadiagonal";  break;
      case Dense:          s = "Dense";          break;
      case Striped:        s = "Striped";        break;
   }

   return s;
}


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
      case Striped:
         fillPentadiagonal();
         break;
      default:
         fillDense();
         break;
   }
}


void VMatrix::bindStripes(std::vector<int> const& stripes)
{
   release();
   m_storage = Striped;
   m_stripes = stripes;

   unsigned nStripes(m_stripes.size());
   unsigned m(std::min(m_nRows,m_nCols));
   m_data = new double[nStripes*m];

   for (int i = 0; i < nStripes*m_nCols; ++i) {
       m_data[i] = 0;
   }

   if (!m_functor) return;

   for (unsigned k = 0; k < nStripes; ++k) {
       int offset(m_stripes[k]);
       if (offset < 0) {
          unsigned max(std::min(m_nRows + offset,m_nCols));
//        std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned j = 0; j < max; ++j) {
/*
               std::cout << "  setting data = " << j + k*m<< " to " 
                         << (*m_functor)(j-offset,j) << std::endl;
*/
              m_data[j + k*m] = (*m_functor)(j-offset,j);
          }
       }else {
          unsigned max(std::min(m_nRows, m_nCols-offset));
 //       std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned i = 0; i < max; ++i) {
/*
               std::cout << "  setting data = " << i + k*m << " to " 
                         << (*m_functor)(i,i+offset) << std::endl;
*/
              m_data[i + k*m] = (*m_functor)(i,i+offset);
          }
       } 
   }

/*
   for (int i = 0; i < nStripes*m; ++i) {
       std::cout << "striped data " << i << " "<< m_data[i]<< std::endl;
   }
*/
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

   //for (unsigned i=0; i < k; ++i) std::cout<< i << " " << m_data[i] << std::endl;

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


void VMatrix::toDense()
{
   double* data = new double[m_nRows*m_nCols];

   unsigned k(0);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j, ++k) {
           data[k] = (*this)(i,j);
       }
   }  

   release();
   m_data = data;
   m_storage = Dense;
}


double VMatrix::operator()(int const i, int const j) const
{
   double value(0.0);

   switch (m_storage) {
      case Zero:
         break;
      case Diagonal:
         if (i==j) {
            value = m_data[i];
         }
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
      case Striped: {
         int stripe(j-i);
         std::vector<int>::const_iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);

         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            unsigned m(std::min(m_nRows,m_nCols));
            unsigned index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            value = m_data[ij + index*m];
         }

         } break;
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
      case Striped:
         //std::cout << "Attempt to set on Striped storage" << std::cout;
         break;
      default:
         break;
   }
}


void VMatrix::print(const char* msg) const
{
   if (msg) {
      std::cout << msg << std::endl;
   }
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

   if (!A.isBound() || !B.isBound() || !C.isBound()) {
      std::cerr << "Unbound matrix encountered in VMatrix::matrix_product" << std::endl;
   }

   StorageT storageA(A.storage()), storageB(B.storage()), storageC(C.storage());

   double* a(A.data());
   double* b(B.data());
   double* c(C.data());

//   std::cerr << "VMatrix multiplication for " << toString(storageC) << " <- " 
//             << toString(storageA) << " x " << toString(storageB) << std::endl;
             
   if (storageA == VMatrix::Zero || 
       storageB == VMatrix::Zero) {
      // Nothing to do

   }else if (storageA == VMatrix::Diagonal && 
             storageB == VMatrix::Diagonal &&
             storageC == VMatrix::Diagonal) {
     unsigned nc(C.nCols());
     for (unsigned i = 0; i < A.nCols(); ++i) {
         c[i] = a[i]*b[i];
     }
 
   }else if (storageA == VMatrix::Diagonal && 
             storageB == VMatrix::Diagonal &&
             storageC == VMatrix::Dense) {
     unsigned nc(C.nCols());
     for (unsigned i = 0; i < A.nCols(); ++i) {
         c[i*nc+i] = a[i]*b[i];
     }
      
   }else if (storageA == VMatrix::Diagonal && 
             storageB == VMatrix::Dense &&
             storageC == VMatrix::Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     for (unsigned i = 0; i < nr; ++i) {
         for (unsigned j = 0; j < nc; ++j) {
             c[i*nc+j] = a[i] * b[i*nc+j];
         }
     }

  }else if (storageA == VMatrix::Dense && 
            storageB == VMatrix::Diagonal &&
            storageC == VMatrix::Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     for (unsigned i = 0; i < nr; ++i) {
         for (unsigned j = 0; j < nc; ++j) {
             c[i*nc+j] = a[i*nc+j] * b[j];
         }
     }

  // Striped multiplication 
  }else if (storageA == VMatrix::Striped && 
            storageB == VMatrix::Dense   &&
            storageC == VMatrix::Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     std::vector<int> const& stripes(A.stripes());
     unsigned nStripes(stripes.size());
     unsigned rowsA(A.nRows());
     unsigned colsA(A.nCols());
     unsigned m(std::min(rowsA,colsA));
  
     for (unsigned s = 0; s < nStripes; ++s) {
         int offset(stripes[s]);

         if (offset < 0) { // lower triangular
            unsigned k(std::min(rowsA + offset,colsA));
//            std::cout << "contraction = " << offset << " running to " << k << std::endl;
            for (unsigned i = 0; i < k; ++i) {
                double x(a[i + s*m]);
                for (unsigned j = 0; j < nc; ++j) {
//std::cout << "C("<<(i-offset) << "," << j <<") = " << x << " x " << b[i*nc+j]<< std::endl;
                    c[(i-offset)*nc+j] += x*b[i*nc+j];
                }
           }
        }else {
          unsigned k(std::min(rowsA, colsA-offset));
//            std::cout << "contraction = " << offset << " running to " << k << std::endl;

          for (unsigned i = 0; i < k; ++i) {
                double x(a[i + s*m]);
                for (unsigned j = 0; j < nc; ++j) {
//std::cout << "C("<<(i-offset) << "," << j <<") = " << x << " x " << b[i*nc+j]<< std::endl;
                    c[i*nc+j] += x*b[(i+offset)*nc+j];
                }
          }
        } 
     }

   }else if (storageA == VMatrix::Dense && 
             storageB == VMatrix::Dense && 
             storageC == VMatrix::Dense) {
      if (C.nCols() == 1) {
         // handle the case when B and C are vectors.
         cblas_dgemv(CblasRowMajor, CblasNoTrans,
            A.nRows(), A.nCols(), 1.0, A.data(), A.nCols(),
            B.data(), 1, 1.0, C.data(), 1);
      } else {
         cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nCols(),
             B.data(), B.nCols(), 1.0, C.data(), C.nCols());
      }

   }else {
      std::cerr << "VMatrix multiplication not defined for " 
                << toString(storageC) << " <- " 
                << toString(storageA) << " x " << toString(storageB) << std::endl;
   }
}
