#include "VMatrix.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
/******************************************************************************
 * 
 *  Virtual matrix class
 * 
 *****************************************************************************/


std::string VMatrix::toString(StorageT storage)
{
   std::string s;
   switch (storage) {
      case Zero:      s = "Zero";      break;
      case Diagonal:  s = "Diagonal";  break;
      case Striped:   s = "Striped";   break;
      case Dense:     s = "Dense";     break;
   }

   return s;
}


VMatrix& VMatrix::init(size_t const nRows, size_t const nCols, StorageT const storage)
{
   release();
   m_nRows   = nRows;
   m_nCols   = nCols;
   m_storage = storage;

   if (m_storage == Striped) {
      std::cerr << "WARNING: Invalid initialization for striped VMatrix" << std::endl;
   }

   return *this;
}


VMatrix& VMatrix::init(size_t const nRows, size_t const nCols, 
   std::vector<int> const& stripes)
{
   release();
   m_nRows   = nRows;
   m_nCols   = nCols;
   m_storage = Striped;
   m_stripes = stripes;

   return *this;
}


void VMatrix::bind()
//std::shared_ptr<double> VMatrix::bind()
{
   release();
   size_t n(0);

   switch (m_storage) {
      case Zero:
         n = 1;
         break;
      case Diagonal:
         n = std::min(m_nRows,m_nCols);
         break;
      case Striped:
         n = std::min(m_nRows,m_nCols) * m_stripes.size();
         break;
      case Dense:
         n = m_nRows*m_nCols;
         break;
   }

   m_data = new double[n];

   //return std::make_shared(m_data);
}


void VMatrix::bind(Functor const& functor)
{
   bind();

   switch (m_storage) {
      case Zero:
         fillZero(functor);
         break;

      case Diagonal:
         fillDiagonal(functor);
         break;

      case Striped:
         fillStriped(functor);
         break;

      case Dense:
         fillDense(functor);
         break;
   }
}


void VMatrix::fillZero(Functor const& functor)
{
   // This represents a zero block matrix where the entries are not
   // explicitly stored.  To initialize a zero block matrix use
   // the appropriate storage type and the ZeroFunctor.
   m_data[0] = 0.0;
}


void VMatrix::fillDense(Functor const& functor)
{
   unsigned k(0);
   for (unsigned i = 0; i < m_nRows; ++i) {
       for (unsigned j = 0; j < m_nCols; ++j, ++k) {
           m_data[k] = functor(i,j); 
       }
   }
}


void VMatrix::fillStriped(Functor const& functor)
{
   unsigned nStripes(m_stripes.size());
   unsigned m(std::min(m_nRows,m_nCols));

   for (unsigned k = 0; k < nStripes; ++k) {
       int offset(m_stripes[k]);
       if (offset < 0) {
          unsigned max(std::min(m_nRows + offset,m_nCols));
//        std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned j = 0; j < max; ++j) {
/*
               std::cout << "  setting data = " << j + k*m<< " to " 
                         << functor(j-offset,j) << std::endl;
*/
              m_data[j + k*m] = functor(j-offset,j);
          }
       }else {
          unsigned max(std::min(m_nRows, m_nCols-offset));
 //       std::cout << "offset = " << offset << " running to " << max << std::endl;
          for (unsigned i = 0; i < max; ++i) {
/*
               std::cout << "  setting data = " << i + k*m << " to " 
                         << functor(i,i+offset) << std::endl;
*/
              m_data[i + k*m] = functor(i,i+offset);
          }
       } 
   }

/*
   for (int i = 0; i < nStripes*m; ++i) {
       std::cout << "striped data " << i << " "<< m_data[i]<< std::endl;
   }
*/
}


void VMatrix::fillDiagonal(Functor const& functor)
{
   unsigned m(std::min(m_nRows,m_nCols));
   for (unsigned i = 0; i < m; ++i) {
       m_data[i] = functor(i,i);
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
   m_stripes.clear();
}



double VMatrix::operator()(unsigned const i, unsigned const j) const
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
   }

   return value;
}


void VMatrix::set(unsigned const i, unsigned const j, double value)
{
   switch (m_storage) {
      case Zero:
         break;

      case Diagonal:
         if (i==j) {
            m_data[i] = value;
         }
         break;

      case Striped: {
         int stripe(j-i);
         std::vector<int>::iterator it;
         it = std::find(m_stripes.begin(), m_stripes.end(), stripe);

         if (it != m_stripes.end()) {
            // We have hit a non-zero element
            unsigned m(std::min(m_nRows,m_nCols));
            unsigned index = std::distance(m_stripes.begin(), it);
            int ij = (stripe < 0) ? j : i;
            m_data[ij + index*m] = value;
         }

         } break;

      case Dense:
         m_data[i*m_nCols + j] = value;
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


/*
// Accumulates the A.B product into C:
//    C += A.B 
void VMatrix::matrix_product(VMatrix& C, VMatrix& A, VMatrix& B)
{
   if (A.nCols() != B.nRows() ||
       B.nCols() != C.nCols() ||
       A.nRows() != C.nRows()) {

       std::cerr << "Barf on the matrix dimensions:" << std::endl;
       std::cerr << A.nCols() << " != " << B.nRows() << " || " << std::endl;
       std::cerr << B.nCols() << " != " << C.nCols() << " || " << std::endl;
       std::cerr << A.nRows() << " != " << C.nRows() << std::endl;
   }

   if (!A.isBound() || !B.isBound() || !C.isBound()) {
      std::cerr << "Unbound matrix encountered in VMatrix::matrix_product" << std::endl;
   }

   StorageT storageA(A.storage());   double* a(A.data());
   StorageT storageB(B.storage());   double* b(B.data());
   StorageT storageC(C.storage());   double* c(C.data());

//   std::cerr << "VMatrix multiplication for " << toString(storageC) << " <- " 
//             << toString(storageA) << " x " << toString(storageB) << std::endl;
             

   //                     ---------- Zero Matrices ----------
   if (storageA == VMatrix::Zero || 
       storageB == VMatrix::Zero) {
      // Nothing to do

   //                     ---------- Diagonal Matrices ----------
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

   //                     ---------- Striped Matrices ----------
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

  }else if (storageA == VMatrix::Dense   && 
            storageB == VMatrix::Striped &&
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
*/
