#ifndef VMATRIX_H
#define VMATRIX_H
/******************************************************************************
 * 
 *  Class declarations for managing block matrices.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero matrices.
 * 
 *  The data in the matrix is not stored explicitly, evaluate() must be 
 *  called to load the give array with  
 *
 *****************************************************************************/

#include <cstddef>
#include <vector>
//#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Functor.h"
#include "Types.h"


// Matrix class representing a tile of the BlockMatrix. 

template <class T>
class VMatrix
{
    public:
       // On construction, the VMatrix does not have any data associated with it.
       // Allocation is only done when calling bind().
       VMatrix(size_t const nRows = 0, size_t const nCols = 0, 
          StorageT const storage = Dense) : m_data(0), m_layout(RowMajor)
       {
          init(nRows, nCols, storage);
       }

       VMatrix(VMatrix<T> const& that) : m_data(0), m_layout(RowMajor)
       {
          copy(that);
       }

       VMatrix<T>& operator=(VMatrix<T> const& that)
       {
          if (this != &that) copy(that);
          return *this;
       }

       VMatrix<T>& init(size_t const nRows, size_t const nCols, StorageT const storage = Dense)
       {
          release();
          m_nRows   = nRows;
          m_nCols   = nCols;
          m_storage = storage;

          if (m_storage == Striped) {
             std::cerr << "WARNING: Invalid initialization for striped VMatrix" << std::endl;
          }else if (m_storage == Banded) {
             std::cerr << "WARNING: Invalid initialization for banded VMatrix" << std::endl;
          }

          return *this;
       }


       VMatrix<T>& init(size_t const nRows, size_t const nCols, 
          std::vector<int> const& stripes)
       {
          release();
          m_nRows   = nRows;
          m_nCols   = nCols;
          m_storage = Striped;
          m_stripes = stripes;

          return *this;
       }

       VMatrix<T>& init(size_t const nRows, size_t const nCols, size_t const lbands, 
          size_t const ubands)
       {
          release();
          m_nRows   = nRows;
          m_nCols   = nCols;
          m_storage = Banded;

          m_stripes.push_back(lbands);
          m_stripes.push_back(ubands);

          return *this;
       }

       ~VMatrix() { 
          release(); 
       }

       void release() 
       {
          if (m_data) {
             delete [] m_data;
          }
       }

       T* data() { return m_data; }
       T const* data() const { return m_data; }

       // Allocates space for the VMatrix.  This should be modifed to
       // return a smart pointer to the data.

       void bind()
       //std::shared_ptr<T> VMatrix::bind()
       {
          release();

          switch (m_storage) {
             case Zero:
                m_nData = 1;
                break;
             case Diagonal:
                m_nData = std::min(m_nRows,m_nCols);
                break;
             case Banded:
                m_nData = m_nRows*(m_stripes[0]+m_stripes[1]+1);  // kl + ku +1
                std::cout << "allocating " << m_nData << " elements for banded matrix (" 
                          << std::setprecision(2) << std::showpoint
                          << 100.0*m_nData/(m_nRows*m_nCols) << " %)" << std::endl;
                break;
             case Striped:
                m_nData = std::min(m_nRows,m_nCols) * m_stripes.size();
                break;
             case Dense:
                m_nData = m_nRows*m_nCols;
                break;
          }

          m_data = new T[m_nData];
          //return std::make_shared(m_data);
       }


       // Allocates space for the VMatrix and computes its elements using the functor.
       void bind(Functor<T> const& functor)
       {
          bind();
          m_layout = RowMajor;

          switch (m_storage) {
             case Zero:
                fillZero(functor);
                break;
             case Diagonal:
                fillDiagonal(functor);
                break;
             case Banded:
                fillBanded(functor);
                break;
             case Striped:
                fillStriped(functor);
                break;
             case Dense:
                fillDense(functor);
                break;
          }
       }

       // Column-major form
       void bindCM(Functor<T> const& functor)
       {
          bind();
          m_layout = ColumnMajor;

          switch (m_storage) {
             case Zero:
                fillZero(functor);
                break;
             case Diagonal:
                fillDiagonal(functor);
                break;
             case Banded:
                fillBandedCM(functor);
                break;
             case Striped:
                std::cerr << "Column-major storage not supported for striped matrices" << std::endl;
                break;
             case Dense:
                fillDenseCM(functor);
                break;
          }
       }

       void invert();

       bool isBound() const { return m_data != 0; }

       size_t nRows() const { return m_nRows; }
       size_t nCols() const { return m_nCols; }

       // Vector containing the indices of the non-zero diagonal stripes 
       // of the matrix.  A value of 0 corresponds to the major diagonal,
       // negative values are offsets below the diagonal and positive
       // values are above the diagonal.
       std::vector<int> const& stripes() const { return m_stripes; }

       void toDense()
       {
          T* data = new T[m_nRows*m_nCols];
          unsigned k(0);

          if (m_layout == RowMajor) {
             for (unsigned i = 0; i < m_nRows; ++i) {
                 for (unsigned j = 0; j < m_nCols; ++j, ++k) {
                     data[k] = (*this)(i,j);
                 }
             }  
          }else {
             for (unsigned j = 0; j < m_nCols; ++j) {
                 for (unsigned i = 0; i < m_nRows; ++i, ++k) {
                     data[k] = (*this)(i,j);
                 }
             }  
          }

          release();
          m_data = data;
          m_storage = Dense;
          m_stripes.clear();
       }

       // These are convenience functions and are very inefficient
       T operator()(unsigned const i, unsigned const j) const
       {
          T value = T();
 
          switch (m_storage) {
             case Zero:
                break;
 
             case Diagonal:
                if (i==j) {
                   value = m_data[i];
                }
                break;
 
             case Banded: {
                int kl(m_stripes[0]);
                int ku(m_stripes[1]);
//             std::cout << "Access element: (" << i << "," << j << ") -> ";
                if (m_layout == RowMajor) {
                   if (std::max(0,(int)i-kl) <= j && j <= std::min((int)m_nCols,(int)i+ku)) {
                      int k(j-i+kl+i*(kl+ku+1));
                      value = m_data[k];
//                    std::cout <<  k << " = ";
                   }
                }else {
                   if (std::max(0,(int)j-ku) <= i && i <= std::min((int)m_nRows,(int)j+kl)) {
                      int k(i-j+ku+j*(kl+ku+1));
                      value = m_data[k];
//                    std::cout <<  k << " = ";
                   }
                }
//              std::cout <<  value << std::endl;
             } break;
 
             case Striped: {
                if (m_layout == RowMajor) {
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
                }else {
                   std::cerr << "Column-major layout for striped matrices not supported" << std::endl;
                }
 
             } break;
 
             case Dense:
                if (m_layout == RowMajor) {
                   value = m_data[i*m_nCols + j];
                }else {
                   value = m_data[i + j*m_nRows];
                }
                break;
          }
 
          return value;
       }
 
        //T& operator() (unsigned const i, unsigned const j);
       void set(unsigned const i, unsigned const j, T const value)
       {
          switch (m_storage) {
             case Zero:
                break;
 
             case Diagonal:
                if (i==j) {
                   m_data[i] = value;
                }
                break;
 
            case Banded:
               std::cerr << "VMatrix::set NYI for Banded matrices" << std::endl;
               break;
 
             case Striped: {
                if (m_layout == RowMajor) {
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
                }else {
                   std::cerr << "VMatrix::set NYI for Striped ColumnMajor matrices" << std::endl;
                }
 
             } break;
 
             case Dense:
                if (m_layout == RowMajor) {
                   m_data[i*m_nCols + j] = value;
                }else {
                   m_data[i + j*m_nRows] = value;
                }
                break;
          }
       }
 
       //This needs to account for the different storage types
       VMatrix<T>& operator+=(VMatrix<T> const& that)
       {
#ifdef DEBUG
          assert(that.m_nCols   == m_nCols);
          assert(that.m_nRows   == m_nRows);
          assert(that.m_nData   == m_n);
          assert(that.m_storage == m_storage);
#endif
          for (unsigned i = 0; i < m_nData; ++i) {
              m_data[i] += that.m_data[i];
          }
          return *this;
       }


       //This needs to account for the different storage types
       VMatrix<T>& operator-=(VMatrix<T> const& that)
       {
#ifdef DEBUG
          assert(that.m_nCols   == m_nCols);
          assert(that.m_nRows   == m_nRows);
          assert(that.m_nData   == m_n);
          assert(that.m_storage == m_storage);
#endif
          for (unsigned i = 0; i < m_nData; ++i) {
              m_data[i] -= that.m_data[i];
          }
          return *this;
       }


       VMatrix<T>& operator-()
       {
          for (unsigned i = 0; i < m_nData; ++i) {
              m_data[i] = -m_data[i];
          }

          return *this;
       }

       double norm2() const;

       void print(const char* msg = 0) const;

       StorageT storage() const { return m_storage; }
       LayoutT  layout() const { return m_layout; }

       bool isZero() const { return m_storage == Zero; }
       bool isDense() const { return m_storage == Dense; }
       bool isStriped() const { return m_storage == Striped; }
       bool isDiagonal() const { return m_storage == Diagonal; }

    protected:
       void fillZero(Functor<T> const& functor)
       {
          // This represents a zero block matrix where the entries are not
          // explicitly stored.  To initialize a zero block matrix use
          // the appropriate storage type and the ZeroFunctor.
          m_data[0] = T();
       }


       void fillDense(Functor<T> const& functor)
       {
          unsigned k(0);
          for (unsigned i = 0; i < m_nRows; ++i) {
              for (unsigned j = 0; j < m_nCols; ++j, ++k) {
                  m_data[k] = functor(i,j); 
              }
          }
       }


       void fillDenseCM(Functor<T> const& functor)
       {
          unsigned k(0);
          for (unsigned j = 0; j < m_nCols; ++j) {
              for (unsigned i = 0; i < m_nRows; ++i, ++k) {
                  m_data[k] = functor(i,j); 
              }
          }
       }


       void fillBanded(Functor<T> const& functor)
       {
          int kl(m_stripes[0]);
          int ku(m_stripes[1]);
          int k;

          for (int i = 0; i < m_nRows; ++i) {
              int jmin = std::max(0,i-kl);
              int jmax = std::min((int)m_nCols,i+ku+1);
//            std::cout << "j range for i = " << i << " -> (" << jmin << "..." << jmax << ")" <<std::endl;
              for (int j = jmin ; j < jmax; ++j) {
                  k = j-i+kl+i*(kl+ku+1);
//                std::cout << "setting "<< k<< " to " << functor(i,j) << std::endl;
                  m_data[k] = functor(i,j);
              }
          }
       }


       void fillBandedCM(Functor<T> const& functor)
       {
          int kl(m_stripes[0]);
          int ku(m_stripes[1]);
          int k;

          for (int j = 0; j < m_nCols; ++j) {
              int imin = std::max(0, j-ku);
              int imax = std::min((int)m_nRows, j+kl+1);
              //std::cout << "i range for j = " << j << " -> (" << imin << "..." << imax << ")" <<std::endl;
              for (int i = imin ; i < imax; ++i) {
                  k = i-j+ku+j*(kl+ku+1);
                  m_data[k] = functor(i,j);
              }
          }
       }


       void fillStriped(Functor<T> const& functor)
       {
          unsigned nStripes(m_stripes.size());
          unsigned m(std::min(m_nRows,m_nCols));

          for (unsigned k = 0; k < nStripes; ++k) {
              int offset(m_stripes[k]);
              if (offset < 0) {
                 unsigned max(std::min(m_nRows + offset,m_nCols));
//               std::cout << "offset = " << offset << " running to " << max << std::endl;
                 for (unsigned j = 0; j < max; ++j) {
//                    std::cout << "  setting data = " << j + k*m<< " to " 
//                              << functor(j-offset,j) << std::endl;
                     m_data[j + k*m] = functor(j-offset,j);
                 }
              }else {
                 unsigned max(std::min(m_nRows, m_nCols-offset));
//               std::cout << "offset = " << offset << " running to " << max << std::endl;
                 for (unsigned i = 0; i < max; ++i) {
//                    std::cout << "  setting data = " << i + k*m << " to " 
//                              << functor(i,i+offset) << std::endl;
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


       void fillDiagonal(Functor<T> const& functor)
       {
          unsigned m(std::min(m_nRows,m_nCols));
          for (unsigned i = 0; i < m; ++i) {
              m_data[i] = functor(i,i);
          }
       }

       size_t m_nRows;       
       size_t m_nCols;       
       size_t m_nData;       

       StorageT m_storage;
       LayoutT  m_layout;
       std::vector<int> m_stripes;
       T* m_data;

   private:
      void copy(VMatrix<T> const& that)
      {
         m_nRows = that.m_nRows;
         m_nCols = that.m_nCols;

         m_storage = that.m_storage;
         m_layout  = that.m_layout;
         m_stripes = that.m_stripes;
         m_data    = 0;
   
         if (that.isBound()) {
            m_nData = that.m_nData;
            m_data = new T[m_nData];
            for (unsigned i = 0; i < m_nData; ++i) {    
                m_data[i] = that.m_data[i];
            }
         }
      }
};

#endif
