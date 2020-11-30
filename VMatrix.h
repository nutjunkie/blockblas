#ifndef VMATRIX_H
#define VMATRIX_H
/******************************************************************************
 * 
 *  Class declarations for managing block matrices.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero matrices.
 * 
 *****************************************************************************/

#include <cstddef>
#include <vector>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Functor.h"
#include "Types.h"


// Matrix class representing a tile of the BlockMatrix. 

template <class T, LayoutT L>
class VMatrix
{
    public:
       // On construction, the VMatrix does not have any data associated with it.
       // Allocation is only done when calling bind().
       VMatrix(size_t const nRows = 0, size_t const nCols = 0, 
          StorageT const storage = Dense) : m_data(0)
       {
          init(nRows, nCols, storage);
       }

       VMatrix(VMatrix<T,L> const& that) : m_data(0)
       {
          copy(that);
       }

       VMatrix<T,L>& operator=(VMatrix<T,L> const& that)
       {
          if (this != &that) copy(that);
          return *this;
       }

       VMatrix<T,L>& fromDouble(VMatrix<double,L> const& that);

       void info(const char* msg = 0) const
       {
          if (msg) {
             std::cout << msg << std::endl;
          }
    
          std::cout << "Storage:    " << toString(m_storage) << std::endl;
          std::cout << "Num data:   " << m_nData << std::endl;
          std::cout << "Type size:  " << sizeof(T) << std::endl;
          std::cout << "Dimensions: " << m_nRows << "x" << m_nCols << std::endl;
       }

       VMatrix<T,L>& init(size_t const nRows, size_t const nCols, 
          StorageT const storage = Dense)
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


       VMatrix<T,L>& init(size_t const nRows, size_t const nCols, 
          std::vector<int> const& stripes)
       {
          release();
          m_nRows   = nRows;
          m_nCols   = nCols;
          m_storage = Striped;
          m_stripes = stripes;

          return *this;
       }

       VMatrix<T,L>& init(size_t const nRows, size_t const nCols, size_t const lbands, 
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
          m_data = 0;
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

       void bind(T const* data)
       {
           bind();
           memcpy(m_data, data, m_nData*sizeof(T));
       }

       void unbind(T* data)
       {
           memcpy(data, m_data, m_nData*sizeof(T));
       }

       void bind(T const* data, unsigned const ld);
       void unbind(T* data, unsigned const ld);

       // Allocates space for the VMatrix and computes its elements using the functor.
       void bind(Functor<T> const& functor)
       {
          bind();

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

       void invert();
       bool isBound() const { return m_data != 0; }

       size_t nRows() const { return m_nRows; }
       size_t nCols() const { return m_nCols; }

       // Vector containing the indices of the non-zero diagonal stripes 
       // of the matrix.  A value of 0 corresponds to the major diagonal,
       // negative values are offsets below the diagonal and positive
       // values are above the diagonal.
       std::vector<int> const& stripes() const { return m_stripes; }

       void toDense();
       
       // These are convenience functions and are very inefficient
       T operator()(unsigned const i, unsigned const j) const;
       void set(unsigned const i, unsigned const j, T const value);

       VMatrix<T,L>& operator+=(VMatrix<T,L> const& that)
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

       VMatrix<T,L>& operator-=(VMatrix<T,L> const& that)
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

/* This is dangerous and should be deprecated */
       VMatrix<T,L>& operator-()
       {
          for (unsigned i = 0; i < m_nData; ++i) {
              m_data[i] = -m_data[i];
          }

          return *this;
       }

       // Adds the value t to the diagonals
       VMatrix<T,L>& operator+=(T const t)
       {
          unsigned k(std::min(m_nCols, m_nRows));
          T val;
          for (unsigned i = 0; i < k; ++i) {
              val = (*this)(i,i);
              set(i,i,val+t);
          }
          return *this;
       }

       double norm2() const;

       void print(const char* msg = 0) const
       {
          if (msg) {
             std::cout << msg << std::endl;
          }   
          std::cout << std::fixed << std::showpoint << std::setprecision(2);
          for (unsigned i = 0; i < m_nRows; ++i) {
              for (unsigned j = 0; j < m_nCols; ++j) {
                  std::cout << std::setw(5) << (*this)(i,j) << " ";
              }   
              std::cout << std::endl;
          }   
          std::cout << std::endl;
       }

       StorageT storage() const { return m_storage; }
       LayoutT layout() const;


       bool isZero() const { return m_storage == Zero; }
       bool isDense() const { return m_storage == Dense; }
       bool isStriped() const { return m_storage == Striped; }
       bool isDiagonal() const { return m_storage == Diagonal; }

// This needs cleaning up 
   public:
//    protected:
         void fillZero(Functor<T> const& functor)
         {
             // This represents a zero block matrix where the entries are not
             // explicitly stored.  To initialize a matrix with zeros use
             // the appropriate storage type and the ZeroFunctor.
             m_data[0] = T();
          }


       void fillDense(Functor<T> const& functor);
       void fillBanded(Functor<T> const& functor);
       void fillStriped(Functor<T> const& functor);

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
       std::vector<int> m_stripes;
       T* m_data;

   private:
      void copy(VMatrix<T,L> const& that)
      {
         m_nRows   = that.m_nRows;
         m_nCols   = that.m_nCols;
         m_storage = that.m_storage;
         m_stripes = that.m_stripes;

         release();
   
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
