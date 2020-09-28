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

#include <string>
#include <cstddef>
#include <vector>
#include <memory>
#include "Functor.h"


// Matrix class representing a tile of the BlockMatrix. 

class VMatrix
{
    public:
       enum StorageT { Zero, Diagonal, Banded, Striped, Dense };
       enum LayoutT { RowMajor, ColumnMajor };

       static std::string toString(StorageT);

       // On construction, the VMatrix does not have any data associated with it.
       // Allocation is only done when calling bind().
       VMatrix(size_t const nRows = 0, size_t const nCols = 0, 
          StorageT const storage = Dense) : m_data(0), m_layout(RowMajor)
       {
          init(nRows, nCols, storage);
       }

       VMatrix(VMatrix const& that) : m_data(0), m_layout(RowMajor)
       {
          copy(that);
       }

       VMatrix& operator=(VMatrix const& that)
       {
          if (this != &that) copy(that);
          return *this;
       }


       VMatrix& init(size_t const nRows, size_t const nCols, 
          StorageT const storage = Dense);

       VMatrix& init(size_t const nRows, size_t const nCols, 
          std::vector<int> const& stripes);

       VMatrix& init(size_t const nRows, size_t const nCols, 
          size_t const lbands, size_t const ubands);

       ~VMatrix() { 
          release(); 
       }

       void release() 
       {
          if (m_data) {
             delete [] m_data;
          }
       }

       // Allocates space for the VMatrix.  This should be modifed to
       // return a smart pointer to the data.
       //std::shared_ptr<double> bind();
       void bind();
       double* data() { return m_data; }
       double const* data() const { return m_data; }

       // Allocates space for the VMatrix and computes its elements
       // using the functor.
       void bind(Functor const& functor);
       // Column-major form
       void bindCM(Functor const& functor);

       bool isBound() const { return m_data != 0; }

       size_t nRows() const { return m_nRows; }
       size_t nCols() const { return m_nCols; }

       // Vector containing the indices of the non-zero diagonal stripes 
       // of the matrix.  A value of 0 corresponds to the major diagonal,
       // negative values are offsets below the diagonal and positive
       // values are above the diagonal.
       std::vector<int> const& stripes() const { return m_stripes; }

       void toDense();

       void invert();

       // These are convenience functions and are very inefficient
       double operator() (unsigned const i, unsigned const j) const;
       //double& operator() (unsigned const i, unsigned const j);
       void set(unsigned const i, unsigned const j, double value);
 
       VMatrix& operator+=(VMatrix const& that);
       VMatrix& operator-=(VMatrix const& that);
       VMatrix& operator-();

       double norm2() const;

       void print(const char* = 0) const;

       StorageT storage() const { return m_storage; }
       LayoutT  layout() const { return m_layout; }

       bool isZero() const { return m_storage == Zero; }
       bool isDense() const { return m_storage == Dense; }
       bool isStriped() const { return m_storage == Striped; }
       bool isDiagonal() const { return m_storage == Diagonal; }

    protected:
       void fillZero(Functor const&);
       void fillDense(Functor const&);
       void fillDenseCM(Functor const&);
       void fillBanded(Functor const&);
       void fillBandedCM(Functor const&);
       void fillStriped(Functor const&);
       void fillDiagonal(Functor const&);

       size_t m_nRows;       
       size_t m_nCols;       
       size_t m_nData;       

       StorageT m_storage;
       LayoutT  m_layout;
       std::vector<int> m_stripes;
       double* m_data;

   private:
      void copy(VMatrix const&);
       
};

#endif
