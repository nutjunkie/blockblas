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
       enum StorageT { Zero, Diagonal, Striped, Dense };

       static std::string toString(StorageT);

       // On construction, the VMatrix does not have any data associated with it.
       // Allocation is only done when calling bind().
       VMatrix(size_t const nRows = 0, size_t const nCols = 0, 
          StorageT const storage = Dense) : m_data(0)
       {
          init(nRows, nCols, storage);
       }

       VMatrix& init(size_t const nRows, size_t const nCols, 
          StorageT const storage = Dense);

       VMatrix& init(size_t const nRows, size_t const nCols, 
          std::vector<int> const& stripes);

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

       // Allocates space for the VMatrix and computes its elements
       // using the functor.
       void bind(Functor const& functor);

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
       double operator() (unsigned const i, unsigned const j) const;
       //double& operator() (unsigned const i, unsigned const j);
       void set(unsigned const i, unsigned const j, double value);

       VMatrix& operator*=(VMatrix const& rhs);

       void print(const char* = 0) const;

       StorageT storage() const { return m_storage; }
       bool isZero() const { return m_storage == Zero; }
       bool isDense() const { return m_storage == Dense; }
       bool isStriped() const { return m_storage == Striped; }
       bool isDiagonal() const { return m_storage == Diagonal; }

    private:
       void fillZero(Functor const&);
       void fillDense(Functor const&);
       void fillStriped(Functor const&);
       void fillDiagonal(Functor const&);

       size_t m_nRows;       
       size_t m_nCols;       

       StorageT m_storage;
       std::vector<int> m_stripes;
       double* m_data;
};

#endif
