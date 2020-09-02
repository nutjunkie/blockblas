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
#include "Functor.h"


// Matrix class representing a tile of the BlockMatrix.  These tiles are "aware" of
// how to compute the elements via the Functor argument

class VMatrix
{
    public:
       enum StorageT { Zero, Diagonal, Tridiagonal, Pentadiagonal, Dense, Striped};

       static std::string toString(StorageT);

       // If no functor is given, the data is left uninitialized on bind()
       VMatrix(unsigned const nRows = 0, unsigned const nCols = 0, 
           StorageT const storage = Dense, Functor const* functor = 0) :
          m_nRows(m_nRows), m_nCols(m_nCols),
          m_storage(storage), m_functor(functor), m_data(0), m_stripes(0)
       { 
       }

       // The functor could be redefined to operate on double* directly, if required
       VMatrix& init(unsigned const nRows, unsigned const nCols, 
          StorageT const storage = Dense, Functor const* functor = 0)
       {
          m_nRows   = nRows;
          m_nCols   = nCols;
          m_storage = storage;
          m_functor = functor;
          return *this;
       }

       ~VMatrix() { 
          release(); 
       }

       unsigned nRows() const { return m_nRows; }
       unsigned nCols() const { return m_nCols; }

       void toDense();

       
       std::vector<int> const&  stripes() { return m_stripes; }

       void bind();

       void bindStripes(std::vector<int> const& stripes);

       void release() 
       {
          if (m_data) {
             delete[] m_data;
             m_data = 0;
          }
       } 


       bool isBound() const { return m_data != 0; }

       // Use this to fill the raw data array
       double* data() { return m_data; }

       double operator() (int const i, int const j) const;

       void set(int const i, int const j, double value);

       VMatrix& operator*=(VMatrix const& rhs);

       void print(const char* = 0) const;
       
       static void matrix_product(VMatrix& c, VMatrix& A, VMatrix& b);

       StorageT storage() const { return m_storage; }
       bool isZero() const { return m_storage == Zero; }
       bool isDense() const { return m_storage == Dense; }
       //bool isBand() const { return m_storage < Striped; }
       bool isStriped() const { return m_storage == Striped; }

    private:
       void fillZero();
       void fillDense();
       void fillDiagonal();
       void fillTridiagonal();
       void fillPentadiagonal();

       unsigned m_nRows;       
       unsigned m_nCols;       

       StorageT m_storage;
       Functor const* m_functor; 
       std::vector<int> m_stripes;
       double* m_data;
};

#endif
