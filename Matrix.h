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
#include "Functor.h"


// Matrix class representing a tile of the BlockMatrix.  These tiles are "aware" of
// how to compute the elements via the Functor argument

class Matrix
{
    public:
       enum Storage { Zero, Diagonal, Tridiagonal, Pentadiagonal, Dense, Symmetric };

       Matrix(Functor const* functor = &s_zero, unsigned const nRows = 0, 
          unsigned const nCols = 0, Storage const storage = Zero) : 
          m_functor(functor), m_nRows(m_nRows), m_nCols(m_nCols),
          m_storage(storage), m_data(0)
       { }

       // The functor could be redefined to operate on double* directly, if required
       void init(Functor const* functor, unsigned const nRows, unsigned const nCols, 
          Storage const storage)
       {
          m_functor = functor;
          m_nRows   = nRows;
          m_nCols   = nCols;
          m_storage = storage;
       }

       ~Matrix() { release(); }

       unsigned nRows() const { return m_nRows; }
       unsigned nCols() const { return m_nCols; }

       void bind();

       void release() 
       {
          if (m_data) {
             delete[] m_data;
             m_data = 0;
          }
       } 

       // This is an extremely inefficient way of accessing the elements
       // and is included for debug purposes only.
       double operator() (int const i, int const j) const;

       Storage storage() const { return m_storage; }
       bool isZero() const { return m_storage == Zero; }
       bool isDense() const { return m_storage == Dense; }
       bool isBand() const { return m_storage < Dense; }

    private:
       static ZeroFunctor s_zero;

       void fillZero();
       void fillDense();
       void fillDiagonal();
       void fillTridiagonal();
       void fillPentadiagonal();

       Functor const* m_functor; 
       Storage  m_storage;
       unsigned m_nRows;       
       unsigned m_nCols;       
       double*  m_data;
};
