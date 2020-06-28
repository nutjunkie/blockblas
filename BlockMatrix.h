/******************************************************************************
 * 
 *  Class declarations for managing block matrices.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero matrices.
 * 
 *  The data in the matrix is not stored explicitly, evaluate() must be 
 *  called to load the give array with  
 *
 *****************************************************************************/

#include "Matrix.h"

class BlockMatrix
{
    public:
       BlockMatrix(unsigned const nRowBlocks, unsigned const nColBlocks) :
          m_nRowBlocks(nRowBlocks), m_nColBlocks(nColBlocks)
       {
          m_blocks = new Matrix[m_nRowBlocks*m_nColBlocks];
       }
       
       ~BlockMatrix()
       {
          delete [] m_blocks;
       }

       Matrix& operator()(unsigned const row, unsigned const col)
       {
           return m_blocks[row*m_nColBlocks + col];
       }

       Matrix const& operator()(unsigned const row, unsigned const col) const
       {
           return m_blocks[row*m_nColBlocks + col];
       }

       unsigned nRowBlocks() const { return m_nRowBlocks; }
       unsigned nColBlocks() const { return m_nColBlocks; }
 
       bool consistent() const;

       Matrix toDense() const;

       void info() const;
       void print() const;
       
    private:
       unsigned nRows() const;
       unsigned nCols() const;

       unsigned m_nRowBlocks;
       unsigned m_nColBlocks;
       Matrix*  m_blocks;
};
