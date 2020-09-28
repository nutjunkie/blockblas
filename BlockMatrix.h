#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H
/******************************************************************************
 * 
 *  Class declarations for managing block matrices.  The interface is desiged
 *  to homogenize the handling of both dense, banded and zero matrices.
 * 
 *  The data in the matrix is not stored explicitly, evaluate() must be 
 *  called to load the give array with  
 *
 *****************************************************************************/

#include "VMatrix.h"

class Functor;

class BlockMatrix
{
    public:
       BlockMatrix(unsigned const nRowBlocks, unsigned const nColBlocks) :
          m_nRowBlocks(nRowBlocks), m_nColBlocks(nColBlocks)
       {
          m_blocks = new VMatrix[m_nRowBlocks*m_nColBlocks];
       }

       BlockMatrix(BlockMatrix const& that);

/*
       BlockMatrix& operator=(BlockMatrix const& that)
       {
          if (this != &that) copy(that);
          return *this;
       }
*/

       BlockMatrix& operator+=(BlockMatrix const&);

       BlockMatrix& operator-=(BlockMatrix const&);
       BlockMatrix& operator-();
       
       ~BlockMatrix()
       {
          delete [] m_blocks;
       }

       VMatrix& operator()(unsigned const row, unsigned const col = 0)
       {
           return m_blocks[row*m_nColBlocks + col];
       }

       VMatrix const& operator()(unsigned const row, unsigned const col = 0) const
       {
           return m_blocks[row*m_nColBlocks + col];
       }

       void bind(Functor&);

       unsigned nRowBlocks() const { return m_nRowBlocks; }
       unsigned nColBlocks() const { return m_nColBlocks; }

       unsigned nRows() const;
       unsigned nCols() const;
 
       bool consistent() const;

       void toDense(VMatrix*) const;

       unsigned rowOffset(unsigned bi) const;
       unsigned colOffset(unsigned bj) const;

       static void matrix_product(BlockMatrix& C, BlockMatrix& A, BlockMatrix& B);

       void info(const char* = 0) const;
       void print(const char* = 0) const;
       
    private:
       unsigned m_nRowBlocks;
       unsigned m_nColBlocks;
       VMatrix*  m_blocks;
};
#endif
