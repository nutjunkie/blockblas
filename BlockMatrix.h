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
#include <algorithm>
#include <iostream>
#include <iomanip>

class Functor;

class BlockMatrix
{
    public:
       BlockMatrix(unsigned const nRowBlocks, unsigned const nColBlocks) :
          m_nRowBlocks(nRowBlocks), m_nColBlocks(nColBlocks)
       {
          m_blocks = new VMatrix[m_nRowBlocks*m_nColBlocks];
       }

       BlockMatrix(BlockMatrix const& that)
       {
          m_nRowBlocks = that.m_nRowBlocks;
          m_nColBlocks = that.m_nColBlocks;
          m_blocks = new VMatrix[m_nRowBlocks*m_nColBlocks];

          for (unsigned row = 0; row < m_nRowBlocks; ++row) {
              for (unsigned col = 0; col < m_nColBlocks; ++col) {
                  m_blocks[row*m_nColBlocks + col] = that(row,col);
             }
          }
       }

/*
       BlockMatrix& operator=(BlockMatrix const& that)
       {
          if (this != &that) copy(that);
          return *this;
       }
*/

       BlockMatrix& BlockMatrix::operator+=(BlockMatrix const& that)
       {
#ifdef DEBUG
          assert(that.nRowBlocks() == m_nRowBlocks);
          assert(that.nColBlocks() == m_nColBlocks);
#endif

          for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
              for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
                  (*this)(bi,bj) += that(bi,bj);
              }
          }
       }

       BlockMatrix& BlockMatrix::operator-=(BlockMatrix const& that)
       {
#ifdef DEBUG
          assert(that.nRowBlocks() == m_nRowBlocks);
          assert(that.nColBlocks() == m_nColBlocks);
#endif

          for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
              for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
                  (*this)(bi,bj) -= that(bi,bj);
              }
          }
       }

       BlockMatrix& BlockMatrix::operator-()
       {
          for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
              for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
                  -(*this)(bi,bj);
              }
          }
       }
       
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

       void bind(Functor& functor) 
       {
          for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
              for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
                  (*this)(bi,bj).bind(functor);
              }
          }
       }

       unsigned nRowBlocks() const { return m_nRowBlocks; }
       unsigned nColBlocks() const { return m_nColBlocks; }

       unsigned nRows() const
       {
          if (m_nRowBlocks == 0 || m_nColBlocks == 0) return 0;

          unsigned n(0);
          for (unsigned i = 0; i < m_nRowBlocks; ++i) {
              n += (*this)(i,0).nRows();
          }
          return n;
       }

       unsigned nCols() const
       {
          if (m_nRowBlocks == 0 || m_nColBlocks == 0) return 0;

          unsigned n(0);
          for (unsigned i = 0; i < m_nColBlocks; ++i) {
              n += (*this)(0,i).nCols();
          }
          return n;
       }
 
       bool consistent() const;
       bool consistent() const
       {
          if (m_nRowBlocks == 0 || m_nColBlocks == 0) return true;

          // Check row heights
          for (unsigned i = 0; i < m_nRowBlocks; ++i) {
              unsigned n = (*this)(i,0).nRows();

              for (unsigned j = 1; j < m_nColBlocks; ++j) {
                  if (n != (*this)(i,j).nRows()) return false;
              }
          }

          // Check column widths
          for (unsigned i = 0; i < m_nColBlocks; ++i) {
              unsigned n = (*this)(0,i).nCols();

              for (unsigned j = 1; j < m_nRowBlocks; ++j) {
                  if (n != (*this)(j,i).nCols()) return false;
              }
          }

          return true;
       }

       void toDense(VMatrix*) const;
       void BlockMatrix::toDense(VMatrix* vm) const
       {
          vm->init(nRows(), nCols(), VMatrix::Dense).bind();

          for (unsigned bi = 0; bi < nRowBlocks(); ++bi) {
              unsigned iOff(rowOffset(bi));
              for (unsigned bj = 0; bj < nColBlocks(); ++bj) {
                  unsigned jOff(colOffset(bj));
//              std::cout << "VMatrix("<< bi << "," << bj << ") with offset (" 
//                        << iOff << "," << jOff << ")" << std::endl;
                  VMatrix const& m((*this)(bi,bj));
                  for (unsigned i = 0; i < m.nRows(); ++i) {
                      for (unsigned j = 0; j < m.nCols(); ++j) {
                          vm->set(iOff+i, jOff+j, m(i,j));
                      }
                  }
              }
          }
       }

       // Inefficient implmentation targeted for debugging.
       unsigned rowOffset(unsigned bi) const
       {
          unsigned offset(0);
          for (unsigned i = 0; i < bi; ++i) {
              offset += (*this)(i,0).nRows();
          }
          return offset;
       }

       // Inefficient implmentation targeted for debugging.
       unsigned colOffset(unsigned bj) const
       {
          unsigned offset(0);
          for (unsigned j = 0; j < bj; ++j) {
              offset += (*this)(0,j).nCols();
          }
          return offset;
       }

       static void matrix_product(BlockMatrix& C, BlockMatrix& A, BlockMatrix& B);

       void info(const char* msg = 0) const
       {
          if (msg) {
             std::cout << msg << std::endl;
          }
          std::cout << "Number of row blocks: :    " << nRowBlocks() << std::endl;
          std::cout << "Number of column blocks:   " << nColBlocks() << std::endl;
          std::cout << "Total number of rows:      " << nRows() << std::endl;
          std::cout << "Total number of columns:   " << nCols() << std::endl;
          std::cout << "BlockMatrix is consistent: " << (consistent() ? "true" : "false") << std::endl;
          std::cout << std::endl;

          std::cout << "Tile sizes:" << std::endl;
          for (unsigned row = 0; row < m_nRowBlocks; ++row) {
              for (unsigned col = 0; col < m_nColBlocks; ++col) {
                  VMatrix const& m((*this)(row, col));
                  std::cout << "(" << m.nRows() <<","<< m.nCols() << ")  ";
              }   
              std::cout << std::endl;
          } 

          std::cout << std::endl;
          std::cout << "Data arrangement:" << std::endl;

          for (unsigned row = 0; row < m_nRowBlocks; ++row) {
              for (unsigned col = 0; col < m_nColBlocks; ++col) {
                  VMatrix::StorageT s((*this)(row, col).storage());
                  switch (s) {
                     case VMatrix::Zero:     std::cout << '.';  break;
                     case VMatrix::Diagonal: std::cout << '\\'; break;
                     case VMatrix::Banded:   std::cout << 'b';  break;
                     case VMatrix::Dense:    std::cout << 'X';  break;

                     case VMatrix::Striped: {
                        unsigned n((*this)(row, col).stripes().size());
                        std::cout << n;
                     } break;
                  }
                  std::cout << "   ";
              }   
              std::cout << std::endl;
          }
       }

       void print(const char* msg = 0) const
       {
          if (msg) {
             std::cout << msg << std::endl;
          }
          for (unsigned bi = 0; bi < m_nRowBlocks; ++bi) {
              unsigned nr = (*this)(bi,0).nRows();
              for (unsigned i = 0; i < nr; ++i) {
                  for (unsigned bj = 0; bj < m_nColBlocks; ++bj) {
                      VMatrix const& m((*this)(bi, bj));
                      unsigned nc = m.nCols();
                      for (unsigned j = 0; j < nc; ++j) {
                          if (m(i,j) == 0) {
                             std::cout << "    . ";
                          }else {
                             std::cout << std::setw(5) << m(i,j) << " ";
                          }
                      }
                      std::cout << " ";
                 }
                 std::cout << std::endl;
              }
             std::cout << std::endl;
          }
       }
       
    private:
       unsigned m_nRowBlocks;
       unsigned m_nColBlocks;
       VMatrix*  m_blocks;
};
#endif
