#include "BlockMatrix.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
/******************************************************************************
 * 
 *  Class for managing block matrices.
 * 
 *****************************************************************************/


void BlockMatrix::info() const
{
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

   char const* chars = ".\\TPXS";

   for (unsigned row = 0; row < m_nRowBlocks; ++row) {
       for (unsigned col = 0; col < m_nColBlocks; ++col) {
           VMatrix::StorageT s((*this)(row, col).storage());
           std::cout << chars[s] << "  ";
       }   
       std::cout << std::endl;
   }
}


unsigned BlockMatrix::nRows() const
{
   if (m_nRowBlocks == 0 || m_nColBlocks == 0) return 0;

   unsigned n(0);
   for (unsigned i = 0; i < m_nRowBlocks; ++i) {
       n += (*this)(i,0).nRows();
   }
   return n;
}


unsigned BlockMatrix::nCols() const
{
   if (m_nRowBlocks == 0 || m_nColBlocks == 0) return 0;

   unsigned n(0);
   for (unsigned i = 0; i < m_nColBlocks; ++i) {
       n += (*this)(0,i).nCols();
   }
   return n;
}


bool BlockMatrix::consistent() const
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


void BlockMatrix::print(const char* msg) const
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


// Inefficient implmentation targeted for debugging.
unsigned BlockMatrix::rowOffset(unsigned bi) const
{
   unsigned offset(0);
   for (unsigned i = 0; i < bi; ++i) {
       offset += (*this)(i,0).nRows();
   }
   return offset;
}


// Inefficient implmentation targeted for debugging.
unsigned BlockMatrix::colOffset(unsigned bj) const
{
   unsigned offset(0);
   for (unsigned j = 0; j < bj; ++j) {
       offset += (*this)(0,j).nCols();
   }
   return offset;
}


void BlockMatrix::toDense(VMatrix* vm) const
{
   vm->init(nRows(), nCols(), VMatrix::Dense).bind();

   for (unsigned bi = 0; bi < nRowBlocks(); ++bi) {
       unsigned iOff(rowOffset(bi));
       for (unsigned bj = 0; bj < nColBlocks(); ++bj) {
           unsigned jOff(colOffset(bj));
//       std::cout << "VMatrix("<< bi << "," << bj << ") with offset (" 
//                 << iOff << "," << jOff << ")" << std::endl;
           VMatrix const& m((*this)(bi,bj));
           for (unsigned i = 0; i < m.nRows(); ++i) {
               for (unsigned j = 0; j < m.nCols(); ++j) {
                   vm->set(iOff+i, jOff+j, m(i,j));
               }
           }
       }
   }
}


void BlockMatrix::matrix_product(BlockMatrix& C, BlockMatrix& A, BlockMatrix& B)
{
   // Should check matrix dimensions
   for (unsigned bi = 0; bi < A.nRowBlocks(); ++bi) {
       for (unsigned bj = 0; bj < B.nColBlocks(); ++bj) {
           for (unsigned k = 0; k < A.nColBlocks(); ++k) {
//           std::cout << "Multiplying block: C(" << bi << "," << bj << ") <- A(" 
//                     << bi << "," << k << ") x B(" << k << "," << bj << ")" << std::endl;
               VMatrix::matrix_product(C(bi,bj), A(bi,k), B(k,bj));
           }
       }
   }
}
