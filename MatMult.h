#ifndef MATMULT_H
#define MATMULT_H

#include "Types.h"
#include <iostream>
#include "BlockMatrix.h"

// Accumulates the AxB product into C:
//    C += A x B 


template<class T, LayoutT L>
void gemm(VMatrix<T,L> const& A, VMatrix<T,L> const& B, VMatrix<T,L>& C); 

template<class T, LayoutT L>
void gbmv(VMatrix<T,L> const& A, VMatrix<T,L> const& B, VMatrix<T,L>& C, int const kl,
   int const ku, int const band, int const vec); 
   

// Accumulates the A.B product into C:
//    C += A.B 
template <class T, LayoutT L>
void matrix_product(VMatrix<T,L>& C, VMatrix<T,L> const& A, VMatrix<T,L> const& B)
{
   if (A.nCols() != B.nRows() ||
       B.nCols() != C.nCols() ||
       A.nRows() != C.nRows()) {

       std::cerr << "Barf on the matrix dimensions:" << std::endl;
       std::cerr << A.nCols() << " != " << B.nRows() << " || " << std::endl;
       std::cerr << B.nCols() << " != " << C.nCols() << " || " << std::endl;
       std::cerr << A.nRows() << " != " << C.nRows() << std::endl;
   }

   if (!A.isBound() || !B.isBound() || !C.isBound()) {
      std::cerr << "Unbound matrix encountered in VMatrix::matrix_product" << std::endl;
   }

   StorageT storageA(A.storage());   T const* a(A.data());
   StorageT storageB(B.storage());   T const* b(B.data());
   StorageT storageC(C.storage());   T* c(C.data());

//   std::cerr << "VMatrix multiplication for " << VMatrix::toString(storageC) << " <- " 
//             << VMatrix::toString(storageA) << " x " << VMatrix::toString(storageB) << std::endl;
             

   //                     ---------- Zero Matrices ----------
   if (storageA == Zero || 
       storageB == Zero) {
      // Nothing to do

   //                     ---------- Diagonal Matrices ----------
   }else if (storageA == Diagonal && 
             storageB == Diagonal &&
             storageC == Diagonal) {
     unsigned nc(C.nCols());
     for (unsigned i = 0; i < A.nCols(); ++i) {
         c[i] += a[i]*b[i];
     }
 
   }else if (storageA == Diagonal && 
             storageB == Diagonal &&
             storageC == Dense) {
     unsigned nc(C.nCols());
     for (unsigned i = 0; i < A.nCols(); ++i) {
         c[i*nc+i] += a[i]*b[i];
     }
      
   }else if (storageA == Diagonal && 
             storageB == Dense &&
             storageC == Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     for (unsigned i = 0; i < nr; ++i) {
         for (unsigned j = 0; j < nc; ++j) {
             c[i*nc+j] += a[i] * b[i*nc+j];
         }
     }

  }else if (storageA == Dense && 
            storageB == Diagonal &&
            storageC == Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     for (unsigned i = 0; i < nr; ++i) {
         for (unsigned j = 0; j < nc; ++j) {
             c[i*nc+j] += a[i*nc+j] * b[j];
         }
     }

   }else if (storageA == Banded && 
             storageB == Dense  &&
             storageC == Dense) {

      unsigned nvec(B.nCols());
      std::vector<int> stripes(A.stripes());
      int kl(stripes[0]);
      int ku(stripes[1]);
      int band(kl+ku+1);

      for (unsigned vec = 0; vec < nvec; ++vec) {
          gbmv(A, B, C, kl, ku, band, vec);
      }

  //  ---------- Striped Matrices ----------
  }else if (storageA == Striped && 
            storageB == Dense   &&
            storageC == Dense) {

      unsigned rowsC(C.nRows());
      unsigned colsC(C.nCols());

      std::vector<int> const& stripes(A.stripes());
      unsigned rowsA(A.nRows());
      unsigned colsA(A.nCols());
      unsigned m(std::min(rowsA,colsA));

      for (unsigned s = 0; s < stripes.size(); ++s) {
          int offset(stripes[s]);
          int offC(std::max(0,-offset));
          int offB(std::max(0, offset));
          // Contraction length
          int len = (offset < 0) ? std::min(rowsA+offset, colsA)
                                 : std::min(rowsA, colsA-offset);

//        std::cout << "offset: " << offset << " offC: " << offC << " offB: " << offB 
//                  << " contraction: " <<  len  << std::endl;;

#pragma omp parallel for
          for (unsigned i = 0; i < len; ++i) {
              T        a0(a[i + s*m]);
              T const* b0(&b[(offB+i)*colsC]);
              T*       c0(&c[(offC+i)*colsC]);
//            cblas_daxpy(colsC, a0, b0, 1, c0, 1);
              for (unsigned j = 0; j < colsC; ++j) {
                  c0[j] += a0*b0[j];
              }
          }
      }

   }else if (storageA == Dense   && 
             storageB == Striped &&
             storageC == Dense) {

      unsigned rowsC(C.nRows());
      unsigned colsC(C.nCols());

      std::vector<int> const& stripes(B.stripes());
      unsigned rowsB(B.nRows());
      unsigned colsB(B.nCols());
      unsigned m(std::min(rowsB,colsB));
  
      for (unsigned s = 0; s < stripes.size(); ++s) {
          int offset(stripes[s]);
          int offC(std::max(0,offset));
          int offA(std::min(0,offset));
          // Contraction length
          int len = (offset < 0) ? std::min(rowsB + offset,colsB)
                                 : std::min(rowsB, colsB-offset);

          for (unsigned i = 0; i < rowsC; ++i) {
              for (unsigned j = 0; j < len; ++j) {
                  c[i*colsC+j+offC] += a[i*rowsB+j-offA] * b[s*m+j];
              }
          }
      }


   }else if (storageA == Striped && 
             storageB == Striped &&
             storageC == Dense) {

      std::vector<int> const& stripesA(A.stripes());
      std::vector<int> const& stripesB(B.stripes());

      unsigned ma(std::min(A.nRows(),A.nCols()));
      unsigned mb(std::min(B.nRows(),B.nCols()));

      for (unsigned sa = 0; sa < stripesA.size(); ++sa) {
          for (unsigned sb = 0; sb < stripesB.size(); ++sb) {
              // column of A needs to match row of B
              // element goes into row of A and col of B
              int colA = std::max(0, stripesA[sa]);
              int rowB = std::max(0,-stripesA[sa]);
              if (colA <= rowB) {
                 colA += rowB;
                 //c[] += a[sa*ma + j] * b[sb*mb + j]
              }
              
              std::cout << "NYI" << std::endl;
              std::cout << "ColA: " << colA << " RowB: " << rowB << std::endl;
          }
      }

   //                     ---------- Dense Matrices ----------
   }else if (storageA == Dense && 
             storageB == Dense && 
             storageC == Dense) {
      gemm(A,B,C);

   }else {
      std::cerr << "VMatrix multiplication not defined for " 
                << toString(storageC) << " <- " 
                << toString(storageA) << " x " 
                << toString(storageB) << std::endl;
   }
}


template <class T, LayoutT L>
void matrix_product(BlockMatrix<T,L>& C, BlockMatrix<T,L> const& A, BlockMatrix<T,L> const& B)
{
   // Should check matrix dimensions
   for (unsigned bi = 0; bi < A.nRowBlocks(); ++bi) {
       for (unsigned bj = 0; bj < B.nColBlocks(); ++bj) {
           for (unsigned k = 0; k < A.nColBlocks(); ++k) {
//           std::cout << "Multiplying block: C(" << bi << "," << bj << ") <- A(" 
//                     << bi << "," << k << ") x B(" << k << "," << bj << ")" << std::endl;
               matrix_product(C(bi,bj), A(bi,k), B(k,bj));
           }
       }
   }
}


template <class T, LayoutT L>
void matrix_product_sans_diagonal(BlockMatrix<T,L>& C, BlockMatrix<T,L> const& A, 
   BlockMatrix<T,L> const& B)
{
   // Should check matrix dimensions
   for (unsigned bi = 0; bi < A.nRowBlocks(); ++bi) {
       for (unsigned bj = 0; bj < B.nColBlocks(); ++bj) {
           for (unsigned k = 0; k < A.nColBlocks(); ++k) {
               if (bi != k ) {
                  matrix_product(C(bi,bj), A(bi,k), B(k,bj));
              }
           }
       }
   }
}

#endif
