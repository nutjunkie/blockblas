#include "VMatrix.h"
#include "BlockMatrix.h"
#include <iostream>
#include <veclib/veclib.h>


// Accumulates the A.B product into C:
//    C += A.B 
void matrix_product(VMatrix& C, VMatrix& A, VMatrix& B)
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

   VMatrix::StorageT storageA(A.storage());   double* a(A.data());
   VMatrix::StorageT storageB(B.storage());   double* b(B.data());
   VMatrix::StorageT storageC(C.storage());   double* c(C.data());

//   std::cerr << "VMatrix multiplication for " << VMatrix::toString(storageC) << " <- " 
//             << VMatrix::toString(storageA) << " x " << VMatrix::toString(storageB) << std::endl;
             

   //                     ---------- Zero Matrices ----------
   if (storageA == VMatrix::Zero || 
       storageB == VMatrix::Zero) {
      // Nothing to do

   //                     ---------- Diagonal Matrices ----------
   }else if (storageA == VMatrix::Diagonal && 
             storageB == VMatrix::Diagonal &&
             storageC == VMatrix::Diagonal) {
     unsigned nc(C.nCols());
     for (unsigned i = 0; i < A.nCols(); ++i) {
         c[i] = a[i]*b[i];
     }
 
   }else if (storageA == VMatrix::Diagonal && 
             storageB == VMatrix::Diagonal &&
             storageC == VMatrix::Dense) {
     unsigned nc(C.nCols());
     for (unsigned i = 0; i < A.nCols(); ++i) {
         c[i*nc+i] = a[i]*b[i];
     }
      
   }else if (storageA == VMatrix::Diagonal && 
             storageB == VMatrix::Dense &&
             storageC == VMatrix::Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     for (unsigned i = 0; i < nr; ++i) {
         for (unsigned j = 0; j < nc; ++j) {
             c[i*nc+j] = a[i] * b[i*nc+j];
         }
     }

  }else if (storageA == VMatrix::Dense && 
            storageB == VMatrix::Diagonal &&
            storageC == VMatrix::Dense) {
     unsigned nc(B.nCols());
     unsigned nr(A.nRows());

     for (unsigned i = 0; i < nr; ++i) {
         for (unsigned j = 0; j < nc; ++j) {
             c[i*nc+j] = a[i*nc+j] * b[j];
         }
     }

   //                     ---------- Striped Matrices ----------
  }else if (storageA == VMatrix::Striped && 
            storageB == VMatrix::Dense   &&
            storageC == VMatrix::Dense) {

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
          int len = (offset < 0) ? std::min(rowsA + offset,colsA)
                                 : std::min(rowsA, colsA-offset);

//        std::cout << "offset: " << offset << " offC: " << offC << " offB: " << offB 
//                  << " contraction: " <<  len  << std::endl;;

          offC *= colsC;
          offB *= colsC;
          
          for (unsigned i = 0; i < len; ++i) {
              double x(a[i + s*m]);
              for (unsigned j = 0; j < colsC; ++j) {
                  c[offC+j] += x*b[offB+j];
//                std::cout << "C("<<(i+offC) << "," << j <<") = " << x << " x " 
//                   << b[(i+offB)*colsC+j] << " = " << c[(i+offC)*colsC+j] << std::endl;
              }
              offC += colsC;
              offB += colsC;
          }
      }

   }else if (storageA == VMatrix::Dense   && 
             storageB == VMatrix::Striped &&
             storageC == VMatrix::Dense) {

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

   }else if (storageA == VMatrix::Striped && 
             storageB == VMatrix::Striped &&
             storageC == VMatrix::Dense) {

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
              
              std::cout << "ColA: " << colA << " RowB: " << rowB << std::endl;
          }
      }

   //                     ---------- Dense Matrices ----------
   }else if (storageA == VMatrix::Dense && 
             storageB == VMatrix::Dense && 
             storageC == VMatrix::Dense) {
      if (C.nCols() == 1) {
         // handle the case when B and C are vectors.
         cblas_dgemv(CblasRowMajor, CblasNoTrans,
            A.nRows(), A.nCols(), 1.0, A.data(), A.nCols(),
            B.data(), 1, 1.0, C.data(), 1);
      } else {
         cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            A.nRows(), B.nCols(), A.nCols(), 1.0, A.data(), A.nCols(),
             B.data(), B.nCols(), 1.0, C.data(), C.nCols());
      }


   }else {
      std::cerr << "VMatrix multiplication not defined for " 
                << VMatrix::toString(storageC) << " <- " 
                << VMatrix::toString(storageA) << " x " 
                << VMatrix::toString(storageB) << std::endl;
   }
}


void matrix_product(BlockMatrix& C, BlockMatrix& A, BlockMatrix& B)
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


