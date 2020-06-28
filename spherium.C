#include "BlockMatrix.h"
#include <iostream>

int main(void) 
{
   unsigned nColBlocks(3);
   unsigned nRowBlocks(4);

   BlockMatrix blockMatrix(nRowBlocks, nColBlocks);

   DebugFunctor debugFunctor;
   ZeroFunctor  zeroFunctor;
   DiagonalFunctor diagonalFunctor;

/*
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           blockMatrix(row,col).init(&debugFunctor,2*(row+1), 2*(col+1)+1,Matrix::Dense);
       }
   }


   blockMatrix(0,0).init(&diagonalFunctor,2,3,Matrix::Diagonal);
   blockMatrix(0,2).init(&zeroFunctor,2,7,Matrix::Zero);

   blockMatrix(1,2).init(&debugFunctor,4,7,Matrix::Tridiagonal);

   blockMatrix(3,0).init(&zeroFunctor,8,3,Matrix::Zero);
   blockMatrix(3,2).init(&diagonalFunctor,8,7,Matrix::Diagonal);

   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           blockMatrix(row,col).bind();
       }
   }



   blockMatrix.info();
   blockMatrix.print();
*/

   Matrix::Storage storage;
   storage = Matrix::Tridiagonal;
   storage = Matrix::Pentadiagonal;
   BlockMatrix bm(1,1);

   bm(0,0).init(&debugFunctor,9,7, storage);
   bm(0,0).bind();
   bm.print();

   bm(0,0).init(&debugFunctor,9,8, storage);
   bm(0,0).bind();
   bm.print();

   bm(0,0).init(&debugFunctor,9,9, storage);
   bm(0,0).bind();
   std::cout << "Rows == Cols" << std::endl;
   bm.print();




/*
   bm(0,0).init(&debugFunctor,9,7, storage);
   bm(0,0).bind();
   std::cout << "Rows > Cols" << std::endl;
   bm.print();
*/



   return 0;
}
