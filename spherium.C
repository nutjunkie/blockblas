#include "BlockMatrix.h"
#include <iostream>

int test_1();
int test_2();
int test_3();

int main(void) 
{
   unsigned nColBlocks(3);
   unsigned nRowBlocks(4);

   BlockMatrix blockMatrix(nRowBlocks, nColBlocks);

   DebugFunctor    debugFunctor;
   ZeroFunctor     zeroFunctor;
   DiagonalFunctor diagonalFunctor;


   VMatrix::StorageT storage(VMatrix::Dense);
   TestFunctor testFunctor(0,0);
   BlockMatrix bm(1,1);

   bm(0,0).init(9,9, storage, &testFunctor);
   bm(0,0).bind();
   bm.print();


   BlockMatrix bm2(2,2);

   TestFunctor fun1(0,0);
   bm2(0,0).init(4,4, storage, &fun1).bind();

   TestFunctor fun2(4,0);
   bm2(1,0).init(5,4, storage, &fun2).bind();

   TestFunctor fun3(0,4);
   bm2(0,1).init(4,5, storage, &fun3).bind();

   TestFunctor fun4(4,4);
   bm2(1,1).init(5,5, storage, &fun4).bind();

   bm2.info();
   bm2.print();

/*
   bm(0,0).init(&debugFunctor,9,8, storage);
   bm(0,0).bind();
   bm.print();

   bm(0,0).init(&debugFunctor,9,9, storage);
   bm(0,0).bind();
   std::cout << "Rows == Cols" << std::endl;
   bm.print();

   bm(0,0).init(&debugFunctor,9,7, storage);
   bm(0,0).bind();
   std::cout << "Rows > Cols" << std::endl;
   bm.print();
*/



   return 0;
}
