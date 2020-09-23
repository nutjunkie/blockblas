#include "BlockMatrix.h"
#include "MatMult.h"
#include "Timer.h"
#include "util.h"


ZeroFunctor     zeroFunctor;
DebugFunctor    debugFunctor;
TestFunctor     testFunctor;
DiagonalFunctor diagonalFunctor;


// Tests the construction and debug printing of the BlockMatrix class for
// dense matrices
int test_1()
{
   print_header(1, "Ctor and debug printing");

   unsigned nRowBlocks(4);
   unsigned nColBlocks(3);

   BlockMatrix bm(nRowBlocks,nColBlocks);
  
   // Initialize the structure of our BlockMatrix.  Note this does not
   // allocate any memory.
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           // These are the dimensions of the current tile
           unsigned nRows(2*(row+1));
           unsigned nCols(2*(col+1)+1);
           bm(row,col).init(nRows,nCols,VMatrix::Dense);
       }
   }

   // Bind our data, this allocates the memory and sets the entries
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           bm(row,col).bind(debugFunctor);
       }
   }

   bm.info();
   bm.print();

   return 0;
}


// Tests the different storage formats for the tiles.
int test_2()
{
   print_header(2, "Block formats");

   unsigned nRowBlocks(4);
   unsigned nColBlocks(3);

   BlockMatrix bm(nRowBlocks,nColBlocks);

   // Initialize the structure of our BlockMatrix.  Note in this case
   // we bind the data at the same time.
   bm(0,0).init(2,3,VMatrix::Diagonal).bind(diagonalFunctor);
   bm(0,1).init(2,7,VMatrix::Zero    ).bind(zeroFunctor);
   bm(0,2).init(2,5,VMatrix::Dense   ).bind(debugFunctor);

   bm(1,0).init(6,3,VMatrix::Diagonal).bind(diagonalFunctor);
   bm(1,1).init(6,7,VMatrix::Diagonal).bind(debugFunctor);
   bm(1,2).init(6,5,VMatrix::Diagonal).bind(diagonalFunctor);

   bm(2,0).init(2,3,VMatrix::Diagonal).bind(diagonalFunctor);
   bm(2,1).init(2,7,VMatrix::Zero    ).bind(zeroFunctor);
   bm(2,2).init(2,5,VMatrix::Diagonal).bind(diagonalFunctor);

   std::vector<int> stripes{-4,-2,-1,1};
   std::vector<int> stripes2{-1,1};

   bm(3,0).init(9,3,VMatrix::Zero    ).bind(debugFunctor);
   bm(3,1).init(9,7,stripes          ).bind(debugFunctor);
   bm(3,2).init(9,5,stripes2         ).bind(debugFunctor);

   bm.info();
   bm.print();

   return 0;
}


int test_3()
{
   print_header(3, "Dense <- Dense x Dense");

   VMatrix a, b, c;

   a.init( 9,10).bind(testFunctor);
   b.init(10, 8).bind(testFunctor);
   c.init( 9, 8).bind(zeroFunctor);
   matrix_product(c, a, b); 

   return matrix_residue(c, "test_3.dat");
}


// Test the block matrix multiplication C = A.A
int test_4()
{
   print_header(4, "Block matrix multiplication");

   BlockMatrix A(3,3);
   BlockMatrix C(3,3);

   for (unsigned i = 0; i < 3; ++i) {
       for (unsigned j = 0; j < 3; ++j) {
           A(i,j).init(4,4);
           C(i,j).init(4,4).bind(zeroFunctor);
       }
   }

   // Need to init before binding A to use the row/colOffset functions
   for (unsigned i = 0; i < 3; ++i) {
       for (unsigned j = 0; j < 3; ++j) {
           A(i,j).bind(TestFunctor(A.rowOffset(i),A.rowOffset(j)));
       }
   }

   // Compute the matrix product via blocks
   matrix_product(C, A, A);

   VMatrix a, b, c;
   a.init(12,12, VMatrix::Dense).bind(testFunctor);
   c.init(12,12, VMatrix::Dense).bind(zeroFunctor);

   // Compute the matrix product all as one
   matrix_product(c,a,a);

   C.toDense(&b);

   return matrix_residue(b,c);
}


//Timing test
int test_5(unsigned n)
{
   print_header(5, "BlockMatrix timing test");

   const unsigned dim(512);
   const unsigned blocks(8);

   BlockMatrix A(blocks,blocks);
   BlockMatrix C(blocks,blocks);

   for (unsigned bi = 0; bi < blocks; ++bi) {
       for (unsigned bj = 0; bj < blocks; ++bj) {
           if (zeroTest(bi,bj)) {
              A(bi,bj).init(dim,dim, VMatrix::Dense).bind(testFunctor);
           }else {
              A(bi,bj).init(dim,dim, VMatrix::Zero ).bind(testFunctor);
           }
           C(bi,bj).init(dim,dim, VMatrix::Dense).bind(zeroFunctor);
       }
   }

   A.info("A Matrix:");
   
   VMatrix a, c;
   C.toDense(&c);
   A.toDense(&a);

   Timer timer;

   std::cout << "Performing " << n << "iterations for time-averaging" << std::endl;
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(C,A,A);
   }
   double elapsed(timer.stop());
   std::cout << " Average BlockMatrix time: " << timer.format(elapsed/n) << std::endl;

   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(c,a,a);
   }
   elapsed = timer.stop();
   std::cout << " Average VMatrix time:     " << timer.format(elapsed/n) << std::endl;

   std::cout << "-----------------" << std::endl;

   VMatrix b;
   C.toDense(&b);

   return matrix_residue(b,c);
}



int test_6()
{
   print_header(6, "dense <- diag x dense");

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   b.init(10,10, VMatrix::Dense   ).bind(testFunctor);
   c.init(10,10, VMatrix::Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

   return matrix_residue(c,"test_6.dat");
}


int test_7()
{
   print_header(7, "dense <- dense x diag");

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Dense   ).bind(testFunctor);
   b.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   c.init(10,10, VMatrix::Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

   return matrix_residue(c,"test_7.dat");
}


int test_8()
{
   print_header(8, "dense <- diag x diag");

   VMatrix a, b, c;

   a.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   b.init(10,10, VMatrix::Diagonal).bind(testFunctor);
   c.init(10,10, VMatrix::Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 
   c.print();

   return matrix_residue(c,"test_8.dat");
}


int test_9()
{
   print_header(9, "dense <- striped x dense");

   VMatrix a, b, c, d;
   std::vector<int> stripes = {-8,-2,2,4};

   a.init(15,16, stripes).bind(debugFunctor);
   b.init(16,10, VMatrix::Dense).bind(debugFunctor);
   c.init(15,10, VMatrix::Dense).bind(zeroFunctor);
   d.init(15,10, VMatrix::Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 

   //a.print("A with stripes");
   //b.print("B is dense");
   //c.print("Product of striped x dense");

   // Do a dense version for checking
   a.toDense();
   matrix_product(d, a, b); 

   return  matrix_residue(c,d);
}


int test_10()
{
   print_header(10, "dense <- dense x striped");

   VMatrix a, b, c, d;
   std::vector<int> stripes = {-8,-2,2,4};

   a.init(15,16, VMatrix::Dense).bind(debugFunctor);
   b.init(16,10, stripes).bind(debugFunctor);
   c.init(15,10, VMatrix::Dense).bind(zeroFunctor);
   d.init(15,10, VMatrix::Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 
   //a.print();
   //b.print();
   //c.print();

   // Do a dense version for checking
   b.toDense();
   matrix_product(d, a, b); 

   return matrix_residue(c,d);
}


int test_11()
{
   print_header(11, "dense <- striped x striped");
   
   VMatrix a, b, c, d;
   std::vector<int> stripesA = {-2};
   std::vector<int> stripesB = {-2};

   a.init(6,7, stripesA).bind(debugFunctor);
   b.init(7,8, stripesB).bind(debugFunctor);
   c.init(6,8, VMatrix::Dense).bind(zeroFunctor);
   d.init(6,8, VMatrix::Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 
   a.print();
   b.print();
   c.print();

   // Do a dense version for checking
   a.toDense();
   b.toDense();
   matrix_product(d, a, b); 

   return matrix_residue(d,c);
}


int test_12()
{
   print_header(12, "Banded matrices");

   VMatrix a, b, c, d;
   unsigned dim(5);
   unsigned nvec(3);

   a.init(dim,dim,1,1).bind(debugFunctor);

   b.init(dim,nvec).bind(debugFunctor);
   c.init(dim,nvec).bind(zeroFunctor);
   d.init(dim,nvec).bind(zeroFunctor);

   a.print("banded matrix print:");
   b.print("vector print:");

   matrix_product(c, a, b); 
   c.print("product print:");
   a.toDense();
   matrix_product(d, a, b); 

   return matrix_residue(d,c);
}

int test_13()
{
   print_header(13, "Colum-major matrices");

   VMatrix a, b, c, d;
   unsigned dim(5);
   unsigned nvec(3);

   a.init(dim,dim,1,1).bindCM(debugFunctor);
   b.init(dim,nvec).bindCM(debugFunctor);

   c.init(dim,nvec).bindCM(zeroFunctor);
   d.init(dim,nvec).bindCM(zeroFunctor);

   a.print("banded matrix print:");
   b.print("vector print:");

   matrix_product(c, a, b); 
   c.print("product print:");
   a.toDense();
   matrix_product(d, a, b); 

   return matrix_residue(d,c);
}

int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok(0);
   ok = ok 
/*
*/
      + test_1()
      + test_2()
      + test_3()
      + test_4()
      + test_6()
      + test_7()
      + test_8()
      + test_9()
      + test_10()
//    + test_11()
      + test_12()
      + test_13()
      ;

    //test_5(5);

   std::cout << std::endl;
   std::cout << "Test run finished: ";

   if (ok) {
      std::cout << "FAIL" << std::endl;
      std::cout << "With " << ok << " failure(s)" << std::endl;
   }else {
      std::cout << "ALL PASS" << std::endl;
   }


   return ok;
}
