#include "BlockMatrix.h"
#include "MatMult.h"
#include "JacobiSolver.h"
#include "Timer.h"
#include "util.h"

#define COLUMNMAJOR

#ifdef COLUMNMAJOR
#define LAYOUT ColumnMajor
#else
#define LAYOUT RowMajor
#endif


ZeroFunctor<double>      zeroFunctor;
DebugFunctor             debugFunctor;
TestFunctor              testFunctor;
DiagonalFunctor<double>  diagonalFunctor;


// Tests the construction and debug printing of the BlockMatrix<double> class for
// dense matrices
int test_1()
{
   print_header(1, "Ctor and debug printing");

   unsigned nRowBlocks(4);
   unsigned nColBlocks(3);

   BlockMatrix<double, LAYOUT> bm(nRowBlocks,nColBlocks);
  
   // Initialize the structure of our BlockMatrix<double>.  Note this does not
   // allocate any memory.
   for (unsigned row = 0; row < nRowBlocks; ++row) {
       for (unsigned col = 0; col < nColBlocks; ++col) {
           // These are the dimensions of the current tile
           unsigned nRows(2*(row+1));
           unsigned nCols(2*(col+1)+1);
           bm(row,col).init(nRows,nCols,Dense);
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

   BlockMatrix<double, LAYOUT> bm(nRowBlocks,nColBlocks);

   // Initialize the structure of our BlockMatrix<double>.  Note in this case
   // we bind the data at the same time.
   bm(0,0).init(2,3,Diagonal).bind(diagonalFunctor);
   bm(0,1).init(2,7,Zero    ).bind(zeroFunctor);
   bm(0,2).init(2,5,Dense   ).bind(debugFunctor);

   bm(1,0).init(6,3,Diagonal).bind(diagonalFunctor);
   bm(1,1).init(6,7,Diagonal).bind(debugFunctor);
   bm(1,2).init(6,5,Diagonal).bind(diagonalFunctor);

   bm(2,0).init(2,3,Diagonal).bind(diagonalFunctor);
   bm(2,1).init(2,7,Zero    ).bind(zeroFunctor);
   bm(2,2).init(2,5,Diagonal).bind(diagonalFunctor);

   std::vector<int> stripes{-4,-2,-1,1};
   std::vector<int> stripes2{-1,1};

   bm(3,0).init(9,3,Zero    ).bind(debugFunctor);
   bm(3,1).init(9,7,stripes          ).bind(debugFunctor);
   bm(3,2).init(9,5,stripes2         ).bind(debugFunctor);

   bm.info();
   bm.print();

   BlockMatrix<double, LAYOUT> bm2(bm);
   bm2.info();
   bm2.print();

   return 0;
}


int test_3()
{
   print_header(3, "Dense <- Dense x Dense");
#ifdef COLUMNMAJOR
   std::cout << "Skipping test, ColumnMajor implementation NYI" << std::endl;
   return 0;
#else
   VMatrix<double, LAYOUT> a, b, c;

   a.init( 9,10).bind(testFunctor);
   b.init(10, 8).bind(testFunctor);
   c.init( 9, 8).bind(zeroFunctor);
   matrix_product(c, a, b); 

   return matrix_residue(c, "test_3_rm.dat");
#endif
}


// Test the block matrix multiplication C = A.A
int test_4()
{
   print_header(4, "Block matrix multiplication");

   BlockMatrix<double, LAYOUT> A(3,3);
   BlockMatrix<double, LAYOUT> C(3,3);

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

   VMatrix<double, LAYOUT> a, b, c;
   a.init(12,12, Dense).bind(testFunctor);
   c.init(12,12, Dense).bind(zeroFunctor);

   // Compute the matrix product all as one
   matrix_product(c,a,a);

   C.toDense(&b);

   return matrix_residue(b,c);
}


//Timing test
int test_5(unsigned n)
{
   print_header(5, "BlockMatrix<double, LAYOUT> timing test");

   const unsigned dim(512);
   const unsigned blocks(8);

   BlockMatrix<double, LAYOUT> A(blocks,blocks);
   BlockMatrix<double, LAYOUT> C(blocks,blocks);

   for (unsigned bi = 0; bi < blocks; ++bi) {
       for (unsigned bj = 0; bj < blocks; ++bj) {
           if (zeroTest(bi,bj)) {
              A(bi,bj).init(dim,dim, Dense).bind(testFunctor);
           }else {
              A(bi,bj).init(dim,dim, Zero ).bind(testFunctor);
           }
           C(bi,bj).init(dim,dim, Dense).bind(zeroFunctor);
       }
   }

   A.info("A Matrix:");
   
   VMatrix<double, LAYOUT> a, c;
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

   VMatrix<double, LAYOUT> b;
   C.toDense(&b);

   return matrix_residue(b,c);
}



int test_6()
{
   print_header(6, "dense <- diag x dense");
#ifdef COLUMNMAJOR 
   std::cout << "Skipping test, ColumnMajor implementation NYI" << std::endl;
   return 0;
#else
   VMatrix<double, LAYOUT> a, b, c;

   a.init(10,10, Diagonal).bind(testFunctor);
   b.init(10,10, Dense   ).bind(testFunctor);
   c.init(10,10, Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

   return matrix_residue(c,"test_6.dat");
#endif
}


int test_7()
{
   print_header(7, "dense <- dense x diag");
#ifdef COLUMNMAJOR 
   std::cout << "Skipping test, ColumnMajor implementation NYI" << std::endl;
   return 0;
#else
   VMatrix<double, LAYOUT> a, b, c;

   a.init(10,10, Dense   ).bind(testFunctor);
   b.init(10,10, Diagonal).bind(testFunctor);
   c.init(10,10, Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 

   return matrix_residue(c,"test_7.dat");
#endif
}


int test_8()
{
   print_header(8, "dense <- diag x diag");
#ifdef COLUMNMAJOR 
   std::cout << "Skipping test, ColumnMajor implementation NYI" << std::endl;
   return 0;
#else
   VMatrix<double, LAYOUT> a, b, c;

   a.init(10,10, Diagonal).bind(testFunctor);
   b.init(10,10, Diagonal).bind(testFunctor);
   c.init(10,10, Dense   ).bind(zeroFunctor);

   matrix_product(c, a, b); 
   c.print();

   return matrix_residue(c,"test_8.dat");
#endif
}


int test_9()
{
   print_header(9, "dense <- striped x dense");

   VMatrix<double, LAYOUT> a, b, c, d;
   std::vector<int> stripes = {-8,-2,2,4};

   a.init(15,16, stripes).bind(debugFunctor);
   b.init(16,10, Dense).bind(debugFunctor);
   c.init(15,10, Dense).bind(zeroFunctor);
   d.init(15,10, Dense).bind(zeroFunctor);

   matrix_product(c, a, b); 

   c.print("Product of striped x dense");

   // Do a dense version for checking
   a.toDense();
   a.info();
   a.print("A with stripes");
   b.print("B is dense");
   matrix_product(d, a, b); 
   d.print("Product dense ");

   return  matrix_residue(c,d);
}


int test_10()
{
   print_header(10, "dense <- dense x striped");

   VMatrix<double, LAYOUT> a, b, c, d;
   std::vector<int> stripes = {-8,-2,2,4};

   a.init(15,16, Dense).bind(debugFunctor);
   b.init(16,10, stripes).bind(debugFunctor);
   c.init(15,10, Dense).bind(zeroFunctor);
   d.init(15,10, Dense).bind(zeroFunctor);

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
   
   VMatrix<double, LAYOUT> a, b, c, d;
   std::vector<int> stripesA = {-2};
   std::vector<int> stripesB = {-2};

   a.init(6,7, stripesA).bind(debugFunctor);
   b.init(7,8, stripesB).bind(debugFunctor);
   c.init(6,8, Dense).bind(zeroFunctor);
   d.init(6,8, Dense).bind(zeroFunctor);

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

   VMatrix<double, LAYOUT> a, b, c, d;
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
   print_header(13, "BlockMatrix::bind");

   VMatrix<double, ColumnMajor> a;
   
   a.init(12,16, Dense).bind(debugFunctor);
   a.print("VMatrix a");

   double* data(a.data());

   unsigned nRowBlocks(2);
   unsigned nColBlocks(4);
   BlockMatrix<double, ColumnMajor> A(nRowBlocks,nColBlocks);

   for (unsigned bi = 0; bi < nRowBlocks; ++bi) {
       for (unsigned bj = 0; bj < nColBlocks; ++bj) {
           A(bi,bj).init(6,4);
       }
   }

   A.bind(data);
   A.print("BlockMatrix A");

   double junk[200];

   //weird bug, if this is set to data it don't work
   //double* d = data;
   double* d = junk;

   A.unbind(d);

   for (int i = 0; i < 20; ++i ){
       std::cout << d[i] << std::endl;
   }

   a.init(12,16).bind(d);
   a.print("VMatrix a");

   a.unbind(d);
   for (int i = 0; i < 20; ++i ){
       std::cout << d[i] << std::endl;
   }



   return 0;
}

int test_14()
{
   print_header(14, "Matrix inversion");

   VMatrix<double, LAYOUT> a, c;
   a.init(12,12).bind(testFunctor);
   c.init(12,12).bind(zeroFunctor);

   VMatrix<double, LAYOUT> b(a);

   a.invert();
   matrix_product(c, a, b); 

   c.print("Should be identity:");

   return 0;
}

int test_16()
{
   print_header(14, "Complex Matrix inversion");

   VMatrix<complex, LAYOUT> a, c;
   a.init(12,12).bind(ComplexTestFunctor());
   c.init(12,12).bind(ZeroFunctor<complex>());

   VMatrix<complex, LAYOUT> b(a);
   a.print("Original A matrix");

   a.invert();
   a.print("Inverse A matrix");
   matrix_product(c, a, b); 

   c.print("Should be identity:");

   return 0;
}

int test_15()
{
   print_header(15, "Jacobi Solver");

   unsigned nBlocks(5);
   unsigned blockDim(5);

   BlockMatrix<double, LAYOUT> A(nBlocks,nBlocks);
   BlockMatrix<double, LAYOUT> b(nBlocks,1);
   BlockMatrix<double, LAYOUT> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).init(blockDim,1,Dense);
       x(row,0).init(blockDim,1,Dense);
       for (unsigned col = 0; col < nBlocks; ++col) {
           // These are the dimensions of the current tile
           A(row,col).init(blockDim,blockDim,Dense);
       }
   }


   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).bind(debugFunctor);
       x(row,0).bind(zeroFunctor);
       for (unsigned col = 0; col < nBlocks; ++col) {
           A(row,col).bind(TestFunctor(A.rowOffset(row),A.rowOffset(col)));
       }
   }

   jacobi_solver(x,A,b);

   return 0;
}

int test_17()
{
   print_header(17, "Complex Jacobi Solver");

   unsigned nBlocks(5);
   unsigned blockDim(5);

   BlockMatrix<complex, LAYOUT> A(nBlocks,nBlocks);
   BlockMatrix<complex, LAYOUT> b(nBlocks,1);
   BlockMatrix<complex, LAYOUT> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).init(blockDim,1,Dense);
       x(row,0).init(blockDim,1,Dense);
       for (unsigned col = 0; col < nBlocks; ++col) {
           // These are the dimensions of the current tile
           A(row,col).init(blockDim,blockDim,Dense);
       }
   }


   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).bind(ComplexDebugFunctor());
       x(row,0).bind(ZeroFunctor<complex>());
       for (unsigned col = 0; col < nBlocks; ++col) {
           A(row,col).bind(ComplexTestFunctor(A.rowOffset(row),A.rowOffset(col)));
       }
   }

   jacobi_solver(x,A,b);

   return 0;
}

int test_18()
{
   print_header(18, "VMatrix -> Block");

   VMatrix<double, LAYOUT> a;
   std::vector<int> stripes{-3,-2,-1,0,1,2,3};
   a.init(11,11, stripes).bind(StencilFunctor());
   a.print("VMatrix a - sparse");

   VMatrix<complex, LAYOUT> z;
   z.fromDouble(a);
   z.info();
   z.print("VMatrix z - sparse/complex");
   

   z.toDense();
   z.info();
   z.print("VMatrix z - complex");

   BlockMatrix<complex, LAYOUT> Z(z);
   
   Z.info();
   Z.print("BlockMatrix Z");

   return 0;
}


int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok(0);
   test_13();
   test_15();
   test_17();

   return 0;
   ok = ok 
      + test_1()
      + test_2()
      + test_3()
      + test_4()
      + test_6()
      + test_7()
      + test_8()
      + test_9()
      + test_10()
      + test_11()
      + test_12()
      + test_13()
      + test_14()
      + test_15()
      + test_16()
      + test_17()
      + test_18()
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
