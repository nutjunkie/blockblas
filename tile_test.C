#include "ZeroTile.h"
#include "DiagonalTile.h"
#include "StripedTile.h"
#include "CMTile.h"
#include "Functor.h"
#include "JacobiSolver.h"

#include "TileArray.h"



void print_header(unsigned n, char const* header)
{
   std::string s(" test_");
   s += std::to_string(n) + ": " + std::string(header);

   unsigned len(s.length());
   std::cout << std::string(len+1, '=') << std::endl;
   std::cout << s << std::endl;
   std::cout << std::string(len+1, '=') << std::endl;
   
}



int test_1()
{
   print_header(1, "ZeroTile");

   ZeroTile<double> td(10,10);
   td.alloc(); 
   td.info();
   td.print();

   ZeroTile<complex> tc(9,8);
   tc.alloc(); 
   tc.info();
   tc.print();

   return 0;
}


int test_2()
{
   print_header(2, "DiagonalTile");

   DebugFunctor f;
   
   DiagonalTile<double> t(5,6);
   t.info();
   t.print();

   t.alloc();
   t.fill(f);
   t.print("After binding");
   t.info();

   return 0;
}



int test_3()
{
   print_header(3, "StripedTile");

   DebugFunctor df;
   StencilFunctor sf;

   std::vector<int> stripes = {-1,1};
   
   StripedTile<double> t(5,6,stripes);
   t.info();
   t.print();

   t.alloc();
   t.fill(df);
   t.print("After binding");
   t.info();

   return 0;
}


int test_4()
{
   print_header(4, "CMTile");

   DebugFunctor df;
   StencilFunctor sf;

/*
   CMTile<double> t(5,6);
   t.alloc();
   t.fill(df);

   {
      CMTile<double> s(5,6);
      s.bind(t.data());
      s.info();
      s.set(1,1,100.0);
      s.print("S before going out of scope");
      double* t = s.release();
      std::cout << "release = " << t << std::endl;
      s.dealloc();
   }
   t.print("After binding");
*/

   double m[25];

   CMTile<double> v(5,5);
   CMTile<double> u(3,3);
   v.bind(m);   
   v.fill();

   u.bind(m,5);   
   u.info();
   u.fill(df);
   u.print("u from array");

   v.info();   
   v.print("v from array");


   for (int i = 0; i < 25; ++i) {
       std::cout << "m[" << i << "] = " << m[i] << std::endl;
   }

   return 0;
}


int test_5()
{
   print_header(5, "CMTile matrix multiplication");

   const size_t n(80);

   double a[n], b[n], c[n];
   memset(a,0,80*sizeof(double));
   memset(b,0,80*sizeof(double));
   memset(c,0,80*sizeof(double));

   CMTile<double> A(10,6), B(6,8), C(10,8);

   DebugFunctor df;

   A.bind(a);
   A.fill(df);

   B.bind(b,10);  // give b a different leading dimension
   B.fill(df);

   C.bind(c);

   A.print("Matrix A");
   B.print("Matrix B");

   tile_product(A,B,0.0,C);
   C.print("Product Matrix C");

   for (int i = 0; i < n; ++i) {
       std::cout << "bufs[" << i << "]:  " << a[i] << "  " << b[i] << "  " << c[i] << std::endl;
   }

}


int test_6()
{
   print_header(6, "Striped-CMTile matrix multiplication");

   std::vector<int> stripes{-3,-1,0,1,3};
   StripedTile<double> A(10,6,stripes);
   A.fill(StencilFunctor());

   CMTile<double> B(6,8), C(10,8);
   B.fill(DebugFunctor());
   C.fill();

   A.print("Matrix A");
   CMTile<complex> Z;
   Z.info();
   Z.from(A);
   Z.info();
   Z.print("Converted A -> Z");

   B.print("Matrix B");

   tile_product(A,B,0.0,C);
   C.print("Product Matrix C");

   CMTile<double> D(10,6);

   D.from(A);
   tile_product(D,B,0.0,C);
   C.print("Product Matrix C A->D");
}


int test_7()
{
   print_header(7,"CMTile += DiagonalTile");

   CMTile<double> A(8,12);
   double a[120];
   A.bind(a,10);
   A.fill(DebugFunctor());

   DiagonalTile<double> D(8,12);
   D.fill(DebugFunctor());

   A += D;
   A.print("After adding diagonal");

   CMTile<double> C(A);
   A += C;
   A.print("After adding C");

   A.scale(-1.0);

   CMTile<double> B(10,12);
   B.bind(a);
   B.print("Buffer contents ");
}


int test_8()
{
   print_header(8,"TileArray");

   unsigned nRows(5);
   unsigned nCols(4);

   unsigned nR(3);
   unsigned nC(2);

   TileArray<double> TA(nRows,nCols);

   for (unsigned bi = 0; bi < nRows; ++bi) {
       for (unsigned bj = 0; bj < nCols; ++bj) {
           TA.set(bi,bj, new CMTile<double>(nR,nC));
           TA(bi,bj).fill(DebugFunctor());
       }
   }

   TA.info("TileArray print");
   TA.print("TileArray print");

   double m[200];
   for (unsigned i = 0; i < 200; ++i) {
       m[i] = 1.0*(i+1);
   }

   TA.bind(m);
   TA.print("bound TileArray from data");

   TileArray<double> TB(TA);
   TB.print("TileArray from copy");

   TA += TB;

   TA.print("TileArray after add ");
   std::cout << "Printing array:" << std::endl;
   for (unsigned i = 0; i < 200; ++i) {
       std::cout << "m[" << i << "] = " << m[i] << std::endl;
   }
}


int test_9()
{
   print_header(9, "Jacobi Solver");

   unsigned nBlocks(5);
   unsigned blockDim(5);

   TileArray<double> A(nBlocks,nBlocks);
   TileArray<double> b(nBlocks,1);
   TileArray<double> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b.set(row, 0, new CMTile<double>(blockDim,1));
       x.set(row, 0, new CMTile<double>(blockDim,1));

       for (unsigned col = 0; col < nBlocks; ++col) {
           // These are the dimensions of the current tile
           A.set(row, col, new CMTile<double>(blockDim,blockDim));
       }
   }


   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).fill(DebugFunctor());
       x(row,0).fill(ZeroFunctor<double>());
       for (unsigned col = 0; col < nBlocks; ++col) {
           A(row,col).fill(TestFunctor(A.rowOffset(row),A.rowOffset(col)));
       }
   }

   jacobi_solver(A,b,x);

   return 0;
}


int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok(0);

   ok = ok 
      + test_1()
      + test_2()
      + test_3()
      + test_4()
      + test_5()
      + test_6()
      + test_7()
      + test_8()
      + test_9()
   ;

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
