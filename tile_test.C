#include "util.h"
#include "ZeroTile.h"
#include "DiagonalTile.h"
#include "StripedTile.h"
#include "CMTile.h"
#include "TileProduct.h"
#include "TileArray.h"


#include "JacobiSolver.h"




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
   std::cout << "Trigger a warning on an unbound tile:" << std::endl;
   t.print();

   t.alloc();
   t.fill(f);
   t.info();
   t.print("After binding");

   return 0;
}



int test_3()
{
   print_header(3, "StripedTile");

   DebugFunctor df;
   StencilFunctor sf;

   std::vector<int> stripes = {-1,1,4};
   StripedTile<double> t(5,6,stripes);

   t.alloc().fill(df);
   t.info();
   t.print("After binding");

   CMTile<complex> Z;
   Z.from(t);
   Z.info();
   Z.print("Converted t to complex");

   return 0;
}


int test_4()
{
   print_header(4, "CMTile");

   DebugFunctor df;
   StencilFunctor sf;

   CMTile<double> t(5,6);
   t.alloc().fill(df);

   {
      CMTile<double> s(5,6);
      s.bind(t.data());
      s.info();
      s.set(1,1,10.0);
      s.print("S before going out of scope");
      s.dealloc();
   }
   t.print("After modification via s");


   double m_[] = { 1.01, 2.01, 3.01, 0.00, 0.00, 1.02, 2.02, 3.02, 0.00,
			       0.00, 1.03, 2.03, 3.03, 0.00, 0.00, 0.00, 0.00, 0.00,
				   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

   double m[25];

   CMTile<double> v(5,5);
   CMTile<double> u(3,3);
   v.bind(m);   
   v.fill();

   u.bind(m,5);   
   u.info();
   u.fill(df);
   u.print("u (3x3) from array");

   v.info();   
   v.print("v (5x5) from array");


   for (int i = 0; i < 25; ++i) {
       if (std::abs(m[i]-m_[i]) > 1e-8) return 1;
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

   double c_[] = { 
        21.972100, 43.032100, 64.092100, 85.152100, 106.21210, 127.27210, 
        148.33210, 169.39210, 190.45210, 211.51210, 22.034200, 43.154200, 
        64.274200, 85.394200, 106.51420, 127.63420, 148.75420, 169.87420, 
        190.99420, 212.11420, 22.096300, 43.276300, 64.456300, 85.636300, 
        106.81630, 127.99630, 149.17630, 170.35630, 191.53630, 
        212.71630, 22.158400, 43.398400, 64.638400, 85.878400, 107.11840, 
        128.35840, 149.59840, 170.83840, 192.07840, 213.31840, 22.220500, 
        43.520500, 64.820500, 86.120500, 107.42050, 128.72050, 150.02050, 
        171.32050, 192.62050, 213.92050, 22.282600, 43.642600, 65.002600, 
        86.362600, 107.72260, 129.08260, 150.44260, 171.80260, 193.16260, 
        214.52260, 22.344700, 43.764700, 65.184700, 86.604700, 108.02470, 
        129.44470, 150.86470, 172.28470, 193.70470, 215.12470, 22.406800, 
        43.886800, 65.366800, 86.846800, 108.32680, 129.80680, 151.28680, 
        172.76680, 194.24680, 215.72680};

   for (int i = 0; i < n; ++i) {
       if (std::abs(c[i]-c_[i]) > 1e-8) {
          std::cout << i << ": " << std::defaultfloat << std::abs(c[i]-c_[i]) 
                    << " -> " << std::defaultfloat << c[i] << "  " << std::defaultfloat<< c_[i] << std::endl;
          std::cout << "FAILED" << std::endl;
          return 1;
       }
   }

   return 0;
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
   B.print("Matrix B");

   tile_product(A,B,0.0,C);
   C.print("Product Matrix C");

   CMTile<double> D(10,6);
   CMTile<double> Cc(C);

   D.from(A);
   tile_product(D,B,0.0,Cc);
   C -= Cc;
   if (std::abs(C.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << C.norm2() << std::endl;
      return 1;
   }

   return 0;
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
   TB.info();
   TB.print("TileArray from copy");

   TA += TB;

   TA.print("TileArray after add ");
   //std::cout << "Printing array:" << std::endl;
   //for (unsigned i = 0; i < 200; ++i) {
   //    std::cout << "m[" << i << "] = " << m[i] << std::endl;
   //}
}


int test_9()
{
   print_header(9, "Jacobi Solver");

   unsigned nBlocks(5);
   unsigned blockDim(5);
   unsigned nEigen(3);

   TileArray<double> A(nBlocks,nBlocks);
   TileArray<double> b(nBlocks,1);
   TileArray<double> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b.set(row, 0, new CMTile<double>(blockDim,nEigen));
       x.set(row, 0, new CMTile<double>(blockDim,nEigen));

       for (unsigned col = 0; col < nBlocks; ++col) {
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


   for (unsigned ev = 0; ev < nEigen; ++ev) {
       x(0,0).set(ev,ev,1.0);
   }

   int rc = jacobi_solverLU(A,x,b);

   if (rc < 0) {
      std::cout << "ERROR: JacobiSolver failed to  converged in " << -rc << " cycles" << std::endl;
   }else {
      TileArray<double> result(x);
      result.fill();
      product(A, x, result);
      result -= b;
      result.print("Residual vector(s)");
   }

   return 0;
}


int test_10()
{
   print_header(10, "Complex Jacobi Solver");

   unsigned nBlocks(5);
   unsigned blockDim(5);
   unsigned nEigen(3);

   TileArray<double> A(nBlocks,nBlocks);
   TileArray<double> b(nBlocks,1);
   TileArray<double> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b.set(row, 0, new CMTile<double>(blockDim,nEigen));
       x.set(row, 0, new CMTile<double>(blockDim,nEigen));

       for (unsigned col = 0; col < nBlocks; ++col) {
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


   for (unsigned ev = 0; ev < nEigen; ++ev) {
       x(0,0).set(ev,ev,1.0);
   }

   int rc = jacobi_solverLU(A,x,b);

   if (rc < 0) {
      std::cout << "ERROR: JacobiSolver failed to  converged in " << -rc << " cycles" << std::endl;
   }else {
      TileArray<double> result(x);
      result.fill();
      product(A, x, result);
      result -= b;
      result.print("Residual vector(s)");
   }

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
