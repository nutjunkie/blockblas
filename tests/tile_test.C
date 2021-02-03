#include "util.h"
#include "ZeroTile.h"
#include "DiagonalTile.h"
#include "StripedTile.h"
#include "CMTile.h"
#include "TileProduct.h"
#include "TileArray.h"
#include "SymmetricTileArray.h"
#include "EigenSolver.h"

#include "JacobiSolver.h"
#include "ConjugateSolver.h"


DebugFunctor<double> df;


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

   DiagonalTile<double> t(5,6);
   t.info();
   std::cout << "Trigger a warning on an unbound tile:" << std::endl;
   t.print();

   t.alloc();
   t.fill(df);
   t.info();
   t.print("After binding");

   return 0;
}



int test_3()
{
   print_header(3, "StripedTile");

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
   v.fill0();

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
   A.fill(StencilFunctor<double>());

   CMTile<double> B(6,8), C(10,8), Cc(10,8);
   B.fill(df);
   C.fill0();
   Cc.fill0();

   B.print("Matrix B");

   CMTile<double> D(10,6);
   D.from(A);

   A.print("Matrix A");
   tile_product(A,B,0.0,C);
   C.print("Product Matrix A.B");

   D.print("Dense Matrix D");
   tile_product(D,B,0.0,Cc);
   Cc.print("Product Matrix D.B");

   C -= Cc;
   //C.print("Difference :Matrix C");
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
   A.fill(df);

   DiagonalTile<double> D(8,12);
   D.fill(df);

   A += D;
   A.print("After adding diagonal");

   CMTile<double> C(A);
   A += C;
   A.print("After adding C");

   A.scale(-1.0);

   CMTile<double> B(10,12);
   B.bind(a);
   B.print("Buffer contents ");

   return 0;
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
           TA(bi,bj).fill(df);
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
       b(row,0).fill(df);
       x(row,0).fill(ZeroFunctor<double>());
       for (unsigned col = 0; col < nBlocks; ++col) {
           A(row,col).fill(TestFunctor(A.rowOffset(row),A.colOffset(col)));
       }
   }


   for (unsigned ev = 0; ev < nEigen; ++ev) {
       x(0,0).set(ev,ev,1.0);
   }

   int rc = jacobi_solverLU(A,x,b);

   if (rc < 0) {
      std::cout << "ERROR: JacobiSolver failed to  converged in " << -rc << " cycles" << std::endl;
      return 1;
   }else {
      return 0;
      TileArray<double> result(x);
      result.fill();
      product(A, x, result);
      result -= b;
      result.print("Residual vector(s)");
   }

   return 0;
}


template <class T>
int test_10()
{
   print_header(10, "Complex Jacobi Solver");

   unsigned nBlocks(3);
   unsigned blockDim(4);
   unsigned nEigen(3);

   CMTile<T> pA(nBlocks*blockDim,nBlocks*blockDim);
   pA.fill(StencilFunctor<T>(1.0));

   TileArray<T> A(nBlocks,nBlocks);
   TileArray<T> b(nBlocks,1);
   TileArray<T> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b.set(row, 0, new CMTile<T>(blockDim,nEigen));
       x.set(row, 0, new CMTile<T>(blockDim,nEigen));

       for (unsigned col = 0; col < nBlocks; ++col) {
           A.set(row, col, new CMTile<T>(blockDim,blockDim));
       }
   }


   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).fill(DebugFunctor<T>());
       x(row,0).fill(ZeroFunctor<T>());
   }

   A.bind(pA.data());

   A.addToDiag(10.0);


   for (unsigned ev = 0; ev < nEigen; ++ev) {
       x(0,0).set(ev,ev,1.0);
   }

   A.print("A matrix");
   b.print("b vector");

   int rc = jacobi_solverLU(A,x,b);

   if (rc < 0) {
      std::cout << "ERROR: JacobiSolver failed to  converged in " << -rc << " cycles" << std::endl;
      return 1;
   }else {
      return 0;
      TileArray<T> result(x);
      result.fill();
      product(A, x, result);
      result -= b;
      result.print("Residual vector(s)");
   }

   return 0;
}



int test_11()
{
   print_header(11, "TileArray -> CMTile conversion");

   unsigned nBlocks(3);
   unsigned blockDim(4);

   TileArray<double> A(nBlocks,nBlocks);

   for (unsigned row = 0; row < nBlocks; ++row) {
       for (unsigned col = 0; col < nBlocks; ++col) {
           A.set(row, col, new CMTile<double>(blockDim,blockDim));
       }
   }

   CMTile<double> pA(A.nRows(), A.nCols());
   pA.fill(DebugFunctor<double>());
   A.bind(pA.data());

   A.print("TileArray A:");

   CMTile<double> D(A);
   D.print("CMTile from A:");

   D -= pA;

   return  (std::abs(D.norm2()) < 1e-8) ? 0 : 1;
}


int test_12()
{
   print_header(12, "LAPACK eigenvalues");

   unsigned nBlocks(3);
   unsigned blockDim(4);
   unsigned nEigen(3);

   CMTile<double> pA(nBlocks*blockDim,nBlocks*blockDim);
   pA.fill(StencilFunctor<double>(1.0));

   TileArray<double> A(nBlocks,nBlocks);

   for (unsigned row = 0; row < nBlocks; ++row) {
       for (unsigned col = 0; col < nBlocks; ++col) {
           A.set(row, col, new CMTile<double>(blockDim,blockDim));
       }
   }

   A.bind(pA.data());
   eigenvalues(A);

   return 0;
}


int test_13()
{
   print_header(13, "LAPACK eigenvalues");

   unsigned nBlocks(3);
   unsigned blockDim(3);

   CMTile<complex> pA(nBlocks*blockDim,nBlocks*blockDim);
   pA.fill(StencilFunctor<complex>(1.0));

   TileArray<complex> A(nBlocks,nBlocks);

   for (unsigned row = 0; row < nBlocks; ++row) {
       for (unsigned col = 0; col < nBlocks; ++col) {
           A.set(row, col, new CMTile<complex>(blockDim,blockDim));
       }
   }

   A.bind(pA.data());

   return 0;
}


int test_14()
{
   print_header(14, "TileProduct test");

   unsigned nBlocks(4);
   unsigned blockDim(4);
   TileArray<double> A(nBlocks,nBlocks);
   TileArray<double> B(nBlocks,nBlocks);
   TileArray<double> C(nBlocks,nBlocks);

   for (unsigned row = 0; row < nBlocks; ++row) {
       for (unsigned col = 0; col < nBlocks; ++col) {
           A.set(row, col, new CMTile<double>(blockDim,blockDim));
           B.set(row, col, new CMTile<double>(blockDim,blockDim));
           C.set(row, col, new CMTile<double>(blockDim,blockDim));
           A(row, col).fill(TestFunctor());
           B(row, col).fill(TestFunctor());
           C(row, col).fill0();
       }
   }

   product(A, B, C);
   C.print("Dense C");
}




template <class T>
int test_15()
{
   print_header(15, "Conjugate solver");

   unsigned nBlocks(3);
   unsigned blockDim(4);
   unsigned nEigen(3);

   CMTile<T> pA(nBlocks*blockDim,nBlocks*blockDim);
   pA.fill(StencilFunctor<T>(1.0));

   TileArray<T> A(nBlocks,nBlocks);
   TileArray<T> b(nBlocks,1);
   TileArray<T> x(nBlocks,1);

   for (unsigned row = 0; row < nBlocks; ++row) {
       b.set(row, 0, new CMTile<T>(blockDim,nEigen));
       x.set(row, 0, new CMTile<T>(blockDim,nEigen));

       for (unsigned col = 0; col < nBlocks; ++col) {
           A.set(row, col, new CMTile<T>(blockDim,blockDim));
       }
   }


   for (unsigned row = 0; row < nBlocks; ++row) {
       b(row,0).fill(DebugFunctor<T>());
       x(row,0).fill(ZeroFunctor<T>());
   }

   A.bind(pA.data());



   for (unsigned ev = 0; ev < nEigen; ++ev) {
       x(0,0).set(ev,ev,1.0);
   }

   A.print("A matrix");
   b.print("b vector");
   x.print("x vector");

   int rc = conjugate_gradient(A,x,b,10.0);

   if (rc < 0) {
      std::cout << "ERROR: ConjugateSolver failed to  converged in " << -rc << " cycles" << std::endl;
      return 1;
   }else {
      return 0;
      TileArray<T> result(x);
      result.fill();
      product(A, x, result);
      result -= b;
      result.print("Residual vector(s)");
   }

   return 0;
}


int test_16()
{
   print_header(16, "Trans CMTile");

   const size_t n(80);

   double a[n], b[n], c[n];
   memset(a,0,80*sizeof(double));
   memset(b,0,80*sizeof(double));
   memset(c,0,80*sizeof(double));

   CMTile<double> A(6,9), B(6,8), C(9,8);

   A.bind(a);
   A.fill(df);

   B.bind(b,10);  // give b a different leading dimension
   B.fill(df);

   C.bind(c,10);

   A.print("Matrix A");
   B.print("Matrix B");

   tile_product(A,B,0.0,C,CblasTrans);

   C.print("Product Matrix C");

   return 0;
}


int test_17()
{
   print_header(17, "Diagonal Tile product");

   const size_t n(80);

   double a[n], b[n], c[n];
   memset(a,0,80*sizeof(double));
   memset(b,0,80*sizeof(double));
   memset(c,0,80*sizeof(double));

   DiagonalTile<double> A(6,6);
   CMTile<double> B(6,8), C(6,8);

   A.bind(a);
   A.fill(df);

   B.bind(b,10);  // give b a different leading dimension
   B.fill(df);

   C.bind(c,10);

   A.print("Matrix A");
   B.print("Matrix B");

   tile_product(A,B,0.0,C);
   C.print("Product Matrix C");

   CMTile<double> D(A);
   tile_product(D,B,-1.0,C);

   if (std::abs(C.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << C.norm2() << std::endl;
      return 1;
   }

   return 0;
}


int test_18()
{
   print_header(18, "Real/Imaginary parts");

   CMTile<complex> A(10,5);
   A.fill(DebugFunctor<complex>());
   A.scale(complex(2.0,1.5));

   CMTile<double>  B(10,5);
   B.fill(DebugFunctor<double>());

   B.scale(2.0);

   CMTile<double>  Ar;

   A.print("Complex Matrix");

   A.getReal(Ar);
   Ar.print("Real part");
   B -= Ar;

   if (std::abs(B.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << B.norm2() << std::endl;
      return 1;
   }

   B.fill(DebugFunctor<double>());
   B.scale(1.5);

   A.getImag(Ar);
   Ar.print("Imag part");

   B -= Ar;
   if (std::abs(B.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << B.norm2() << std::endl;
      return 1;
   }

   return 0;
}


int test_19()
{
   print_header(19, "Real/Imaginary parts");

   CMTile<complex> A(10,5);
   A.fill(DebugFunctor<complex>());
   A.scale(complex(2.0,1.5));

   CMTile<double>  B(10,5);
   B.fill(DebugFunctor<double>());
   B.scale(-2.0);

   A.addReal(B);

   B.fill(DebugFunctor<double>());
   B.scale(-1.5);

   A.addImag(B);
   A.print("zeroed matrix");

   if (std::abs(A.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << A.norm2() << std::endl;
      return 1;
   }

   return 0;
}


int test_20()
{
   print_header(20, "Real*Imaginary product");

   CMTile<double>  A(10,5);
   CMTile<complex> B(5, 8);
   CMTile<complex> C1(10,8);
   CMTile<complex> C2(10,8);
   CMTile<complex> D(10,5);

   complex zero(0.0);
   complex scale(1.0,1.5);

   A.fill(DebugFunctor<double>());
   B.fill(DebugFunctor<complex>());
   B.scale(scale);
   D.fill0();
   D.addReal(A);
   C1.alloc();
   C2.alloc();


   tile_product(D, B, zero, C2);
   C2.print("Complex product");

   tile_product(A, B, zero, C1);
   C1.print("Real product");

   C1 -= C2;

   if (std::abs(C1.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << C1.norm2() << std::endl;
      return 1;
   }

   return 0;
}


int test_21()
{
   print_header(21, "Transpose striped product ");

   std::vector<int> stripes = {-1,1,4};
   StripedTile<double> st(6,10,stripes);
   st.fill(df);
   st.print("Striped matrix");

   CMTile<double> A(6,8);
   A.fill(DebugFunctor<double>());
   A.print("Dense matrix");

   CMTile<double> B(10,8);
   B.alloc();

   tile_product(st, A, 0.0, B,CblasTrans);

   B.print("product");

   return 0;
}


int test_22()
{
   print_header(22, "SymmetricTileArray");

   size_t nTiles(4);
   size_t nRows(5);
   size_t nCols(3);
   SymmetricTileArray<double> STA(nTiles);
   TileArray<double> TA(nTiles, nTiles);
   TileArray<double> TB(nTiles, nTiles);
   TileArray<double> TC(nTiles, nTiles);

   int dim[] = {1,2,3,4};

   // This is not good, the TestFunctor can't be used here becase it relies on 
   // the existence of the last colum of tiles.  So we need to allocate before
   // filling.
   for (unsigned bj = 0; bj < nTiles; ++bj) {
       for (unsigned bi = 0; bi < nTiles; ++bi) {
           TA.set(bi,bj, new CMTile<double>(dim[bi],dim[bj]));
           TB.set(bi,bj, new CMTile<double>(dim[bi],dim[bj]));
           TC.set(bi,bj, new CMTile<double>(dim[bi],dim[bj]));
           if (bi <= bj) {
              STA.set(bi,bj, new CMTile<double>(dim[bi],dim[bj]));
           }
       }
   }

   for (unsigned bj = 0; bj < nTiles; ++bj) {
       for (unsigned bi = 0; bi < nTiles; ++bi) {
           TA(bi,bj).fill(TestFunctor(STA.rowOffset(bi),STA.colOffset(bj)));
           TB(bi,bj).fill0();
           TC(bi,bj).fill0();
           if (bi <= bj) {
              STA(bi,bj).fill(TestFunctor(STA.rowOffset(bi),STA.colOffset(bj)));
           }
       }
   }

   TA.print("TileArray");
   STA.info("SymmetricTileArray");


   product(TA,TA,TB);
   TB.print("Tile product");

   product(STA,TA,TC);
   TC.print("SymmetricTileArray product");

   TC -= TB;

   TC.print("SymmetricTileArray difference");
   
   if (std::abs(TC.norm2()) > 1e-8) {
      std::cout << "FAILED: norm2 =" << TC.norm2() << std::endl;
      return 1;
   }

   return 0;
}


int test_23()
{
   print_header(23, "Tile sorting");

   size_t nTiles(4);
   SymmetricTileArray<double> STA(nTiles);

   int dim[] = {1,2,3,4};

   std::vector<int> stripes1 = {-1};
   std::vector<int> stripes2 = {-1,1,2};

   // This is not good, the TestFunctor can't be used here becase it relies on 
   // the existence of the last colum of tiles.  So we need to allocate before
   // filling.
   
   STA.set(0,0, new       CMTile<double>(dim[0],dim[0]));
   STA.set(0,1, new DiagonalTile<double>(dim[0],dim[1]));
   STA.set(0,2, new     ZeroTile<double>(dim[0],dim[2]));
   STA.set(0,3, new DiagonalTile<double>(dim[0],dim[3]));

   STA.set(1,1, new       CMTile<double>(dim[1],dim[1]));
   STA.set(1,2, new  StripedTile<double>(dim[1],dim[2],stripes1));
   STA.set(1,3, new DiagonalTile<double>(dim[1],dim[3]));

   STA.set(2,2, new       CMTile<double>(dim[2],dim[2]));
   STA.set(2,3, new  StripedTile<double>(dim[2],dim[3],stripes1));

   STA.set(3,3, new  StripedTile<double>(dim[3],dim[3],stripes2));

   for (unsigned bi = 0; bi < nTiles; ++bi) {
       for (unsigned bj = bi; bj < nTiles; ++bj) {
           STA(bi,bj).fill(TestFunctor(STA.rowOffset(bi),STA.colOffset(bj)));
       }
   }

   STA.info("Symmetric Tile");
   STA.print("Symmetric Tile");

   std::vector<TileIndex> indices(STA.sort());
   std::vector<TileIndex>::iterator iter;
   for (iter = indices.begin(); iter  != indices.end(); ++iter) {
       size_t ai(std::min(iter->first, iter->second));
       size_t aj(std::max(iter->first, iter->second));
       std::cout << "(" << iter->first << "," << iter->second << ")" << "   " << STA(ai,aj).numData() <<  std::endl;
   }


   return 0;
}





int main()
{
   std::cout << "Running tests:" << std::endl;
   int ok(0);

   ok = ok 
/*
      + test_1()
      + test_2()
      + test_3()
      + test_4()
      + test_5()
      + test_6()
      + test_7()
      + test_8()
      + test_9()
      + test_10<double>()
      + test_11()
      + test_12()
      + test_14()
      + test_16()
//      + test_15<double>()
      + test_17()
      + test_18()
      + test_19()
      + test_20()
      + test_21()
*/
      + test_22()
      + test_23()
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
