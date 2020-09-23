#include <iostream>

#include "BlockMatrix.h"
#include "MatMult.h"
#include "Timer.h"
#include "util.h"


ZeroFunctor zeroFunctor;
TestFunctor testFunctor;

int test_1(unsigned blocks, unsigned dim)
{
   std::cout << "==============================" << std::endl;
   std::cout << " test_1: Dense block multiply " << std::endl;
   std::cout << "==============================" << std::endl;

   BlockMatrix A(blocks,blocks);
   BlockMatrix C(blocks,blocks);

   makeDense(A, dim, DebugFunctor());
   makeDense(C, dim, ZeroFunctor());

   VMatrix a, c;
   C.toDense(&c);
   A.toDense(&a);

   //A.info("A Matrix:");
   std::cout << "Blocks: " << blocks << " x " << blocks << std::endl;
   std::cout << "Tiles:  " << dim << " x " << dim << std::endl;
   std::cout << "Size:   " << A.nRows() << " x " << A.nCols() << std::endl;
   
   Timer timer;
   
   timer.start();
   matrix_product(C,A,A);
   timer.stop();
   std::cout << "BlockMatrix time: " << timer.format() << std::endl;;

   //  C.print("Product from BlockMatrix multiplication");
   timer.start();
   matrix_product(c,a,a);
   timer.stop();
   std::cout << "VMatrix time:     " << timer.format() << std::endl;

   VMatrix b;
   C.toDense(&b);

   double* data_c(c.data());
   double* data_b(b.data());
   double  res(0);
   unsigned n(c.nCols()*c.nRows());

   for (unsigned i = 0; i < n; ++i) {
       res = std::max(res, std::abs( data_c[i] - data_b[i]));
   }

   std::cout << "Max deviation:    " << res << std::endl;

   std::cout << "-----------------" << std::endl;
   if (res > 1e-8) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   } 

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


int test_2(unsigned dim)
{
   print_header(2, "Striped - Dense multiply");

   VMatrix a, b, c, d;
   std::vector<int> stripes = {-3,-4, 0 ,4 ,3 };
   
   a.init(dim, dim, stripes).bind(TestFunctor());
   b.init(dim, dim, VMatrix::Dense).bind(TestFunctor());
   c.init(dim, dim, VMatrix::Dense).bind(ZeroFunctor());
   d.init(dim, dim, VMatrix::Dense).bind(ZeroFunctor());

   std::cout << "Size:   " << a.nRows() << " x " << a.nCols() << std::endl;
   std::cout << "Press Enter to begin" << std::endl;
   int ch = getchar();
   

   Timer timer;

   timer.start();
   matrix_product(c,a,b);
   double elapsed(timer.stop());

   std::cout << "Striped Matrix time: " << elapsed << std::endl;

   a.toDense();
   timer.start();
   matrix_product(d,a,b);
   std::cout << "Dense Matrix time:   " << timer.stop() << std::endl;

   return matrix_residue(d,c);
}


int test_3(unsigned n)
{
   print_header(5, "Striped x dense timing test");

   unsigned dim;
   dim = 16384;
   dim = 4096;
   dim = 8192;
   unsigned hdim(dim);

   VMatrix as, bs, cs;
   std::vector<int> stripes = {-4,-2,-1,0,1,2,4};
   as.init(dim, dim, stripes).bind(TestFunctor());
   bs.init(dim, hdim, VMatrix::Dense).bind(TestFunctor());
   cs.init(dim, hdim, VMatrix::Dense).bind(ZeroFunctor());

   VMatrix ab, bb, cb;
   ab.init(dim, dim, 4,4).bindCM(TestFunctor());
   bb.init(dim, hdim, VMatrix::Dense).bindCM(TestFunctor());
   cb.init(dim, hdim, VMatrix::Dense).bindCM(ZeroFunctor());

   std::cout << "Size:   " << as.nRows() << " x " << as.nCols() << std::endl;
   std::cout << "Performing " << n << " iterations for time-averaging" << std::endl;
   std::cout << "Press Enter to begin:" << std::endl;
   int ch = getchar();

   Timer timer;
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(cs,as,bs);
   }   
   double elapsed(timer.stop());
   std::cout << " Average striped time: " << timer.format(elapsed/n) << std::endl;

   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(cb, ab, bb);
   }   
   elapsed = timer.stop();
   std::cout << " Average banded  time: " << timer.format(elapsed/n) << std::endl;


   VMatrix ds;
   ds.init(dim, hdim, VMatrix::Dense).bind(ZeroFunctor());
   as.toDense();
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(ds,as,bs);
   }   
   elapsed = timer.stop();
   std::cout << " Average dense   time: " << timer.format(elapsed/n) << std::endl;

   std::cout << "-----------------" << std::endl;

   // This will fail because of the different data orientation
   return matrix_residue(ds,cs); // + matrix_residue(ds,cb);
}


int test_4(unsigned n)
{
   print_header(5, "Banded x dense timing test");

   unsigned dim;
   dim = 4096;
   dim = 16384;
   dim = 8192;
   unsigned hdim(dim);

   VMatrix a, b, c, d, e, f;
   std::vector<int> stripes = {-4,-1,0,1,4};
   
   a.init(dim, dim, stripes).bind(TestFunctor());

   b.init(dim, hdim, VMatrix::Dense).bindCM(TestFunctor());
   c.init(dim, hdim, VMatrix::Dense).bindCM(ZeroFunctor());
   d.init(dim, hdim, VMatrix::Dense).bindCM(ZeroFunctor());
   e.init(dim, dim, 4, 4).bindCM(TestFunctor());
   f.init(dim, hdim, VMatrix::Dense).bindCM(ZeroFunctor());

   std::cout << "Size:   " << a.nRows() << " x " << a.nCols() << std::endl;
   std::cout << "Press Enter to begin" << std::endl;
   int ch = getchar();
 
   std::cout << "Performing " << n << " iterations for time-averaging" << std::endl;

   Timer timer;
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(c,a,b);
   }   
   double elapsed(timer.stop());
   std::cout << " Average striped time: " << timer.format(elapsed/n) << std::endl;

   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(f,e,b);
   }   
   elapsed = timer.stop();
   std::cout << " Average banded  time: " << timer.format(elapsed/n) << std::endl;

   a.toDense();
   timer.start();
   for (unsigned i = 0; i < n; ++i) {
       std::cout << "." << std::flush;
       matrix_product(d,a,b);
   }   
   elapsed = timer.stop();
   std::cout << " Average dense   time: " << timer.format(elapsed/n) << std::endl;

   std::cout << "-----------------" << std::endl;
   c.print("cmat");
   d.print("dmat");
   f.print("fmat");

   return matrix_residue(d,c) + matrix_residue(d,f);
}


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


int main()
{
   std::cout << "Running tests:" << std::endl;

   test_3(3);

   return 0;

   unsigned total(std::pow(2,10));
   int ok(0);

   for (unsigned p = 14; p < 15; ++p) {
       ok += test_2(std::pow(2,p));
   }

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
