#include "ZeroTile.h"
#include "DiagonalTile.h"
#include "StripedTile.h"
#include "CMTile.h"
#include "Functor.h"



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

   ZeroTile<int> ti(5,5);
   ti.alloc(); 
   ti.info();
   ti.print();


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
   v.init();

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
   C.init();

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
