#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "Tile.h"


namespace Log
{
   void warn(std::string const& msg) 
   {
      std::cerr << "Warning: " << msg << std::endl;
   }

   void error(std::string const& msg) 
   {
      std::cerr << "ERROR: " << msg << std::endl;
   }

}


bool zeroTest(unsigned i, unsigned j)
{  
   // diagonal
   //return i == j;
   // checkerboard => 50% sparsity
   return  ((i+j) % 2 == 0);
   // dense 
   //return  true;
}


void print_header(unsigned n, char const* header)
{
   std::string s(" test_");
   s += std::to_string(n) + ": " + std::string(header);

   unsigned len(s.length());
   std::cout << std::string(len+1, '=') << std::endl;
   std::cout << s << std::endl;
   std::cout << std::string(len+1, '=') << std::endl;
   
}


template <class T>
int matrix_residue(T const* a, T const* b, unsigned const n)
{
   double  res(0);
   double  max(0);
   for (unsigned i = 0; i < n; ++i) {
       max = std::max(max, std::abs(a[i]));
       res = std::max(res, std::abs(a[i] - b[i]));
   }

   std::cout << std::scientific;
   std::cout << "Max matrix deviation:   " << res << std::endl;
   std::cout << "Max matrix entry:       " << max << std::endl;
   if (res > 1e-8) {
      std::cout << "FAIL" << std::endl << std::endl;
      return 1;
   }

   std::cout << "PASS" << std::endl << std::endl;
   return 0;
}


template <class T>
int matrix_residue(Tile<T> const& a, Tile<T> const& b)
{
   unsigned n(a.nCols()*a.nRows());
   unsigned m(a.nCols()*a.nRows());

   if (n!=m) {
      std::cerr << "Mismatched matrix dimensions in matrix_residue" << std::endl;
      return 1;
   }

   return matrix_residue(a.data(), b.data(), n);
}


#endif
