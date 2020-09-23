#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "VMatrix.h"
#include "BlockMatrix.h"


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


int matrix_residue(double const* a, double const* b, unsigned const n)
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


int matrix_residue(VMatrix const& a, VMatrix const& b)
{
   unsigned n(a.nCols()*a.nRows());
   unsigned m(a.nCols()*a.nRows());

   if (n!=m) {
      std::cerr << "Mismatched matrix dimensions in matrix_residue" << std::endl;
      return 1;
   }

   return matrix_residue(a.data(), b.data(), n);
}


int matrix_residue(VMatrix const& a, std::string const& fname)
{
   std::ifstream ifs(fname.c_str(), std::ios::in);

   if (!ifs.is_open()) {
      std::cerr << "Failed to open flie " << fname << std::endl;
      return 1;
   }

   double x;
   std::vector<double> vec;
   while (ifs >> x) { vec.push_back(x); }

   unsigned n(a.nCols()*a.nRows());
   if (n != vec.size()) {
      std::cerr << "Mismatched file data in " << fname << std::endl;
      return 1;
   }

   return matrix_residue(a.data(), vec.data(), n);
}


void makeDense(BlockMatrix& bm, unsigned dim, Functor const& functor)
{ 
   unsigned nRows(bm.nRowBlocks());
   unsigned nCols(bm.nColBlocks());

   for (unsigned bi = 0; bi < nRows; ++bi) {
       for (unsigned bj = 0; bj < nCols; ++bj) {
           bm(bi,bj).init(dim,dim, VMatrix::Dense).bind(functor);
       }
   }
}


void makeDiagonal(BlockMatrix& bm, unsigned dim, Functor const& functor)
{ 
   unsigned nRows(bm.nRowBlocks());
   unsigned nCols(bm.nColBlocks());

   for (unsigned bi = 0; bi < nRows; ++bi) {
       for (unsigned bj = 0; bj < nCols; ++bj) {
           if (bi ==  bj) {
              bm(bi,bj).init(dim,dim, VMatrix::Dense).bind(functor);
           }else {
           }
       }
   }
}

#endif
