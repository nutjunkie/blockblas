#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>

#include "TileArray.h"



std::vector<double> readFile(std::string const& filename);


void readMatrix(std::string const& fname, TileArray<double>& TA);


bool zeroTest(unsigned i, unsigned j);


void print_header(unsigned n, char const* header);


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
