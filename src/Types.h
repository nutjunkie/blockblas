#ifndef TYPES_H
#define TYPES_H

#include <complex>

#ifdef __INTEL_COMPILER
#define MKL_Complex16 complex
#endif
typedef std::complex<double> complex;


#include <string>
#include <utility>

#define MAX_ITER 50

enum StorageT { Zero, Diagonal, Banded, Striped, CMDense, RMDense, Dense }; 
enum LayoutT { RowMajor, ColumnMajor };

inline std::string toString(StorageT storage)
{
   std::string s;
   switch (storage) {
      case Zero:      s = "Zero";        break;
      case Diagonal:  s = "Diagonal";    break;
      case Banded:    s = "Banded";      break;
      case Striped:   s = "Striped";     break;
      case CMDense:   s = "Dense (CM)";  break;
      case RMDense:   s = "Dense (RM)";  break;
   }

   return s;
}

typedef std::pair<size_t, size_t> TileIndex;

#endif
