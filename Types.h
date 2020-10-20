#ifndef TYPES_H
#define TYPES_H

#include <complex>
#include <string>

#ifdef __INTEL_COMPILER
#define MKL_Complex16 std::complex<double>
typedef std::complex<double> complex;
#else
typedef std::complex<double> complex;
#endif

#ifdef __INTEL_COMPILER
#include <mkl.h>
#else
#include <veclib/veclib.h>
#endif

enum StorageT { Zero, Diagonal, Banded, Striped, Dense };
enum LayoutT { RowMajor, ColumnMajor };

inline std::string toString(StorageT storage)
{
   std::string s;
   switch (storage) {
      case Zero:      s = "Zero";      break;
      case Diagonal:  s = "Diagonal";  break;
      case Banded:    s = "Banded";    break;
      case Striped:   s = "Striped";   break;
      case Dense:     s = "Dense";     break;
   }

   return s;
}

#endif
