#ifndef TYPES_H
#define TYPES_H

#include <complex>

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

#include <string>

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

#endif
