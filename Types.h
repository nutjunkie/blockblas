#ifndef TYPES_H
#define TYPES_H

#include <complex>

#ifdef __INTEL_COMPILER
#define MKL_Complex16 std::complex<double>
typedef std::complex<double> complex;
#else
typedef std::complex<double> complex;
#endif

#ifdef __MAC_OS__
#include <veclib/veclib.h>
#else
#include <mkl.h>
#endif

#include <string>

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

#endif
