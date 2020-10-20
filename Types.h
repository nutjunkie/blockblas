#ifndef TYPES_H
#define TYPES_H

#include <complex>
#include <string>

typedef std::complex<double> complex;

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
