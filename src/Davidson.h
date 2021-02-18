#ifndef DAVIDSON_H
#define DAVIDSON_H

#include "TileArray.h"

// A    - CI matrix in SymmetricTileArray form
// vec0 - Initial guess eigenvector
// lam0 - Initial guess eigenvalue

template <template<class> class TT>
int DavidsonMethod(TT<double> const& A, TileArray<double> const& vec0, double lam0);

#endif
