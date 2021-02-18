#ifndef DAVIDSON_H
#define DAVIDSON_H

#include "TileArray.h"

int DavidsonIteration(TileArray<double> const& A, TileArray<double> const& vec0, double lam0);

template <template<class> class TT>
int DavidsonMethod(TT<double> const& A, TileArray<double> const& vec0, double lam0);

#endif
