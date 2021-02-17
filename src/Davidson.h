#ifndef DAVIDSON_H
#define DAVIDSON_H

#include "TileArray.h"

int DavidsonIteration(TileArray<double> const& A, TileArray<double> const& vec0, double lam0);

int DavidsonMethod(TileArray<double> const& A, TileArray<double> const& vec0, double lam0);

#endif
