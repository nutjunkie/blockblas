#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H

#include "BlockMatrix.h"

void jacobi_solver(BlockMatrix<double>& x,  BlockMatrix<double> const& A, BlockMatrix<double> const& b);

#endif
