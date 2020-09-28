#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H

class BlockMatrix;

void jacobi_solver(BlockMatrix& x,  BlockMatrix const& A, BlockMatrix const& b);

#endif
