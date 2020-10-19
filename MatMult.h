#include "BlockMatrix.h"

// Accumulates the AxB product into C:
//    C += A x B 

void matrix_product(BlockMatrix<double>& C, BlockMatrix<double> const& A, BlockMatrix<double> const& B);
void matrix_product_sans_diagonal(BlockMatrix<double>& C, BlockMatrix<double> const& A, BlockMatrix<double> const& B);
void matrix_product(VMatrix<double>& C, VMatrix<double> const& A, VMatrix<double> const& B);
