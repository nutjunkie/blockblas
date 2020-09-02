class BlockMatrix;
class VMatrix;

// Accumulates the AxB product into C:
//    C += A x B 

void matrix_product(BlockMatrix& C, BlockMatrix& A, BlockMatrix& B);

void matrix_product(VMatrix& C, VMatrix& A, VMatrix& B);
