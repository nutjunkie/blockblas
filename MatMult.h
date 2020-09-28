class BlockMatrix;
class VMatrix;

// Accumulates the AxB product into C:
//    C += A x B 

void matrix_product(BlockMatrix& C, BlockMatrix const& A, BlockMatrix const& B);
void matrix_product_sans_diagonal(BlockMatrix& C, BlockMatrix const& A, BlockMatrix const& B);

void matrix_product(VMatrix& C, VMatrix const& A, VMatrix const& B);
