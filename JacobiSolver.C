#include "JacobiSolver.h"


// Solves A.x = b where A has been LU factorized and b is overwritten with the solution x.
template <>
void lu_solve(Tile<double> const& At, Tile<double>& bt, int* ipiv)
{
   if (At.storage() != CMDense || bt.storage() != CMDense) {
      Log::error("Invalid tiles passed to lu_solve");
      return;
   }

   CMTile<double> const& A = dynamic_cast< CMTile<double> const&>(At);
   CMTile<double>& b = dynamic_cast< CMTile<double>&>(bt);

   if (A.nRows() != A.nCols() || A.nCols() != b.nRows() || !A.isBound() || !b.isBound() ) {
      std::stringstream ss("lu_solve() called on invalid Tile combination: ");
      ss << "(" << A.nRows() << "," << A.nCols() << ") x "
         << "(" << b.nRows() << "," << b.nCols() << ") " << std::endl;
      ss << "isBound " << A.isBound() << " " << b.isBound() << std::endl;
      Log::error(ss.str());
      return;
   }

   int n(A.nRows());
   int nrhs(b.nCols());
   int lda(A.leadingDim());
   int ldb(b.leadingDim());

   int info = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'U', n, nrhs, A.data(), lda, ipiv, b.data(), ldb); 

   if (info != 0) {
      std::string s("Bad return from dsytrs: ");
      s += std::to_string(info);
      Log::error(s);
   }
}


