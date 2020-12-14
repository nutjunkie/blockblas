#include "EigenSolver.h"
#include "CMTile.h"
#include "Log.h"

template <>
void eigenvalues(TileArray<double> const& A)
{

   CMTile<double> a(A);

   int nRows(a.nRows());
   int nCols(a.nCols());
   int lda(a.leadingDim());

   if (nRows != nCols) {
      std::stringstream ss("eigenvalues called on invalid Tile (");
      ss << nRows << "," << nCols << ")";
      Log::error(ss.str());
      return;
   }

   double* work = new double[nRows];

   int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', nRows, a.data(), lda, work);

   for (unsigned i = 0; i < std::min(nRows, 10); ++i) {
       std::cout << std::fixed << std::showpoint << std::setprecision(10);
       std::cout << "Eigenvalue: " << i  << "   " << work[i] << std::endl;
   }

   if (info !=0) {
      std::string s("Bad return from dsyev: ");
      s += std::to_string(info);
      Log::error(s);
   }

   delete [] work;
}
