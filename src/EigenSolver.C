#include "EigenSolver.h"
#include "CMTile.h"
#include "Log.h"
#include "Timer.h"
#include <mkl.h>

template <>
void eigenvalues(TileArray<double> const& A)
{

   CMTile<double> a(A);
   CMTile<double> b(A);

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

   Timer timer;
   timer.start();
   int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', nRows, a.data(), lda, work);
   timer.stop();
   std::cout << "DSYEV time: " << timer.format() << std::endl;

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



   unsigned subspace(3);

   timer.start();
   double vl(0.0), vu(0.0), abstol(-1.0);
   int il(1),iu(subspace),m;

   double* w = new double[subspace];
   double* z = new double[subspace*nRows];
   int* isuppz  = new int[nRows];

   /*
    * The following gives a good description of the arguments.
    * https://www.ibm.com/support/knowledgecenter/SSFHY8_6.2/reference/am5gr_hsspevx.html#am5gr_hsspevx__am5gr_orsspevx
    */

   info = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'I', 'U', nRows, b.data(), lda, 
      vl, vu, il, iu, abstol, &m, w, z, nRows, isuppz);

   timer.stop();
   std::cout << "DSYEVR time: " << timer.format() << std::endl;

   for (unsigned i = 0; i < std::min(nRows, m); ++i) {
       std::cout << std::fixed << std::showpoint << std::setprecision(10);
       std::cout << "Eigenvalue: " << i  << "   " << w[i] << std::endl;
   }

   if (info !=0) {
      std::string s("Bad return from dsyev: ");
      s += std::to_string(info);
      Log::error(s);
   }

   delete [] z;
   delete [] w;
}
