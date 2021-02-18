#include "TileArray.h"
#include "JacobiSolver.h"
#include "ConjugateSolver.h"
#include "TileProduct.h"
#include "EigenSolver.h"
#include "Diagonalize.h"
#include "Davidson.h"
#include "Timer.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include "mkl.h"
#ifdef MYMPI
#include <mpi.h>
#endif 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>



int simonSays();
int simonSaysLess();
int spherium90_512();



int main(int argc, char **argv)
{
    int rank(0), numprocs(1), rv(0);

#ifdef MYMPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
 
    if (rank == 0)  {
       std::cout << "Running on " << numprocs << " procesors" << std::endl;
    }
    //rv = simonSaysLess();
    rv = spherium90_512();

#ifdef MYMPI
    MPI_Finalize();
#endif

    return rv;
}


int spherium90_512()
{
   Timer timer;
   timer.start();

    int rank(0);
#ifdef MYMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

   unsigned nBlocks  = 11;
   //unsigned blocks[] = {1, 248, 248, 247, 247, 246, 247, 248, 247, 246, 245};
   //std::string fname("mat250");

   unsigned blocks[] =  {1, 510, 510, 509, 509, 508, 509, 510, 509, 508, 507};
   std::string fname("mat512");

   TileArray<double> Gv(nBlocks,1);
   TileArray<double> TA(nBlocks,nBlocks);
   for (unsigned bi = 0; bi < TA.nRowTiles(); ++bi) {
       Gv.set(bi,0, new CMTile<double>(blocks[bi],1));
       for (unsigned bj = 0; bj < TA.nColTiles(); ++bj) {
           TA.set(bi,bj, new CMTile<double>(blocks[bi],blocks[bj]));
       }
   }

   timer.stop();
   if (rank == 0) std::cout << "Setup time:     " << timer.format() << std::endl;

   timer.start();
   std::vector<double> mat(readFile(fname));
   if (rank == 0) std::cout << "Number of matrix entries read: " << mat.size() << std::endl;
   TA.bind(&mat[0]);
   timer.stop();
   if (rank == 0) std::cout << "Read time:      " << timer.format() << std::endl;

/*
   timer.start();
   eigenvalues(TA);
   timer.stop();
   std::cout << "LAPACK time:    " << timer.format() << std::endl;
*/

   unsigned subspace(2);
   double Emin(46.0);
   double Emax(47.0);
   int rv(0);

   subspace = 5;
   Emin = 46.0;
   Emax = 55.0;

   timer.start();
   TA.reduce();
   timer.stop();
   if (rank == 0) {
      std::cout << "Reduction time: " << timer.format() << std::endl;
      TA.info("Spherium (4,0) Matrix");
   }

   SymmetricTileArray<double> STA(TA);
   if (rank == 0) {
      STA.info("Symmetric Spherium (4,0) Matrix");
   }

   timer.start();
   rv = diagonalize(STA, subspace, Emin, Emax);
   timer.stop();
   if (rank == 0) std::cout << "FEAST time:     " << timer.format() << std::endl;

   unsigned N(TA.nRows());
   double* gv = new double[N];
   memset(gv, 0, N*sizeof(double));
   gv[0] = 1.0;
   Gv.bind(gv);

   timer.start();
   rv = DavidsonMethod(STA, Gv, 46.51);
   timer.stop();
   if (rank == 0) std::cout << "Davidson time:     " << timer.format() << std::endl;
   delete [] gv;

   return rv;
}



int simonSays()
{
   TileArray<double> TA;

   readMatrix("matrix.bin", TA);

   CMTile<double>* hf(new CMTile<double>(1,1));
   hf->alloc();
   hf->set(0,0,0.0);

   TA.set(0,0, hf);
   TA.addToDiag(1.00);
   TA.info("TA info");
   TA.print("Full Matrix");

   Timer timer;
   timer.start();
   eigenvalues(TA);
   timer.stop();
   std::cout << "LAPACK time: " << timer.format() << std::endl;

   int rv(0);
   unsigned subspace(2);
   double const Emin(0.0);
   double const Emax(1.0);

   timer.start();
   rv = diagonalize(TA, subspace, Emin, Emax);
   timer.stop();

   std::cout << "FEAST time: " << timer.format() << std::endl;

   return rv;
}


int simonSaysLess()
{
   TileArray<double> TA;
   readMatrix("matrix.bin", TA);
   CMTile<double>* hf(new CMTile<double>(1,1));
   hf->alloc();
   hf->set(0,0,0.0);
   TA.set(0,0, hf);
   TA.addToDiag(1.00);
   TA.info("TA info");
   TA.print("Full Matrix");

   SymmetricTileArray<double> STA;
   readMatrix("matrix.bin", STA);
   hf = new CMTile<double>(1,1);
   hf->alloc();
   hf->set(0,0,0.0);
   STA.set(0,0, hf);
   STA.addToDiag(1.00);
   STA.info("STA info");
   STA.print("Full Matrix");

   Timer timer;
   timer.start();
   eigenvalues(TA);
   timer.stop();
   std::cout << "LAPACK time: " << timer.format() << std::endl;

   int rv(0);
   unsigned subspace(2);
   double const Emin(0.0);
   double const Emax(1.0);

   timer.start();
   rv = diagonalize(STA, subspace, Emin, Emax);
   timer.stop();

   std::cout << "FEAST time: " << timer.format() << std::endl;

   return rv;
}

