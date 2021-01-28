#include "TileArray.h"
#include "JacobiSolver.h"
#include "ConjugateSolver.h"
#include "TileProduct.h"
#include "EigenSolver.h"
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

#define max(a, b) (a) < (b) ? (b): (a)



void readMatrix(std::string const& fname, TileArray<double>& TA);
int simonSays();
int spherium90_10();
int spherium90_512();
int diagonalize(TileArray<double>& A, unsigned const subspace, const double Emin, const double Emax);




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
    rv = simonSays();
    //rv = spherium90_10();
    //rv = spherium90_512();

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

   TileArray<double> TA(nBlocks,nBlocks);
   for (unsigned bi = 0; bi < TA.nRowTiles(); ++bi) {
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

   unsigned subspace(5);
   double Emin(46.0);
   double Emax(55.0);

   timer.start();
   TA.reduce();
   timer.stop();
   if (rank == 0) {
      std::cout << "Reduction time: " << timer.format() << std::endl;
      TA.info("Spherium (4,0) Matrix");
   }

   timer.start();
   int rv = diagonalize(TA, subspace, Emin, Emax);
   timer.stop();

   if (rank == 0) std::cout << "FEAST time:     " << timer.format() << std::endl;

   return rv;
}



int spherium90_10()
{
   unsigned nBlocks  = 11;
   unsigned blocks[] = { 1, 6, 6, 5, 5, 4, 5, 6, 5, 4, 3};

   TileArray<double> TA(nBlocks,nBlocks);
   for (unsigned bi = 0; bi < TA.nRowTiles(); ++bi) {
       for (unsigned bj = 0; bj < TA.nColTiles(); ++bj) {
           TA.set(bi,bj, new CMTile<double>(blocks[bi],blocks[bj]));
       }
   }

   std::vector<double> mat(readFile("mat10"));
   std::cout << "Number of matrix entries read: " << mat.size() << std::endl;
   TA.bind(&mat[0]);

   Timer timer;
   timer.start();
   eigenvalues(TA);
   timer.stop();
   std::cout << "LAPACK time: " << timer.format() << std::endl;

   unsigned subspace(5);
   double const Emin(46.0);
   double const Emax(55.0);

   // Reduce the dense matrices to bespoke storage
   // TA.reduce();

   timer.start();
   int rv = diagonalize(TA, subspace, Emin, Emax);
   timer.stop();

   std::cout << "FEAST time: " << timer.format() << std::endl;

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

