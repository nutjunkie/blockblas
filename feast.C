/*
! Example program for using Intel MKL Extended Eigensolvers (sparse format).
!
! Consider the matrix A:
!
!                 |  5   2   1   1   0   0   0   0   0   0   0  |
!                 |  2   6   3   1   1   0   0   0   0   0   0  |
!                 |  1   3   6   3   1   1   0   0   0   0   0  |
!                 |  1   1   3   6   3   1   1   0   0   0   0  |
!                 |  0   1   1   3   6   3   1   1   0   0   0  |
!    A    =       |  0   0   1   1   3   6   3   1   1   0   0  |,
!                 |  0   0   0   1   1   3   6   3   1   1   0  |
!                 |  0   0   0   0   1   1   3   6   3   1   1  |
!                 |  0   0   0   0   0   1   1   3   6   3   1  |
!                 |  0   0   0   0   0   0   1   1   3   6   2  |
!                 |  0   0   0   0   0   0   0   1   1   2   5  |
!
! stored as sparse matrix  matrix (DOUBLE PRECISION version).
! B is a unit matrix:
!
!                 |  1   0   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   0   1  |
!
!*******************************************************************************/

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



int diagonalize(TileArray<double>&, unsigned const subspace, const double Emin, const double Emax);
int sample();
int stephen();
int stephen100();

std::vector<double> readFile(std::string const& filename)
{

   std::vector<double> mat;
   std::ifstream ifs(filename.c_str(), std::ios::in);
   if (!ifs.is_open()) {
      std::cerr << "Failed to open flie " << filename << std::endl;
      return mat;
   }

   double x;
   while (ifs >> x) { mat.push_back(x); }

   ifs.close();

   return mat;
}



bool run_sample  = false;
bool run_stephen = true;

int main(int argc, char **argv)
{
    int rank(1),numprocs;
#ifdef MYMPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
 
    //if (run_sample) return sample();
    if (run_stephen && rank==1) {
        int rv = stephen100();
#ifdef MYMPI
        MPI_Finalize();
#endif
        return rv;
    }

    int rv(0);
    unsigned const nBlocks(1);
    unsigned const blockSize(21);
    unsigned const matsize(nBlocks*blockSize);
    unsigned const subspace(8);
    
    std::vector<int> stripes{-3,-2,-1,0,1,2,3};


        // Striped test matrix
        StripedTile<double> A(matsize,matsize, stripes);
        A.fill(StencilFunctor<double>(0.1));
        //A.fill(TestFunctor());

        A.info("A Matrix");
        A.print("A Matrix");

    
        TileArray<double> bA(nBlocks,nBlocks);
        for (unsigned bi = 0; bi < bA.nRowTiles(); ++bi) {
            for (unsigned bj = 0; bj < bA.nColTiles(); ++bj) {
                bA.set(bi,bj, new CMTile<double>(blockSize,blockSize));
            }
        }

        CMTile<double> dense(A);
        bA.bind(dense.data());


    TileArray<double> bB(nBlocks,nBlocks);
    for (unsigned bi = 0; bi < bB.nRowTiles(); ++bi) {
        for (unsigned bj = 0; bj < bB.nColTiles(); ++bj) {
            if (bi == bj) {
               bB.set(bi,bj, new CMTile<double>(blockSize,blockSize));
               bB(bi,bj).fill(TestFunctor());
               bB(bi,bj).scale(1.0/(bi+bj+1));
            }else {
               //bB.set(bi,bj, new StripedTile<double>(blockSize,blockSize,stripes));
               //bB(bi,bj).fill(StencilFunctor<double>(0.1));
               bB.set(bi,bj, new ZeroTile<double>(blockSize,blockSize));
               bB(bi,bj).fill0();
            }
        }
    }

   // bB.print("bB Matrix");

    double const Emin(5.00);
    double const Emax(5.5);

    Timer timer;
    bB.info();

    timer.start();
    eigenvalues(bA);
    timer.stop();
    std::cout << "LAPACK time: " << timer.format() << std::endl;


    timer.start();
    rv = diagonalize(bA, subspace, Emin, Emax);
    timer.stop();

    std::cout << "FEAST time: " << timer.format() << std::endl;

    return rv;
}



int diagonalize(TileArray<double>& A, TileArray<complex>& bmAc, unsigned const subspace, const double Emin, const double Emax)
{
    const MKL_INT N = A.nCols();

    double* res  = new double[N];
    double* E    = new double[N];
    double* X    = new double[N*subspace];
    double* work = new double[N*subspace];
    double* Aq   = new double[subspace*subspace];
    double* Sq   = new double[subspace*subspace];

    complex* workc = new complex[N*subspace];
    
    double  epsout(0.0);
    MKL_Complex16 Ze;
    MKL_INT M0(subspace);
    MKL_INT M = M0;
    MKL_INT loop(0);
    MKL_INT info(0);
    MKL_INT ijob(-1);
    MKL_INT fpm[128];

    std::cout << "Subspace size: " << subspace << std::endl;

    feastinit(fpm);

    fpm[0] = 1; /* Generate runtime messages */
    fpm[5] = 1; /* Second stopping criteria  */

    while ( ijob != 0 )
    {
        //  dfeast_srci(&ijob,&N,&Ze,work,workc,Aq,Sq,fpm,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);
        //
        //  void dfeast_srci (
        //   - MKL_INT* ijob,        job indicator
        //   - const MKL_INT* n,     the size of the problem
        //     MKL_Complex16* ze,    the coordinate of the complex countour
        //   - double* work,         workspace dimension n by m0
        //   - MKL_Complex16* workc, workspace dimension n by m0
        //   - double* aq,           workspace dimension m0 by m0
        //   - double* sq,           workspace dimension m0 by m0
        //   - MKL_INT* fpm,         parameters
        //     double* epsout,       relative error
        //     MKL_INT* loop,        number of refinement loops executed
        //   - const double* emin,   lower bound for interval
        //   - const double* emax,   upper rbound for interval
        //   - MKL_INT* m0,          guess for subspace dimension
        //     double* lambda,       eigenvalues
        //   - double* q,            n x m basis for subspace
        //     MKL_INT* m,           number of eigenvalues found in interval
        //     double* res,          relative residual vectors
        //     MKL_INT* info         output
        //   );

        dfeast_srci(&ijob,&N,&Ze,work,workc,Aq,Sq,fpm,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);

        if ( info != 0 ) {
            printf("DFEAST_SRCI info %i \n", info);
            return 1;
        }

        switch ( ijob )
        {
            case -2:
               // New loop
               break;
            case 0:
               // End
               break;
            case 10:
               // Preconditioner, we do nothing here
               break;
            case 11: {
               std::cout<< std::setprecision(7);

//             std::cout << "Complex root: "<< Ze << std::endl;
//             std::cout << "Num RHS:      "<< fpm[23-1] << std::endl;
//             zC = Ze * zB + zA
//
               // Solve (ZeB-A) caux = workc[0:N-1][0:M0-1]
               // and put result into  workc

               bmAc.addToDiag(Ze); 

               unsigned nTiles(bmAc.nRowTiles());
               TileArray<complex> bmBc(nTiles,1);
               TileArray<complex> bmQc(nTiles,1);

               for (unsigned bi = 0; bi < nTiles; ++bi) {
                   bmBc.set(bi,0, new CMTile<complex>(bmAc(bi,0).nRows(),M0));
                   bmQc.set(bi,0, new CMTile<complex>(bmAc(bi,0).nRows(),M0));
               }

               CMTile<complex> vmQc(N,M0);
               DiagonalFunctor<complex> diag(complex(1.0,0.0));
               vmQc.fill(diag);

               bmQc.bind(vmQc.data());
               bmBc.bind(workc);
               //bmBc.print("Bound work director");

               //int rc = jacobi_solver(bmAc, bmQc, bmBc);
               //int rc = conjugate_gradient(bmAc, bmQc, bmBc);
               int rc = conjugate_gradientPC(bmAc, bmQc, bmBc);

               if (rc < 0) {
                  //Log::error("Jacobi failed to converge");
               }
               
               memcpy(workc, vmQc.data(), N*M0*sizeof(complex));
               //for (unsigned i = 0; i < 27; ++i) {
               //     std::cout << "Workc[" << i << "] = " << workc[i] << std::endl;
               //}
               bmAc.addToDiag(-Ze); 

            } break;

            case 30: {
               // Perform multiplication A x[0:N-1][i:j]
               // and put result into   work[0:N-1][i:j]
               // where i = fpm[23]-1, j = fpm[23]+fpm[24]-2
               MKL_INT colsX = fpm[24];
               MKL_INT imem = N*(fpm[23]-1);

               unsigned nTiles(A.nRowTiles());

               TileArray<double> taX(nTiles,1);
               TileArray<double> taW(nTiles,1); 

               for (unsigned bi = 0; bi < nTiles; ++bi) {
                   taX.set(bi,0, new CMTile<double>(A(bi,0).nRows(),colsX));
                   taW.set(bi,0, new CMTile<double>(A(bi,0).nRows(),colsX));
               }

               memset(work+imem, 0, N*colsX*sizeof(double));

               taX.bind(X   +imem);
               taW.bind(work+imem);
               product(A, taX, taW);


            } break;
               

            case 40: {
                // Perform multiplication B x[0:N-1][i:j]
                // and put result into   work[0:N-1][i:j]
                // where i = fpm[23]-1, j = fpm[23]+fpm[24]-2
                MKL_INT colsX = fpm[24];
                MKL_INT imem  = N*(fpm[23]-1);

                // B is the identity, so we just do a copy work <- X
                memcpy(work+imem, X+imem, colsX*N*sizeof(double));

            } break;
 
            default:
                printf("Wrong ijob %i", ijob); fflush(0);
                return 1;
        }
    }


   
    printf("\n");
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n\n");
    printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
    printf("# Modes found/subspace: %d %d     Iterations: %d\n",M,M0,loop);

    double trace(0.0);
    for (int i = 0; i < M; i++) trace = trace+E[i];
    printf("Trace %.15e \n", trace);
    printf("Relative error on the Trace %.15e\n\n",epsout );

    if (run_sample) {
       // This is the print out for the sample problem
       printf("   Computed    |    Expected  \n");
       printf("   Eigenvalues |    Eigenvalues \n");

       //!!!!!!!!!!!!!!! Exact eigenvalues in range (3.0, 7.0) !!!!!!!!!!!!!!!!!!!!!!
       double  Eig[11];
       for (int i=0; i<N; i++ )
           Eig[i] = (double)0.0;

       Eig[0] = (double)3.1715728752538100;
       Eig[1] = (double)4.0000000000000000;
       Eig[2] = (double)4.0000000000000000;
       Eig[3] = (double)4.1292484841890931;
       Eig[4] = (double)4.4066499006731521;
       Eig[5] = (double)6.0000000000000000;

       double eigabs(0.0);
       double r;

       for (int i=0; i<M; i++ )
       {
           r = fabs(E[i]-Eig[i]);
           eigabs = max(eigabs, r);
           printf("%.15e %.15e \n", E[i], Eig[i]);
       }
       printf(" Max value of computed eigenvalue - expected eigenvalues %.15e \n\n", eigabs);
    }else {
       printf("   Computed  Eigenvalues  \n");
       std::cout << std::fixed << std::showpoint << std::setprecision(10);
       for (int i=0; i<M; i++ ) {
           std::cout << "Eigenvalue: " << i  << "   " << E[i] << std::endl;
       }
    }


    delete [] res;
    delete [] E;
    delete [] X;
    delete [] work;
    delete [] Sq;
    delete [] Aq;
    delete [] workc;

    return 0;
}


int diagonalize(TileArray<double>& A, unsigned const subspace, const double Emin, const double Emax)
{
    TileArray<complex> bmAc;
    bmAc.from(A);
    bmAc.scale(-1.0);
    bmAc.info("Complex copied info");
    return diagonalize(A, bmAc, subspace, Emin, Emax);
}



int stephen100()
{
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

   std::vector<double> mat(readFile(fname));
   std::cout << "Number of matrix entries read: " << mat.size() << std::endl;
   TA.bind(&mat[0]);
   TA.info("Stephen Matrix");
   Timer timer;

   timer.start();
   eigenvalues(TA);
   timer.stop();
   std::cout << "LAPACK time: " << timer.format() << std::endl;
/*
*/

   unsigned subspace(5);
   double Emin(46.0);
   double Emax(55.0);

   //Emin = 46.0;
   //Emax = 47.0;

   TA.reduce();

    TileArray<complex> bmAc;
    bmAc.from(TA);
    bmAc.scale(-1.0);
    bmAc.info("Complex copied info");

   timer.start();
   int rv = diagonalize(TA, bmAc, subspace, Emin, Emax);
   timer.stop();

   std::cout << "FEAST time: " << timer.format() << std::endl;

   return rv;
}



int stephen()
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

   //TA.reduce();

   timer.start();
   int rv = diagonalize(TA, subspace, Emin, Emax);
   timer.stop();

   std::cout << "FEAST time: " << timer.format() << std::endl;

   return rv;
}



int sample()
{   

    unsigned const matsize(11);
    unsigned const subspace(8);

    std::vector<int> stripes{-3,-2,-1,0,1,2,3};
    StripedTile<double> A(matsize,matsize,stripes);
    A.fill(StencilFunctor<double>());

    A.set( 0, 0, 5.0);
    A.set( 1, 0, 2.0);
    A.set( 0, 1, 2.0);
    A.set( 9,10, 2.0);
    A.set(10, 9, 2.0);
    A.set(10,10, 5.0);


    double const Emin(3.0);
    double const Emax(7.0);

    CMTile<double> dense(A);
    dense.print("A Matrix");

    TileArray<double> bA(dense);

    int rv = diagonalize(bA, subspace, Emin, Emax);
    return rv;
}


