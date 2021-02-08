#ifndef DIAGONALIZE_H
#define DIAGONALIZE_H

#include "TileArray.h"
#include "SymmetricTileArray.h"
#include "JacobiSolver.h"
#include "ConjugateSolver.h"
#include "TileProduct.h"
#include "EigenSolver.h"
#include "Timer.h"
#include "Log.h"

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




template <class TT>
int diagonalize(TT& A, unsigned const subspace, const double Emin, const double Emax)
//int diagonalize(TileArray<double>& A, unsigned const subspace, const double Emin, const double Emax)
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
               // Solve (A-ZeB) caux = -workc[0:N-1][0:M0-1]
               // and put result into  workc

               unsigned nTiles(A.nRowTiles());
               TileArray<complex> bmBc(nTiles,1);
               TileArray<complex> bmQc(nTiles,1);

               for (unsigned bi = 0; bi < nTiles; ++bi) {
                   bmBc.set(bi,0, new CMTile<complex>(A(bi,nTiles-1).nRows(),M0));
                   bmQc.set(bi,0, new CMTile<complex>(A(bi,nTiles-1).nRows(),M0));
               }

               CMTile<complex> vmQc(N,M0);
               DiagonalFunctor<complex> diag(complex(1.0,0.0));
               vmQc.Tile<complex>::fill(diag);

               bmQc.bind(vmQc.data());
               bmBc.bind(workc);
               bmBc.scale(-1.0);
               //bmBc.print("Bound work director");

               complex zero(0.0);
               //int rc = jacobi_solver(bmAc, bmQc, bmBc, Ze);
               //int rc = conjugate_gradient(bmAc, bmQc, bmBc, -Ze);
               int rc = conjugate_gradient(A, bmQc, bmBc, -Ze);

               if (rc < 0) {
                  Log::error("Solver failed to converge");
               }
               
               memcpy(workc, vmQc.data(), N*M0*sizeof(complex));
               //for (unsigned i = 0; i < 27; ++i) {
               //     std::cout << "Workc[" << i << "] = " << workc[i] << std::endl;
               //}

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
                   taX.set(bi,0, new CMTile<double>(A(bi,nTiles-1).nRows(),colsX));
                   taW.set(bi,0, new CMTile<double>(A(bi,nTiles-1).nRows(),colsX));
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


    int rank(0);
#ifdef MYMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    if (rank == 0) {
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

#endif
