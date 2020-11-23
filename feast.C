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

#include "BlockMatrix.h"
#include "JacobiSolver.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include "mkl.h"

#include <iostream>
#include <iomanip>

#define max(a, b) (a) < (b) ? (b): (a)



int diagonalize(VMatrix<double,ColumnMajor>&, unsigned const subspace, const double Emin, const double Emax);
int sample();



int main()
{
    int rv(0);
    unsigned const matsize(21);
    unsigned const subspace(11);

    std::vector<int> stripes{-3,-2,-1,0,1,2,3};
    VMatrix<double,ColumnMajor> A;
    A.init(matsize,matsize, stripes).bind(StencilFunctor());

    A.set( 0, 0, 5.0);
    A.set( 1, 0, 2.0);
    A.set( 0, 1, 2.0);
    A.set( 9,10, 2.0);
    A.set(10, 9, 2.0);
    A.set(10,10, 5.0);

    //A.print("A Matrix");

    double const Emin(3.0);
    double const Emax(7.0);

    rv = diagonalize(A, subspace, Emin, Emax);
    return rv;
}



int diagonalize(VMatrix<double,ColumnMajor>& A, unsigned const subspace, const double Emin, const double Emax)
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
        if ( info != 0 )
        {
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
                // Preconditioner, we can do nothing here
                break;
            case 11: {
                std::cout<< std::setprecision(7);

//             std::cout << "Complex root: "<< Ze << std::endl;
//             std::cout << "Num RHS:      "<< fpm[23-1] << std::endl;
//                 
//             zC = Ze * zB + zA
               VMatrix<complex,ColumnMajor> vmAc;
               vmAc.fromDouble(A);
               -vmAc;
               vmAc += Ze;   // overloaded operator only adjusts the diagonal
               vmAc.toDense();
//             vmAc.print("A Matrix");

               VMatrix<complex,ColumnMajor> vmBc;
               vmBc.init(N,M0,Dense).bind(workc);
//             vmBc.print("complex RHS");

               DiagonalFunctor<complex> diag(complex(1.0,0.0));
               VMatrix<complex,ColumnMajor> vmQc;
               vmQc.init(N,M0,Dense).bind(diag);

               // Solve (ZeB-A) caux = workc[0:N-1][0:M0-1]
               // and put result into  workc

               BlockMatrix<complex,ColumnMajor> bmAc(vmAc);
               BlockMatrix<complex,ColumnMajor> bmBc(vmBc);
               BlockMatrix<complex,ColumnMajor> bmQc(vmQc);

               jacobi_solver(bmQc, bmAc, bmBc);
               bmQc(0,0).unbind(workc);

            } break;

            case 30: {
                // Perform multiplication A x[0:N-1][i:j]
                // and put result into   work[0:N-1][i:j]
                // where i = fpm[23]-1, j = fpm[23]+fpm[24]-2
                MKL_INT colsX = fpm[24];
                MKL_INT imem = N*(fpm[23]-1);

                VMatrix<double,ColumnMajor> vmX, vmW;
                vmX.init(N,colsX).bind(X+imem);
                vmW.init(N,colsX).bind(ZeroFunctor<double>());
                A.toDense();
                matrix_product(vmW, A, vmX);
                vmW.unbind(work+imem);

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

    if (false) {
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
       for (int i=0; i<M; i++ ) printf("%.15e  \n", E[i] );
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



int sample()
{   
    unsigned const matsize(11);
    unsigned const subspace(8);

    std::vector<int> stripes{-3,-2,-1,0,1,2,3};
    VMatrix<double,ColumnMajor> A;
    A.init(matsize,matsize, stripes).bind(StencilFunctor());

    A.set( 0, 0, 5.0);
    A.set( 1, 0, 2.0);
    A.set( 0, 1, 2.0);
    A.set( 9,10, 2.0);
    A.set(10, 9, 2.0);
    A.set(10,10, 5.0);

    A.print("A Matrix");

    double const Emin(3.0);
    double const Emax(7.0);

    int rv = diagonalize(A, subspace, Emin, Emax);
    return rv;
}


