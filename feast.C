/*******************************************************************************
* Copyright 2005-2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Extended Eigensolvers
!             C example
!
!*******************************************************************************
!
! Example program for using Intel MKL Extended Eigensolvers (sparse format).
!
! The following routines are used in the example:
!          DGEMM MKL_SPARSE_D_MM MKL_SPARSE_Z_ADD DFEAST_SRCI PARDISO.
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
!  In what follows the symbol ' represents a transpose operation.
!
!  The test performs the following operations :
!
!       1. The code calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!       2. The  code solves  the generalized eigenvalue problem  Ax=eBx using
!          DFEAST_SRCI.
!
!       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!           are the expected eigenvalues  and E(i) are eigenvalues computed
!           with the help of DFEAST_SRCI().
!
!       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I
!          where X is the matrix of eigenvectors computed with the help of
!          DFEAST_SRCI. DGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
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

int diagonalize();


int main()
{
    int rv(0);

    rv= diagonalize();
    return rv;
}



int diagonalize()
{
    std::vector<int> stripes{-3,-2,-1,0,1,2,3};
    VMatrix<double> vmA, vmb;
    vmA.init(11,11, stripes).bind(StencilFunctor());
    vmb.init(11,1,  Dense).bind(DebugFunctor());

    vmA.set(0,0, 5.0);
    vmA.set(1,0, 2.0);
    vmA.set(0,1, 2.0);
    vmA.set( 9,10, 2.0);
    vmA.set(10, 9, 2.0);
    vmA.set(10,10, 5.0);

    //vmA.print("A Matrix");


    //!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const MKL_INT N = 11;
    double   val[38];  //! upper triangle of A
    double   valb[11]; //! vector b

    MKL_Complex16 cval[38], cvalb[11], *cvalz, caux[8*11];                   // cval is complex version of A

    //!!!!!!!!!!!!!!! Declaration of Spblas variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!! A, B - Handles containing a sparse matrix in internal data structure!!!!!!
    //!!! descrA - Structure specifying sparse matrix properties              !!!!!!
    //!!! indexing - Indicates how input arrays are indexed                   !!!!!!
    //!!! layout - Describes the storage scheme for the dense matrix in mm    !!!!!!
    //!!! opeartion - Specifies operation op() on input matrix in mm          !!!!!!

    //!!!!!!!!!!!!!!! Declaration of FEAST variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
    MKL_INT       fpm[128];
    double        Emin, Emax;
    double        epsout;
    MKL_INT       loop;
    MKL_INT       ijob;
    MKL_Complex16 Ze;
    double        work[8*11];
    MKL_Complex16 workc[8*11];
    double        Aq[8*8], Sq[8*8];
    MKL_INT       L = 8;
    MKL_INT       M0, M, info;
    double        E[11];
    double        X[121];
    double        res[11];

    //!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!! Eig - array for storing exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!
    double        Eig[11];
    double        R[11];
    double        Y[121];
    MKL_INT       i,j;
    MKL_INT       ldx = 11, ldy = 11;
    double        trace;
    double        smax, eigabs;
    char          DGEMMC = 'T', DGEMMN = 'N';
    double        one, zero;
    MKL_INT       colsX, imem;

    printf("\n    FEAST DFEAST_SRCI EXAMPLE PROGRAM\n");
    one  = (double)1.0;
    zero = (double)0.0;

    //!!!!!!!!!!!!!!! Exact eigenvalues in range (3.0, 7.0) !!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<N; i++ )
        Eig[i] = (double)0.0;

    Eig[0] = (double)3.1715728752538100;
    Eig[1] = (double)4.0000000000000000;
    Eig[2] = (double)4.0000000000000000;
    Eig[3] = (double)4.1292484841890931;
    Eig[4] = (double)4.4066499006731521;
    Eig[5] = (double)6.0000000000000000;

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!     Initialize matrices     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for ( i=0; i<N*N; i++ )
    {
        X[i] = zero;
    }

    //!!!!!!!  Initialize upper part of symmetric matrix A   !!!!!!!!!!!!!!!!!!!!!!!
    val[0] = (double)5.0;
    val[1] = (double)2.0;
    val[2] = (double)1.0;
    val[3] = (double)1.0;
    for ( i=1; i<9; i++ )
    {
        val[i*4] = (double)6.0;
        val[i*4+1] = (double)3.0;
        val[i*4+2] = (double)1.0;
        val[i*4+3] = (double)1.0;
    }
    val[35] = (double)6.0;
    val[36] = (double)2.0;
    val[37] = (double)5.0;

    //!!!!!!!!!!!!!!!!!!    Initialize unit matrix B   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<N; i++ )
    {
        valb[i]  = (double)1.0;
    }


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!  Initialize upper part of complex symmetric matrix -A !!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for ( i=0; i<38; i++ )
    {
        cval[i] = complex(-val[i], 0.0);
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!   Print matrix dimension          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    printf("Sparse matrix size  %i \n", N);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!      Testing FEAST sparse format drivers !!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!! Search interval [Emin,Emax] including M eigenpairs !!!!!!!!!!!!
    Emin = (double)3.0;
    Emax = (double)7.0;
    M0 = L;
    M = M0;
    printf(" Search interval [ %.15e, %.15e ]  \n", Emin, Emax);
    epsout = (double)0.0;
    loop = 0;
    info = 0;
    //!!!!!!!!!!!!     Initialize PARDISO                   !!!!!!!!!!!!!!!!

    //         Copies all real array valb to complex array cvalb
    for ( i=0; i<N; i++ )
    {
        cvalb[i] = complex(valb[i], 0.0);
    }

    //         Task 1. Initialize FEAST to define the default values for the input
    //         FEAST parameters.
    //
    ijob = -1;
    info = 0;
    feastinit(fpm);
    fpm[0] = 1; /* Generate runtime messages */
    fpm[5] = 1; /* Second stopping criteria  */
    //
    //         Task 2. The  code solves  the generalized eigenvalue problem  Ax=eBx using
    //         DFEAST_SRCI.
    //
    //
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
               VMatrix<complex> vmAc;
               vmAc.fromDouble(vmA);
               -vmAc;
               vmAc += Ze;   // overloaded operator only adjusts the diagonal
               vmAc.toDense();
//             vmAc.print("A Matrix");

               VMatrix<complex> vmBc;
               vmBc.init(N,M0,Dense).bindCM(workc);
//             vmBc.print("complex RHS");

               DiagonalFunctor<complex> diag(complex(1.0,0.0));
               VMatrix<complex> vmQc;
               vmQc.init(N,M0,Dense).bind(diag);

               //!!!!!!!!!!!!!!! Solve (ZeB-A) caux = workc[0:N-1][0:M0-1] !!!!!!!!
               //!!!!!!!!!!!!!!! and put result into workc                 !!!!!!!!

               BlockMatrix<complex> bmAc(vmAc);
               BlockMatrix<complex> bmBc(vmBc);
               BlockMatrix<complex> bmQc(vmQc);

               jacobi_solver(bmQc, bmAc, bmBc);
               bmQc(0,0).unbindCM(workc);

            } break;

            case 30: {
                //!!!!!!!!!!!!! Perform multiplication A x[0:N-1][i:j]      !!!!!!!!
                //!!!!!!!!!!!!! and put result into work[0:N-1][i:j]        !!!!!!!!
                //!!!!!!!!!!!!! where i = fpm[23]-1, j = fpm[23]+fpm[24]-2  !!!!!!!!
                colsX = fpm[24];
                imem = N*(fpm[23]-1);

                VMatrix<double> vmX, vmW;
                vmX.init(N,colsX).bindCM(X+imem);
                vmW.init(N,colsX).bindCM(ZeroFunctor<double>());
                matrix_product(vmW, vmA, vmX);

                vmW.unbindCM(work+imem);
             
/*
                vmX.init(N,colsX).bindCM(work+imem);
                vmX.print("Case 30: work+imem CM");

                vmX.init(N,colsX).bindCM(work+imem);
                vmX.print("Case 30: work+imem");
*/


            } break;

               

            case 40: {
        
                //!!!!!!!!!!!!! Perform multiplication B x[0:N-1][i:j]      !!!!!!!!
                //!!!!!!!!!!!!! and put result into work[0:N-1][i:j]        !!!!!!!!
                //!!!!!!!!!!!!! where i = fpm[23]-1, j = fpm[23]+fpm[24]-2  !!!!!!!!
                colsX = fpm[24];
                imem = N*(fpm[23]-1);

                VMatrix<double> vmX, vmW;
                vmX.init(N,colsX).bindCM(X+imem);
                vmW.init(N,colsX).bindCM(work+imem);

                // B is the identity, so we just do a copy work <- X
                memcpy(work+imem, X+imem, colsX*N*sizeof(double));
                
                //vmX.print("Case 40: vmX");
                //vmW.print("Case 40: vmW");

                } break;
 
            default:
                printf("Wrong ijob %i", ijob); fflush(0);
                return 1;
        }
    }
    //!!!!!!!!!!!!!!! Release memory                              !!!!!!!!!!!!!!!!!!!!
    //
    //         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    //         are the expected eigenvalues  and E(i) are eigenvalues computed
    //         with the help of  DFEAST_SRCI().
    //
    printf("  Number of eigenvalues found %d \n", M);
    printf("   Computed    |    Expected  \n");
    printf("   Eigenvalues |    Eigenvalues \n");
    eigabs = (double)0.0;
    for ( i=0; i<M; i++ )
    {
        R[i] = fabs(E[i]-Eig[i]);
        eigabs = max(eigabs, R[i]);
        printf("%.15e %.15e \n", E[i], Eig[i]);
    }
    printf(" Max value of | computed eigenvalue -expected eigenvalues | %.15e \n", eigabs);
    //
    //         Task 4.  Compute the maximum absolute value of the matrix
    //         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
    //         the help of DFEAST_SRCI.
    //         Call DGEMM (BLAS Level 3 Routine) to compute (X')*X.
    //
    dgemm(&DGEMMC,&DGEMMN,&M,&M,&N,&one,X,&ldx,X,&ldx,&zero,Y,&ldy);

    //          Compute Y=Y-I.
    for ( i=0; i<M; i++ )
        Y[i*M + i] = Y[i*M + i]-(double)1.0;

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace = (double)0.0;
    for ( i=0; i<M; i++ )
    {
        trace = trace+E[i];
    }
    printf("TRACE %.15e \n", trace);
    printf("Relative error on the Trace %.15e \n",epsout );
    printf("Eigenvalues/Residuals\n");
    for ( i=0; i<M; i++ )
    {
        printf("   %d  %.15e %.15e \n",i, E[i], res[i]);
    }
    smax = (double)0.0;
    for ( i=0; i<M; i++ )
    {
        for ( j=0; j<M; j++ )
        {
            smax = max(smax, fabs(Y[i*M + j]));
        }
    }
    printf( "Max of (transposed of X)*X-I %.15e \n", smax);

    return 0;
}
