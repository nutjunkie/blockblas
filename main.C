#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
 
void Multiply(int n, double** a, double** b, double** c)
{
    int bi=0;
    int bj=0;
    int bk=0;
    int i=0;
    int j=0;
    int k=0;
    int blockSize=10; 
    
    for(bi=0; bi<n; bi+=blockSize)
        for(bj=0; bj<n; bj+=blockSize)
            for(bk=0; bk<n; bk+=blockSize)
                for(i=0; i<blockSize; i++)
                    for(j=0; j<blockSize; j++)
                        for(k=0; k<blockSize; k++)
                            c[bi+i][bj+j] += a[bi+i][bk+k]*b[bk+k][bj+j];
}
 
int test(void)
{
    int n;
    double** A;
    double** B;
    double** C;
    int numreps = 5;
    int i=0;
    int j=0;
    struct timeval tv1, tv2;
    struct timezone tz;
    double elapsed;
    printf ("Please enter matrix dimension n : ");
    scanf("%d", &n);
    // allocate memory for the matrices
     
    ///////////////////// Matrix A //////////////////////////
     
    A =(double **)malloc(n*sizeof(double *));
    A[0] = (double *)malloc(n*n*sizeof(double));
    if(!A)
    {
        printf("memory failed \n");
        exit(1);
    }
    for(i=1; i<n; i++)
    {
        A[i] = A[0]+i*n;
        if (!A[i])
        {
            printf("memory failed \n");
            exit(1);
        }
    }
 
    ///////////////////// Matrix B //////////////////////////
 
    B =(double **)malloc(n*sizeof(double *));
    B[0] = (double *)malloc(n*n*sizeof(double));
    if(!B)
    {
        printf("memory failed \n");
        exit(1);
    }
    for(i=1; i<n; i++)
    {
        B[i] = B[0]+i*n;
        if (!B[i])
        {
            printf("memory failed \n");
            exit(1);
        }
    }
 
    ///////////////////// Matrix C //////////////////////////
    C =(double **)malloc(n*sizeof(double *));
    C[0] = (double *)malloc(n*n*sizeof(double));
    if(!C)
    {
        printf("memory failed \n");
        exit(1);
    }
    for(i=1; i<n; i++)
    {
        C[i] = C[0]+i*n;
        if (!C[i])
        {
            printf("memory failed \n");
            exit(1);
        }
    }
 
    // initialize the matrices
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            A[i][j] = 1;
            B[i][j] = 2;
        }
    }
 
    //multiply matrices
 
    printf("Multiply matrices %d times...\n", numreps);
    for (i=0; i<numreps; i++)
    {
        gettimeofday(&tv1, &tz);
        Multiply(n,A,B,C);
        gettimeofday(&tv2, &tz);
        elapsed += (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
    }
    printf("Time = %lf \n",elapsed);
     
    //deallocate memory
     
    free(A[0]);
    free(A);
    free(B[0]);
    free(B);
    free(C[0]);
    free(C);
    return 0;
}


int main(void
