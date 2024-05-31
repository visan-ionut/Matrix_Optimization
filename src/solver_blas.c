/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>
#include <string.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double* A, double* B) {

    // Allocate memory for matrices
    double* B_copy = (double*)malloc(N * N * sizeof(double));
    double* ATB = (double*)malloc(N * N * sizeof(double));
    double* BXA = (double*)malloc(N * N * sizeof(double));
    double* result = (double*)malloc(N * N * sizeof(double));

    // Make a copy of B since B will be modified
    cblas_dcopy(N * N, B, 1, B_copy, 1);

    // Compute AT X B using the copy of B
    // A is upper triangular, so we need to use CblasTrans to multiply with AT
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, 1.0, A, N, B_copy, N);
    cblas_dcopy(N * N, B_copy, 1, ATB, 1);

    // Compute B X A using the original B
    cblas_dcopy(N * N, B, 1, BXA, 1);
    cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1.0, A, N, BXA, N);

    // Compute AT X B + B X A
    cblas_daxpy(N * N, 1.0, BXA, 1, ATB, 1);

    // Compute (AT X B + B X A) X B
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, ATB, N, B, N, 0.0, result, N);

    // Free allocated memory
    free(B_copy);
    free(ATB);
    free(BXA);

    return result;
}
