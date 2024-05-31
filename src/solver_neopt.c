/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <string.h>

 /*
  * Add your unoptimized implementation here
  */
double* my_solver(int N, double* A, double* B) {

	// A transpose
	// Allocate memory for the transpose of matrix A
	double* AT = (double*)malloc(N * N * sizeof(double));
	memset(AT, 0, N * N * sizeof(double));
	// Compute the transpose of the upper triangular matrix A
	for (int i = 0; i < N; i++) {
		// Since A is upper triangular, j starts from i to N
		for (int j = i; j < N; j++) {
			AT[j * N + i] = A[i * N + j];
		}
	}

	// B transpose
	// Allocate memory for the transpose of matrix B
	double* BT = (double*)malloc(N * N * sizeof(double));
	memset(BT, 0, N * N * sizeof(double));
	// Compute the transpose of matrix B
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			BT[j * N + i] = B[i * N + j];
		}
	}

	// AT X B
	// Allocate memory for the result of AT * B
	double* ATXB = (double*)malloc(N * N * sizeof(double));
	memset(ATXB, 0, N * N * sizeof(double));
	// Compute the matrix multiplication AT * B
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			// Since AT is lower triangular, k goes only to i
			for (int k = 0; k <= i; k++) {
				ATXB[i * N + j] += AT[i * N + k] * B[k * N + j];
			}
		}
	}

	// B X A
	// Allocate memory for the result of B * A
	double* BXA = (double*)malloc(N * N * sizeof(double));
	memset(BXA, 0, N * N * sizeof(double));
	// Compute the matrix multiplication B * A
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			// Since A is upper triangular, k goes only to j
			for (int k = 0; k <= j; k++) {
				BXA[i * N + j] += B[i * N + k] * A[k * N + j];
			}
		}
	}

	// AT X B + B X A
	// Allocate memory for the sum of ATXB and BXA
	double* sum = (double*)malloc(N * N * sizeof(double));
	memset(sum, 0, N * N * sizeof(double));
	// Compute the element-wise sum of ATXB and BXA
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sum[i * N + j] = ATXB[i * N + j] + BXA[i * N + j];
		}
	}

	// sum X BT
	// Allocate memory for the final result
	double* result = (double*)malloc(N * N * sizeof(double));
	memset(result, 0, N * N * sizeof(double));
	// Compute the matrix multiplication sum * BT
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				result[i * N + j] += sum[i * N + k] * BT[k * N + j];
			}
		}
	}

	// Free allocated memory
	free(AT);
	free(BT);
	free(ATXB);
	free(BXA);
	free(sum);

	return result;
}
