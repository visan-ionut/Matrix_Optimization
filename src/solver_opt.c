/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <string.h>

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double* A, double* B) {

	// A transpose
	// Allocate memory for the transpose of matrix A
	double* AT = (double*)malloc(N * N * sizeof(double));
	memset(AT, 0, N * N * sizeof(double));
	// Compute the transpose of the upper triangular matrix A
	for (int i = 0; i < N; i++) {
		// Since A is upper triangular, j starts from i
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
			register double* a = &AT[i * N];
			register double* b = &B[j];

			register double sum = 0.0;
			// Since AT is lower triangular, k goes only to i
			for (int k = 0; k <= i; k++) {
				sum = sum + *a * *b;
				a++;
				b += N;
			}
			ATXB[i * N + j] = sum;
		}
	}

	// B X A
	// Allocate memory for the result of B * A
	double* BXA = (double*)malloc(N * N * sizeof(double));
	memset(BXA, 0, N * N * sizeof(double));
	// Compute the matrix multiplication B * A
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			register double* a = &B[i * N];
			register double* b = &A[j];

			register double sum = 0.0;
			// Since A is upper triangular, k goes only to j
			for (int k = 0; k <= j; k++) {
				sum = sum + *a * *b;
				a++;
				b += N;
			}
			BXA[i * N + j] = sum;
		}
	}

	// AT X B + B X A
	// Allocate memory for the sum of ATXB and BXA
	double* sum_up = (double*)malloc(N * N * sizeof(double));
	memset(sum_up, 0, N * N * sizeof(double));
	// Compute the element-wise sum of ATXB and BXA
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sum_up[i * N + j] = ATXB[i * N + j] + BXA[i * N + j];
		}
	}

	// sum_up X BT
	// Allocate memory for the final result
	double* result = (double*)malloc(N * N * sizeof(double));
	memset(result, 0, N * N * sizeof(double));
	// Compute the matrix multiplication sum_up * BT
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			register double* a = &sum_up[i * N];
			register double* b = &BT[j];

			register double sum = 0.0;
			for (int k = 0; k < N; k++) {
				sum = sum + *a * *b;
				a++;
				b += N;
			}
			result[i * N + j] = sum;
		}
	}

	// Free allocated memory
	free(AT);
	free(BT);
	free(ATXB);
	free(BXA);
	free(sum_up);

	return result;
}
