Implementation

tema3_neopt:

This is the basic implementation without optimizations. The purpose of
this version is to directly implement matrix multiplication and addition
operations without using optimized functions from external libraries.

Implementation Description:

Matrix A Transposition:
To compute AT (the transpose of matrix A), we allocate a vector of size
N×N to store the transpose. Since A is an upper triangular matrix, we only
traverse the elements above the main diagonal (inclusive) and copy them to AT
at the corresponding positions. The copy formula used is AT[j * N + i] = A[i * N + j],
where i and j are the indices of the matrix rows and columns.

Matrix B Transposition:
Similarly to A, we allocate a vector of size N×N to store the transpose of matrix B.
We traverse all elements in B and copy them to BT using the formula BT[j * N + i] = B[i * N + j].

Computing the Product AT×B:
We allocate a vector for the result of the product AT×B. For each element of the resulting
matrix, we calculate the sum of the corresponding elements' products from AT and B.
Since AT is a lower triangular matrix (the transpose of an upper triangular matrix),
we only traverse the necessary elements from AT up to the main diagonal.

Computing the Product B×A:
Similar to the previous step, we allocate a vector for the result of the product B×A.
For each element of the resulting matrix, we calculate the sum of the corresponding elements'
products from B and A. Since A is an upper triangular matrix, we only traverse the necessary
elements from A up to the main diagonal.

Adding the Two Products:
We allocate a vector to store the sum of the matrices resulting from AT×B and B×A.
We traverse all elements and compute their sum element by element.

Computing the Final Product (AT×B+B×A)×BT:
We allocate a vector for the final result. For each element of the resulting matrix, we calculate
the sum of the corresponding elements' products from the summed matrix and BT.

Memory Management:
All intermediate and final matrices are dynamically allocated using malloc. At the end of the
computation, we free the memory allocated for the intermediate matrices to avoid memory leaks.

tema3_blas:

This variant uses functions from the BLAS library to optimize matrix operations. BLAS (Basic
Linear Algebra Subprograms) provides optimized routines for basic operations on vectors and
matrices, which can significantly reduce execution time.

Implementation Description:

Memory Allocation:
Four matrices are allocated to store copies and intermediate results: B_copy, ATB, BXA, and result.

Copying Matrix B:
Since B will be modified during operations, we make a copy of it using cblas_dcopy.

Computing the Product AT×B:
We use cblas_dtrmm to multiply AT with B's copy. This function is optimized for multiplying a
triangular matrix (in this case, A is upper triangular) and allows specifying transposition.

Computing the Product B×A:
Similarly to the previous step, we use cblas_dtrmm to multiply B with A, specifying that A is
triangular and does not need to be transposed.

Adding the Two Products:
We use cblas_daxpy to add the results of the two intermediate products element by element.

Computing the Final Product (AT×B+B×A)×BT:
We use cblas_dgemm to multiply the intermediate sum with BT. This function is optimized for
multiplying two matrices and allows specifying the transposition of the input matrices.

Memory Management:
We free the memory allocated for intermediate matrices using free to avoid memory leaks.

tema3_opt_m:

This is a manually optimized version of the brute-force (neopt) implementation. The aim is to
reduce execution time through code-level optimizations, such as using register variables for
quick memory access and efficient iteration using pointers.

Implementation Description:

Transposing Matrices:
A transposed:
Memory is allocated for the transposed matrix AT. Since A is an upper triangular matrix, we
only traverse elements above the main diagonal (inclusive) and copy them to AT at the
corresponding positions.
B transposed:
Similarly, memory is allocated for BT, and all elements from B are copied into BT.

Computing the Product AT×B:
Optimization:
Using register variables for partial sums and iterating through the use of pointers.
Implementation:
Memory is allocated for the result ATXB.
For each element of the resulting matrix, we use pointers (a and b) to traverse the matrices AT and B.
We use register variables to store partial sums, which allows for quick data access and improved performance.
Since AT is a lower triangular matrix, we iterate up to the diagonal element (inclusive).

Computing the Product B×A:
Optimization:
Similar to the previous step, we use pointers and register variables to traverse matrices
B and A, ensuring efficient memory access.
Implementation:
Memory is allocated for the result BXA.
We use pointers and register variables to traverse matrices B and A.

Adding the Two Products:
Optimization:
The sum calculation element by element is optimized by simply adding the two resulting matrices.
Implementation:
Memory is allocated for the sum of matrices ATXB and BXA.
We calculate the sum element by element and store the result in sum_up.

Computing the Final Product (AT×B+B×A)×BT:
Optimization:
Using pointers and register variables optimizes the calculation of the final product.
Implementation:
Memory is allocated for the final result result.
For each element of the resulting matrix, we use pointers and register variables to traverse
the matrices sum_up and BT.

Memory Management:
We ensure the release of memory allocated for intermediate matrices immediately after they are
no longer needed, to minimize memory usage.

Advantages of Optimization:
Using register variables:
Storing partial sums in a register variable allows the processor to quickly access and update the
value, as register variables are stored in the processor's registers, which are much faster than
main memory.

Efficient Iteration with Pointers:
Using pointers to iterate through the matrix allows direct access to matrix elements without needing to compute the index each time. This reduces computational overhead and improves performance.
