package com.kwakyecodes;

final public class Matrix {
    private final int M; // Number of rows of matrix
    private final int N; // Number of columns of matrix
    private final double[][] entries; // Elements of Matrix
    public static double DELTA = 1e-9;

    /* Zero matrix constructor */
    public Matrix(int M, int N) {
        if (M < 1 || N < 1)
            throw new RuntimeException("Invalid matrix size");

        this.M = M;
        this.N = N;
        entries = new double[M][N];
    }

    /* 2d-array Matrix Constructor */
    public Matrix(double[][] entries) {
        this(entries.length, entries[0].length);

        // Assign elements to Matrix data
        for (int i = 0; i < M; i++) {
            // Check if 2d-array is a valid matrix
            if (entries[i].length != N)
                throw new RuntimeException("2d-array is not a valid matrix.");
            
            for (int j = 0; j < N; j++) {
                this.entries[i][j] = entries[i][j];
            }
        }
    }

    /* Identity matrix */
    public static Matrix identity(int N) {
        Matrix A = new Matrix(N, N);
        for (int i = 0; i < N; i++) {
            A.entries[i][i] = 1;
        }
        return A;
    }

    /* Random matrix */
    public static Matrix random(int M, int N) {
        Matrix A = new Matrix(M, N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++){
                A.entries[i][j] = Math.random();
            }
        }
        return A;
    }

    /* Method that returns the number of rows of the Matrix */
    public int getRowSize() {
        return M;
    }

    /* Method that returns the number of columns of the Matrix */
    public int getColumnSize() {
        return N;
    }

    /* Method that computes the sum of two matrices */
    public Matrix plus(Matrix A) {
        if (M != A.M || N != A.N) { // Check if size of matrices are equal
            throw new RuntimeException("The size of the matrices must be equal to perform this operation.");
        }
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                C.entries[i][j] = entries[i][j] + A.entries[i][j];
            }
        }
        return C;
    }

    /* Method that computes the difference of two matrices */
    public Matrix minus(Matrix A) {
        if (M != A.M || N != A.N) { // Check if size of matrices are equal
            throw new RuntimeException("The size of the matrices must be equal to perform this operation.");
        }
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++) {
                C.entries[i][j] = entries[i][j] - A.entries[i][j];
            }
        }
        return C;
    }

    /* Method that checks if two matrices are equal */
    public boolean isEqualTo(Matrix A) {
        if (M != A.M || N != A.N)
            return false;
        for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++) {
                if (!compareDoubles(entries[i][j], A.entries[i][j]))
                    return false;
            }
        }
        return true;
    }

    /* Method that computes the transpose of a matrix */
    public Matrix transpose() {
        Matrix C = new Matrix(N, M);
        for (int i = 0; i < C.M; i++) {
            for (int j = 0; j < C.N; j++)
                C.entries[i][j] = entries[j][i];
        }
        return C;
    }

    /* Method that computes the product of two matrices */
    public Matrix times(Matrix A) {
        if (N != A.M) 
            throw new RuntimeException("Operation cannot be performed due to column-row size mismatch.");

        Matrix C = new Matrix(M, A.N), aT = A.transpose();
        for (int i = 0; i < C.M; i++){
            for (int j = 0; j < C.N; j++) {
                C.entries[i][j] = vectorDotProduct(entries[i], aT.entries[j]);
            }
        }
        return C;
    }

    /* Method that computes the product of a matrix and a scalar */
    public Matrix times(double scalar) {
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++) {
                C.entries[i][j] = entries[i][j] * scalar;
            }
        }
        return C;
    }

    /* Method that computes the product of a matrix and a vector array */
    public Matrix times(double[] vector) {
        if (M != vector.length)
            throw new RuntimeException("Operation cannot be performed due to column-row size mismatch.");

        Matrix C = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            C.entries[i][0] = vectorDotProduct(entries[i], vector);
        }
        return C;
    }
 
    /* Method that returns a string representation of the matrix object */
    public String toString() {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < this.M; i++) {
            for (int j = 0; j < this.N; j++) {
                builder.append(String.format("%15f ", this.entries[i][j]));
            }
            builder.append("\n"); // Add a new line after each row
        }
        return builder.toString();
    }


    /* Method to compute the determinant of the matrix object */
    public double determinant() {
        if (N != M) // Check if the matrix is a squre matrix
            throw new RuntimeException("The determinant of non-square matrix is not defined.");
        
        Matrix temp  = new Matrix(entries);
        double D = 1;

        // Perform row operations to put matrix in triangular form
        // And multiply diagonals to calculate the determinant
        for (int i = 0; i < temp.M; i++){ // Iterate columns
            for (int j = i; j < temp.M; j++) { // Iterate rows
                if (i == j){
                    if (!compareDoubles(temp.entries[j][i], 0)) {
                        D *= temp.entries[j][i];
                        continue;
                    }
                    else {
                        for (int k = j+1; k < temp.M; k++) {
                            if (!compareDoubles(temp.entries[k][i], 0)) { // Swap this row with the current row
                                temp.swapRows(j, k);
                                break;
                            }
                        }
                        // If current row is not pivot row, matrix is not full rank, and so the determinant is zero
                        if (compareDoubles(temp.entries[j][i], 0)) return 0;
                        D *= -temp.entries[j][i];
                    }
                }
                else {
                    if (compareDoubles(temp.entries[j][i], 0)) continue;
                    else {
                        temp.addToRow(j, i, -temp.entries[j][i]/temp.entries[i][i]);
                    }
                }
            }
        }
        return D;
    }

    /* Method to computes and returns the echelon form of a matrix */
    public Matrix rowEchelonForm() {
        Matrix A = new Matrix(entries);
        A.rowOperationsForward();
        return A;
    }

    /* Method to computes and returns the reduced echelon form of a matrix */
    public Matrix reducedRowEchelonForm() {
        Matrix A = new Matrix(entries);
        A.rowOperationsForward();

        for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++) {
                if (!compareDoubles(A.entries[i][j], 0)) {
                    A.multiplyRow(i, 1/A.entries[i][j]);
                    break;
                }
            }
        }

        A.rowOperationsBackward();
        return A;
    }

    /* Method that checks if the vectors in a matrix are linearly independent or not */
    public boolean areLinearlyIndependent() {
        // If there are more columns than rows, then vectors are linearly dependent
        if (N > M) return false;

        // If the number of pivot rows/columns is less than the number of columns, 
        // then vectors are linearly dependent
        if (numberOfPivotRows()[0] < N) return false;

        return true;
    }

    /* Method that computes the rank of a matrix */
    public int rank() {
        return numberOfPivotRows()[0];
    }
    
    /* Method that computes the nullity of the matrix */
    public int nullity() {
        return N - numberOfPivotRows()[0];
    }

    /* Method that checks of a vector is in the null space of a matrix */
    public boolean isInNullSpace(double[] vector) {
        if (vector.length != N)
            return false; 

        Matrix res = times(vector);
        return isZeroVector(res);
    }

    /* Method that checks of a vector is in the null space of a matrix */
    public boolean isInNullSpace(Matrix vector) {
        if (vector.M != N || vector.N != 1)
            return false;

        Matrix res = times(vector);
        return isZeroVector(res);
    }

    /* Method that checks if a given vector is in the span of the matrix. */
    public boolean spans(double[] vector) {
        if (vector.length != M)
            return false;

        Matrix A = getAugmentedMatrix(vector);
        return A.isConsistent();
    }

    /* Method that checks if a given vector is in the span of the matrix. */
    public boolean spans(Matrix vector) {
        if (vector.M != M || vector.N != 1)
            return false;

        Matrix A = getAugmentedMatrix(vector);
        return A.isConsistent();
    }

    /* Method that checks if a given matrix is invertible or not */
    public boolean isInvertible() {
        if (M != N) // Non-square matrices are non-invertible
            return false;

        return !compareDoubles(determinant(), 0);
    }

    /* Method that computes the inverse of a given matrix */
    public Matrix inverse() {
        if (!isInvertible())
            throw new RuntimeException("Matrix is non-invertible.");

        Matrix A = new Matrix(entries); // Make a copy of current matrix
        Matrix augMat = A.getAugmentedMatrix(identity(N));
        Matrix augMatReduced = augMat.reducedRowEchelonForm();
        Matrix res = new Matrix(N, N);
        for (int i = 0; i < res.M; i++) {
            for (int j = 0; j < res.N; j++) {
                res.entries[i][j] = augMatReduced.entries[i][N + j];
            }
        }
        return res;
    }

    /* Method that computes the solution x of Ax = b */
    public Matrix solve(double[] rhs) {
        if (!spans(rhs)) {
            System.out.println("There are no solutions.");
            System.exit(0);
        }
        
        if (!isInvertible()) { // Return the general solution in this case
            System.out.println("There are infinitely many solutions.");
            return getAugmentedMatrix(rhs).reducedRowEchelonForm();
        }

        return inverse().times(rhs);
    }

    /* Method that computes the solution x of Ax = b */
    public Matrix solve(Matrix rhs) {
        if (!spans(rhs)) {
            System.out.println("There are no solutions.");
            System.exit(0);
        }
        
        if (!isInvertible()) { // Return the general solution in this case
            System.out.println("There are infinitely many solutions.");
            return getAugmentedMatrix(rhs).reducedRowEchelonForm();
        }

        return inverse().times(rhs);
    }

    /* Method to check if a given scalar is an eigenvalue of the matrix */
    public boolean isEigenValue(double scalar) {
        if (M != N)
            throw new RuntimeException("Non-square matrices do not have eigenvalues.");

        // Compute A - lambda*I
        Matrix C = minus(identity(N).times(scalar));

        // scalar is eigen value if determinant of C is 0
        return compareDoubles(C.determinant(), 0);
    }

    /* Method that checks if a given vector is an eigenvector of the matrix */
    public boolean isEigenVector(double[] vector) {
        if (M != N)
            throw new RuntimeException("Non-square matrices do not have eigenvectors.");

        if (vector.length != N)
            return false;

        if (isZeroVector(vector))
            return false;

        // Compute A*v
        Matrix C = times(vector);
        if (isZeroVector(C)) // If A*v is a zero vector, then eigenvalue is 0
            return true;

        double lambda = getMultipleFactor(vector, C);
        return !compareDoubles(lambda, 0);
    }

    /* Method that checks if a given vector is an eigenvector of the matrix */
    public boolean isEigenVector(Matrix vector) {
        if (M != N)
        throw new RuntimeException("Non-square matrices do not have eigenvectors.");

        if (vector.M != M && vector.N != 1)
            return false;
        
        if (isZeroVector(vector))
            return false;

        double lambda = getMultipleFactor(vector, times(vector));
        return !compareDoubles(lambda, 0);
    }

    /* Method that checks if a vector is a zero vector */
    private static boolean isZeroVector(double[] vector) {
        for (double ele: vector){
            if (!compareDoubles(ele, 0)) 
                return false;
        }
        return true;
    }

    /* Method that checks if a vector is a zero vector */
    private static boolean isZeroVector(Matrix vector) {
        for (int i = 0; i < vector.M; i++){
            if (!compareDoubles(vector.entries[i][0], 0))
                return false;
        }
        return true;
    }

    /* Method that computes the factor of two vectors if they are multiples.
     * Otherwise the function just returns 0.0 */
    private static double getMultipleFactor(double[] vec1, Matrix vec2) {
        double factor = 0;
        for (int i = 0; i < vec1.length; i++){
            if (compareDoubles(vec2.entries[i][0], 0) && compareDoubles(vec1[i], 0))
                continue;
            if (compareDoubles(vec2.entries[i][0], 0) || compareDoubles(vec1[i], 0))
                return 0;
            
            if (compareDoubles(factor, 0))
                factor = vec2.entries[i][0] / vec1[i];
            else {
                if (!compareDoubles(vec2.entries[i][0]/vec1[i], factor))
                    return 0;
            }
        }
        return factor;
    }

    /* Method that computes the factor of two vectors if they are multiples.
     * Otherwise the function just returns 0.0 */
    private static double getMultipleFactor(Matrix vec1, Matrix vec2) {
        double factor = 0;
        for (int i = 0; i < vec1.M; i++){
            if (compareDoubles(vec2.entries[i][0], 0) && compareDoubles(vec1.entries[i][0], 0))
                continue;
            if (compareDoubles(vec2.entries[i][0], 0) || compareDoubles(vec1.entries[i][0], 0))
                return 0;
            
            if (compareDoubles(factor, 0))
                factor = vec2.entries[i][0] / vec1.entries[i][0];
            else {
                if (!compareDoubles(vec2.entries[i][0]/vec1.entries[i][0], factor))
                    return 0;
            }
        }
        return factor;
    }

    /* Method that returns the augmented matrix of a matrix and a vector */
    private Matrix getAugmentedMatrix(double[] vector) {
        Matrix A = new Matrix(M, N+1);
        for (int i = 0; i < A.M; i++){
            for (int j = 0; j < A.N; j++) {
                A.entries[i][j] = (j < N)? entries[i][j]:vector[i];
            }
        }
        return A;
    }

    /* Method that returns the augmented matrix of two matrices */
    private Matrix getAugmentedMatrix(Matrix C) {
        Matrix A = new Matrix(M, N + C.N);
        for (int i = 0; i < A.M; i++) {
            for (int j = 0; j < A.N; j++) {
                A.entries[i][j] = (j < N)? entries[i][j]:C.entries[i][j - N];
            }
        }
        return A;
    }

    /* Method that checks if a matrix is consistent or inconsistent */
    private boolean isConsistent() {
        Matrix A = new Matrix(entries);
        A.rowOperationsForward();

        for (int i = A.M-1; i >= 0; i--) {
            for (int j = 0; j < A.N; j++) {
                if (j < A.N-1 && !compareDoubles(A.entries[i][j], 0))
                    break;
                // Return false if last column is a pivot column
                if (j == A.N-1 && !compareDoubles(A.entries[i][j], 0))
                    return false;
            }
        }
        return true;
    }

    /* Method that finds the number of pivots rows in a matrix 
     * In addition to that, it returns the index of the last pivot columns plus 1 */
    private int[] numberOfPivotRows() {
        Matrix A = new Matrix(entries);
        A.rowOperationsForward();

        int[] pivotRows = {0, 0};
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (!compareDoubles(A.entries[i][j], 0)) {
                    pivotRows[0]++;
                    pivotRows[1] = j;
                    break;
                }
            }
        }
        return pivotRows;
    }

    /* Method that performs row operations on matrix to put it in row echelon form */
    private void rowOperationsForward() {
        for (int i = 0; i < N; i++){ // Iterate columns
            for (int j = i; j < M; j++) { // Iterate rows
                if (i == j){
                    if (compareDoubles(entries[j][i], 0)) { 
                        // Find a row which has a pivot element in this column
                        for (int k = j+1; k < M; k++) {
                            if (!compareDoubles(entries[k][i], 0)) { // Swap this row with the current row
                                swapRows(j, k);
                                break;
                            }
                        }
                    }
                }
                else {
                    if (!compareDoubles(entries[j][i], 0))
                        addToRow(j, i, -entries[j][i]/entries[i][i]);
                }
            }
        }
    }

    /* Method that performs row operations to reduce a matrix that is already in row echelon form */
    private void rowOperationsBackward() {
        int[] pivotsCount = numberOfPivotRows();
        for (int i = pivotsCount[1]; i >= 0; i--) { // Iterate columns
            for (int j = pivotsCount[0]-1; j >= 0; j--) { // Iterate rows
                if (j == pivotsCount[0]-1 ){
                    if (!compareDoubles(entries[j][i], 1)) break; // Skip non-pivot columns
                    else continue;
                }
                if (!compareDoubles(entries[j][i], 0))
                    addToRow(j, pivotsCount[0]-1, -entries[j][i]/entries[pivotsCount[0]-1][i]);
                if (j == 0) pivotsCount[0]--;
            }
        }
    }

    /* Method that swaps two rows */
    private void swapRows(int row1, int row2) {
        double[] temp = entries[row1];
        entries[row1] = entries[row2];
        entries[row2] = temp;
    }

    /* Method that multiplies a given row by a constant */
    private void multiplyRow(int row, double k) {
        for (int i = 0; i < N; i++){
            entries[row][i] = k * entries[row][i];
        }
    }

    /* Method that adds a multiple of row2 to row1 to another */
    private void addToRow(int row1, int row2, double k) {
        for (int i = 0; i < N; i++){
            entries[row1][i] = entries[row1][i] + (k * entries[row2][i]);
        }
    }

    /* Method that computes the dot product of two vectors/arrays */
    private static double vectorDotProduct(double[] A, double[] B) {
        int size = A.length;
        double res = 0;
        for (int i = 0; i < size; i++)
            res += A[i] * B[i];
        return res;
    }

    /* Method that compares the values of two decimals */
    private static boolean compareDoubles(double d1, double d2) {
        return Math.abs(d1 - d2) < DELTA;
    }
}