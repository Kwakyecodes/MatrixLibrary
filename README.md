# Matrix Library by Kwakyecodes

The Matrix Library is a comprehensive Java library designed to facilitate the creation and manipulation of matrix objects. Offering a wide range of functionalities from basic matrix creation to advanced operations, this library aims to be a powerful tool for developers working on mathematical, engineering, and scientific applications.

## Features

- Creation of zero, identity, and random matrices.
- Support for basic matrix operations (to be detailed further).
- Customizable precision through a delta variable for floating-point comparisons.

## Getting Started

To use the Matrix Library in your project, include the `com.kwakyecodes` package. The library provides several constructors and static methods to create and initialize matrices in various ways.

### Prerequisites

Ensure you have Java installed on your system to use this library.

## Installation

To integrate the Matrix library into your project, you can use Maven. 

### Step 1: Add the JitPack Repository

First, add the JitPack repository to your `pom.xml` file to enable Maven to fetch the library from your GitHub repository. Insert the following snippet into the `<repositories>` section of your `pom.xml`:

```xml
<repositories>
    <repository>
        <id>jitpack.io</id>
        <url>https://jitpack.io</url>
    </repository>
</repositories>
```

### Step 2: Add the Dependency

Next, add the Matrix library as a dependency. Replace `Tag` with the version of the library you wish to use (e.g., the release tag on GitHub like `v1.0`). Alternatively, you can use the commit hash or `main-SNAPSHOT` to use the latest commit on the `main` branch.

```xml
<dependencies>
    <dependency>
        <groupId>com.github.Kwakyecodes</groupId>
        <artifactId>MatrixLibrary</artifactId>
        <version>Tag</version>
    </dependency>
</dependencies>
```

### Step 3: Build Your Project

With the repository and dependency added to your `pom.xml`, you can now build your project. Maven will fetch the Matrix library from your GitHub repository and include it in your project's build path.

Run the following command in your project's root directory:

```shell
mvn clean install
```

This command compiles your project and downloads the necessary dependencies, including the Matrix library.

## Usage

Here are some examples of how to create and work with matrix objects using the Matrix Library.

### Creating a Zero Matrix

```java
Matrix zeroMatrix = new Matrix(3, 4); // Creates a 3x4 matrix filled with zeros
```

### Creating a Matrix from a 2D Array

```java
double[][] data = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
Matrix matrixFrom2DArray = new Matrix(data); // Creates a matrix with the specified entries
```

### Generating an Identity Matrix

```java
Matrix identityMatrix = Matrix.identity(4); // Creates a 4x4 identity matrix
```

### Generating a Random Matrix

```java
Matrix randomMatrix = Matrix.random(3, 3); // Creates a 3x3 matrix with random entries between 0 and 1
```

## Constructors

- `Matrix(int M, int N)`: Constructs a zero matrix with `M` rows and `N` columns.
- `Matrix(double[][] entries)`: Constructs a matrix from a given 2D array, throwing a `RuntimeException` if the array is not rectangular.

## Static Methods

- `public static Matrix identity(int N)`: Returns an `N`x`N` identity matrix.
- `public static Matrix random(int M, int N)`: Returns an `M`x`N` matrix with randomly generated entries.

In addition to creating matrices, the Matrix Library supports various operations to manipulate and compare matrix objects. Here are some key methods and how to use them.

### Retrieving Matrix Dimensions

```java
Matrix matrix = new Matrix(3, 5);
int rows = matrix.getRowSize(); // Returns 3
int columns = matrix.getColumnSize(); // Returns 5
```

### Adding Two Matrices

```java
Matrix A = new Matrix(new double[][]{{1, 2}, {3, 4}});
Matrix B = new Matrix(new double[][]{{5, 6}, {7, 8}});
Matrix sum = A.plus(B); // Computes the sum of A and B
```

### Subtracting Two Matrices

```java
Matrix C = A.minus(B); // Computes the difference between A and B
```

### Checking Matrix Equality

```java
boolean areEqual = A.isEqualTo(B); // Checks if matrices A and B are equal
```

### Computing the Transpose of a Matrix

```java
Matrix transpose = A.transpose(); // Computes the transpose of matrix A
```

### Precision in Comparisons

When checking for equality between matrix entries, a precision delta (`DELTA`) is used to account for the inaccuracy of floating-point representations. This is an adjustable static field that can be set according to the desired precision level.

```java
Matrix.DELTA = 1e-5; // Adjusting the precision delta
```

### Multiplying Two Matrices

```java
Matrix A = new Matrix(new double[][]{{1, 2}, {3, 4}});
Matrix B = new Matrix(new double[][]{{2, 0}, {1, 3}});
Matrix productAB = A.times(B); // Computes the product of matrices A and B
```

### Multiplying a Matrix by a Scalar

```java
double scalar = 5;
Matrix scaledA = A.times(scalar); // Scales matrix A by the scalar value
```

### Multiplying a Matrix by a Vector

```java
double[] vector = {1, 2};
Matrix productVector = A.times(vector); // Computes the product of matrix A and vector
```

### Computing the Determinant of a Matrix

```java
double determinant = A.determinant(); // Computes the determinant of matrix A
```

### Echelon Forms

#### Row Echelon Form

```java
Matrix A = new Matrix(new double[][]{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
Matrix ref = A.rowEchelonForm(); // Computes the row echelon form of matrix A
```

#### Reduced Row Echelon Form

```java
Matrix rref = A.reducedRowEchelonForm(); // Computes the reduced row echelon form of matrix A
```

### Linear Independence

```java
boolean areIndependent = A.areLinearlyIndependent(); // Checks if the vectors of matrix A are linearly independent
```

### Matrix Rank

```java
int rank = A.rank(); // Computes the rank of matrix A
```

### Matrix Nullity

```java
int nullity = A.nullity(); // Computes the nullity of matrix A
```

### Null Space and Span

#### Checking for Null Space

```java
double[] vector = {1, 2, 3};
boolean isInNullSpace = A.isInNullSpace(vector); // Checks if a vector is in the null space of matrix A
```

#### Checking for Span

```java
boolean isSpanned = A.spans(vector); // Checks if a given vector is in the span of the matrix A
```

## Invertibility and Solving Equations

### Checking Invertibility

Determines if a matrix is invertible, which is only true for square matrices with a non-zero determinant.

```java
boolean isInv = A.isInvertible(); // Checks if matrix A is invertible
```

### Inverse of a Matrix

Computes the inverse of an invertible matrix, utilizing the augmented matrix method.

```java
Matrix invA = A.inverse(); // Computes the inverse of matrix A
```

### Solving Linear Equations

Solves the equation \(Ax = b\) for \(x\), supporting both singular and non-singular cases.

```java
double[] b = {1, 2, 3};
Matrix x = A.solve(b); // Finds x in Ax = b
```

## Eigenvalues and Eigenvectors

### Checking for Eigenvalues

Determines if a given scalar is an eigenvalue of the matrix (applicable only to square matrices).

```java
boolean isEigenVal = A.isEigenValue(5); // Checks if 5 is an eigenvalue of A
```

### Checking for Eigenvectors

Verifies if a given vector is an eigenvector of the matrix (applicable only to square matrices).

```java
double[] vector = {1, 2, 3};
boolean isEigenVec = A.isEigenVector(vector); // Checks if vector is an eigenvector of A
```

## Method Descriptions

- `public int getRowSize()`: Returns the number of rows in the matrix.
- `public int getColumnSize()`: Returns the number of columns in the matrix.
- `public Matrix plus(Matrix A)`: Returns the sum of the current matrix and matrix `A`, throwing a `RuntimeException` if the matrices are not of the same size.
- `public Matrix minus(Matrix A)`: Returns the difference between the current matrix and matrix `A`, throwing a `RuntimeException` if the matrices are not of the same size.
- `public boolean isEqualTo(Matrix A)`: Checks if the current matrix is equal to matrix `A` in dimensions and all entries, considering a precision delta for floating-point comparisons.
- `public Matrix transpose()`: Returns the transpose of the current matrix.
- `public Matrix times(Matrix A)`: Returns the product of the current matrix with another matrix `A`, ensuring the number of columns in the current matrix matches the number of rows in `A`.
- `public Matrix times(double scalar)`: Returns a new matrix resulting from the multiplication of every entry of the current matrix by a scalar.
- `public Matrix times(double[] vector)`: Multiplies the matrix by a vector, returning the resultant matrix as a single-column matrix.
- `public double determinant()`: Calculates and returns the determinant of the matrix, throwing a `RuntimeException` if the matrix is not square.
- `public Matrix rowEchelonForm()`: Returns the row echelon form of the matrix.
- `public Matrix reducedRowEchelonForm()`: Returns the reduced row echelon form of the matrix.
- `public boolean areLinearlyIndependent()`: Checks if the vectors within the matrix are linearly independent.
- `public int rank()`: Computes and returns the rank of the matrix.
- `public int nullity()`: Computes and returns the nullity of the matrix.
- `public boolean isInNullSpace(double[] vector)`: Checks if a vector is in the null space of the matrix.
- `public boolean spans(double[] vector)`: Checks if a given vector is in the span of the matrix.
- `public boolean isInvertible()`: Checks if the matrix is invertible.
- `public Matrix inverse()`: Computes the inverse of the matrix if it's invertible.
- `public Matrix solve(double[] rhs)`: Solves \(Ax = b\) for \(x\) using the matrix and a right-hand side vector.
- `public Matrix solve(Matrix rhs)`: Solves \(Ax = b\) for \(x\) using the matrix and a right-hand side matrix (for vector \(b\)).
- `public boolean isEigenValue(double scalar)`: Checks if a given scalar is an eigenvalue of the matrix.
- `public boolean isEigenVector(double[] vector)`: Checks if a given vector is an eigenvector of the matrix.
- `public boolean isEigenVector(Matrix vector)`: Checks if a given matrix (vector) is an eigenvector of the matrix.

## Important Considerations

- **Matrix Multiplication**: The method for matrix multiplication utilizes the transpose of matrix `A` for efficiency, leveraging a `vectorDotProduct` utility method for calculating the dot product of matrix rows and columns. Ensure the column size of the first matrix matches the row size of the second.
- **Scalar and Vector Multiplication**: These operations allow for flexible manipulation of matrices, useful in various computational scenarios.
- **Determinant Calculation**: The determinant is calculated using row operations to convert the matrix into an upper triangular form, then multiplying its diagonal entries. This method only applies to square matrices and employs precision checks for zero values to handle floating-point arithmetic.
- **Echelon Forms**: Echelon forms of a matrix are useful for various matrix analyses, including solving linear systems and determining rank or independence.
- **Linear Independence**: This is determined based on the matrix's echelon form, evaluating if the vectors form a set of linearly independent vectors.
- **Rank and Nullity**: These fundamental concepts are closely tied to the dimensions of the vector space spanned by the matrix's columns (rank) and the dimension of the solution space of the homogeneous system (nullity).
- **Vector Relationships**: Understanding whether vectors are in the null space or span of a matrix is crucial for solving equations and analyzing spaces defined by matrices.
- **Invertibility and Solving Equations**: These functionalities are crucial for many applications, including systems of linear equations and matrix theory. The library provides mechanisms to handle both singular (non-invertible) and non-singular matrices.
- **Eigenvalues and Eigenvectors**: Identifying eigenvalues and eigenvectors is essential in various fields, such as differential equations, quantum mechanics, and vibration analysis. The library supports evaluations based on the determinant and algebraic properties of matrices.