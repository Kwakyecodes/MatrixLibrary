package com.kwakyecodes;

import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;

/**
 * Unit test for the Matrix class.
 */
public class MatrixTest {
    
    private Matrix zeroMat;
    private Matrix matrix;
    private Matrix identityMat;
    private Matrix smallIdenMatrix;
    private Matrix randomMat;
    private Matrix randomVector;
    private double[] eightZeros = {0, 0, 0, 0, 0, 0, 0, 0};
    private double[] smallVector = {4, 5};
    private Matrix eightZerosMat;
    private Matrix hundredZerosMat;

    @Before
    public void setUp() {
        zeroMat = new Matrix(500, 100);
        double[][] entries = {
            {12, 3, 0, 0, 0, 0, 0, 0},
            {-5, 1, 0, 0, 0, 0, 0, 0},
            {0, 0, 4, 5, 2, 0, 0, 0},
            {0, 0, -1, 0, -2, 0, 0, 0},
            {0, 0, 3, 5, -1, 0, 0, 0},
            {0, 0, 0, 0, 0, 1, 2, 3},
            {0, 0, 0, 0, 0, -2, -1, 5},
            {0, 0, 0, 0, 0, 0, -2, 0}
        }; 
        matrix = new Matrix(entries);
        identityMat = Matrix.identity(100);
        randomMat = Matrix.random(500, 100);
        eightZerosMat = new Matrix(8, 1);
        randomVector = Matrix.random(100, 1);
        System.out.println(randomVector.toString());
        hundredZerosMat = new Matrix(100, 1);
        smallIdenMatrix = Matrix.identity(2);
    }

    @Test
    public void testValidMatrixCreation() {
        // Test creating a zero matrix with given dimensions
        assertTrue(zeroMat instanceof Matrix);
        assertEquals(zeroMat.getRowSize(), 500);
        assertEquals(zeroMat.getColumnSize(), 100);

        // Test creating a matrix with 2d array of doubles
        assertTrue(matrix instanceof Matrix);
        assertEquals(matrix.getRowSize(), 8);
        assertEquals(matrix.getColumnSize(), 8);

        // Test creating an identity matrix with given dimensions
        assertTrue(identityMat instanceof Matrix);
        assertEquals(identityMat.getRowSize(), 100);
        assertEquals(identityMat.getColumnSize(), 100);

        // Test creating a random matrix with given dimensions
        assertTrue(randomMat instanceof Matrix);
        assertEquals(randomMat.getRowSize(), 500);
        assertEquals(randomMat.getColumnSize(), 100);
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidMatrixCreation1() {
        @SuppressWarnings("unused")
        Matrix temp = new Matrix(-5, 9);
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidMatrixCreation2() {
        @SuppressWarnings("unused")
        Matrix temp = new Matrix(1000, 0);
    }

    @Test
    public void testValidMatrixAddition() {
        assertTrue(randomMat.plus(randomMat).isEqualTo(randomMat.times(2)));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidMatrixAddition() {
        identityMat.plus(zeroMat);
    }

    @Test
    public void testValidMatrixSubtraction() {
        assertTrue(randomMat.minus(randomMat).isEqualTo(zeroMat));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidMatrixSubtraction() {
        identityMat.minus(zeroMat);
    }

    @Test
    public void testTranspose() {
        assertTrue(identityMat.transpose().isEqualTo(identityMat));
    }

    @Test
    public void testValidMatrixMatrixMultiplication() {
        assertTrue(randomMat.times(identityMat).isEqualTo(randomMat));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidMatrixMatrixMultiplication() {
        identityMat.times(randomMat);
    }

    @Test
    public void testValidMatrixVectorMultiplication() {
        assertTrue(matrix.times(eightZeros).isEqualTo(eightZerosMat));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidMatrixVectorMultiplication() {
        randomMat.times(eightZeros);
    }

    @Test
    public void testValidPrint() {
        assertEquals(eightZerosMat.toString(), "       0.000000 \n       0.000000 \n       0.000000 \n       0.000000 \n       0.000000 \n       0.000000 \n       0.000000 \n       0.000000 \n");
    }

    @Test
    public void testValidDeterminant() {
        assertEquals(matrix.determinant(), -2970, Matrix.DELTA);
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidDeterminant() {
        randomMat.determinant();
    }

    @Test
    public void testValidRowEchelonForm() {
        assertTrue(matrix.reducedRowEchelonForm().isEqualTo(Matrix.identity(8)));
    }

    @Test
    public void testAreLinearlyIndependent() {
        assertTrue(matrix.areLinearlyIndependent());
        assertFalse(zeroMat.areLinearlyIndependent());
    }

    @Test
    public void testRank() {
        assertEquals(matrix.rank(), 8);
        assertEquals(identityMat.rank(), 100);
    }

    @Test
    public void testNullity() {
        assertEquals(matrix.nullity(), 0);
        assertEquals(identityMat.nullity(), 0);
    }

    @Test
    public void testIsInNullSpace() {
        assertTrue(matrix.isInNullSpace(eightZerosMat));
        assertTrue(matrix.isInNullSpace(eightZeros));
        assertFalse(matrix.isInNullSpace(randomMat));
    }

    @Test
    public void testSpans() {
        assertTrue(matrix.spans(Matrix.random(8, 1)));
        assertTrue(matrix.spans(eightZeros));
        assertFalse(matrix.spans(randomMat));
    }

    @Test
    public void testIsInvertible() {
        assertTrue(matrix.isInvertible());
        assertFalse(randomMat.isInvertible());
    }

    @Test
    public void testValidinverse() {
        assertTrue(identityMat.inverse().isEqualTo(identityMat));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidInverse() {
        randomMat.inverse();
    }

    @Test
    public void testValidSolve() {
        assertTrue(matrix.solve(eightZeros).isEqualTo(eightZerosMat));
        assertTrue(matrix.solve(eightZerosMat).isEqualTo(eightZerosMat));
    }

    @Test
    public void testValidIsEigenValue() {
        assertTrue(identityMat.isEigenValue(1));
        assertFalse(identityMat.isEigenValue(2));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidIsEigenValue() {
        randomMat.isEigenValue(5);
    }

    @Test
    public void testValidIsEigenVector() {
        assertFalse(identityMat.isEigenVector(eightZerosMat));
        assertFalse(identityMat.isEigenVector(hundredZerosMat));
        assertTrue(identityMat.isEigenVector(randomVector));
        assertTrue(smallIdenMatrix.isEigenVector(smallVector));
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidIsEigenVector() {
        randomMat.isEigenVector(hundredZerosMat);
    }
}