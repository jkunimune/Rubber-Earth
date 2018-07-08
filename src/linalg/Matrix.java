/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package linalg;

import java.util.Arrays;

/**
 * A two-dimensional array of numbers.
 * 
 * @author Justin Kunimune
 */
public class Matrix {
	
	protected final double[][] values;
	
	
	/**
	 * Instantiate an nxm matrix with all zeroes.
	 * @param n - The height
	 * @param m - The width
	 */
	public Matrix(int n, int m) {
		this.values = new double[n][m];
	}
	
	/**
	 * Instantiate a matrix based on an existing 2D array
	 * @param ds
	 */
	public Matrix(double[][] values) {
		this.values = values;
	}
	
	/**
	 * Instantiate a matrix based on a 1D array given dimensiosn
	 * @param n - The height
	 * @param m - The width
	 * @param values - The values from left to right then top to bottom
	 */
	public Matrix(int n, int m, double... values) {
		if (values.length != n*m)
			throw new IllegalArgumentException(values.length+" values do not fit in a "+n+"x"+m+" matrix.");
		this.values = new double[n][m];
		for (int i = 0; i < values.length; i ++)
			this.set(i/m, i%m, values[i]);
	}
	
	/**
	 * Instantiate a zero Matrix of size nxm.
	 * Alias of Matrix(int, int)
	 * @param n - The height
	 * @param m - The width
	 * @return I_n
	 */
	public static Matrix zeroes(int n, int m) {
		return new Matrix(n, m);
	}
	
	/**
	 * Instantiate an identity Matrix of size nxn.
	 * @param n - The dimension
	 * @return I_n
	 */
	public static Matrix identity(int n) {
		Matrix I = new Matrix(n, n);
		for (int i = 0; i < n; i ++)
			I.set(i, i, 1);
		return I;
	}
	
	/**
	 * Vertically concatenate two Matrices of similar widths.
	 * @param matrices - The matrices to be stacked
	 * @return the concatenation
	 */
	public static Matrix vertcat(Matrix... matrices) {
		int n = 0;
		for (Matrix mat: matrices) {
			if (mat.getM() != matrices[0].getM())
				throw new IllegalArgumentException("Cannot compute this determinant. Matrices "+Arrays.toString(matrices)+" do not have matching widths.");
			n += mat.getN();
		}
		Matrix cat = new Matrix(n, matrices[0].getM());
		n = 0;
		for (Matrix mat: matrices) {
			for (int i = 0; i < mat.getN(); i ++)
				for (int j = 0; j < mat.getM(); j ++)
					cat.set(i+n, j, mat.get(i, j));
			n += mat.getN();
		}
		return cat;
	}
	
	/**
	 * Horizontally concatenate two Matrices of similar heights.
	 * @param matrices - The matrices to be concatenated
	 * @return the concatenation
	 */
	public static Matrix horzcat(Matrix... matrices) {
		int m = 0;
		for (Matrix mat: matrices) {
			if (mat.getN() != matrices[0].getN())
				throw new IllegalArgumentException("Cannot compute this determinant. Matrices "+Arrays.toString(matrices)+" do not have matching heights.");
			m += mat.getM();
		}
		Matrix cat = new Matrix(matrices[0].getN(), m);
		m = 0;
		for (Matrix mat: matrices) {
			for (int i = 0; i < mat.getN(); i ++)
				for (int j = 0; j < mat.getM(); j ++)
					cat.set(i, j+m, mat.get(i, j));
			m += mat.getM();
		}
		return cat;
	}
	
	
	/**
	 * Compute the transpose.
	 * @return the transpose of this
	 */
	public Matrix T() {
		Matrix tp = new Matrix(this.getN(), this.getM());
		for (int i = 0; i < tp.getN(); i ++)
			for (int j = 0; j < tp.getM(); j ++)
				tp.set(i, j, this.get(j, i));
		return tp;
	}
	
	/**
	 * Compute a version of this where each column has magnitude one
	 * @return
	 */
	public Matrix norm() {
		Matrix hat = new Matrix(this.getN(), this.getM());
		for (int j = 0; j < this.getM(); j ++) {
			double colSum = 0;
			for (int i = 0; i < this.getN(); i ++)
				colSum += this.get(i, j)*this.get(i, j);
			if (colSum != 0)
				for (int i = 0; i < this.getN(); i ++)
					hat.set(i, j, this.get(i, j)/Math.sqrt(colSum));
			else
				for (int i = 0; i < this.getN(); i ++) // for zero columns,
					hat.set(i, j, i == 0 ? 1 : 0); // default to rightward
		}
		return hat;
	}
	
	/**
	 * Compute the determinant.
	 * @return the determinant of this
	 */
	public double det() {
		if (this.getN() != this.getM())
			throw new IllegalArgumentException("Cannot compute this determinant. Dimensions "+getN()+"x"+getM()+" are not square.");
		if (this.getN() == 2)
			return this.get(0, 0)*this.get(1, 1) - this.get(0, 1)*this.get(1, 0);
		else
			throw new IllegalArgumentException("I haven't told it how to do those determinants yet.");
	}
	
	/**
	 * Compute the trace.
	 * @return the trace of this
	 */
	public double tr() {
		double tr = 0;
		for (int i = 0; i < Math.min(this.getN(), this.getM()); i ++)
			tr += this.get(i, i);
		return tr;
	}
	
	/**
	 * Compute the square of the magnitude of this matrix.
	 * @return the sum of all the squares
	 */
	public double sqr() {
		double mag2 = 0;
		for (int i = 0; i < this.getN(); i ++)
			for (int j = 0; j < this.getM(); j ++)
				mag2 += this.get(i, j)*this.get(i, j);
		return mag2;
	}
	
	/**
	 * Compute the magnitude of this matrix.
	 * @return the square root of the sum of all the squares
	 */
	public double mag() {
		return Math.sqrt(this.sqr());
	}
	
	/**
	 * Multiply this Matrix by a scalar.
	 * @param a - The factor
	 * @return the product
	 */
	public Matrix times(Matrix that) {
		if (this.getM() != that.getN())
			throw new IllegalArgumentException("Cannot multiply these matrices. Dimensions "+this.getN()+"x"+this.getM()+" and "+that.getN()+"x"+that.getM()+" do not agree.");
		Matrix prod = new Matrix(this.getN(), that.getM());
		for (int i = 0; i < this.getN(); i ++)
			for (int j = 0; j < that.getM(); j ++)
				for (int k = 0; k < this.getM(); k ++)
					prod.set(i,j,prod.get(i,j) + this.get(i, k)*that.get(k, j));
		return prod;
	}
	
	/**
	 * Multiply this Matrix by a scalar.
	 * @param a - The factor
	 * @return the product
	 */
	public Matrix times(double a) {
		Matrix product = new Matrix(this.getN(), this.getM());
		for (int i = 0; i < this.getN(); i ++)
			for (int j = 0; j < this.getM(); j ++)
				product.set(i, j, this.get(i, j) * a);
		return product;
	}
	
	/**
	 * Divide this Matrix by a scalar.
	 * @param a - The divisor
	 * @return the quotient
	 */
	public Matrix over(double a) {
		return this.times(1./a);
	}
	
	/**
	 * Add two Matrices.
	 * @param that - The addend
	 * @return the sum of this and that
	 */
	public Matrix plus(Matrix that) {
		if (this.getN() != that.getN() || this.getM() != that.getM())
			throw new IllegalArgumentException("Cannot add these matrices. Dimensions "+this.getN()+"x"+this.getM()+" and "+that.getN()+"x"+that.getM()+" do not match.");
		Matrix sum = new Matrix(this.getN(), this.getM());
		for (int i = 0; i < this.getN(); i ++)
			for (int j = 0; j < this.getM(); j ++)
				sum.set(i, j, this.get(i, j) + that.get(i, j));
		return sum;
	}
	
	/**
	 * Subtract two Matrices.
	 * @param that - The subtractee
	 * @return the difference of this and that
	 */
	public Matrix minus(Matrix that) {
		return this.plus(that.times(-1));
	}
	
	/**
	 * Am I NaN?
	 * @return whether any of the arguments are NaN
	 */
	public boolean isNaN() {
		for (int i = 0; i < this.getN(); i ++)
			for (int j = 0; j < this.getM(); j ++)
				if (Double.isNaN(this.get(i, j)))
					return true;
		return false;
	}
	
	/**
	 * Compute the height.
	 * @return the height of this
	 */
	public int getN() {
		return this.values.length;
	}
	
	/**
	 * Compute the width.
	 * @return the width of this
	 */
	public int getM() {
		if (this.values.length > 0)
			return this.values[0].length;
		else
			return 0;
	}
	
	/**
	 * Extract a single scalar value.
	 * @param i
	 * @param j
	 * @return the value this_{i,j}
	 */
	public double get(int i, int j) {
		return this.values[i][j];
	}
	
	/**
	 * Set a single scalar value.
	 * @param i
	 * @param j
	 * @param a
	 */
	public void set(int i, int j, double a) {
		this.values[i][j] = a;
	}
	
	/**
	 * Add to a single scalar value.
	 * @param i
	 * @param j
	 * @param a
	 */
	public void add(int i, int j, double a) {
		this.values[i][j] += a;
	}
	
	public String toString() {
		String str = "Matrix("+this.getN()+", "+this.getM()+",\n  ";
		for (int i = 0; i < this.getN(); i ++) {
			for (int j = 0; j < this.getM(); j ++)
				str += String.format("%- 6.4f, ", this.get(i, j));
			str += "\n  ";
		}
		return str.substring(0, str.length()-5) + ")";
	}
}
