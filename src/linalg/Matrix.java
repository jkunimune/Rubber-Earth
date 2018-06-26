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
	 * Add two Matrices.
	 * @param that - The addend
	 * @return the sum of this and that
	 */
	public Matrix plus(Matrix that) {
		if (this.getN() != that.getN() || this.getM() != that.getM())
			throw new IllegalArgumentException("Cannot sum these matrices; the dimensions "+this.getN()+"x"+this.getM()+" and "+that.getN()+"x"+that.getM()+" do not match.");
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
		return this.values[0].length;
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
}
