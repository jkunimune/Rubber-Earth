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
 * A matrix that happens to be square.
 * 
 * @author Justin Kunimune
 */
public class SquareMatrix extends Matrix {
	
	/**
	 * Instantiate an nxn Matrix with all zeroes.
	 * @param n - The dimension
	 */
	public SquareMatrix(int n) {
		super(n, n);
	}
	
	/**
	 * Convert a square Matrix to a SquareMatrix
	 * @param mat - The nxn Matrix
	 */
	public SquareMatrix(Matrix mat) {
		super(mat.getN(), mat.getN());
		if (mat.getN() != mat.getM())
			throw new IllegalArgumentException("Matrix of dimensions "+mat.getN()+"x"+mat.getM()+" cannot be cast to a SquareMatrix.");
		for (int i = 0; i < mat.getN(); i ++)
			for (int j = 0; j < mat.getM(); j ++)
				this.set(i, j, mat.get(i, j));
	}
	
	/**
	 * Instantiate an identity Matrix of size nxn.
	 * @param n - The dimension
	 * @return I_n
	 */
	public static SquareMatrix identity(int n) {
		SquareMatrix I = new SquareMatrix(n);
		for (int i = 0; i < n; i ++)
			I.set(i, i, 1);
		return I;
	}
	
	/**
	 * Compute the determinant.
	 * @return the determinant of this
	 */
	public double det() {
		//TODO
		return 0;
	}
	
	/**
	 * Compute the trace.
	 * @return the trace of this
	 */
	public double tr() {
		//TODO
		return 0;
	}
	
	@Override
	public SquareMatrix plus(Matrix that) {
		return new SquareMatrix(super.plus(that));
	}
	
	@Override
	public SquareMatrix minus(Matrix that) {
		return new SquareMatrix(super.minus(that));
	}
	
	@Override
	public SquareMatrix times(double a) {
		return new SquareMatrix(super.times(a));
	}
	
	@Override
	public SquareMatrix over(double a) {
		return new SquareMatrix(super.over(a));
	}
}
