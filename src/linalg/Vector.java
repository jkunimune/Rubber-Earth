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
 * A vertical one-dimensional array of numbers.
 * 
 * @author Justin Kunimune
 */
public class Vector extends Matrix {
	
	/**
	 * Instantiate an nx1 matrix with all zeroes.
	 * @param n - The dimension
	 */
	public Vector(int n) {
		super(n, 1);
	}
	
	
	/**
	 * Extract a single scalar value.
	 * @param i
	 * @return the value this_{i}
	 */
	public double get(int i) {
		return this.values[i][0];
	}
	
	/**
	 * Set a single scalar value.
	 * @param i
	 * @param a
	 */
	public void set(int i, double a) {
		this.values[i][0] = a;
	}
}
