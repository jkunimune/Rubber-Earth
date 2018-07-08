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
package model;


/**
 * A simple interface to make these generics easier with which to deal.
 * Represents all of the mesh vertices from the same point on the globe.
 * 
 * @author Justin Kunimune
 */
public class Cell {
	
	public static final int NORTHEAST = 0, NORTHWEST = 1, SOUTHWEST = 2, SOUTHEAST = 3;
	
	private final Vertex[] corners = new Vertex[4]; //which vertex is attached to each sector (there must be exactly one)
	private final double lambda, mu; // the elastic properties
	private final double delP, delL; // the latitudinal and longitudinal spans
	private double energy; // the stored elastic potential energy
	
	
	
	public Cell(double lambda, double mu, double delP, double delL, Vertex ne, Vertex nw, Vertex sw, Vertex se) {
		this.lambda = lambda;
		this.mu = mu;
		this.delP = delP;
		this.delL = delL;
		corners[NORTHEAST] = ne;
		corners[NORTHWEST] = nw;
		corners[SOUTHWEST] = sw;
		corners[SOUTHEAST] = se;
	}
	
	
	/**
	 * Compute the strain energy density times volume for this cell, save it and return it.
	 */
	double computeEnergy() {
		this.energy = 0*delP*delL;
		return energy;
	}
	
	
	public Vertex getCorner(int direction) {
		return corners[direction];
	}
	
	
	@Override
	public String toString() {
		String s = "Cell(";
		for (Vertex v: this.corners)
			s += v+", ";
		return s.substring(0,s.length()-2) + ")";
	}
}
