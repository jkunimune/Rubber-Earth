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
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	public static final int ENE = 0, NNE = 1, NNW = 2, WNW = 3;
	public static final int WSW = 4, SSW = 5, SSE = 6, ESE = 7;
	
	private final Cell[] neighbors = new Cell[4]; // the eight neighbors (might be null)
//	private final VertexSet sisters; // any vertices that occupy the same spot on the globe
	private double x, y; // the current planar coordinates
	private double forceX, forceY;
	
	
	public Vertex(double x, double y) {
		this.x = x;
		this.y = y;
	}
	
	
	void setForce(double fx, double fy) {
		this.forceX = fx;
		this.forceY = fy;
	}
	
	
	void descend(double timestep) {
		this.x += timestep*this.forceX;
		this.y += timestep*this.forceY;
	}
	
	
	public Cell getNeighbor(int direction) {
		return neighbors[direction];
	}
	
	
	public double getX() {
		return this.x;
	}
	
	
	public double getY() {
		return this.y;
	}
	
	
	@Override
	public String toString() {
		return "Vertex("+getX()+", "+getY()+")";
	}
}
