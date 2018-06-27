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

import linalg.Matrix;

/**
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	public static final int ESE = 0, ENE = 1, NNE = 2, NNW = 3;
	public static final int WNW = 4, WSW = 5, SSW = 6, SSE = 7;
	
	private final Vertex[] neighbors = new Vertex[8]; //the eight neighbors (might be null)
//	private final VertexSet sisters; //any vertices that occupy the same spot on the globe
	private final double delP, delL; //the latitudinal and longitudinal spans
	private double mass; //its inertia
	private double x, y; //the current planar coordinates
	private Matrix netForceDensity;
	
	
	public Vertex(double delP, double delL, double mass, double x, double y) {
		this.delP = delP;
		this.delL = delL;
		this.mass = mass;
		this.x = x;
		this.y = y;
//		this.sisters = new VertexSet(this);
	}
	
	
	public Vertex(Vertex sister) {
		this.delP = sister.delP;
		this.delL = sister.delL;
		sister.mass /= 2; //TODO: not always half
		this.mass = sister.mass;
		this.x = sister.x;
		this.y = sister.y;
//		this.sisters = sister.sisters;
//		this.sisters.add(this);
	}
	
	
	/**
	 * Compute the force on this vertex, the stress divergence, and save it.
	 */
	void computenetForceDensity() {
		this.netForceDensity = new Matrix(2, 1);
		netForceDensity.set(0, 0, Math.random()*2-1);
		netForceDensity.set(1, 0, Math.random()*2-1);
	}
	
	
	Matrix getnetForceDensity() {
		return this.netForceDensity;
	}
	
	
	void setnetForceDensity(Matrix netForceDensity) {
		this.netForceDensity = netForceDensity;
	}
	
	
	void descend(double timestep) {
		this.x += timestep*netForceDensity.get(0, 0);
		this.y += timestep*netForceDensity.get(1, 0);
	}
	
	
	void connectTo(int direction, Vertex neighbor) {
		assert this != neighbor;
		this.neighbors[direction] = neighbor; // set that as this neighbor
		if (direction%2 == 0) // if the direction is clockwise of a cardinal
			neighbor.neighbors[(direction+5)%8] = this; // the reverse is an advance of 5
		else // if it is widdershins of a cardinal
			neighbor.neighbors[(direction+3)%8] = this; // the reverse is an advance of 3
	}
	
	
	void disconnectFrom(int direction) {
		if (direction%2 == 0) // if the direction is clockwise of a cardinal
			this.neighbors[direction].neighbors[(direction+5)%8] = null; // the reverse is an advance of 5
		else // if it is widdershins of a cardinal
			this.neighbors[direction].neighbors[(direction+3)%8] = null; // the reverse is an advance of 3
		this.neighbors[direction] = null; // now forget it
	}
	
	
	public double getX() {
		return x;
	}
	
	
	public double getY() {
		return y;
	}
	
	
	public Vertex getNeighbor(int direction) {
		return neighbors[direction];
	}
}
