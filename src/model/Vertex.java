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

import java.util.Arrays;

import linalg.Matrix;

/**
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	public static final int ENE = 0, NNE = 1, NNW = 2, WNW = 3;
	public static final int WSW = 4, SSW = 5, SSE = 6, ESE = 7;
	
	private final Vertex[] neighbors = new Vertex[8]; // the eight neighbors (might be null)
//	private final VertexSet sisters; // any vertices that occupy the same spot on the globe
	private final double lambda, mu; // the elastic properties
	private final double delP, delL; // the latitudinal and longitudinal spans
	private double mass; // its inertia
	private double x, y; // the current planar coordinates
	private Matrix netForce;
	
	
	public Vertex(double lambda, double mu, double delP, double delL, double mass, double x, double y) {
		this.lambda = lambda;
		this.mu = mu;
		this.delP = delP;
		this.delL = delL;
		this.mass = mass;
		this.x = x;
		this.y = y;
//		this.sisters = new VertexSet(this);
	}
	
	
	public Vertex(Vertex sister) {
		this.lambda = sister.lambda;
		this.mu = sister.mu;
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
	void computeNetForce() {
		int regime = 0;
		Matrix[] stresses = new Matrix[4];
		for (int i = ENE; i < ESE; i += 2) {
			if (this.neighbors[i] != null) {
				if (this.neighbors[i+1] == null)
					System.out.println(Arrays.toString(neighbors));
				stresses[i/2] = cauchyStress(this, neighbors[i], neighbors[i+1], i/2);
				regime += 1;
			}
		}
		if (delL == 0)
			regime = 0; // use this to represent poles
		
		Matrix forceM; // the net force in material coordinates
		switch (regime) {
		case 0: // pole
			forceM = new Matrix(2, 1);
			for (int i = 0; i < 4; i ++) { // look only at the latitudinal connections
				int k = i/2*4 + i%2 + NNE;
				if (neighbors[k] != null) {
					Matrix dl = new Matrix(2, 1,
							neighbors[k].getX() - this.getX(),
							neighbors[k].getY() - this.getY());
					double stretch = dl.mag()/delP;
					double pressure = mu*(1 - 1/(stretch*stretch)) + lambda*(stretch - 1); // approximate poles as uniaxial extension
					int sign = (k < 4) ? 1 : -1;
					forceM.add(1, 0, sign*pressure*neighbors[k].delL/4); // /2 because averaging between the polar and nonpolar and /2 again because only on one side
				}
			}
			break;
		
		case 1: // corner
			forceM = new Matrix(2, 1);
			break;
			
		case 2: // edge
			Matrix sig = new Matrix(2, 2);
			for (Matrix sigi: stresses)
				if (sigi != null)
					sig = sig.plus(sigi);
			sig = sig.over(2);
			
			Matrix nHat;
			double surfArea;
			if (neighbors[VertexSet.NORTHEAST*2] == null) { // if the edge is to the N or E
				if (neighbors[VertexSet.NORTHWEST*2] == null) // if the edge is to the N or W
					nHat = new Matrix(2, 1, 0., -1.); // it's to the N
				else
					nHat = new Matrix(2, 1, -1., 0.); // it's to the E
			}
			else {
				if (neighbors[VertexSet.NORTHWEST*2] == null) // if the edge is to the N or W
					nHat = new Matrix(2, 1, 1., 0.); // it's to the W
				else
					nHat = new Matrix(2, 1, 0., 1.); // it's to the S
			}
			if (nHat.get(0, 0) == 0) // if the edge is to the N or S
				surfArea = delL;
			else
				surfArea = delP;
			if (surfArea == 0) 	surfArea = delP; //XXX: this is a workaround; I need to work with the poles later
			
			forceM = sig.times(nHat).times(surfArea);
			break;
			
		case 3: // inner corner
			forceM = Matrix.zeroes(2, 1);
			break; //TODO
			
		case 4: // completely surrounded
			Matrix divSig = new Matrix(2, 1);
			for (int i = 0; i < 2; i ++) { // i is the element of the final vector
				for (int j = 0; j < 2; j ++) { // j is the side on which we measure the derivatives
					divSig.add(i, 0,
							(stresses[(j+3)%4].get(0, i) - stresses[(j+1)%4].get(0, i))/delL +
							(stresses[(j+0)%4].get(1, i) - stresses[(j+2)%4].get(1, i))/delP);
				}
			}
			forceM = divSig.times(delL*delP);
			break;
			
		default:
			assert false;
			return;
		}
		assert !forceM.isNaN() : "Nonnumeric force calculated by "+this+"!";
		
		Matrix bases = new Matrix(2, 2);
		for (int i = 0; i < 2; i ++) {
			for (int j = 0; j < 4; j ++) {
				int k = (j/2*4 + j%2 + i*2 + 3)%8;
				if (neighbors[k] != null) {
					int s = (j < 2) ? -1 : 1; // protip: switch the sign on this for some wicked art
					bases.add(0, i, s*(neighbors[k].getX()-this.getX()));
					bases.add(1, i, s*(neighbors[k].getY()-this.getY()));
				}
			}
		}
		this.netForce = bases.norm().times(forceM);
	}
	
	
	Matrix getNetForce() {
		return this.netForce;
	}
	
	
	void setNetForce(Matrix netForceDensity) {
		this.netForce = netForceDensity;
	}
	
	
	double getSpeed2() {
		return this.netForce.sqr()/(this.mass*this.mass);
	}
	
	
	void descend(double timestep) {
		this.x += timestep*netForce.get(0, 0)/mass;
		this.y += timestep*netForce.get(1, 0)/mass;
	}
	
	
	void connectTo(int direction, Vertex neighbor) {
		assert this != neighbor;
		if (neighbor == null) 	return;
		this.neighbors[direction] = neighbor; // set that as this neighbor
		if (direction%2 == 0) // if the direction is widdershins of a cardinal
			neighbor.neighbors[(direction+3)%8] = this; // the reverse is an advance of 3
		else // if it is clockwise of a cardinal
			neighbor.neighbors[(direction+5)%8] = this; // the reverse is an advance of 5
	}
	
	
	void disconnectFrom(int direction) {
		if (direction%2 == 0) // if the direction is widdershins of a cardinal
			this.neighbors[direction].neighbors[(direction+3)%8] = null; // the reverse is an advance of 3
		else // if it is clockwise of a cardinal
			this.neighbors[direction].neighbors[(direction+5)%8] = null; // the reverse is an advance of 5
		this.neighbors[direction] = null; // now forget it
	}
	
	
	double distanceTo(Vertex that) {
		return Math.hypot(this.getX()-that.getX(), this.getY()-that.getY());
	}
	
	
	/**
	 * Compute the Cauchy stress tensor around these points
	 * using a Neo-Hookean model
	 * @return the true stress, in the material coordinates
	 */
	private Matrix cauchyStress(Vertex oo, Vertex io, Vertex oi, int direction) {
		Matrix I = Matrix.identity(2);
		double delx = (direction%2 == 0) ? (oo.delL + io.delL)/2 : (oo.delP + io.delP)/2; // the distance from oo to io
		double dely = (direction%2 == 0) ? (oo.delP + oi.delP)/2 : (oo.delL + oi.delL)/2; // the distance from oo to oi
		Matrix F = new Matrix(new double[][] {
			{(io.getX()-oo.getX())/delx, (io.getY()-oo.getY())/delx},
			{(oi.getX()-oo.getX())/dely, (oi.getY()-oo.getY())/dely}});
		Matrix B = F.times(F.T());
		double J = F.det();
		double i1 = B.tr();
		Matrix s = I.times(lambda*(J-1)).plus(B.times(mu/J/J)).minus(I.times(mu/2/J/J*i1));
		if (direction%2 == 1)
			s = new Matrix(new double[][] { // flip x and y if we've been in a weird direction
				{ s.get(1, 1),-s.get(0, 1)},
				{-s.get(0, 1), s.get(0, 0)}});
		return s;
	}
	
	
	public double getX() {
		return this.x;
	}
	
	
	public double getY() {
		return this.y;
	}
	
	
	public double getMass() {
		return this.mass;
	}
	
	
	public Vertex getNeighbor(int direction) {
		return this.neighbors[direction];
	}
	
	
	public String toString() {
		return "Vertex("+getX()+", "+getY()+")";
	}
}
