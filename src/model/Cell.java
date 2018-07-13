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
import java.util.Collections;
import java.util.List;

import linalg.Matrix;

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
	private double defaultEnergy; // the stored elastic potential energy
	
	
	
	public Cell(double lambda, double mu, double delP, double delL,
			Vertex ne, Vertex nw, Vertex sw, Vertex se) {
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
	 * @return the energy felt by this cell
	 */
	double computeAndSaveEnergy() {
		this.defaultEnergy = getCurrentEnergy();
		return this.defaultEnergy;
	}
	
	/**
	 * Compute the difference in energy in this cell between this state and the default state.
	 * @return the increase in energy
	 */
	double computeDeltaEnergy() {
		return getCurrentEnergy() - this.defaultEnergy;
	}
	
	/**
	 * Compute the total energy in the cell in this current configuration.
	 * @return the current energy.
	 */
	private double getCurrentEnergy() {
		Vertex  ne = corners[NORTHEAST], nw = corners[NORTHWEST],
				sw = corners[SOUTHWEST], se = corners[SOUTHEAST];
		Matrix F = new Matrix(new double[][] {
			{
				((ne.getX()+se.getX())/2 - (nw.getX()+sw.getX())/2)/delL,
				((ne.getY()+se.getY())/2 - (nw.getY()+sw.getY())/2)/delL
			}, {
				((ne.getX()+nw.getX())/2 - (se.getX()+sw.getX())/2)/delP,
				((ne.getY()+nw.getY())/2 - (se.getY()+sw.getY())/2)/delP
			}
		});
		Matrix B = F.times(F.T());
		double J = F.det();
		double i1 = B.tr();
		return mu/2*(i1 - 2 - 2*Math.log(J)) + lambda/2*Math.pow(Math.log(J), 2) * delP*delL; // is this volume term supposed to be undeformed or deformed volume? I can't find a good answer on the internet.
	}
	
	public double getVolume() {
		return this.delP * this.delL;
	}
	
	
	public double getCX() {
		double x = 0;
		for (Vertex v: corners)
			x += v.getX();
		return x/4;
	}
	
	public double getCY() {
		double y = 0;
		for (Vertex v: corners)
			y += v.getY();
		return y/4;
	}
	
	
	public Vertex getCorner(int direction) {
		return this.corners[direction];
	}
	
	public void setCorner(int direction, Vertex corner) {
		this.corners[direction] = corner;
	}
	
	public List<Vertex> getCornersUnmodifiable() {
		return Collections.unmodifiableList(Arrays.asList(this.corners));
	}
	
	
	@Override
	public String toString() {
		String s = "Cell(";
		for (Vertex v: this.corners)
			s += v+", ";
		return s.substring(0,s.length()-2) + ")";
	}
}
