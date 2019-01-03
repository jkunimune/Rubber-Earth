/**
 * This is free and unencumbered software released into the public domain.
 * 
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 * 
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 * 
 * For more information, please refer to <http://unlicense.org>
 */
package model;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import utils.Matrix;

/**
 * A simple interface to make these generics easier with which to deal.
 * Represents all of the mesh vertices from the same point on the globe.
 * 
 * @author Justin Kunimune
 */
public class Cell {
	
	public static final int NE = 0, NW = 1, SW = 2, SE = 3;
	private static final int CELLS_IN_CELL = 2; // the higher this number is, the better the approximation
	
	private final Vertex[] corners = new Vertex[4]; //which vertex is attached to each sector (there must be exactly one)
	private final double strength; // the max stress before tearing
	private final double lambda, mu; // the elastic properties
	private final double delP; // the latitudinal span
	private final double delLN, delLS; // the longitudinal spans at various latitudes (lengths, not angles)
	private double defaultEnergy; // the stored elastic potential energy
	
	
	
	public Cell(double strength, double lambda, double mu, double scale,
			Vertex nw, Vertex ne, Vertex sw, Vertex se) {
		this.strength = strength;
		this.lambda = lambda;
		this.mu = mu;
		this.delP = scale;
		this.delLN = scale*Math.cos(ne.getLat()); // we need multiple delL values for good approximations
		this.delLS = scale*Math.cos(se.getLat());
		corners[NE] = ne;
		corners[NW] = nw;
		corners[SW] = sw;
		corners[SE] = se;
		
		for (Vertex v: corners) // make sure the reference goes both ways
			v.addNeighbor(this);
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
		Vertex  ne = corners[NE], nw = corners[NW],
				sw = corners[SW], se = corners[SE];
		double energy = 0;
		for (int x = 0; x < CELLS_IN_CELL; x ++) { // iterate over the minicells in the cell
			for (int y = 0; y < CELLS_IN_CELL; y ++) {
				double cw = (x+.5)/CELLS_IN_CELL, ce = 1 - cw; // use linear interpolation to measure energy in part of the cell
				double cs = (y+.5)/CELLS_IN_CELL, cn = 1 - cs; // (these are the interpolation coefficients)
				double delXP = cw*(nw.getX()-sw.getX()) + ce*(ne.getX()-se.getX()); // compute the deformed side lengths
				double delYP = cw*(nw.getY()-sw.getY()) + ce*(ne.getY()-se.getY());
				double delXL = cs*(se.getX()-sw.getX()) + cn*(ne.getX()-nw.getX());
				double delYL = cs*(se.getY()-sw.getY()) + cn*(ne.getY()-nw.getY());
				double delL = cs*delLS + cn*delLN; // and the undeformed side lengths
				Matrix F = new Matrix(2, 2,
						delXL/delL, delYL/delL, // then get the deformation gradient from that
						delXP/delP, delYP/delP);
				Matrix B = F.times(F.T()); // the rest is fancy Neo-Hookean stuff
				double J = F.det();
				double i1 = B.tr();
				energy += (mu/2*(i1 - 2 - 2*Math.log(J)) + lambda/2*Math.pow(Math.log(J), 2)) *
						delP*delL/(CELLS_IN_CELL*CELLS_IN_CELL); // don't forget to multiply energy density by unformed volume
			}
		}
		
		return energy;
	}
	
	public double getVolume() {
		return this.delP * (this.delLN+this.delLS)/2.;
	}
	
	public double getStrength() {
		return this.strength;
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
	
	
	public boolean isAdjacentTo(Cell that) { // kitty-corner cells don't count
		Set<Vertex> sharedVertices = new HashSet<Vertex>();
		for (Vertex v: this.getCornersUnmodifiable())
			if (that.getCornersUnmodifiable().contains(v))
				sharedVertices.add(v); // it's kind of annoying that I need a Hash and can't use an int,
		return sharedVertices.size() >= 2; // but it doesn't count as two shared corners if they're the same Vertex
	}
	
	public boolean isAdjacentTo(Vertex v) {
		return this.getCornersUnmodifiable().contains(v);
	}
	
	
	@Override
	public String toString() {
		String s = "Cell(";
		for (Vertex v: this.corners)
			s += v+", ";
		return s.substring(0,s.length()-2) + ")";
	}
}
