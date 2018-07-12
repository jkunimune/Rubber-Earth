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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

/**
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	private final double lat, lon; // the spherical coordinates
	private double x, y; // the current planar coordinates
	private double velX, velY;
	private final Map<Cell, double[]> forces; // the attached cells and the forces they exert
	private Vertex clockwise, widershin; // the next vertices along the edge
	
	
	public Vertex(double lat, double lon) {
		this.lat = lat;
		this.lon = lon;
		this.forces = new HashMap<Cell, double[]>();
		this.clockwise = null;
		this.widershin = null; // null, null, and NaN are the defaults for non-edges
	}
	
	public Vertex(double lat, double lon, double x, double y) {
		this(lat, lon);
		this.x = x;
		this.y = y;
	}
	
	public Vertex(double lat, double lon, Function<double[], double[]> projection) {
		this(lat, lon);
		double[] coordinates = projection.apply(new double[] {lat, lon});
		this.x = coordinates[0];
		this.y = coordinates[1];
	}
	
	public Vertex(Vertex that) {
		this(that.getLat(), that.getLon());
		this.x = that.getX();
		this.y = that.getY();
	}
	
	
	void setForce(Cell exerter, double forceX, double forceY) {
		this.forces.get(exerter)[0] = forceX;
		this.forces.get(exerter)[1] = forceY;
	}
	
	double getForceX(Cell exerter) {
		return this.forces.get(exerter)[0];
	}
	
	double getForceY(Cell exerter) {
		return this.forces.get(exerter)[1];
	}
	
	void setVel(double velX, double velY) {
		this.velX = velX;
		this.velY = velY;
	}
	
	void descend(double timestep) {
		this.x += timestep*this.velX;
		this.y += timestep*this.velY;
	}
	
	boolean isEdge() {
		return this.clockwise != null;
	}
	
	double[] getEdgeDirection() {
		double tX = this.widershin.getX()-this.clockwise.getX();
		double tY = this.widershin.getY()-this.clockwise.getY();
		double t = Math.hypot(tX, tY);
		return new double[] {tX/t, tY/t};
	}
	
	void stepX(double step) {
		this.x += step;
	}
	
	void stepY(double step) {
		this.y += step;
	}
	
	public double getX() {
		return this.x;
	}
	
	public double getY() {
		return this.y;
	}
	
	public double getLat() {
		return this.lat;
	}
	
	public double getLon() {
		return this.lon;
	}
	
	void setClockwiseNeighbor(Vertex neighbor) {
		this.clockwise = neighbor;
		neighbor.widershin = this;
	}
	
	void setWidershinNeighbor(Vertex neighbor) {
		this.widershin = neighbor;
		neighbor.clockwise = this;
	}
	
	public Set<Cell> getNeighborsUnmodifiable() {
		return Collections.unmodifiableSet(this.forces.keySet());
	}
	
	public Set<Cell> getNeighborsUnmodifiable(boolean clone) {
		if (clone)
			return new HashSet<Cell>(getNeighborsUnmodifiable());
		else
			return getNeighborsUnmodifiable();
	}
	
	public void addNeighbor(Cell neighbor) {
		this.forces.put(neighbor, new double[] {0.,0.});
	}
	
	public void transferNeighbor(Cell neighbor, Vertex repl) {
		this.forces.remove(neighbor);
		repl.forces.put(neighbor, new double[] {0., 0.});
		for (int i = 0; i < 4; i ++)
			if (neighbor.getCorner(i) == this)
				neighbor.setCorner(i, repl);
	}
	
	@Override
	public String toString() {
		return "Vertex("+getX()+", "+getY()+")";
	}
}
