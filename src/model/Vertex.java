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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

/**
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	public static final int NOT_CONNECTED = 0;
	public static final int CLOCKWISE = 1;
	public static final int WIDERSHIN = 2;
	
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
	
	double getVelX() {
		return this.velX;
	}
	
	double getVelY() {
		return this.velY;
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
	
	public double getWeight() {
		double w = 0;
		for (Cell c: this.getNeighborsUnmodifiable())
			w += c.getWeight()/this.getNeighborsUnmodifiable().size();
		return w;
	}
	
	void setClockwiseNeighbor(Vertex neighbor) {
		this.clockwise = neighbor;
		neighbor.widershin = this;
	}
	
	void setWidershinNeighbor(Vertex neighbor) {
		this.widershin = neighbor;
		neighbor.clockwise = this;
	}
	
	public Vertex getClockwiseNeighbor() {
		return this.clockwise;
	}
	
	public Vertex getWidershinNeighbor() {
		return this.widershin;
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
	
	public List<Cell> getNeighborsInOrder() { // from the widdershins neighbour to the clockwise neighbour
		LinkedList<Cell> out = new LinkedList<Cell>();
		for (Cell c: this.getNeighborsUnmodifiable()) { // start with the unordered set
			if (c.getCornersUnmodifiable().contains(this.getWidershinNeighbor())) { // find the one adjacent to the widdershins neighbour
				out.addFirst(c);
				break;
			}
		}
		while (!out.getLast().isAdjacentTo(this.getClockwiseNeighbor())) { // then until you are adjacent to the clockwise neighbour
			for (Cell c: this.getNeighborsUnmodifiable()) { // look for the next cell
				if (!out.contains(c) && c.isAdjacentTo(out.getLast())) { // that is not already in the list and adjacent to the last one
					out.addLast(c); // and add it
					break;
				}
			}
		}
		return out;
	}
	
	void addNeighbor(Cell neighbor) {
		this.forces.put(neighbor, new double[] {0.,0.});
	}
	
	void transferNeighbor(Cell neighbor, Vertex repl) {
		this.forces.remove(neighbor);
		repl.forces.put(neighbor, new double[] {0., 0.});
		for (int i = 0; i < 4; i ++)
			if (neighbor.getCorner(i) == this)
				neighbor.setCorner(i, repl);
	}
	
	public Collection<Vertex> getLinks() { // only vertices with which we have an edge
		Set<Vertex> out = new HashSet<Vertex>();
		for (Cell c: this.getNeighborsUnmodifiable()) {
			for (int i = 0; i < 4; i ++) {
				if (this == c.getCorner(i)) {
					out.add(c.getCorner((i+1)%4));
					out.add(c.getCorner((i+3)%4));
				}
			}
		}
		out.remove(this);
		return out;
	}
	
	public int directionTo(Vertex that) {
		for (Cell c: this.getNeighborsUnmodifiable()) {
			for (int i = 0; i < 4; i ++) {
				if (c.getCorner(i) == this && c.getCorner((i+1)%4) == that)
					return WIDERSHIN;
				else if (c.getCorner(i) == that && c.getCorner((i+1)%4) == this)
					return CLOCKWISE;
			}
		}
		return NOT_CONNECTED;
	}
	
	public double distanceTo(Vertex that) {
		return Math.hypot(this.getX()-that.getX(), this.getY()-that.getY());
	}
	
	@Override
	public String toString() {
		return "Vertex("+getX()+", "+getY()+")";
	}
}
