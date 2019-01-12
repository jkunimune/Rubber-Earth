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
	private double X, Y; // the current planar coordinates
	private double velX, velY;
	private final Map<Element, double[]> forces; // the attached cells and the forces they exert
	private Vertex clockwise, widershin; // the next vertices along the edge
	
	
	public Vertex(double lat, double lon) {
		this.lat = lat;
		this.lon = lon;
		this.forces = new HashMap<Element, double[]>();
		this.clockwise = null;
		this.widershin = null; // null, null, and NaN are the defaults for non-edges
	}
	
	public Vertex(double lat, double lon, double x, double y) {
		this(lat, lon);
		this.X = x;
		this.Y = y;
	}
	
	public Vertex(double lat, double lon, Function<double[], double[]> projection) {
		this(lat, lon);
		double[] coordinates = projection.apply(new double[] {lat, lon});
		this.X = coordinates[0];
		this.Y = coordinates[1];
	}
	
	public Vertex(Vertex that) {
		this(that.lat, that.lon, that.X, that.Y);
	}
	
	
	void setForce(Element exerter, double forceX, double forceY) {
		this.forces.get(exerter)[0] = forceX;
		this.forces.get(exerter)[1] = forceY;
	}
	
	double getForceX(Element exerter) {
		return this.forces.get(exerter)[0];
	}
	
	double getForceY(Element exerter) {
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
		this.X += timestep*this.velX;
		this.Y += timestep*this.velY;
	}
	
	boolean isEdge() {
		return this.clockwise != null;
	}
	
	double[] getEdgeDirection() { // TODO Do I need this?
		double tX = this.widershin.getX()-this.clockwise.getX();
		double tY = this.widershin.getY()-this.clockwise.getY();
		double t = Math.hypot(tX, tY);
		return new double[] {tX/t, tY/t};
	}
	
	void stepX(double step) {
		this.X += step;
	}
	
	void stepY(double step) {
		this.Y += step;
	}
	
	public double getX() {
		return this.X;
	}
	
	public double getY() {
		return this.Y;
	}
	
	public double getR() {
		return Math.hypot(this.getX(), this.getY());
	}
	
	void setPos(double X, double Y) {
		this.X = X;
		this.Y = Y;
	}
	
	public double getTransformedX(double... transform) { // this seems as good a place to put this method as any
		double cx = transform[0], cy = transform[1], th = transform[2];
		return   (this.getX()-cx)*Math.cos(th) + (this.getY()-cy)*Math.sin(th);
	}
	
	public double getTransformedY(double... transform) {
		double cx = transform[0], cy = transform[1], th = transform[2];
		return - (this.getX()-cx)*Math.sin(th) + (this.getY()-cy)*Math.cos(th);
	}
	
	public double getLat() {
		return this.lat;
	}
	
	public double getLon() {
		return this.lon;
	}
	
	public double getStrength() {
		double w = 0;
		for (Element c: this.getNeighborsUnmodifiable())
			w += c.getStrength()/this.getNeighborsUnmodifiable().size();
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
	
	public Set<Element> getNeighborsUnmodifiable() {
		return Collections.unmodifiableSet(this.forces.keySet());
	}
	
	public Set<Element> getNeighborsUnmodifiable(boolean clone) {
		if (clone)
			return new HashSet<Element>(getNeighborsUnmodifiable());
		else
			return getNeighborsUnmodifiable();
	}
	
	public List<Element> getNeighborsInOrder() { // from the widdershins neighbour to the clockwise neighbour
		LinkedList<Element> out = new LinkedList<Element>();
		for (Element c: this.getNeighborsUnmodifiable()) { // start with the unordered set
			if (c.getVerticesUnmodifiable().contains(this.getWidershinNeighbor())) { // find the one adjacent to the widdershins neighbour
				out.addFirst(c);
				break;
			}
		}
		untilComplete:
		while (!out.getLast().isAdjacentTo(this.getClockwiseNeighbor())) { // then until you are adjacent to the clockwise neighbour
			for (Element c: this.getNeighborsUnmodifiable()) { // look for the next cell
				if (!out.contains(c) && c.isAdjacentTo(out.getLast())) { // that is not already in the list and adjacent to the last one
					out.addLast(c); // and add it
					continue untilComplete;
				}
			}
			assert false: "There was nothing to add.";
		}
		assert out.size() == this.getNeighborsUnmodifiable().size();
		return out;
	}
	
	void addNeighbor(Element neighbor) {
		this.forces.put(neighbor, new double[] {0.,0.});
	}
	
	void transferNeighbor(Element neighbor, Vertex repl) {
		this.forces.remove(neighbor);
		repl.forces.put(neighbor, new double[] {0., 0.});
		for (int i = 0; i < 3; i ++)
			if (neighbor.getVertex(i) == this)
				neighbor.setVertex(i, repl);
	}
	
	public Collection<Vertex> getLinks() { // all vertices with which we have an edge
		Set<Vertex> out = new HashSet<Vertex>();
		for (Element e: this.getNeighborsUnmodifiable())
			out.addAll(e.getVerticesUnmodifiable());
		return out;
	}
	
	public double undeformedDistanceTo(Vertex that) {
		for (Element e: this.getNeighborsUnmodifiable()) { // look for an Element they share
			if (e.isAdjacentTo(that)) { // (there should be two; it doesn't matter which)
				double[] Ra = e.getUndeformedPos(this);
				double[] Rb = e.getUndeformedPos(that);
				return Math.hypot(Ra[0] - Rb[0], Ra[1] - Rb[1]);
			}
		}
		throw new IllegalArgumentException("Vertices must be adjacent.");
	}
	
	public double distanceTo(Vertex that) {
		return Math.hypot(this.getX()-that.getX(), this.getY()-that.getY());
	}
	
	public boolean isLeftOf(Vertex v0, Vertex v2) {
		return (this.getX()-v0.getX()) * (v2.getY()-this.getY())
				- (this.getY()-v0.getY()) * (v2.getX()-this.getX()) < 0;
	}
	
	@Override
	public String toString() {
		return "Vertex("+getX()+", "+getY()+")";
	}
}
