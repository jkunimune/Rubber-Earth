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

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Function;

/**
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	private final Set<Cell> neighbors; // the eight neighbors (might be null)
	private final double lat, lon; // the spherical coordinates
	private double x, y; // the current planar coordinates
	private double forceX, forceY;
	
	
	public Vertex(double lat, double lon) {
		this.lat = lat;
		this.lon = lon;
		this.neighbors = new HashSet<Cell>();
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
	
	
	void setForce(double fx, double fy) {
		this.forceX = fx;
		this.forceY = fy;
	}
	
	void descend(double timestep) {
		this.x += timestep*this.forceX;
		this.y += timestep*this.forceY;
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
	
	public Set<Cell> getNeighborsUnmodifiable() {
		return Collections.unmodifiableSet(this.neighbors);
	}
	
	public void addNeighbor(Cell neighbor) {
		this.neighbors.add(neighbor);
	}
	
	@Override
	public String toString() {
		return "Vertex("+getX()+", "+getY()+")";
	}
}
