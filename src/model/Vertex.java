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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * A single point in the rubber mesh.
 * 
 * @author Justin Kunimune
 */
public class Vertex {
	
	private final Set<Cell> neighbors; // the eight neighbors (might be null)
//	private final VertexSet sisters; // any vertices that occupy the same spot on the globe
	private double x, y; // the current planar coordinates
	private double forceX, forceY;
	
	
	public Vertex(double x, double y) {
		this.x = x;
		this.y = y;
		this.neighbors = new HashSet<Cell>();
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
