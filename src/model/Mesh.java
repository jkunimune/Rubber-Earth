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

import java.util.Iterator;

/**
 * An array of points that represents the Earth
 * 
 * @author Justin Kunimune
 */
public class Mesh implements Iterable<Vertex> {
	
	private final VertexSet[][] vertices;
	private boolean done = false;
	
	
	
	public Mesh(int resolution, InitialConfiguration init) {
		this.vertices = new VertexSet[2*resolution+1][4*resolution];
		for (int i = 0; i < vertices.length; i ++)
			for (int j = 0; j < vertices[i].length; j ++)
				vertices[i][j] = init.initialVertexSet(i, j, resolution);
	}
	
	
	
	/**
	 * Move all vertices to a slightly more favourable position
	 */
	public void update() {
		// TODO: Implement this
		
	}
	
	
	/**
	 * Should we stop?
	 * @return whether it is done
	 */
	public boolean isDone() {
		return this.done;
	}
	
	
	
	public Iterator<Vertex> iterator() {
		return new Iterator<Vertex>() {
			private int i = 0; //the row of the current VertexSet
			private int j = 0; //the col of the current VertexSet
			private Iterator<Vertex> currentSet = vertices[0][0].iterator();
			
			public boolean hasNext() {
				return	currentSet.hasNext() ||
						j+1 < vertices[i].length ||
						i+1 < vertices.length;
			}
			
			public Vertex next() {
				if (currentSet == null || !currentSet.hasNext()) {
					j ++;
					if (j >= vertices[i].length) {
						i ++;
						if (i >= vertices.length) {
							return null;
						}
						j = 0;
					}
					currentSet = vertices[i][j].iterator();
				}
				return currentSet.next();
			}
		};
	}
	
	
	
	/**
	 * Determines how the thing will start out.
	 * 
	 * @author Justin Kunimune
	 */
	public enum InitialConfiguration {
		SINUSOIDAL {
			public VertexSet initialVertexSet(int i, int j, int res) {
				double phi = Math.PI/2 * (res - i)/res;
				double lam = Math.PI/2 * (j - 2*res)/res;
				double delP = Math.PI/2 / res;
				double delL = delP * Math.cos(phi);
				double x = lam * Math.cos(phi);
				double y = phi;
				if (j != 0) // for the majority of things
					return new VertexSet(new Vertex(delP, delL, x, y));
				else { //for the prime meridian
					Vertex eastern = new Vertex(delP, delL, x, y);
					Vertex western = new Vertex(delP, delL, -x, y);
					return new VertexSet(eastern, western, western, eastern);
				}
			}
		},
		
		SINUSOIDAL_FLORENCE {
			public VertexSet initialVertexSet(int i, int j, int res) {
				// TODO: Implement this
				return null;
			}
		},
		
		AZIMUTHAL {
			public VertexSet initialVertexSet(int i, int j, int res) {
				// TODO: Implement this
				return null;
			}
		};
		
		
		public abstract VertexSet initialVertexSet(int i, int j, int res);
	}
}
