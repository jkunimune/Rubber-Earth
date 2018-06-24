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

/**
 * An array of points that represents the Earth
 * 
 * @author Justin Kunimune
 */
public class Mesh {
	
	private final VertexSet[][] vertices;
	private boolean done = false;
	
	
	
	public Mesh(int resolution, InitialConfiguration init) {
		this.vertices = new VertexSet[2*resolution][4*resolution];
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
	
	
	
	/**
	 * Determines how the thing will start out.
	 * 
	 * @author Justin Kunimune
	 */
	public enum InitialConfiguration {
		SINUSOIDAL {
			public VertexSet initialVertexSet(int i, int j, int res) {
				// TODO: Implement this
				return null;
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
