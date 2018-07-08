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

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;

/**
 * An array of points that represents the Earth
 * 
 * @author Justin Kunimune
 */
public class Mesh {
	
	private final Collection<Cell> cells;
	private final Collection<Vertex> vertices;
	private final double precision;
	private boolean done = false;
	
	
	
	public Mesh(int resolution, InitialConfiguration init,
			double lambda, double mu, double precision) {
		this.cells = new LinkedList<Cell>();
		this.vertices = new LinkedList<Vertex>();
		this.precision = precision;
		for (int i = 0; i < 2*resolution; i ++) {
			for (int j = 0; j < 4*resolution; j ++) {
				init.spawnCell(i, j, resolution, lambda, mu, cells, vertices);
			}
		}
	}
	
	
	
	/**
	 * Move all vertices to a slightly more favourable position
	 */
	public void update() {
		double totEnergy = 0;
		for (Cell c: cells)
			totEnergy += c.computeEnergy(); // compute the elastic potential energy
		
		for (Vertex v: vertices)
			v.setForce(Math.random()-.5, Math.random()-.5);
		
		for (Vertex v: vertices)
			v.descend(0.0003);
	}
	
	
	/**
	 * Should we stop?
	 * @return whether it is done
	 */
	public boolean isDone() {
		return this.done;
	}
	
	
	public Collection<Cell> getCellsUnmodifiable() {
		return Collections.unmodifiableCollection(this.cells);
	}
	
	
	public Collection<Vertex> getVerticesUnmodifiable() {
		return Collections.unmodifiableCollection(this.vertices);
	}
	
	
	
	/**
	 * Does the tricky setup stuff. Generates the mesh in its initial position.
	 * 
	 * @author Justin Kunimune
	 */
	public enum InitialConfiguration {
		SINUSOIDAL {
			private Vertex[][] vertexArray = null;
			
			public void spawnCell(
					int i, int j, int res, double lambda, double mu,
					Collection<Cell> cells, Collection<Vertex> vertices) {
				if (vertexArray == null) 	vertexArray = new Vertex[2*res+1][4*res+1];
				
				double phiN = Math.PI/2 * (res - i)/res; // compute some coordinates
				double phiS = Math.PI/2 * (res - i-1)/res;
				double lamW = Math.PI/2 * (j - 2*res)/res;
				double lamE = Math.PI/2 * (j+1 - 2*res)/res;
				
				if (i == 0 && j == 0) // create the upper left hand corner
					vertexArray[i][j] = new Vertex(lamW*Math.cos(phiN), phiN);
				if (j == 0) // create the left prime meridian, if necessary
					vertexArray[i+1][j] = new Vertex(lamW*Math.cos(phiS), phiS);
				if (i == 0) // make sure the North Pole is one vertex
					vertexArray[i][j+1] = vertexArray[i][j];
				if (i == 2*res-1) // make sure the South Pole is one vertex
					vertexArray[i+1][j+1] = vertexArray[i+1][j];
				else // generate new interior Vertices if necessary
					vertexArray[i+1][j+1] = new Vertex(lamE*Math.cos(phiS), phiS);
				Vertex nw = vertexArray[i][j], ne = vertexArray[i][j+1], // reuse all other vertices
						sw = vertexArray[i+1][j], se = vertexArray[i+1][j+1];
				
				double delP = Math.PI/2 / res;
				double delL = delP * Math.cos((phiN+phiS)/2);
				Cell cell = new Cell(lambda, mu, delP, delL, ne, nw, sw, se);
				cells.add(cell);
				
				for (int k = 0; k < 4; k ++) { // look at those vertices
					if (!vertices.contains(cell.getCorner(k))) // if we just created this one
						vertices.add(cell.getCorner(k)); // add it to the Collection
				}
			}
		},
		
		SINUSOIDAL_FLORENCE {
			public void spawnCell(
					int i, int j, int res, double lambda, double mu, Collection<Cell> cells,
					Collection<Vertex> vertices
			) {
				// TODO: Implement this
			}
		},
		
		AZIMUTHAL {

			public void spawnCell(
					int i, int j, int res, double lambda, double mu, Collection<Cell> cells,
					Collection<Vertex> vertices
			) {
				// TODO: Implement this
			}
		};
		
		
		/**
		 * Create a new cell, and vertices if necessary, and add all created Objects to cells and vertices.
		 * @param i - The vertical index of the desired cell, from 0 (north pole) to 2*res (south pole)
		 * @param j - The horizontal index of the desired cell, from 0 (west) to 4*res (east)
		 * @param res - The number of cells in this mesh from equator to pole
		 * @param lambda - A material constant to pass on to new cells.
		 * @param mu - A material constant to pass on to new cells.
		 * @param cells - The Collection to which to add the new Cell
		 * @param vertices - The Collection to which to add any new Vertices
		 */
		public abstract void spawnCell(int i, int j, int res, double lambda, double mu,
				Collection<Cell> cells, Collection<Vertex> vertices);
	}
}
