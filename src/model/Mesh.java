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
	
	private static final double STEP = 1e-8; // an arbitrarily small number
	private static final double ARMIJO_GOLDSTEIN_C = 0.7; // for the backtracking
	private static final double ARMIJO_GOLDSTEIN_T = 0.5; // for the backtracking
	
	private final Collection<Cell> cells;
	private final Collection<Vertex> vertices;
	private final double precision;
	private final double lengthScale;
	private double elasticEnergy;
	private double tearLength;
	
	
	
	public Mesh(int resolution, InitialConfig init,
			double lambda, double mu, double precision) {
		this.cells = new LinkedList<Cell>();
		this.vertices = new LinkedList<Vertex>();
		this.precision = precision;
		this.lengthScale = Math.PI/2 / resolution;
		
		for (int i = 0; i < 2*resolution; i ++)
			for (int j = 0; j < 4*resolution; j ++) // let init populate the mesh
				init.spawnCell(i, j, resolution, lambda, mu, cells, vertices);
		init.cleanup(); // let init finish up
		for (Cell c: cells)
			for (Vertex v: c.getCornersUnmodifiable()) // make sure these relationships are mutual
				v.addNeighbor(c);
		
		this.elasticEnergy = computeTotEnergy();
		this.tearLength = 0;
	}
	
	
	
	/**
	 * Move all vertices to a slightly more favourable position
	 * @param dampFactor - The amount to damp the poles; 0 for not at all and 1 for all the way
	 * @return true if it successfully found a more favourable configuration,
	 * 		false if it thinks it's time to quit.
	 */
	public boolean update(double dampFactor) {
		double Ui = this.elasticEnergy;
		
		double maxVel = 0;
		double gradDotVel = 0;
		for (Vertex v: vertices) {
			double netForceX = 0, netForceY = 0;
			for (Cell c: v.getNeighborsUnmodifiable()) {
				v.stepX(STEP);
				double forceX = -c.computeDeltaEnergy()/STEP; // compute the force by computing the energy gradient
				v.stepX(-STEP);
				v.stepY(STEP);
				double forceY = -c.computeDeltaEnergy()/STEP;
				v.stepY(-STEP);
				v.setForce(c, forceX, forceY); // store the force from each cell individually for later
				netForceX += forceX;
				netForceY += forceY;
			}
			
			double damping = (1-dampFactor) + dampFactor*Math.cos(v.getLat());
			double velX = netForceX*damping; // damp forces nearer the poles
			double velY = netForceY*damping; // because they have smaller length scales
			v.setVel(velX, velY);
			
			assert !Double.isNaN(velX) && !Double.isNaN(velY);
			double vel = Math.hypot(velX, velY);
			if (vel > maxVel)
				maxVel = vel;
			gradDotVel += - netForceX*velX - netForceY*velY;
		}
		
		double timestep = .5*lengthScale/maxVel;
		for (Vertex v: vertices) // the first timestep is whatever makes the fastest one move one half cell-length
			v.descend(timestep);
		
		double Uf = computeTotEnergy();
		while ((Double.isNaN(Uf) || Uf - Ui > ARMIJO_GOLDSTEIN_C*gradDotVel*timestep) && // if the energy didn't decrease enough
				maxVel*timestep >= precision) {
			for (Vertex v: vertices)
				v.descend(-timestep*(1-ARMIJO_GOLDSTEIN_T)); // backstep and try again
			timestep *= ARMIJO_GOLDSTEIN_T;
			Uf = computeTotEnergy();
		}
		
		if (maxVel*timestep < precision) { // if our steps are really small, then we're done
			for (Vertex v: vertices) // just reset to before we started backtracking
				v.descend(-timestep);
			this.elasticEnergy = computeTotEnergy();
			return false;
		}
		
		this.elasticEnergy = Uf;
		return true;
	}
	
	
	/**
	 * Find the vertex with the highest strain, and separate it into two vertices.
	 * @return true if it successfully tore, false if it could find nothing to tear
	 */
	public boolean rupture() {
		double maxStrain = 0;
		Vertex v0 = null;
		for (Vertex v: this.getVerticesUnmodifiable()) { // first we have to choose where to rupture
			if (v.isEdge()) {
				double[] edge = v.getEdgeDirection();
				double strain = 0;
				double volume = 0;
				for (Cell c: v.getNeighborsUnmodifiable()) {
					double forceDotEdge = v.getForceX(c)*edge[0] + v.getForceY(c)*edge[1];
					double crDotEdge = (c.getCX()-v.getX())*edge[0] + (c.getCY()-v.getY())*edge[1];
					strain += Math.signum(crDotEdge)*forceDotEdge;
					
					volume += c.getVolume(); // get the total involved volume (area)
				}
				strain /= Math.sqrt(volume); // use this as an approximation for surface area (length)
				if (strain > maxStrain) {
					maxStrain = strain;
					v0 = v;
				}
			}
		}
		if (v0 == null)
			return false;
		
		Vertex v1 = new Vertex(v0);
		this.vertices.add(v1);
		double[] edge = v0.getEdgeDirection();
		for (Cell c: v0.getNeighborsUnmodifiable(true)) // look at the cells
			if (v0.getForceX(c)*edge[0] + v0.getForceY(c)*edge[1] < 0) // if it is pulling clockwise
				v0.transferNeighbor(c, v1); // detatch it
		
		for (Vertex vM: v0.getLinks()) { // now look at the Vertices that are still connected to v0
			if (vM.getLinks().contains(v1)) { // find one that is also linked to v1 now
				if (v0.directionTo(vM) == Vertex.WIDERSHIN) { // update the edge chain accordingly
					v1.setWidershinNeighbor(v0.getWidershinNeighbor());
					vM.setWidershinNeighbor(v1);
					v0.setWidershinNeighbor(vM);
				}
				else {
					assert (v0.directionTo(vM) == Vertex.CLOCKWISE);
					v1.setClockwiseNeighbor(v0.getClockwiseNeighbor());
					vM.setClockwiseNeighbor(v1);
					v0.setClockwiseNeighbor(vM);
				}
				this.tearLength += v0.distanceTo(vM); // this is the length of the tear
				break;
			}
		}
		
		return true;
	}
	
	
	/**
	 * Compute the total elastic energy in the system and save it as the "default" state.
	 * @return the total elastic energy.
	 */
	private double computeTotEnergy() {
		double U = 0;
		for (Cell c: cells)
			U += c.computeAndSaveEnergy();
		return U;
	}
	
	
	/** Convert spherical coordinates to cartesian coordinates using the current mest configuration.
	 * @param lat - The latitude of the point to map
	 * @param lon - The longitude of the point to map
	 * @return an array of two elements: {x, y}
	 */
	public double[] map(double lat, double lon) {
		return new double[] {lon, lat};
	}
	
	
	public double getTotalTearLength() {
		return this.tearLength;
	}
	
	
	/**
	 * Return the total elastic energy in the system from the last step.
	 * This differs from computeTotEnergy in that it uses a saved field, not an on-the-spot computation.
	 * This means that it is faster and less prone to oscillating.
	 * @return 
	 */
	public double getTotEnergy() {
		return this.elasticEnergy;
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
	public enum InitialConfig {
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
					vertexArray[i][j] = new Vertex(phiN, lamW, InitialConfig::sinusoidalProj);
				if (j == 0) // create the left prime meridian, if necessary
					vertexArray[i+1][j] = new Vertex(phiS, lamW, InitialConfig::sinusoidalProj);
				if (i == 0) // make sure the North Pole is one vertex
					vertexArray[i][j+1] = vertexArray[i][j];
				if (i == 2*res-1) // make sure the South Pole is one vertex
					vertexArray[i+1][j+1] = vertexArray[i+1][j];
				else // generate new interior Vertices if necessary
					vertexArray[i+1][j+1] = new Vertex(phiS, lamE, InitialConfig::sinusoidalProj);
				Vertex  nw = vertexArray[i][j], ne = vertexArray[i][j+1], // reuse all other vertices
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
			
			public void cleanup() {
				for (int i = 0; i < vertexArray.length-1; i ++) {
					vertexArray[i][0].setWidershinNeighbor(vertexArray[i+1][0]);
					vertexArray[i][vertexArray[i].length-1].setClockwiseNeighbor(
							vertexArray[i+1][vertexArray[i].length-1]);
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
		
		/**
		 * Do anything that needs to be done once all the cells are spawned.
		 */
		public void cleanup() {}
		
		
		public static double[] sinusoidalProj(double[] sphereCoords) {
			return new double[] { sphereCoords[1]*Math.cos(sphereCoords[0]), sphereCoords[0] };
		}
	}
}
