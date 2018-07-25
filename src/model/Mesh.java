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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import linalg.Matrix;

/**
 * An array of points that represents the Earth.
 * 
 * Includes methods that seek out lower-energy configurations using L-BFGS-T, described at
 * 
 * Nocedal, Jorge. “Updating Quasi-Newton Matrices with Limited Storage.” Mathematics of
 * 		Computation, vol. 35, no. 151, 1980, pp. 773–782. JSTOR, JSTOR,
 * 		www.jstor.org/stable/2006193.
 * 
 * @author Justin Kunimune
 */
public class Mesh {
	
	private static final double STEP = 1e-8; // an arbitrarily small number
	private static final double ARMIJO_GOLDSTEIN_C = 0.7;
	private static final double BACKSTEP_TAU = 0.5;
	private static final double L_BFGS_M = 6; // the memory size
	
	private final Cell[][] cells;
	private final List<Vertex> vertices;
	private final Vertex edgeStart; // the start Vertex for iterating around the edge
	private final double precision; // determines how far we update before declaring that we have settled
	private double elasticEnergy; // the potential energy currently stored
	private double tearLength; // the length of the edge in radians
	
	private LinkedList<Matrix> sHist;
	private LinkedList<Matrix> yHist;
	private Matrix gkMinus1 = null; // the previous gradient
	
	
	
	public Mesh(int resolution, InitialConfig init,
			double lambda, double mu, double precision) {
		this.cells = new Cell[2*resolution][4*resolution];
		this.vertices = new ArrayList<Vertex>();
		this.precision = precision;
		for (int i = 0; i < 2*resolution; i ++)
			for (int j = 0; j < 4*resolution; j ++) // let init populate the mesh
				init.spawnCell(i, j, resolution, lambda, mu, cells, vertices);
		this.tearLength = init.cleanup(); // let init finish up
		for (Cell c: getCellsUnmodifiable())
			for (Vertex v: c.getCornersUnmodifiable()) // make sure these relationships are mutual
				v.addNeighbor(c);
		this.edgeStart = getAnyEdge();
		
		this.elasticEnergy = computeTotEnergy();
		this.sHist = new LinkedList<Matrix>();
		this.yHist = new LinkedList<Matrix>();
	}
	
	
	
	/**
	 * Move all vertices to a slightly more favourable position, using L-BFGS optimisation.
	 * @param dampFactor - The amount to damp the poles; 0 for not at all and 1 for all the way
	 * @return true if it successfully found a more favourable configuration,
	 * 		false if it thinks it's time to quit.
	 * @throws InterruptedException 
	 */
	public boolean update() throws InterruptedException {
		double Ui = this.elasticEnergy;
		
		Matrix gk = new Matrix(2*vertices.size(), 1); // STEP 1: compute the gradient
		for (int i = 0; i < vertices.size(); i ++) {
			Vertex v = vertices.get(i);
			for (Cell c: v.getNeighborsUnmodifiable()) {
				v.stepX(STEP);
				double gradX = c.computeDeltaEnergy()/STEP; // compute the force by computing the energy gradient
				v.stepX(-STEP);
				v.stepY(STEP);
				double gradY = c.computeDeltaEnergy()/STEP;
				v.stepY(-STEP);
				v.setForce(c, -gradX, -gradY); // store the force from each cell individually for strain calculations later
				gk.add(2*i+0, 0, gradX);
				gk.add(2*i+1, 0, gradY); // and sum the individual gradients to get the total gradient
			}
		}
		
		if (gkMinus1 != null) // STEP 5 (cont.): save historical vector information
			this.yHist.addLast(gk.minus(gkMinus1));
		if (sHist.size() > L_BFGS_M) {
			sHist.removeFirst();
			yHist.removeFirst();
		}
		
		double[] alpha = new double[sHist.size()];
		Matrix dk = gk.times(-1); // STEP 2: choose the step direction
		for (int i = sHist.size()-1; i >= 0; i --) { // this is where it gets complicated
			alpha[i] = sHist.get(i).dot(dk)/yHist.get(i).dot(sHist.get(i)); // see the paper cited at the top, page 779.
			dk = dk.plus(yHist.get(i).times(-alpha[i]));
		}
		if (!sHist.isEmpty())
			dk = dk.times(yHist.getLast().dot(sHist.getLast())/yHist.getLast().dot(yHist.getLast())); // Cholesky factor
		for (int i = 0; i < sHist.size(); i ++) {
			double beta = yHist.get(i).dot(dk)/yHist.get(i).dot(sHist.get(i));
			dk = dk.plus(sHist.get(i).times(alpha[i]-beta));
		}
		for (int i = 0; i < vertices.size(); i ++) // save the chosen step direction in the vertices
			vertices.get(i).setVel(dk.get(2*i+0, 0), dk.get(2*i+1, 0));
		
		double gradDotVel = gk.dot(dk);
		assert gradDotVel < 0;
		double timestep = 1.; // STEP 3: choose the step size
		for (Vertex v: vertices)
			v.descend(timestep);
		double Uf = computeTotEnergy();
		while ((Double.isNaN(Uf) || Uf - Ui > ARMIJO_GOLDSTEIN_C*timestep*gradDotVel)) { // if the energy didn't decrease enough
			for (Vertex v: vertices)// TODO try weak Wolfe
				v.descend(-timestep*(1-BACKSTEP_TAU)); // backstep and try again
			timestep *= BACKSTEP_TAU;
			Uf = computeTotEnergy();
		}
		
		if ((Ui - Uf)/Ui < precision) { // STEP 4: stop condition
			for (Vertex v: vertices) // if the energy isn't really changing, then we're done
				v.descend(-timestep); // just reset to before we started backtracking
			this.elasticEnergy = computeTotEnergy();
			return false;
		}
		
		Matrix sk = new Matrix(2*vertices.size(), 1); // STEP 5: save historical vector information
		for (int i = 0; i < vertices.size(); i ++) {
			sk.set(2*i+0, 0, vertices.get(i).getVelX()*timestep);
			sk.set(2*i+1, 0, vertices.get(i).getVelY()*timestep);
		}
		this.sHist.addLast(sk);
		this.gkMinus1 = gk;
		
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
		
		this.sHist = new LinkedList<Matrix>(); // with a new number of vertices, these are no longer relevant
		this.yHist = new LinkedList<Matrix>(); // erase them.
		this.gkMinus1 = null;
		
		return true;
	}
	
	
	/**
	 * Compute the total elastic energy in the system and save it as the "default" state.
	 * @return the total elastic energy.
	 */
	private double computeTotEnergy() {
		double U = 0;
		for (Cell c: getCellsUnmodifiable())
			U += c.computeAndSaveEnergy();
		return U;
	}
	
	
	/** Convert spherical coordinates to cartesian coordinates using the current mest configuration.
	 * @param lat - The latitude of the point to map
	 * @param lon - The longitude of the point to map
	 * @return an array of two elements: {x, y}
	 */
	public double[] map(double phi, double lam) {
		double i = (.5 - phi/Math.PI)*cells.length;
		double j = (1. + lam/Math.PI)*cells.length;
		int i0 = (int)Math.min(i, cells.length-1), j0 = (int)Math.min(j, cells[0].length-1);
		double di = i - i0, dj = j - j0;
		Cell c = cells[i0][j0]; // the cell is easy to find in the array
		double x =
				(1-di)*(1-dj)*c.getCorner(Cell.NW).getX() + (1-di)*(dj)*c.getCorner(Cell.NE).getX()
				+ (di)*(1-dj)*c.getCorner(Cell.SW).getX() + (di)*(dj)*c.getCorner(Cell.SE).getX();
		double y = // then just do linear interpolation
				(1-di)*(1-dj)*c.getCorner(Cell.NW).getY() + (1-di)*(dj)*c.getCorner(Cell.NE).getY()
				+ (di)*(1-dj)*c.getCorner(Cell.SW).getY() + (di)*(dj)*c.getCorner(Cell.SE).getY();
		return new double[] {x, y};
	}
	
	
	public double getTotalTearLength() {
		return this.tearLength;
	}
	
	
	/**
	 * Return the total elastic energy in the system from the last step.
	 * This differs from computeTotEnergy in that it uses a saved field, not an on-the-spot computation.
	 * This means that it is faster and less prone to oscillating.
	 * @return the saved total energy
	 */
	public double getTotEnergy() {
		return this.elasticEnergy;
	}
	
	
	/**
	 * Iterate over the vertices that compose the edge of this mesh. This method will always return
	 * Vertices in the same order, widdershins, but may add elements when tears are made.
	 * @return the iterable of edge vertices
	 */
	public Iterable<Vertex> getEdge() {
		return () -> new Iterator<Vertex>() {
			Vertex current = null;
			
			public boolean hasNext() {
				return current != edgeStart.getClockwiseNeighbor();
			}
			
			public Vertex next() {
				this.current = (current == null) ? edgeStart : current.getWidershinNeighbor();
				return current;
			}
		};
	}
	
	
	public Collection<Cell> getCellsUnmodifiable() {
		return Collections.unmodifiableCollection(
				Arrays.stream(cells).flatMap(Arrays::stream).collect(Collectors.toList()));
	}
	
	
	public Collection<Vertex> getVerticesUnmodifiable() {
		return Collections.unmodifiableCollection(this.vertices);
	}
	
	
	private Vertex getAnyEdge() {
		for (Vertex v: getVerticesUnmodifiable())
			if (v.isEdge())
				return v;
		assert false : "If you're reading this, init is doing something very wrong. There is no edge!";
		return null;
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
					Cell[][] cells, Collection<Vertex> vertices) {
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
				cells[i][j] = cell;
				
				for (int k = 0; k < 4; k ++) { // look at those vertices
					if (!vertices.contains(cell.getCorner(k))) // if we just created this one
						vertices.add(cell.getCorner(k)); // add it to the Collection
				}
			}
			
			public double cleanup() {
				for (int i = 0; i < vertexArray.length-1; i ++) {
					vertexArray[i][0].setWidershinNeighbor(vertexArray[i+1][0]);
					vertexArray[i][vertexArray[i].length-1].setClockwiseNeighbor(
							vertexArray[i+1][vertexArray[i].length-1]);
				}
				return Math.PI;
			}
		},
		
		SINUSOIDAL_FLORENCE {
			public void spawnCell(
					int i, int j, int res, double lambda, double mu, Cell[][] cells,
					Collection<Vertex> vertices) {
				// TODO: Implement this
			}
		},
		
		AZIMUTHAL {
			public void spawnCell(
					int i, int j, int res, double lambda, double mu, Cell[][] cells,
					Collection<Vertex> vertices) {
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
				Cell[][] cells, Collection<Vertex> vertices);
		
		/**
		 * Do anything that needs to be done once all the cells are spawned.
		 * @return the total tear length associated with this initial configuration
		 */
		public double cleanup() {return 0;}
		
		
		public static double[] sinusoidalProj(double[] sphereCoords) {
			return new double[] { sphereCoords[1]*Math.cos(sphereCoords[0]), sphereCoords[0] };
		}
	}
}
