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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import utils.Matrix;

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
	private final double precision; // determines how far we update before declaring that we have settled
	private final double maxTearLength; // determines when we stop tearing and declare the map done
	private List<Vertex> edge; // the start Vertex for iterating around the edge
	private double elasticEnergy; // the potential energy currently stored
	private double tearLength; // the length of the edge in radians
	
	private LinkedList<Matrix> sHist; // history of $s$ from the L-BFGS algorithm
	private LinkedList<Matrix> yHist; // history of $y$ from the L-BFGS algorithm
	private Matrix gkMinus1 = null; // the previous value of $g$ from the L-BFGS algorithm
	
	
	
	public Mesh(int resolution, InitialConfig init, double lambda, double mu,
			double precision, double maxTearLength, double[][] weights, double[][] scales) {
		this.cells = new Cell[2*resolution][4*resolution];
		this.vertices = new ArrayList<Vertex>();
		this.precision = precision;
		this.maxTearLength = maxTearLength;
		
		init.instantiate(resolution);
		for (int i = 0; i < 2*resolution; i ++)
			for (int j = 0; j < 4*resolution; j ++) // let init populate the mesh
				init.spawnCell(i, j, lambda*weights[i][j], mu*weights[i][j], scales[i][j], cells, vertices);
		this.tearLength = init.cleanup(); // let init finish up
		for (Cell c: getCellsUnmodifiable())
			for (Vertex v: c.getCornersUnmodifiable()) // make sure these relationships are mutual
				v.addNeighbor(c);
		this.edge = traceEdge();
		
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
			for (Vertex v: vertices)
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
		if (tearLength >= maxTearLength)
			return false;
		
		double maxStrain = -1; // maximise this to tear in the right place
		Vertex v0max = null, v1max = null; // the start and end locations of the tear (the parameters to maximise)
//		System.out.print("[");
		for (Vertex v0: this.edge) { // first we have to choose where to rupture
			for (Vertex v1: v0.getLinks()) {
				if (v1.isEdge())	continue; // iterate over all possible start and end points
				
				double length = v0.distanceTo(v1);
				double[] direction = {-(v0.getY()-v1.getY())/length, (v0.getX()-v1.getX())/length}; // this points widdershins, perpendicular to the potential tear
				double force = 0;
				double sign = 1; // this changes halfway through this loop here when we pass the tear
				for (Cell c: v0.getNeighborsInOrder()) { // compute the force pulling it apart
					force += sign*(v0.getForceX(c)*direction[0] + v0.getForceY(c)*direction[1]);
					if (c.getCornersUnmodifiable().contains(v1))
						sign = -1; // flip the sign when the cell is adjacent to the tear
				}
				
				double weight = Math.pow(v0.getWeight() + v1.getWeight(), 1);
				double strain = force/length/weight;
				if (strain > maxStrain) {
					maxStrain = strain;
					v0max = v0;
					v1max = v1;
				}
//				System.out.print("["+v0.getX()+","+v0.getY()+","+v1.getX()+","+v1.getY()+","+strain+"],");
			}
		}
//		System.out.println("]");
		
		if (v0max == null)
			return false;
		
		Vertex v2 = new Vertex(v0max); // split the vertex
		this.vertices.add(v2);
		for (Cell c: v0max.getNeighborsInOrder()) { // look at the cells
			v0max.transferNeighbor(c, v2); // and detach them
			if (c.getCornersUnmodifiable().contains(v1max))
				break; // until you hit the tear, anyway
		}
		
		v2.setWidershinNeighbor(v0max.getWidershinNeighbor()); // finally, update the edge chain
		v1max.setWidershinNeighbor(v2);
		v0max.setWidershinNeighbor(v1max);
		
		this.tearLength += v0max.distanceTo(v1max);
		this.edge = traceEdge(); // update the edge so that the Renderer knows about this
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
	
	
	/**
	 * Find the edge and save it as a list for thread-safe use later
	 * @return the list of Vertices in the edge. The Vertices may be modified, but the list itself
	 * will not.
	 */
	private List<Vertex> traceEdge() {
		LinkedList<Vertex> output = new LinkedList<Vertex>();
		for (Vertex v: getVerticesUnmodifiable()) {
			if (v.isEdge()) {
				output.add(v);
				break;
			}
		}
		if (output.isEmpty())
			throw new IndexOutOfBoundsException("There are no edges in this thing!");
		while (output.getLast().getWidershinNeighbor() != output.getFirst())
			output.add(output.getLast().getWidershinNeighbor());
		return output;
	}
	
	
	/** Convert spherical coordinates to cartesian coordinates using the current mesh configuration.
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
		return this.edge;
	}
	
	
	public Collection<Cell> getCellsUnmodifiable() {
		return Collections.unmodifiableCollection(
				Arrays.stream(cells).flatMap(Arrays::stream).collect(Collectors.toList()));
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
			private double lam0;
			private int res;
			private Vertex[][] vertexArray = null; // the vertex array is indexed from the edge, not the meridian
			
			public void instantiate(int res) {
				this.lam0 = params[0];
				this.res = res;
				
				vertexArray = new Vertex[2*res+1][4*res+1]; // set up the vertex array
				for (int i = 0; i <= 2*res; i ++) {
					for (int j = 0; j <= 4*res; j ++) {
						if ((i == 0 || i == 2*res) && j > 0)
							vertexArray[i][j] = vertexArray[i][0]; // make sure the poles are all one tile
						else
							vertexArray[i][j] = new Vertex( // but other than that make every vertex from scratch
									Math.PI/2 * (res - i)/res,
									Math.PI/2 * (j - 2*res)/res,
									Mesh::sinusoidalProj);
					}
				}
			}
			
			public void spawnCell(
					int i, int j, double lambda, double mu, double scale,
					Cell[][] cells, Collection<Vertex> vertices) {
				double lamC = Math.PI/2/res; // the angle associated with a single cell
				int vi = i; // the indices of the northwest vertex
				int vj = (int)Math.floorMod(j - Math.round(lam0/lamC), 4*res);
				
				Cell cell = new Cell(lambda, mu, Math.PI/2/res*Math.sqrt(scale),
						vertexArray[vi][vj], vertexArray[vi][vj+1],
						vertexArray[vi+1][vj], vertexArray[vi+1][vj+1]);
				cells[i][j] = cell;
				
				for (int k = 0; k < 4; k ++) { // look at those vertices
					if (!vertices.contains(cell.getCorner(k))) // if we just created this one
						vertices.add(cell.getCorner(k)); // add it to the Collection
				}
			}
			
			public double cleanup() {
				for (int i = 0; i < vertexArray.length-1; i ++) { // make the edges neighbours to each other
					vertexArray[i][0].setWidershinNeighbor(vertexArray[i+1][0]);
					vertexArray[i][vertexArray[i].length-1].setClockwiseNeighbor(
							vertexArray[i+1][vertexArray[i].length-1]);
				}
				return Math.PI;
			}
		},
		
		AZIMUTHAL {
			public void spawnCell(
					int i, int j, double lambda, double mu, double scale,
					Cell[][] cells, Collection<Vertex> vertices) {
				throw new IllegalArgumentException("Go away!"); // TODO: this
			}
		};
		
		
		protected double[] params;
		
		
		/**
		 * Prepare oneself to start creating cells
		 * @param res - The number of cells in this mesh from equator to pole
		 */
		public void instantiate(int res) {}
		
		/**
		 * Create a new cell, and vertices if necessary, and add all created Objects to cells and vertices.
		 * @param i - The vertical index of the desired cell, from 0 (north pole) to 2*res (south pole)
		 * @param j - The horizontal index of the desired cell, from 0 (west) to 4*res (east)
		 * @param lambda - A material constant to pass on to new cells.
		 * @param mu - A material constant to pass on to new cells.
		 * @param cells - The Collection to which to add the new Cell
		 * @param vertices - The Collection to which to add any new Vertices
		 */
		public abstract void spawnCell(int i, int j, double lambda, double mu,
				double scale, Cell[][] cells, Collection<Vertex> vertices);
		
		/**
		 * Do anything that needs to be done once all the cells are spawned.
		 * @return the total tear length associated with this initial configuration
		 */
		public double cleanup() {return 0;}
		
		public static InitialConfig fromName(String name) {
			if (name.equals("sinusoidal")) {
				SINUSOIDAL.params = new double[] {0};
				return SINUSOIDAL;
			}
			else if (name.equals("sinusoidal_florence")) {
				SINUSOIDAL.params = new double[] {Math.toRadians(12)};
				return SINUSOIDAL;
			}
			else if (name.equals("azimuthal")) {
				AZIMUTHAL.params = new double[] {Math.toRadians(90), Math.toRadians(0)};
				return AZIMUTHAL;
			}
			else if (name.equals("azimuthal_nemo")) {
				AZIMUTHAL.params = new double[] {Math.toRadians(-49), Math.toRadians(-123)};
				return AZIMUTHAL;
			}
			else
				throw new IllegalArgumentException(name);
		}
	}
	
	
	public static double[] sinusoidalProj(double[] sphereCoords) {
		return new double[] { sphereCoords[1]*Math.cos(sphereCoords[0]), sphereCoords[0] };
	}
}
