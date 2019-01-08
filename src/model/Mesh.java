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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import utils.Math2;
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
	
	private static final double SHEAR_WEIGHT = 0;//.167; // how much strong shear can cause tears compared to strain
	
	private final Cell[][] cells; // all of the Cells, arranged by latitude (North->South) and longitude (West->East)
	private final List<Vertex> vertices; // all of the Vertices. It doesn't matter what order this List is in, but it must have an order so that I can put them in a Vector.
	private final double precision; // determines how far we update before declaring that we have settled
	private final double maxTearLength; // determines when we stop tearing and declare the map done
	private List<Vertex> edge; // the start Vertex for iterating around the edge
	private double elasticEnergy; // the potential energy currently stored
	private double tearLength; // the length of the edge in radians
	
	private LinkedList<Matrix> sHist; // history of $s$ from the L-BFGS algorithm
	private LinkedList<Matrix> yHist; // history of $y$ from the L-BFGS algorithm
	private Matrix gkMinus1 = null; // the previous value of $g$ from the L-BFGS algorithm
	
	
	
	public Mesh(int resolution, String initialCondition, double lambda, double mu,
			double precision, double maxTearLength, double[][] weights, double[][] scales) {
		this.precision = precision;
		this.maxTearLength = maxTearLength;
		
		InitialConfig init = new InitialConfig(initialCondition, weights, scales, lambda, mu, resolution);
		this.vertices = new ArrayList<Vertex>(init.vertices);
		this.cells = init.cells;
		this.tearLength = init.tearLength;
		
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
	public boolean update() {
		double Ui = this.elasticEnergy;
		
		Matrix gk = new Matrix(2*vertices.size(), 1); // STEP 1: compute the gradient
		for (int i = 0; i < vertices.size(); i ++) {
			Vertex v = vertices.get(i);
			for (Element e: v.getNeighborsUnmodifiable()) {
				v.stepX(STEP);
				double gradX = e.computeDeltaEnergy()/STEP; // compute the force by computing the energy gradient
				v.stepX(-STEP);
				v.stepY(STEP);
				double gradY = e.computeDeltaEnergy()/STEP;
				v.stepY(-STEP);
				v.setForce(e, -gradX, -gradY); // store the force from each cell individually for strain calculations later
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
		assert gradDotVel < 0 : gradDotVel;
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
		System.out.println(Uf);
		return true;
	}
	
	
	/**
	 * Find the vertex with the highest strain, and separate it into two vertices.
	 * @return true if it successfully tore, false if it could find nothing to tear
	 */
	public boolean rupture() {
		if (tearLength >= maxTearLength)
			return false;
		
		double maxValue = -1; // maximise this to tear in the right place
		Vertex v0max = null, v1max = null; // the start and end locations of the tear (the parameters to maximise)
//		System.out.print("[");
		for (Vertex v0: this.edge) { // first we have to choose where to rupture
			for (Vertex v1: v0.getLinks()) {
				if (v1.isEdge())	continue; // iterate over all possible start and end points
				
				double length = v0.distanceTo(v1);
				double[] direction = {-(v0.getY()-v1.getY())/length, (v0.getX()-v1.getX())/length}; // this points widdershins, perpendicular to the potential tear
				double strain = 0, shear = 0;
				double sign = 1; // this changes halfway through this loop here when we pass the tear
				for (Element c: v0.getNeighborsInOrder()) { // compute the force pulling it apart
					strain += sign*(v0.getForceX(c)*direction[0] + v0.getForceY(c)*direction[1]);
					shear += sign*(v0.getForceX(c)*direction[1] - v0.getForceY(c)*direction[0]);
					if (c.getVerticesUnmodifiable().contains(v1))
						sign = -1; // flip the sign when the cell is adjacent to the tear
				}
				
				double strength = (v0.getStrength() + v1.getStrength())/2;
				assert strength >= 0 && strength < 1 : strength;
				double stress = (strain + SHEAR_WEIGHT*Math.abs(shear))/length/strength; // divide the pressure by the strength (proportional to the Lamé params) to get deformation
				double tearValue = stress*(1 - strength); // then throw in a factor of 1-strength to prevent tears from going over continents
				if (tearValue > maxValue) {
					maxValue = tearValue;
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
		for (Element c: v0max.getNeighborsInOrder()) { // look at the cells
			v0max.transferNeighbor(c, v2); // and detach them
			if (c.getVerticesUnmodifiable().contains(v1max))
				break; // until you hit the tear, anyway
		}
		
		v2.setWidershinNeighbor(v0max.getWidershinNeighbor()); // finally, update the edge chain
		v1max.setWidershinNeighbor(v2);
		v0max.setWidershinNeighbor(v1max);
		
		this.tearLength += v0max.distanceTo(v1max); // TODO: unstretched length
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
			for (Element e: c.getElementsUnmodifiable())
				U += e.computeAndSaveEnergy();
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
		assert !output.isEmpty();
		while (output.getLast().getWidershinNeighbor() != output.getFirst())
			output.add(output.getLast().getWidershinNeighbor());
		return output;
	}
	
	
	/** Convert spherical coordinates to Cartesian coordinates using the current mesh configuration.
	 * @param lat - The latitude of the point to map
	 * @param lon - The longitude of the point to map
	 * @return an array of two elements: {x, y}
	 */
	public double[] map(double phi, double lam) {
		int i = Math.min((int)((.5 - phi/Math.PI)*cells.length), cells.length-1);
		int j = Math.min((int)((1. + lam/Math.PI)*cells.length), cells[i].length-1);
		return cells[i][j].map(phi, lam);
	}
	
	
	/**
	 * Compute and return the smallest surrounding rectangle of this mesh.
	 * @return { X centre, Y centre, rotation from horizontal, width, height }
	 */
	public double[] getLinearTransform() {
		LinkedList<Vertex> hull = new LinkedList<Vertex>();
		for (Vertex v: getEdge()) { // do a Graham Scan to get the convex hull
			hull.addFirst(v);
			while (hull.size() >= 3 && hull.get(1).isLeftOf(hull.get(2), hull.get(0)))
				hull.remove(1); // it's really easy, since the edge is already an approximation of the hull
		}
		
		double minArea = Double.POSITIVE_INFINITY;
		double[] bestRectangle = null;
		for (int i = 0; i < hull.size(); i ++) { // now for each segment of the hull
			double theta = Math.atan2(
					hull.get(i).getY()-hull.get((i+1)%hull.size()).getY(),
					hull.get(i).getX()-hull.get((i+1)%hull.size()).getX()); // take the angle
			double aMin = Double.POSITIVE_INFINITY, aMax = Double.NEGATIVE_INFINITY;
			double bMin = Double.POSITIVE_INFINITY, bMax = Double.NEGATIVE_INFINITY; // and fit a rectangle about it
			for (Vertex v: hull) {
				double a = v.getTransformedX(0, 0, theta), b = v.getTransformedY(0, 0, theta);
				if (a < aMin)
					aMin = a;
				if (a > aMax)
					aMax = a;
				if (b < bMin)
					bMin = b;
				if (b > bMax)
					bMax = b;
			}
			
			if ((aMax - aMin)*(bMax - bMin) < minArea) { // finally, evaluate it on its area
				minArea = (aMax - aMin)*(bMax - bMin); // if it passes our test,
				double ca = (aMax+aMin)/2, cb = (bMax+bMin)/2; // find the centre
				double cx = ca*Math.cos(theta) - cb*Math.sin(theta);
				double cy = ca*Math.sin(theta) + cb*Math.cos(theta);
				bestRectangle = new double[] {cx, cy, theta, aMax-aMin, bMax-bMin};
			}
		}
		
		if (bestRectangle[3] < bestRectangle[4]) { // rotate it if it's portrait
			bestRectangle[2] += Math.PI/2;
			double temp = bestRectangle[3];
			bestRectangle[3] = bestRectangle[4];
			bestRectangle[4] = temp;
		}
		bestRectangle[2] = Math2.floorMod(Math.PI/2 + bestRectangle[2], Math.PI) - Math.PI/2; // or if it's upside down
		return bestRectangle;
	}
	
	
	/**
	 * Save this mesh to an ASCII print stream in the following format:
	 * 
	 *      The first line is the comma-separated number of vertices l, height of cell table n, and width of cell table m,
	 *      the width of the map w, and the height of the map h.
	 * 
	 *      This is followed by l rows of comma-separated x and y values for each vertex, in order.
	 * 
	 *      This is followed by n*m rows of comma-separated integers, where each row is a cell (going left to right then right
	 *      to left), and the four integers are the four vertices listed counterclockwise from NE.
	 * 
	 *      This is followed by a long comma-separated list of integers, which are the indices of the vertices in the edge.
	 * @param out - the print stream to which to print all this information.
	 */
	public void save(PrintStream out) { // TODO inverse
		double[] transform = getLinearTransform(); // get the transform so you can apply it before you save
		out.printf("%d,%d,%d,%f,%f,\n", vertices.size(), cells.length, cells[0].length, transform[3], transform[4]); // the header
		for (int i = 0; i < vertices.size(); i ++) // the vertex coordinates
			out.printf("%f,%f,\n",
					vertices.get(i).getTransformedX(transform), vertices.get(i).getTransformedY(transform));
//		for (int y = 0; y < cells.length; y ++) { // the cell corners
//			for (int x = 0; x < cells[y].length; x ++) {
//				for (int i = 0; i < 4; i ++)
//					for (int j = 0; j < 3; j ++)
//						out.printf("%d,", vertices.indexOf(cells[y][x].getElement(i).getVertex(j)));
//				out.printf("\n");
//			}
//		} // TODO: figure out a good way to do this with triangular Elements
		for (Vertex v: edge) // the edge
			out.printf("%d,", vertices.indexOf(v));
		out.printf("\n");
		out.close();
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
	public class InitialConfig {
		
		protected double[] params;
		public Cell[][] cells; // array of cells in order
		public Collection<Vertex> vertices; // list of all vertices
		public double tearLength; // initial amount of tear
		
		
		public InitialConfig(String name, double[][] weights, double[][] scales, double lambda, double mu, int res) {
			if (name.equals("sinusoidal"))
				hammerInit(0, weights, scales, lambda, mu, res);
			else if (name.equals("sinusoidal_florence"))
				hammerInit(Math.toRadians(11), weights, scales, lambda, mu, res);
			else if (name.equals("azimuthal_nemo"))
				azimuthalInit(Math.toRadians(-49), Math.toRadians(-123), weights, scales, lambda, mu, res);
			else if (name.equals("azimuthal_epia"))
				azimuthalInit(Math.toRadians(45), Math.toRadians(85), weights, scales, lambda, mu, res);
			else if (name.equals("polar"))
				polarInit(weights, scales, lambda, mu, res);
			else
				throw new IllegalArgumentException(name);
		}
		
		
		/**
		 * Compute the initial values for a simple lenticular map (a Hammer projection)
		 * @param lam0 - The standard parallel.
		 * @param weights - The table of cell importances. Must be 2*res×4*res.
		 * @param scales - The table of cell size scaling factors. Must be 2*res×4*res.
		 * @param lambda - The base value for the first Lamé parameter.
		 * @param mu - The base value for the second Lamé parameter.
		 * @param res - The number of cells between the poles and the equator.
		 */
		private void hammerInit(double lam0,
				double[][] weights, double[][] scales, double lambda, double mu, int res) {
			double lamC = Math.PI/2/res; // the angle associated with a single Cell
			lam0 = Math.round(lam0/lamC)*lamC; // round meridian to nearest cell
			this.tearLength = Math.PI;
			
			Vertex[][] vertexArray = new Vertex[2*res+1][4*res+1]; // set up the vertex array
			for (int i = 0; i <= 2*res; i ++) {
				for (int j = 0; j <= 4*res; j ++) {
					if ((i == 0 || i == 2*res) && j > 0) {
						vertexArray[i][j] = vertexArray[i][0]; // make sure the poles are all one tile
					}
					else {
						double phi = lamC * (res - i);
						double lam = lamC * (j - 2*res);
						double z = Math.sqrt(1+Math.cos(phi)*Math.cos(lam/2));
						double x = Math.sqrt(8)*Math.cos(phi)*Math.sin(lam/2)/z;
						double y = Math.sqrt(2)*Math.sin(phi)/z;
						vertexArray[i][j] = new Vertex(phi, lam+lam0, x, y); // but other than that make every vertex from scratch
					}
				}
			}
			
			this.cells = new Cell[2*res][4*res];
			for (int i = 0; i < 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) { // populate the mesh with cells
					int vi = i; // the indices of the northwest vertex
					int vj = (int)Math.floorMod(Math.round(j - lam0/lamC), 4*res);
					cells[i][j] = new Cell(weights[i][j], lambda*weights[i][j],
							mu*weights[i][j], lamC*Math.sqrt(scales[i][j]),
							vertexArray[vi][vj], vertexArray[vi][vj+1],
							vertexArray[vi+1][vj], vertexArray[vi+1][vj+1]);
				}
			}
			
			this.vertices = new HashSet<Vertex>();
			for (Cell[] row: cells)
				for (Cell c: row)
					for (Element e: c.getElementsUnmodifiable())
						vertices.addAll(e.getVerticesUnmodifiable()); // collect all Vertices in a List
			
			for (int i = 0; i < vertexArray.length-1; i ++) { // make the edges neighbours to each other
				vertexArray[i][0].setWidershinNeighbor(vertexArray[i+1][0]);
				vertexArray[i][vertexArray[i].length-1].setClockwiseNeighbor(
						vertexArray[i+1][vertexArray[i].length-1]);
			}
		}
		
		
		/**
		 * Compute the initial values for a simple lenticular map (a Hammer projection)
		 * @param lam0 - The standard parallel.
		 * @param weights - The table of cell importances. Must be 2*res×4*res.
		 * @param scales - The table of cell size scaling factors. Must be 2*res×4*res.
		 * @param lambda - The base value for the first Lamé parameter.
		 * @param mu - The base value for the second Lamé parameter.
		 * @param res - The number of cells between the poles and the equator.
		 */
		private void azimuthalInit(double phiP, double lamP,
				double[][] weights, double[][] scales, double lambda, double mu, int res) {
			int pi = res - (int)Math.round(phiP/(Math.PI/2/res)); // round to the nearest joint
			int pj = 2*res + (int)Math.round(lamP/(Math.PI/2/res));
			double phi0 = pi*(Math.PI/2/res) - Math.PI/2; // and move the centre to the antipode of the given point
			double lam0 = pj*(Math.PI/2/res);
			this.tearLength = Math.PI/2/res * (2 + 2*Math.cos(phi0));
				
			Vertex[][] vertexArray = new Vertex[2*res+1][4*res]; // set up the vertex array
			for (int i = 0; i <= 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) {
					if ((i == 0 || i == 2*res) && j > 0) {
						vertexArray[i][j] = vertexArray[i][0]; // make sure the poles are all one vertex
					}
					else {
						double phi = Math.PI/2/res * (res - i);
						double lam = Math.PI/2/res * (j - 2*res);
						double phi1 = Math.asin(Math.sin(phi0)*Math.sin(phi) + Math.cos(phi0)*Math.cos(phi)*Math.cos(lam0-lam)); // relative latitude
						double lam1 = Math.acos((Math.cos(phi0)*Math.sin(phi) - Math.sin(phi0)*Math.cos(phi)*Math.cos(lam0-lam))/Math.cos(phi1))-Math.PI; // relative longitude
						if (Double.isNaN(lam1)) {
							if ((Math.cos(lam0-lam) >= 0 && phi < phi0) || (Math.cos(lam0-lam) < 0 && phi < -phi0))
								lam1 = 0;
							else
								lam1 = -Math.PI;
						}
						else if (Math.sin(lam - lam0) > 0) // it's a plus-or-minus arccos.
							lam1 = -lam1;
						double r = Math.tan((Math.PI/2-phi1)/2); // stereographic, to keep lines from crossing badly
						double x = r*Math.sin(lam1);
						double y =-r*Math.cos(lam1);
						vertexArray[i][j] = new Vertex(phi, lam, x, y); // but other than that make every vertex from scratch
					}
				}
			}
			
			Vertex[] pVertices = new Vertex[4];
			double phi = vertexArray[pi][pj].getLat();
			double lam = vertexArray[pi][pj].getLon();
			double R = Math.tan((Math.PI-Math.PI/4/res)/2);
			for (int k = 0; k < 4; k ++) {
				double th = Math.PI/2*(-k-.5);
				pVertices[k] = new Vertex(phi, lam, R*Math.cos(th), R*Math.sin(th)); // fill in the special pole vertices
			}
			vertexArray[pi][pj] = null; // we will henceforth never use this instance and want to avoid doing so accidentally
			
			this.cells = new Cell[2*res][4*res];
			for (int i = 0; i < 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) { // populate the mesh with cells
					Vertex[] vertices = { vertexArray[i][(j+1)%(4*res)], vertexArray[i][j],
							vertexArray[i+1][j], vertexArray[i+1][(j+1)%(4*res)]};
					for (int k = 0; k < 4; k ++)
						if (vertices[k] == null) // if one of these corners is the special pole one
							vertices[k] = pVertices[k];
					cells[i][j] = new Cell(weights[i][j], lambda*weights[i][j],
							mu*weights[i][j], Math.PI/2/res*Math.sqrt(scales[i][j]),
							vertices[1], vertices[0],
							vertices[2], vertices[3]);
				}
			}
			
			this.vertices = new HashSet<Vertex>();
			for (Cell[] row: cells)
				for (Cell c: row)
					for (Element e: c.getElementsUnmodifiable())
						vertices.addAll(e.getVerticesUnmodifiable()); // collect all Vertices in a List
			
			pVertices[0].setWidershinNeighbor(vertexArray[pi][pj-1]); // finally, define the edge
			pVertices[0].setClockwiseNeighbor(vertexArray[pi+1][pj]);
			pVertices[1].setWidershinNeighbor(vertexArray[pi+1][pj]);
			pVertices[1].setClockwiseNeighbor(vertexArray[pi][pj+1]);
			pVertices[2].setWidershinNeighbor(vertexArray[pi][pj+1]);
			pVertices[2].setClockwiseNeighbor(vertexArray[pi-1][pj]);
			pVertices[3].setWidershinNeighbor(vertexArray[pi-1][pj]);
			pVertices[3].setClockwiseNeighbor(vertexArray[pi][pj-1]);
		}
		
		
		/**
		 * Compute the initial values for a simple polar azimuthal equidistant map
		 * @param weights - The table of cell importances. Must be 2*res×4*res.
		 * @param scales - The table of cell size scaling factors. Must be 2*res×4*res.
		 * @param lambda - The base value for the first Lamé parameter.
		 * @param mu - The base value for the second Lamé parameter.
		 * @param res - The number of cells between the poles and the equator.
		 */
		private void polarInit(
				double[][] weights, double[][] scales, double lambda, double mu, int res) {
			this.tearLength = 0;
			
			Vertex[][] vertexArray = new Vertex[2*res+1][4*res]; // set up the vertex array
			for (int i = 0; i <= 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) {
					if (i == 0 && j > 0) {
						vertexArray[i][j] = vertexArray[i][0]; // make sure the North Pole is all one tile
					}
					else {
						double phi = Math.PI/2/res * (res - i);
						double lam = Math.PI/2/res * (j - 2*res);
						vertexArray[i][j] = new Vertex(phi, lam,
								(Math.PI/2-phi)*Math.sin(lam), -(Math.PI/2-phi)*Math.cos(lam)); // but other than that make every vertex from scratch
					}
				}
			}
			
			this.cells = new Cell[2*res][4*res];
			for (int i = 0; i < 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) { // populate the mesh with cells
					cells[i][j] = new Cell(weights[i][j], lambda*weights[i][j],
							mu*weights[i][j], Math.PI/2/res*Math.sqrt(scales[i][j]),
							vertexArray[i][j], vertexArray[i][(j+1)%(4*res)],
							vertexArray[i+1][j], vertexArray[i+1][(j+1)%(4*res)]);
				}
			}
			
			this.vertices = new HashSet<Vertex>();
			for (Cell[] row: cells)
				for (Cell c: row)
					for (Element e: c.getElementsUnmodifiable())
						vertices.addAll(e.getVerticesUnmodifiable()); // collect all Vertices in a List
			
			for (int j = 0; j < 4*res; j ++) { // make the edges neighbours to each other
				vertexArray[vertexArray.length-1][j].setWidershinNeighbor(
						vertexArray[vertexArray.length-1][(j+1)%(4*res)]);
			}
		}
	}
}
