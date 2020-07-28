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
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Set;

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
	
	private static final double ARMIJO_GOLDSTEIN_C = 0.7;
	private static final double BACKSTEP_TAU = 0.5;
	private static final double L_BFGS_M = 12; // the memory size
	
	private static final double SHEAR_WEIGHT = 0;//.167; // how much strong shear can cause tears compared to strain
	
	private final Cell[][] cells; // all of the Cells, arranged by latitude (North->South) and longitude (West->East)
	private final List<Vertex> vertices; // all of the Vertices. It doesn't matter what order this List is in, but it must have an order so that I can put them in a Vector.
	private final double precision; // determines how far we update before declaring that we have settled
	private final double maxTearLength; // determines when we stop tearing and declare the map done
	private final Set<Vertex> stitchHistory; // all of the Vertices we have stitched
	private List<Vertex> edge; // the start Vertex for iterating around the edge
	private double elasticEnergy; // the potential energy currently stored
	private double tearLength; // the length of the edge in radians
	private boolean active; // are we done yet?
	
	private LinkedList<Matrix> sHist; // history of $s$ from the L-BFGS algorithm
	private LinkedList<Matrix> yHist; // history of $y$ from the L-BFGS algorithm
	private Matrix gkMinus1 = null; // the previous value of $g$ from the L-BFGS algorithm
	
	
	
	public Mesh(int resolution, String initialCondition, double lambda, double mu,
			double precision, double maxTearLength, double[][] weights, double[][] scales,
			double eccentricity) {
		this.precision = precision;
		this.maxTearLength = maxTearLength;
		
		InitialConfig init = new InitialConfig(initialCondition, weights, scales, lambda, mu, resolution, eccentricity);
		this.vertices = new ArrayList<Vertex>(init.vertices);
		this.cells = init.cells;
		this.tearLength = init.tearLength;
		
		this.edge = traceEdge();
		this.stitchHistory = new HashSet<Vertex>();
		
		this.active = true;
		
		this.elasticEnergy = getTotEnergy(false);
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
		if (!this.isActive())	throw new IllegalStateException(); // no more updating once we've finalised
		
		double Ui = getTotEnergy(true);
		
		Matrix gk = getEnergyGradient();
		
		if (gkMinus1 != null) // STEP 5 (cont.): save historical vector information
			this.yHist.addLast(gk.minus(gkMinus1));
		if (sHist.size() > L_BFGS_M) {
			sHist.removeFirst();
			yHist.removeFirst();
		}
		
		Matrix q = gk.times(-1); // STEP 2: choose the step direction
		double[] alpha = new double[sHist.size()];
		for (int i = sHist.size()-1; i >= 0; i --) { // this is where it gets complicated
			alpha[i] = sHist.get(i).dot(q)/yHist.get(i).dot(sHist.get(i)); // see the paper cited at the top, page 779.
			q = q.plus(yHist.get(i).times(-alpha[i]));
		}
		double H0;
		if (!sHist.isEmpty())
			H0 = sHist.getLast().dot(yHist.getLast())/yHist.getLast().dot(yHist.getLast()); // this is our very rough estimate of the inverse Hessian
		else
			H0 = 1;
		Matrix dk = q.times(H0);
		for (int i = 0; i < sHist.size(); i ++) {
			double beta = yHist.get(i).dot(dk)/yHist.get(i).dot(sHist.get(i));
			dk = dk.plus(sHist.get(i).times(alpha[i]-beta));
		}
		
		double gradDotVel = gk.dot(dk);
		if (gradDotVel > 0) { // ensure this number is never positive
			System.err.printf("WARN: It tried to step uphill with g_k \\cdot d_k = %f. I don't know what that means.\n", gradDotVel);
			dk = dk.times(-1);
			gradDotVel = gk.dot(dk);
		}
		for (int i = 0; i < vertices.size(); i ++) // save the chosen step direction in the vertices
			vertices.get(i).setVel(dk.get(2*i+0, 0), dk.get(2*i+1, 0));
		
		double timestep = 1.; // STEP 3: choose the step size
		for (Vertex v: vertices)
			v.descend(timestep);
		double Uf = getTotEnergy(false);
		while ((Double.isNaN(Uf) || Uf - Ui > ARMIJO_GOLDSTEIN_C*timestep*gradDotVel)) { // if the energy didn't decrease enough
			for (Vertex v: vertices)
				v.descend(-timestep*(1-BACKSTEP_TAU)); // backstep and try again
			timestep *= BACKSTEP_TAU;
			Uf = getTotEnergy(false);
		}
		
		if ((Ui - Uf)/Ui < precision) { // STEP 4: stop condition
			for (Vertex v: vertices) // if the energy isn't really changing, then we're done
				v.descend(-timestep); // just reset to before we started backtracking
			this.elasticEnergy = getTotEnergy(false);
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
		if (!this.isActive())	throw new IllegalStateException(); // no more updating once we've finalised
		
		if (tearLength >= maxTearLength)
			return false;
		
		getTotEnergy(true);
		getEnergyGradient(); // these lines are necessary to bring the forces up to date
		
		double maxValue = -1; // maximise this to tear in the right place
		Vertex v0max = null, v1max = null; // the start and end locations of the tear (the parameters to maximise)
		for (Vertex v0: this.edge) { // first we have to choose where to rupture
			for (Vertex v1: v0.getLinks()) {
				if (v1.isEdge())	continue; // iterate over all possible start and end points
				
				double length = v0.distanceTo(v1);
				double[] direction = {-(v0.getY()-v1.getY())/length, (v0.getX()-v1.getX())/length}; // this points widdershins, perpendicular to the potential tear
				double strain = 0, shear = 0;
				double sign = 1; // this changes halfway through this loop here when we pass the tear
				for (Element c: v0.getNeighborsUnmodifiableInOrder()) { // compute the force pulling it apart
					strain += sign*(v0.getForceX(c)*direction[0] + v0.getForceY(c)*direction[1]);
					shear += sign*(v0.getForceX(c)*direction[1] - v0.getForceY(c)*direction[0]);
					if (c.getVerticesUnmodifiable().contains(v1))
						sign = -1; // flip the sign when the cell is adjacent to the tear
				}
				
				double strength = 0; // the strength is a bit weird to calculate,
				for (Element e: v0.getNeighborsUnmodifiable()) // because it's a quantity of Elements, not Vertices
					if (e.isAdjacentTo(v1))
						strength += e.getStrength()/2;
				assert strength >= 0 && strength <= 1 : strength;
				
				double stress = (strain + SHEAR_WEIGHT*Math.abs(shear))/length/strength; // divide the pressure by the strength (proportional to the Lamé params) to get deformation
				double tearValue = stress*(1 - strength); // then throw in a factor of 1-strength to prevent tears from going over continents
				if (tearValue > maxValue) {
					maxValue = tearValue;
					v0max = v0;
					v1max = v1;
				}
			}
		}
		
		if (v0max == null)
			return false;
		
		Vertex v2 = new Vertex(v0max); // split the vertex
		this.vertices.add(v2);
		for (Element c: v0max.getNeighborsUnmodifiableInOrder()) { // look at the cells
			v0max.transferNeighbor(c, v2); // and detach them
			if (c.getVerticesUnmodifiable().contains(v1max))
				break; // until you hit the tear
		}
		
		v2.setWidershinNeighbor(v0max.getWidershinNeighbor()); // finally, update the edge chain
		v1max.setWidershinNeighbor(v2);
		v0max.setWidershinNeighbor(v1max);
		
		this.tearLength += v0max.geographicDistanceTo(v1max);
		this.edge = traceEdge(); // update this.edge in a Thread-safe manner so that the Renderer knows about this
		this.sHist = new LinkedList<Matrix>(); // with a new number of vertices, these are no longer relevant
		this.yHist = new LinkedList<Matrix>(); // erase them.
		this.gkMinus1 = null;
		
		return true;
	}
	
	
	/**
	 * Find the tear ending that seems like it would most benefit from being repaired. Then, close it.
	 * @return true if it successfully stitched something that isn't in stitchHistory.
	 */
	public boolean stitch() {
		if (maxTearLength == 0) // if there was no tearing at all
			return false; // don't actually do this
		
		List<Vertex> stitchCandidates = new LinkedList<Vertex>(); // first, choose a vertex to stitch
		for (Vertex v: edge)
			if (v.getWidershinNeighbor().isSiblingOf(v.getClockwiseNeighbor())) // get all the ends of tears
				stitchCandidates.add(v);
		stitchCandidates.sort((va, vb) -> (int)Math.signum(va.getEdgeAngle() - vb.getEdgeAngle())); // and sort them by the angle they make
		stitchCandidates.removeAll(stitchHistory);
		if (stitchCandidates.isEmpty()) // quit if we've exhausted all options
			return false;
		
		Vertex v0 = stitchCandidates.get(0);
		Vertex w1 = v0.getWidershinNeighbor(), c1 = v0.getClockwiseNeighbor(); // now, begin the tear re-stitching process!
		Vertex c2 = c1.getClockwiseNeighbor();
		
		this.vertices.remove(c1); // delete c1
		this.tearLength -= v0.geographicDistanceTo(w1); // delete the tear from the total tear length
		w1.setClockwiseNeighbor(c2); // rewrite the edge chain to cut v0 out
		v0.internalise(); // make sure v0 knows of its new status
		for (Element e: c1.getNeighborsUnmodifiable(true)) // and re-attach all Elements from now nonexistent c1 to its sibling
			c1.transferNeighbor(e, w1);
		this.edge = traceEdge();
		
		this.sHist = new LinkedList<Matrix>(); // with a new number of vertices, these are no longer relevant
		this.yHist = new LinkedList<Matrix>(); // erase them.
		this.gkMinus1 = null;
		
		this.stitchHistory.add(v0); // remember this
		return true;
	}
	
	
	/**
	 * Prevent any further changes to this Mesh, and let anyone who asks know that it is no longer active.
	 */
	public void finalise() {
		this.active = false;
	}
	
	
	/**
	 * Compute the total elastic energy in the system and return it.
	 * @param propGradient - If true, it will also compute energy gradient information and save
	 * it in the Elements. This method must be called with prepGradient set to true before
	 * calling computeEnergyGradient().
	 * @return the total elastic energy.
	 */
	private double getTotEnergy(boolean prepGradient) {
		double U = 0;
		for (Element e: getElementsUnmodifiable()) {
			if (prepGradient)
				e.computeEnergyAndForce();
			else
				e.computeEnergy();
			U += e.getEnergy();
		}
		return U;
	}
	
	
	/**
	 * Compute the gradient vector of the energy with respect to all Vertex coordinates,
	 * and return it.
	 * @return the total energy gradient.
	 */
	private Matrix getEnergyGradient() {
		Matrix grad = new Matrix(2*vertices.size(), 1);
		for (int i = 0; i < vertices.size(); i ++) {
			Vertex v = vertices.get(i);
			for (Element e: v.getNeighborsUnmodifiable()) {
				double[] force = e.getForce(v);
				v.setForce(e, force[0], force[1]);
				grad.add(2*i+0, 0, -force[0]);
				grad.add(2*i+1, 0, -force[1]);
			}
		}
		return grad;
	}
	
	
	/**
	 * Compute some basic distortion criteria for this Mesh.
	 * @param weights - The array of weights to apply to each cell in the averaging
	 * @return {mean scale, StD of scale, RMS of shape} in Np.
	 */
	public double[] getCriteria(double[][] weights) {
		double scale = 0, scale2 = 0;
		double shape2 = 0;
		double total = 0;
		for (int i = 0; i < cells.length; i ++) {
			for (int j = 0; j < cells[i].length; j ++) {
				for (Element e: cells[i][j].getElementsUnmodifiable()) {
					double weight = weights[i][j]*e.getGeographicArea();
					Matrix F = e.computeAndGetDeformationGradient(true);
					Matrix B = F.times(F.T());
					double[] lambda = B.getEigenvalues();
					double a = Math.sqrt(lambda[0]), b = Math.sqrt(lambda[1]);
					scale += weight*Math.log(a*b);
					scale2 += weight*Math.pow(Math.log(a*b), 2);
					shape2 += weight*Math.pow(Math.log(a/b), 2);
					total += weight;
				}
			}
		}
		
		return new double[] {scale/total, Math.sqrt(scale2/total - scale*scale/(total*total)), Math.sqrt(shape2/total)};
	}
	
	
	/**
	 * Find the edge and save it as a list for thread-safe use later
	 * @return the list of Vertices in the edge. The Vertices may be modified, but the list itself
	 * 	will not.
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
	
	
	/**
	 * Is this point inside or outside the edge as currently traced?
	 * @param X - The X of the point to check.
	 * @param Y - The Y of the point to check.
	 * @return true if it's inside, false if it's outside. It's not precise enough to work well on the edge.
	 */
	private boolean insideEdge(double X, double Y) {
		int inside = 0;
		for (int i = 0; i < edge.size(); i ++) {
			double X0 = edge.get(i).getX(), X1 = edge.get((i+1)%edge.size()).getX(); // for each segment of the edge
			double Y0 = edge.get(i).getY(), Y1 = edge.get((i+1)%edge.size()).getY();
			if ((Y0 > Y) != (Y1 > Y)) // if the two points fall on either side of a rightward ray from (X,Y)
				if ((Y-Y0)/(Y1-Y0)*(X1-X0)+X0 > X) // and the line between them intersects our ray right of (X,Y)
					inside += Math.signum(Y1-Y0);
		}
		return inside > 0; // use the nonzero method, not the even-odd method
	}
	
	
	/** Convert spherical coordinates to Cartesian coordinates using the current mesh configuration.
	 * @param lat - The latitude of the point to map.
	 * @param lon - The longitude of the point to map.
	 * @return an array of two elements: {x, y}.
	 */
	public double[] map(double... philam) {
		double phi = philam[0], lam = philam[1];
		int i = Math.min((int)((.5 - phi/Math.PI)*cells.length), cells.length-1);
		int j = Math.min((int)((lam/Math.PI + 1.)*cells.length), cells[i].length-1);
		double phi0 = (.5 - (double)(i+1)/cells.length)*Math.PI;
		double lam0 = ((double)j/cells.length - 1.)*Math.PI;
		return cells[i][j].map(phi - phi0, lam - lam0);
	}
	
	
	/**
	 * Convert Cartesian coordinates to spherical coordinates using the current mesh configuration.
	 * Warning: this method is pretty slow. Don't use it too often.
	 * @param X - The X coordinate of the point to map.
	 * @param Y - The Y coordinate of the point to map.
	 * @return an array of two elements: {phi, lam}.
	 */
	public double[] inverseMap(double... XY) {
		double X = XY[0], Y = XY[1];
		if (insideEdge(X, Y)) { // if it's inside
			for (Element e: this.getElementsUnmodifiable()) // check each Element once
				if (e.containsDeformed(X, Y, false))
					return e.mapDeformedToSpherical(X, Y);
		} // if it wasn't in any of them, just continue
		
		final double d = 1e-4; // an arbitrarily small number, which doesn't actually have to be that small
		
		double[][] planeVecs = new double[edge.size()][2]; // if it's not in any of those,
		double[][] sfereVecs = new double[edge.size()][3]; // turn to linear extrapolation
		for (int i = 0; i < edge.size(); i ++) {
			Vertex v = edge.get(i);
			planeVecs[i] = v.getEdgeDirection(); // define the direction of the edge here
			assert insideEdge(v.getX()-d*planeVecs[i][0], v.getY()-d*planeVecs[i][1]);
			double[] inCoords = inverseMap(v.getX()-d*planeVecs[i][0], v.getY()-d*planeVecs[i][1]); // then use that to map a point slightly inside
			sfereVecs[i][0] = (Math.cos(v.getPhi())*Math.cos(v.getLam()) - Math.cos(inCoords[0])*Math.cos(inCoords[1]))/d;
			sfereVecs[i][1] = (Math.cos(v.getPhi())*Math.sin(v.getLam()) - Math.cos(inCoords[0])*Math.sin(inCoords[1]))/d; // so that we can map that edge vector to sphere space
			sfereVecs[i][2] = (Math.sin(v.getPhi())                      - Math.sin(inCoords[0]))/d; // (I could do this analytically, but I really don't want to)
		}
		
		int bestI = -1;
		double bestW = Double.NaN, bestS = Double.POSITIVE_INFINITY;
		for (int i = 0; i < edge.size(); i ++) { // now take each edge segment,
			double X0 = edge.get(i).getX(), Y0 = edge.get(i).getY(); // try to extrapolate the point on it
			double u0 = planeVecs[i][0], v0 = planeVecs[i][1];
			double X1 = edge.get((i+1)%edge.size()).getX(), Y1 = edge.get((i+1)%edge.size()).getY();
			double u1 = planeVecs[(i+1)%planeVecs.length][0], v1 = planeVecs[(i+1)%planeVecs.length][1];
			double a = (u0-u1)*(Y0-Y1) - (v0-v1)*(X0-X1); // by solving this quadratic
			double b = (v0-v1)*X - v0*X1 - v1*X0 + 2*v1*X1 - (u0-u1)*Y + u0*Y1 + u1*Y0 - 2*u1*Y1;
			double c = v1*(X-X1) - u1*(Y-Y1);
			for (double sign = -1; sign <= 1; sign += 2) {
				double w = (-b +sign* Math.sqrt(b*b - 4*a*c))/(2*a);
				double s = (X - w*X0 - (1-w)*X1)/(w*u0 + (1-w)*u1);
				if (Double.isFinite(w) && w >= 0 && w <= 1) { // if you can,
					if (Math.abs(s) < Math.abs(bestS)) { // and it's closer to this segment than to the last
						bestI = i; // save it!
						bestW = w;
						bestS = s;
					}
				}
			}
		}
		
		if (bestI >= 0) {
			double[][] vec = {sfereVecs[bestI], sfereVecs[(bestI+1)%sfereVecs.length]}; // use those extrapolation constants
			Vertex[] vtx = {edge.get(bestI), edge.get((bestI+1)%edge.size())}; // and vertices
			double[][] org = new double[2][];
			for (int i = 0; i < 2; i ++)
				org[i] = new double[] {
						Math.cos(vtx[i].getPhi())*Math.cos(vtx[i].getLam()),
						Math.cos(vtx[i].getPhi())*Math.sin(vtx[i].getLam()),
						Math.sin(vtx[i].getPhi())};
			double[] coords3d = new double[3]; // to compute its linearly extrapolated space in the 3d space the globe occupies
			for (int i = 0; i < 3; i ++)
				coords3d[i] = bestW*(org[0][i] + bestS*vec[0][i]) + (1-bestW)*(org[1][i] + bestS*vec[1][i]);
			return new double[] { // finally, gnomonically project that to the globe
					Math.atan2(coords3d[2], Math.hypot(coords3d[0], coords3d[1])),
					Math.atan2(coords3d[1], coords3d[0])};
		}
		else
			throw new IllegalArgumentException(String.format("Could not extrapolate %f, %f", X, Y)); // null means it could not be extrapolated
	}
	
	
	/**
	 * Compute and return the parameters of the smallest surrounding rectangle of this mesh.
	 * @param allowRotation - If true, the bounding box may be rotated to fit the Mesh better.
	 * @return { unrotated X centre, unrotated Y centre, rotation from horizontal, width, height }.
	 */
	public double[] getBoundingBox(boolean allowRotation) {
		double bestTheta = 0;
		double[] bestRectangle = null;
		
		if (allowRotation) { // if rotation is allowed...
			LinkedList<Vertex> hull = new LinkedList<Vertex>(); // SIGH then we have to try this with all the different thetas
			for (Vertex v: getEdge()) { // do a Graham Scan to get the convex hull
				hull.addFirst(v);
				while (hull.size() >= 3 && hull.get(1).isLeftOf(hull.get(2), hull.get(0)))
					hull.remove(1); // it's really easy, since the edge is already an approximation of the hull
			}
			
			double minArea = Double.POSITIVE_INFINITY;
			for (int i = 0; i < hull.size(); i ++) { // now for each segment of the hull
				double theta = Math.atan2(
						hull.get(i).getY()-hull.get((i+1)%hull.size()).getY(),
						hull.get(i).getX()-hull.get((i+1)%hull.size()).getX()); // take the angle
				double[] rectangle = getBoundingBox(theta); // and fit a rectangle about it
				
				if (rectangle[2]*rectangle[3] < minArea) { // finally, evaluate it on its area
					bestTheta = theta;
					bestRectangle = rectangle;
					minArea = rectangle[2]*rectangle[3];
				}
			}
		}
		else { // if not...
			bestRectangle = getBoundingBox(0); // we can just call it zero and be done!
		}
		
		double rotatedCX = bestRectangle[0], rotatedCY = bestRectangle[1], // put it in the format we want
				width = bestRectangle[2], height = bestRectangle[3];
		double correctedCX = rotatedCX*Math.cos(bestTheta) - rotatedCY*Math.sin(bestTheta);
		double correctedCY = rotatedCX*Math.sin(bestTheta) + rotatedCY*Math.cos(bestTheta);
		double[] finalRectangle = new double[] {correctedCX, correctedCY, bestTheta, width, height};
		
		if (allowRotation && finalRectangle[3] < finalRectangle[4]) { // rotate it if it's portrait
			finalRectangle[2] += Math.PI/2;
			double temp = finalRectangle[3];
			finalRectangle[3] = finalRectangle[4];
			finalRectangle[4] = temp;
		}
		finalRectangle[2] = Math2.floorMod(Math.PI/2 + finalRectangle[2], Math.PI) - Math.PI/2; // or if it's upside down
		return finalRectangle;
	}
	
	
	/**
	 * Get the parameters of the smallest rectangle that can fit this Mesh, given the
	 * angle the rectangle's base must make with the horizontal.
	 * @param rotation
	 * @return { rotated center x, rotated center y, width, height }
	 */
	public double[] getBoundingBox(double rotation) {
		double aMin = Double.POSITIVE_INFINITY, aMax = Double.NEGATIVE_INFINITY;
		double bMin = Double.POSITIVE_INFINITY, bMax = Double.NEGATIVE_INFINITY; // and fit a rectangle about it
		for (Vertex v: edge) {
			double[] ab = applyTransform(v.getX(), v.getY(), 0, 0, rotation);
			double a = ab[0], b = ab[1];
			if (a < aMin)
				aMin = a;
			if (a > aMax)
				aMax = a;
			if (b < bMin)
				bMin = b;
			if (b > bMax)
				bMax = b;
		}
		return new double[] {(aMax+aMin)/2, (bMax+bMin)/2, aMax-aMin, bMax-bMin};
	}
	
	
	/**
	 * Save this mesh to an ASCII print stream in the following format:
	 * <br>
	 * 	The first line is the comma-separated number of vertices l, height of cell table n, width of cell table m,
	 * 	number of vertices in edge e, height of pixel table o, width of pixel table p, width of map w, and height of
	 * map h.
	 * <br>
	 * 	This is followed by l rows of comma-separated x and y values for each vertex, in order.
	 * <br>
	 * 	This is followed by n*m rows of comma-separated integers, where each row is a cell (going left to right then top
	 * 	to bottom), the first integer is the slope of the Cell (-1 for split into NE and SW, 1 for split into SE and NW,
	 * 	0 for one element), and the other integers are the indices of the Vertices:
	 * 		ne,nnw,wnw,sw,sse,ese for -1;
	 * 		ene,nne,nw,wsw,ssw,sw for 1;
	 * 		ne,nw,sw,se for 0.
	 * <br>
	 * 	This is followed by a long comma-separated list of integers, which are the indices of the vertices in the edge.
	 * <br>
	 * 	This is followed by o*p rows of comma-separated floats, representing the latitude and longitude at each pixel,
	 *  or the word "NULL" if this point is not on the map.
	 * 	(going left to right, then top to bottom).
	 * @param out - the print stream to which to print all this information.
	 */
	public void save(PrintStream out) {
		double[] transform = getBoundingBox(true); // get the transform so you can apply it before you save
		double width = transform[3], height = transform[4];
		int o = (int)(Math.sqrt(height/width)*2*cells.length);
		int p = (int)(Math.sqrt(width/height)*2*cells.length);
		out.printf(Locale.US, "%d,%d,%d,%d,%d,%d,%f,%f\n",
				vertices.size(), cells.length, cells[0].length, edge.size(), o, p, width, height); // the header
		
		for (int i = 0; i < vertices.size(); i ++) { // the vertex coordinates
			double[] coords = applyTransform(vertices.get(i).getX(), vertices.get(i).getY(), transform);
			out.printf(Locale.US, "%f,%f\n", coords[0], coords[1]);
		}
		
		for (int i = 0; i < cells.length; i ++) { // the cell corners
			for (int j = 0; j < cells[i].length; j ++) {
				Cell cell = cells[i][j]; // this part is surprisingly complicated
				Element ew = cell.getElement(0), ee = null; // because I have to account for all the different possible Cell types
				if (cell.getElementsUnmodifiable().size() >= 2)
					ee = cell.getElement(1);
				Vertex[] vs;
				int shape;
				if (i == 0) { // north pole
					vs = new Vertex[] {ew.getVertex(0), ew.getVertex(0), ew.getVertex(1), ew.getVertex(2)};
					shape = 0;
				}
				else if (i == cells.length-1) { // south pole
					vs = new Vertex[] {ew.getVertex(1), ew.getVertex(2), ew.getVertex(0), ew.getVertex(0)};
					shape = 0;
				}
				else if ((i+j)%2 == 0) { // negative slopes
					vs = new Vertex[] {ee.getVertex(0), ee.getVertex(1), ew.getVertex(2), ew.getVertex(0),
							ew.getVertex(1), ee.getVertex(2)};
					shape = -1;
				}
				else { // positive slopes
					vs = new Vertex[] {ee.getVertex(1), ew.getVertex(2), ew.getVertex(0), ew.getVertex(1),
							ee.getVertex(2), ee.getVertex(0)};
					shape = +1;
				}
				String rowStr = String.format("%d,", shape);
				for (Vertex v: vs)
					rowStr += String.format("%d,", vertices.indexOf(v));
				out.printf(rowStr.substring(0, rowStr.length()-1)+"\n");
			}
		}
		
		for (Vertex v: edge) // the edge
			out.printf(Locale.US, "%d\n", vertices.indexOf(v));
		
		for (int i = 0; i < o; i ++) {
			for (int j = 0; j < p; j ++) {
				double Y = height/2 - i*height/(o-1);
				double X = j*width/(p-1) - width/2;
				double[] coords = this.inverseMap(inverseTransform(X, Y, transform));
				out.printf(Locale.US, "%f,%f\n", coords[0], coords[1]);
			}
		}
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
	
	
	public boolean isActive() {
		return this.active;
	}
	
	
	private static double[] applyTransform(double X, double Y, double... trans) {
		double cX = trans[0], cY = trans[1], th = trans[2];
		return new double[] {
				  (X - cX)*Math.cos(th) + (Y - cY)*Math.sin(th),
				- (X - cX)*Math.sin(th) + (Y - cY)*Math.cos(th) };
	}
	
	
	private static double[] inverseTransform(double X, double Y, double... trans) {
		double cX = trans[0], cY = trans[1], th = trans[2];
		return new double[] {
				X*Math.cos(th) - Y*Math.sin(th) + cX,
				X*Math.sin(th) + Y*Math.cos(th) + cY };
	}
	
	
	public Iterable<Element> getElementsUnmodifiable() {
		return new Iterable<Element>() {
			public Iterator<Element> iterator() {
				return new Iterator<Element>() {
					int i = 0, j = 0, k = 0; // indices of the next Element
					
					public boolean hasNext() {
						return i < cells.length;
					}
					
					public Element next() {
						Element next = cells[i][j].getElement(k);
						k ++;
						if (k >= cells[i][j].getElementsUnmodifiable().size()) {
							k = 0;
							j ++;
							if (j >= cells[i].length) {
								j = 0;
								i ++;
							}
						}
						return next;
					}
				};
			}
		};
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
		
		
		public InitialConfig(String name, double[][] weights, double[][] scales,
				double lambda, double mu, int res, double eccentricity) {
			if (name.equals("hammer"))
				hammerInit(0, weights, scales, lambda, mu, res, eccentricity);
			else if (name.equals("hammer_florence"))
				hammerInit(Math.toRadians(11), weights, scales, lambda, mu, res, eccentricity);
			else if (name.equals("azimuthal_nemo"))
				azimuthalInit(Math.toRadians(-49), Math.toRadians(-123), weights, scales, lambda, mu, res, eccentricity);
			else if (name.equals("azimuthal_epia"))
				azimuthalInit(Math.toRadians(45), Math.toRadians(85), weights, scales, lambda, mu, res, eccentricity);
			else if (name.equals("polar"))
				polarInit(weights, scales, lambda, mu, res, eccentricity);
			else
				throw new IllegalArgumentException(name);
		}
		
		
		/**
		 * Compute the initial values for a simple lenticular map (a Hammer projection)
		 * @param lam0 - The standard meridian.
		 * @param weights - The table of cell importances. Must be 2*res×4*res.
		 * @param scales - The table of cell size scaling factors. Must be 2*res×4*res.
		 * @param lambda - The base value for the first Lamé parameter.
		 * @param mu - The base value for the second Lamé parameter.
		 * @param res - The number of cells between the poles and the equator.
		 */
		private void hammerInit(double lam0,
				double[][] weights, double[][] scales, double lambda, double mu, int res,
				double eccentricity) {
			double lamC = Math.PI/2/res; // the angle associated with a single Cell
			lam0 = Math.round(lam0/lamC)*lamC; // round meridian to nearest cell
			System.out.printf("Initializing to hammer with central meridian %.2f°E.\n", Math.toDegrees(lam0));
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
			
			this.vertices = new ArrayList<Vertex>();
			for (Vertex[] row: vertexArray)
				for (Vertex vtx: row)
					if (!vertices.contains(vtx))
						vertices.add(vtx); // collect all Vertices in a List
			
			this.cells = new Cell[2*res][4*res];
			for (int i = 0; i < 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) { // populate the mesh with cells
					int vi = i; // the indices of the northwest vertex
					int vj = (int)Math.floorMod(Math.round(j - lam0/lamC), 4*res);
					int sign;
					if (i == 0 || i == cells.length-1)	sign = 0; // the orientation of the cell
					else								sign = ((i+j)%2 == 0) ? -1 : 1;
					
					cells[i][j] = new Cell(weights[i][j], scales[i][j],
							lambda*weights[i][j], mu*weights[i][j], lamC*Math.sqrt(scales[i][j]),
							vertexArray[vi][vj], vertexArray[vi][vj+1],
							vertexArray[vi+1][vj], vertexArray[vi+1][vj+1], sign,
							eccentricity);
				}
			}
			
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
		private void azimuthalInit(double phi0, double lam0,
				double[][] weights, double[][] scales, double lambda, double mu, int res,
				double eccentricity) {
			double size = Math.PI/2/res; // the basic angular/undeformed Cell size
			int pi = res - (int)Math.round(phi0/size); // round to the nearest joint
			int pj = 2*res + (int)Math.round(lam0/size);
			if ((pi+pj)%2 == 1)	pj --; // make sure the split point happens where there are enough vertices to handle it
			phi0 = pi*size - Math.PI/2; // and move the centre to the antipode of the given point
			lam0 = pj*size;
			System.out.printf("Initializing to azimuthal equidistant with aspect %.2f°N %.2f°E.\n", Math.toDegrees(phi0), Math.toDegrees(lam0));
			this.tearLength = size * (2 + 2*Math.cos(phi0));
				
			Vertex[][] vertexArray = new Vertex[2*res+1][4*res]; // set up the vertex array
			for (int i = 0; i <= 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) {
					if ((i == 0 || i == 2*res) && j > 0) {
						vertexArray[i][j] = vertexArray[i][0]; // make sure the poles are all one vertex
					}
					else {
						double phi = size * (res - i);
						double lam = size * (j - 2*res);
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
						double R = Math.PI/2-phi1;
						vertexArray[i][j] = new Vertex(phi, lam, R*Math.sin(lam1), -R*Math.cos(lam1)); // besides the poles, make every vertex from scratch
					}
				}
			}
			
			Vertex[] pVertices = new Vertex[8]; // fill in the special pole vertices
			double phi = vertexArray[pi][pj].getPhi();
			double lam = vertexArray[pi][pj].getLam();
			double R = Math.PI;
			for (int k = 0; k < 8; k ++) {
				double th = Math.PI - Math.PI/4*(k+.5);
				pVertices[k] = new Vertex(phi, lam, R*Math.cos(th), R*Math.sin(th));
			}
			
			List<Vertex> orderedVertices = new ArrayList<Vertex>(); // store this as a List
			vertices = orderedVertices; // as well as a Collection
			for (Vertex[] row: vertexArray) // now save them in the List for export,
				for (Vertex vtx: row)
					if (!vertices.contains(vtx))
						vertices.add(vtx);
			vertices.remove(vertexArray[pi][pj]); // being sure to use the special pVertices instead of the boring one
			vertices.addAll(Arrays.asList(pVertices)); // we put in the array
			
			this.cells = new Cell[2*res][4*res];
			for (int i = 0; i < 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) { // populate the mesh with Cells
					if ((i == pi || i+1 == pi) && (j == pj || j+1 == pj))
						continue; // those adjacent to the pole will be done separately
					
					int sign;
					if (i == 0 || i == cells.length-1)	sign = 0; // the orientation of the Cell
					else								sign = ((i+j)%2 == 0) ? -1 : 1;
					
					cells[i][j] = new Cell(weights[i][j], scales[i][j],
							lambda*weights[i][j], mu*weights[i][j], size*Math.sqrt(scales[i][j]),
							vertexArray[i][j], vertexArray[i][(j+1)%(4*res)],
							vertexArray[i+1][j], vertexArray[i+1][(j+1)%(4*res)], sign,
							eccentricity);
				}
			}
			
			cells[pi-1][pj] = new Cell(weights[pi-1][pj], scales[pi-1][pj],
					lambda*weights[pi-1][pj], mu*weights[pi-1][pj], size*Math.sqrt(scales[pi-1][pj]),
					vertexArray[pi-1][pj], vertexArray[pi-1][pj+1], vertexArray[pi-1][pj+1],
					pVertices[1], pVertices[0], vertexArray[pi][pj+1], 1,
					eccentricity); // the Cell northeast of the pole
			cells[pi-1][pj-1] = new Cell(weights[pi-1][pj-1], scales[pi-1][pj-1],
					lambda*weights[pi-1][pj-1], mu*weights[pi-1][pj-1], size*Math.sqrt(scales[pi-1][pj-1]),
					vertexArray[pi-1][pj-1], vertexArray[pi-1][pj-1], vertexArray[pi-1][pj],
					vertexArray[pi][pj-1], pVertices[3], pVertices[2], -1,
					eccentricity); // the Cell northwest of the pole
			cells[pi][pj-1] = new Cell(weights[pi][pj-1], scales[pi][pj-1],
					lambda*weights[pi][pj-1], mu*weights[pi][pj-1], size*Math.sqrt(scales[pi][pj-1]),
					vertexArray[pi][pj-1], pVertices[4], pVertices[5],
					vertexArray[pi+1][pj-1], vertexArray[pi+1][pj-1], vertexArray[pi+1][pj], 1,
					eccentricity); // the Cell southwest of the pole
			cells[pi][pj] = new Cell(weights[pi][pj], scales[pi][pj],
					lambda*weights[pi][pj], mu*weights[pi][pj], size*Math.sqrt(scales[pi][pj]),
					pVertices[6], pVertices[7], vertexArray[pi][pj+1],
					vertexArray[pi+1][pj], vertexArray[pi+1][pj+1], vertexArray[pi+1][pj+1], -1,
					eccentricity); // the Cell southeast of the pole
			
			Collections.sort(orderedVertices, (a,b) -> (int)Math.signum(a.getR()-b.getR())); // sort it by radius
			for (Vertex vtx: orderedVertices) { // finally, with the graph topology done, put the Vertices in more reasonable positions
//				double minR = vtx.getR(); // going from the centre out,
//				assert maxR < 1;
				double minR = 0; // find out if we have to move it along its azimuthal position
				for (Element e: vtx.getNeighborsUnmodifiable()) { // in order to give its neighbors positive area
					if (isLaterIn(e, orderedVertices, vtx)) 	continue; // skip any Elements with Vertices that we're going to move later
					if (!e.isInverted()) 						continue; // if e is positive, skip it
					
					Vertex[] v = e.getVerticesUnmodifiable().toArray(new Vertex[0]);
					int a = e.indexOf(vtx); // do this by computing the intersection of A-O with B-C,
					int b = (a+1)%3, c = (a+2)%3; // where A is this Vertex, O is the origin, and B and C are the other two Vertices of e
					double Rint = v[a].getR()*(v[b].getY()*v[c].getX() - v[b].getX()*v[c].getY()) / 
							(v[a].getX()*(v[b].getY()-v[c].getY()) - v[a].getY()*(v[b].getX()-v[c].getX())); // Rint is the radius at which it would correct e
					if (Rint > minR) // this part may be a bit unintuitive...
						minR = Rint; // but we're setting minR so that by keeping R above it, it will cross all Rints.
				}
				if (minR != 0) { // if you found any restrictions at all
					double finalR = minR + size/2; // either move it to minR plus an Element width or leave it where it lies,
					double factor = finalR/vtx.getR(); // whichever makes the map smaller
					vtx.setPos(factor*vtx.getX(), factor*vtx.getY()); // and vwalla! A reasonable start condition!
				}
			}
			
			for (Vertex pVtx: pVertices) { // next, define the edge, taking advantage of the fact that only pVertices have edges,
				Element neighbor = pVtx.getNeighborsUnmodifiable().iterator().next(); // and each of those only has one neighbor
				List<Vertex> tresVértices = neighbor.getVerticesUnmodifiable(); // which has three Vertices:
				int kp = tresVértices.indexOf(pVtx); // the pole,
				pVtx.setWidershinNeighbor(tresVértices.get((kp+1)%3)); // the pole's widdershins neighbor,
				pVtx.setClockwiseNeighbor(tresVértices.get((kp+2)%3)); // and the pole's clockwise neighbor
			}
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
				double[][] weights, double[][] scales, double lambda, double mu, int res,
				double eccentricity) {
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
			
			this.vertices = new ArrayList<Vertex>();
			for (Vertex[] row: vertexArray)
				for (Vertex vtx: row)
					if (!vertices.contains(vtx))
						vertices.add(vtx); // collect all Vertices in a List
			
			this.cells = new Cell[2*res][4*res];
			for (int i = 0; i < 2*res; i ++) {
				for (int j = 0; j < 4*res; j ++) { // populate the mesh with cells
					int sign;
					if (i == 0)	sign = 0; // the orientation of the Cell
					else		sign = ((i+j)%2 == 0) ? -1 : 1;
					
					cells[i][j] = new Cell(weights[i][j], scales[i][j],
							lambda*weights[i][j], mu*weights[i][j], Math.PI/2/res*Math.sqrt(scales[i][j]),
							vertexArray[i][j], vertexArray[i][(j+1)%(4*res)],
							vertexArray[i+1][j], vertexArray[i+1][(j+1)%(4*res)], sign,
							eccentricity);
				}
			}
			
			for (int j = 0; j < 4*res; j ++) { // make the edges neighbours to each other
				vertexArray[vertexArray.length-1][j].setWidershinNeighbor(
						vertexArray[vertexArray.length-1][(j+1)%(4*res)]);
			}
		}
		
		
		private boolean isLaterIn(Element e, List<Vertex> l, Vertex v) {
			for (Vertex u: e.getVerticesUnmodifiable())
				if (l.indexOf(u) > l.indexOf(v))
					return true;
			return false;
		}
	}
}
