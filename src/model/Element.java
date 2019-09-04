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

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import utils.Matrix;

/**
 * An independent trianglular finite element, which has three vertices that
 * may be shared with other Elements, as well as some material properties.
 * 
 * @author Justin Kunimune
 */
public class Element {
	
	private final double strength; // the max stress before tearing
	private final double scale; // the scale factor
	private final double lambda, mu; // the elastic properties
	
	private final Vertex[] vertices; //which vertex is attached to each sector (there must be exactly one)
	private final double[][] undeformedCoords; // the undeformed x-y coordinates of each vertex (each Vertex has a different undeformed coordinate in each Element to which it's connected
	private final double area; // undeformed area
	
	private double energy; // the stored elastic potential energy
	private Matrix forces; // the gradient of energy with respect to all Vertices
	private Matrix deformationGradient; // the deformation gradient
	
	
	
	public Element(double strength, double scaleFactor, double lambda, double mu,
			Vertex[] vertices, double[][] coords) {
		this.strength = strength;
		this.scale = scaleFactor;
		this.lambda = lambda;
		this.mu = mu;
		
		this.vertices = vertices;
		this.undeformedCoords = coords;
		this.area = (coords[0][0]*coords[1][1] + coords[1][0]*coords[2][1] + coords[2][0]*coords[0][1] -
				(coords[1][0]*coords[0][1] + coords[2][0]*coords[1][1] + coords[0][0]*coords[2][1])) / 2.;
		assert area >= 0;
		
		assert vertices.length == 3 : vertices.length;
		for (Vertex v: vertices) // make sure this reference goes both ways
			v.addNeighbor(this);
	}
	
	
	
	public void computeEnergy() {
		computeEnergyOrForces(true, false, false);
	}
	
	public double computeAndGetEnergy() {
		computeEnergy();
		return getEnergy();
	}
	
	public void computeForces() {
		computeEnergyOrForces(false, true, false);
	}
	
	public double[][] computeAndGetForces() {
		computeForces();
		return getForces();
	}
	
	public void computeEnergyAndForce() {
		computeEnergyOrForces(true, true, false);
	}
	
	public void computeDeformationGradient(boolean geographic) {
		computeEnergyOrForces(false, false, geographic);
	}
	
	public Matrix computeAndGetDeformationGradient(boolean geographic) {
		computeDeformationGradient(geographic);
		return getDeformationGradient();
	}
	
	/**
	 * Compute the potential energy in this Element and the force that this Element applies to each
	 * of its Vertices in this configuration. Save both values for later use.
	 * @param energy - whether computing the potential energy is worth my time
	 * @param force - whether computing the gradient is worth my time
	 * @param geographic - whether I should manually remove the scale factor
	 */
	private void computeEnergyOrForces(boolean energy, boolean force, boolean geographic) {
		Matrix F = new Matrix(2, 2); // this is the deformation gradient
		for (int i = 0; i < 3; i ++) { // it has a lot of terms
			double XA = vertices[i].getX(), YA = vertices[i].getY();
			double xB = undeformedCoords[(i+1)%3][0], yB = undeformedCoords[(i+1)%3][1];
			double xC = undeformedCoords[(i+2)%3][0], yC = undeformedCoords[(i+2)%3][1];
			F = F.plus(new Matrix(2, 2, // so I populate it in this for-loop
					XA*(xB - xC), XA*(yB - yC),
					YA*(xB - xC), YA*(yB - yC)).over(2*area));
		}
		if (geographic)
			F = F.times(scale);
		
		Matrix gradF = new Matrix(2, 3); // this is the derivative of that with respect to the Vertex coordinates
		for (int i = 0; i < 2; i ++) {
			for (int j = 0; j < 3; j ++)
				gradF.set(i, j, (undeformedCoords[(j+1)%3][i] - undeformedCoords[(j+2)%3][i])/(2*area));
		}
		
		double J = F.det();
		if (energy) {
			Matrix B = F.times(F.T()); // the rest is fancy Neo-Hookean stuff
			double i1 = B.tr();
			this.energy = (mu/2*(i1 - 2 - 2*Math.log(J)) + lambda/2*Math.pow(Math.log(J), 2)) * area; // don't forget to multiply energy density by undeformed volume
		}
		if (force) {
			Matrix Ⅎ = new Matrix(0.,1.,-1.,0.).times(F).times(new Matrix(0.,-1.,1.,0.));
			this.forces = ((F.minus(Ⅎ.over(J)).times(mu)).plus(Ⅎ.times(Math.log(J)/J*lambda))).times(gradF).times(-area); // don't forget the negative sign, since F = - gradU
		}
		this.deformationGradient = F;
	}
	
	
	public double getEnergy() {
		return this.energy;
	}
	
	public double[][] getForces() {
		return this.forces.T().getValues();
	}
	
	public double[] getForce(Vertex v) {
		int i = this.getVerticesUnmodifiable().indexOf(v);
		return new double[] {this.forces.get(0, i), this.forces.get(1, i)};
	}
	
	public Matrix getDeformationGradient() {
		return this.deformationGradient;
	}
	
	/**
	 * Does this element contain these undeformed coordinates?
	 * @param x - the undeformed x to test
	 * @param y - the undeformed y to test
	 * @return true if these undeformed coordinates fall inside the undeformed element, false otherwise
	 */
	public boolean containsUndeformed(double x, double y) {
		return containsUndeformed(x, y, -1);
	}
	
	/**
	 * Does this element contain these undeformed coordinates if we open up one side and
	 * make it an intersection of two hemiplanes?
	 * @param x - the undeformed x to test
	 * @param y - the undeformed y to test
	 * @param openSide - the index of the Vertex across from the open side
	 * @return true if these undeformed coordinates fall inside the open undeformed element, false otherwise
	 */
	public boolean containsUndeformed(double x, double y, int openSide) {
		double xa = undeformedCoords[0][0], ya = undeformedCoords[0][1];
		double xb = undeformedCoords[1][0], yb = undeformedCoords[1][1];
		double xc = undeformedCoords[2][0], yc = undeformedCoords[2][1];
		double[] w = new double[3];
		double denom = (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc);
		w[0] = ((yb-yc)*(x-xc) - (xb-xc)*(y-yc)) / denom; // simple barycentric coordinates
		w[1] = ((yc-ya)*(x-xc) - (xc-xa)*(y-yc)) / denom;
		w[2] = 1 - w[0] - w[1];
		
		for (int i = 0; i < 3; i ++) // if any of the barycentric coords are less than 0
			if (i != openSide && w[i] < 0)
				return false;
		return true;
	}
	
	/**
	 * Does this element contain these undeformed coordinates?
	 * @param x - the undeformed x to test.
	 * @param y - the undeformed y to test.
	 * @param edgeBleed - If two vertices are edges, ignore the edge between them.
	 * @return true if these undeformed coordinates fall inside the open undeformed element, false otherwise
	 */
	public boolean containsDeformed(double x, double y, boolean edgeBleed) {
		int openSide = -1;
		if (edgeBleed) {
			for (int i = 0; i < 3; i ++) {
				if (vertices[i].isEdge() && vertices[(i+1)%3].isEdge())
					openSide = (i+2)%3;
			}
		}
		double Xa = vertices[0].getX(), Ya = vertices[0].getY();
		double Xb = vertices[1].getX(), Yb = vertices[1].getY();
		double Xc = vertices[2].getX(), Yc = vertices[2].getY();
		double[] w = new double[3];
		double denom = (Yb-Yc)*(Xa-Xc) - (Xb-Xc)*(Ya-Yc);
		w[0] = ((Yb-Yc)*(x-Xc) - (Xb-Xc)*(y-Yc)) / denom; // simple barycentric coordinates
		w[1] = ((Yc-Ya)*(x-Xc) - (Xc-Xa)*(y-Yc)) / denom;
		w[2] = 1 - w[0] - w[1];
		
		for (int i = 0; i < 3; i ++) // if any of the barycentric coords are less than 0
			if (i != openSide && w[i] < -1e-12) // (accounting for roundoff)
				return false;
		return true;
	}
	
	/** Linearly interpolate a pair of X-Y coordinates from the given x-y coordinates.
	 * @param x - the input x
	 * @param y - the input y
	 * @return {output X, output Y}
	 */
	public double[] mapUndeformedToDeformed(double x, double y) {
		double xa = undeformedCoords[0][0], ya = undeformedCoords[0][1], Xa = vertices[0].getX(), Ya = vertices[0].getY();
		double xb = undeformedCoords[1][0], yb = undeformedCoords[1][1], Xb = vertices[1].getX(), Yb = vertices[1].getY();
		double xc = undeformedCoords[2][0], yc = undeformedCoords[2][1], Xc = vertices[2].getX(), Yc = vertices[2].getY();
		
		double denom = (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc); // simple barycentric coordinates
		double wa = ((yb-yc)*(x-xc) - (xb-xc)*(y-yc)) / denom;
		double wb = ((yc-ya)*(x-xc) - (xc-xa)*(y-yc)) / denom;
		double wc = 1 - wa - wb;
		
		return new double[] {wa*Xa + wb*Xb + wc*Xc, wa*Ya + wb*Yb + wc*Yc};
	}
	
	/** Linearly interpolate a pair of x-y from the given X-Y coordinates.
	 * @param X - the input X
	 * @param Y - the input Y
	 * @return {output x, output y}
	 */
	public double[] mapDeformedToSpherical(double X, double Y) {
		double pa = vertices[0].getPhi(), la = vertices[0].getLam(), Xa = vertices[0].getX(), Ya = vertices[0].getY();
		double pb = vertices[1].getPhi(), lb = vertices[1].getLam(), Xb = vertices[1].getX(), Yb = vertices[1].getY();
		double pc = vertices[2].getPhi(), lc = vertices[2].getLam(), Xc = vertices[2].getX(), Yc = vertices[2].getY();
		double xa = Math.cos(pa)*Math.cos(la), ya = Math.cos(pa)*Math.sin(la), za = Math.sin(pa);
		double xb = Math.cos(pb)*Math.cos(lb), yb = Math.cos(pb)*Math.sin(lb), zb = Math.sin(pb);
		double xc = Math.cos(pc)*Math.cos(lc), yc = Math.cos(pc)*Math.sin(lc), zc = Math.sin(pc);
		
		double denom = (Yb-Yc)*(Xa-Xc) - (Xb-Xc)*(Ya-Yc); // simple barycentric coordinates
		double wa = ((Yb-Yc)*(X-Xc) - (Xb-Xc)*(Y-Yc)) / denom;
		double wb = ((Yc-Ya)*(X-Xc) - (Xc-Xa)*(Y-Yc)) / denom;
		double wc = 1 - wa - wb;
		
		double x = wa*xa + wb*xb + wc*xc, y = wa*ya + wb*yb + wc*yc, z = wa*za + wb*zb + wc*zc; // now we have 3D spherical
		return new double[] {Math.atan2(z, Math.hypot(x, y)), Math.atan2(y, x)};
	}
	
	
	public boolean isDegenerate() {
		return vertices[0] == vertices[1] || vertices[1] == vertices[2] || vertices[2] == vertices[0];
	}
	
	
	public boolean isInverted() {
		return    vertices[0].getX()*vertices[1].getY()
				+ vertices[1].getX()*vertices[2].getY()
				+ vertices[2].getX()*vertices[0].getY()
				- vertices[1].getX()*vertices[0].getY()
				- vertices[2].getX()*vertices[1].getY()
				- vertices[0].getX()*vertices[2].getY() <= 0;
	}
	
	
	public double getStrength() {
		return this.strength;
	}
	
	
	public double getCX() {
		double x = 0;
		for (Vertex v: vertices)
			x += v.getX();
		return x/4;
	}
	
	public double getCY() {
		double y = 0;
		for (Vertex v: vertices)
			y += v.getY();
		return y/4;
	}
	
	public double getScale() {
		return this.scale;
	}
	
	public double getUndeformedArea() {
		return this.area;
	}
	
	public double getGeographicArea() {
		return this.area/this.scale;
	}
	
	
	public Vertex getVertex(int i) {
		return this.vertices[i];
	}
	
	public void setVertex(int i, Vertex corner) {
		this.vertices[i] = corner;
	}
	
	public List<Vertex> getVerticesUnmodifiable() {
		return Collections.unmodifiableList(Arrays.asList(this.vertices));
	}
	
	public int indexOf(Vertex v) {
		return this.getVerticesUnmodifiable().indexOf(v);
	}
	
	
	public double[] getUndeformedPos(Vertex v) {
		return this.undeformedCoords[this.indexOf(v)];
	}
	
	
	public boolean isAdjacentTo(Element that) { // kitty-corner cells don't count
		Set<Vertex> sharedVertices = new HashSet<Vertex>();
		for (Vertex v: this.getVerticesUnmodifiable())
			if (that.getVerticesUnmodifiable().contains(v))
				sharedVertices.add(v); // it's kind of annoying that I need a Hash and can't use an int,
		return sharedVertices.size() >= 2; // but it doesn't count as two shared corners if they're the same Vertex
	}
	
	public boolean isAdjacentTo(Vertex v) {
		return this.getVerticesUnmodifiable().contains(v);
	}
	
	
	@Override
	public String toString() {
		String s = "Element(";
		for (Vertex v: this.vertices)
			s += v+", ";
		return s.substring(0,s.length()-2) + ")";
	}
}
