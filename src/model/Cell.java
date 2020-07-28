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

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * A quadrilateral representing a section of the globe within a certain
 * latitude range and a certain longitude range, and comprising and managing
 * four triangular Elements.
 * 
 * @author Justin Kunimune
 */
public class Cell {
	
	private List<Element> elements; // the four elements in order
	
	private final int sign;
	private final double yN, yS, xN, xS; // some undeformed dimensions
	private final double phiSpan;
	
	
	
	/**
	 * Create a new Cell given the four corners.
	 * @param strength - The structural strength of this cell.
	 * @param scale - The scale factor.
	 * @param lambda - The first Lamé parameter.
	 * @param mu - The second Lamé parameter.
	 * @param size - The height of this cell in undeformed coordinates.
	 * @param nw - The northwest corner.
	 * @param ne - The northeast corner.
	 * @param sw - The southwest corner.
	 * @param se - The southeast corner.
	 * @param sign - +1 to divide it into NW and SE Elements, -1 to divide it into NE and SW,
	 * 		0 to not divide it at all.
	 */
	public Cell(double strength, double scale, double lambda, double mu, double size,
			Vertex nw, Vertex ne, Vertex sw, Vertex se, int sign, double eccentricity) {
		this(strength, scale, lambda, mu, size,
				nw, (sign < 0) ? nw : ne, ne,
				sw, (sign > 0) ? sw : se, se, sign, eccentricity);
	}
	
	
	/**
	 * Create a new Cell given the six possible Vertices (there are six because of the cut down the middle).
	 * @param strength - The structural strength of this cell.
	 * @param scale - The scale factor.
	 * @param lambda - The first Lamé parameter.
	 * @param mu - The second Lamé parameter.
	 * @param size - The height of this cell in undeformed coordinates.
	 * @param nw - The northwest corner.
	 * @param ne - The northeast corner.
	 * @param sw - The southwest corner.
	 * @param se - The southeast corner.
	 * @param sign - +1 to divide it into NW and SE Elements, -1 to divide it into NE and SW,
	 * 		0 to not divide it at all.
	 * @param eccentricity - The eccentricity of the globe.
	 */
	public Cell(double strength, double scale, double lambda, double mu, double size,
			Vertex nw, Vertex n, Vertex ne, Vertex sw, Vertex s, Vertex se, int sign,
			double eccentricity) {
		this.sign = sign;
		this.yN = size/2*(1 - Math.pow(eccentricity, 2))*Math.pow(1 - Math.pow(eccentricity*Math.sin(ne.getPhi()), 2), -3/2.);
		this.yS = size/2*(1 - Math.pow(eccentricity, 2))*Math.pow(1 - Math.pow(eccentricity*Math.sin(se.getPhi()), 2), -3/2.);
		this.xN = size/2*Math.cos(ne.getPhi())*Math.pow(1 - Math.pow(eccentricity*Math.sin(ne.getPhi()), 2), -1/2.);
		this.xS = size/2*Math.cos(se.getPhi());
		this.phiSpan = ne.getPhi()-se.getPhi(); // the angular size of this Cell
		
		this.elements = new LinkedList<Element>();
		if (sign > 0) {
			this.elements.add(new Element(strength, scale, lambda, mu,
					new Vertex[] {nw, sw, n}, new double[][] {{-xN,yN}, {-xS,-yS}, {xN,yN}})); // northwest Element
			this.elements.add(new Element(strength, scale, lambda, mu,
					new Vertex[] {se, ne, s}, new double[][] {{xS,-yS}, {xN,yN}, {-xS,-yS}})); // southeast Element
		}
		else if (sign < 0) {
			this.elements.add(new Element(strength, scale, lambda, mu,
					new Vertex[] {sw, s, nw}, new double[][] {{-xS,-yS}, {xS,-yS}, {-xN,yN}})); // southwest Element
			this.elements.add(new Element(strength, scale, lambda, mu,
					new Vertex[] {ne, n, se}, new double[][] {{xN,yN}, {-xN,yN}, {xS,-yS}})); // northeast Element
		}
		else if (ne == nw)
			this.elements.add(new Element(strength, scale, lambda, mu,
					new Vertex[] {nw, sw, se}, new double[][] {{0,yN}, {-xS,-yS}, {xS,-yS}})); // sole element
		else if (se == sw)
			this.elements.add(new Element(strength, scale, lambda, mu,
					new Vertex[] {se, ne, nw}, new double[][] {{0,-yS}, {xN,yN}, {-xN,yN}})); // sole element
		else
			throw new IllegalArgumentException(nw+","+ne+","+sw+","+se+", "+sign);
	}
	
	
	
	public Collection<Element> getElementsUnmodifiable() {
		return Collections.unmodifiableList(this.elements);
	}
	
	
	public Element getElement(int i) {
		return this.elements.get(i);
	}
	
	
	public double[] map(double delPhi, double delLam) { // get the X-Y coordinates of a point in this Cell given its relative latitude and longitude
		double y = delPhi/phiSpan*(yN + yS) - yS; // y is a simple linear mapping
		double c = delPhi/phiSpan; // (like y, but goes from 0 to 1)
		double x = (delLam/(phiSpan/2) - 1)*(c*xN + (1-c)*xS); // with x we need to account for sphericalness
		if (sign == 0 || x <= sign*(-xS + (xN+xS)/(yN+yS) * (y+yS))) // if it's in the eastern Element (or there's only one Element)
			return elements.get(0).mapUndeformedToDeformed(x, y);
		else // otherwise
			return elements.get(1).mapUndeformedToDeformed(x, y);
	}
}
