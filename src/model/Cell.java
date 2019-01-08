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
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import utils.Math2;

/**
 * A quadrilateral representing a section of the globe within a certain
 * latitude range and a certain longitude range, and comprising and managing
 * four triangular Elements.
 * 
 * @author Justin Kunimune
 */
public class Cell {
	
	private List<Element> elements; // the four elements in order
	
	private double phiC, lamC; // the coordinates of the centre vertex
	private final double yM, xN, xS; // some undeformed dimensions
	
	
	
	public Cell(double strength, double lambda, double mu, double scale,
			Vertex nw, Vertex ne, Vertex c, Vertex sw, Vertex se) {
		this.yM = scale/2;
		this.xN = scale/2*Math.cos(ne.getLat());
		this.xS = scale/2*Math.cos(se.getLat());
		
		this.elements = new LinkedList<Element>();
		if (se != ne)
			this.elements.add(new Element(strength, lambda, mu,
					new Vertex[] {c, se, ne}, new double[][] {{0,0}, {xS,-yM}, {xN,yM}})); // east Element
		if (ne != nw)
			this.elements.add(new Element(strength, lambda, mu,
					new Vertex[] {c, ne, nw}, new double[][] {{0,0}, {xN,yM}, {-xN,yM}})); // north Element
		if (nw != sw)
			this.elements.add(new Element(strength, lambda, mu,
					new Vertex[] {c, nw, sw}, new double[][] {{0,0}, {-xN,yM}, {-xS,-yM}})); // west Element
		if (sw != se)
			this.elements.add(new Element(strength, lambda, mu,
					new Vertex[] {c, sw, se}, new double[][] {{0,0}, {-xS,-yM}, {xS,-yM}})); // south Element
		
		this.phiC = c.getLat();
		this.lamC = c.getLon();
	}
	
	
	public Cell(double strength, double lambda, double mu, double scale,
			Vertex nw, Vertex ne, Vertex sw, Vertex se) {
		this(strength, lambda, mu, scale,
				nw, ne,
				new Vertex(nw, ne, sw, se),
				sw, se);
	}
	
	
	
	public Collection<Element> getElementsUnmodifiable() {
		return Collections.unmodifiableList(this.elements);
	}
	
	
	public Element getElement(int i) {
		return this.elements.get(i);
	}
	
	
	public double[] map(double phi, double lam) { // get the X-Y coordinates of a point in this Cell given its latitude and longitude
		double y = (phi - phiC); // is a simple one-to-one mapping
		double c = (y/yM + 1)/2.; // like y, but goes from 0 to 1
		double delLam = Math2.floorMod(lam - lamC + Math.PI, 2*Math.PI) - Math.PI;
		double x = delLam*(c*xN + (1-c)*xS)/yM; // with x we need to account for sphericalness
		for (Element e: elements)
			if (e.containsUndeformed(x, y, 0))
				return e.mapUndeformedToDeformed(x, y);
		System.out.println("stuff");
		throw new IllegalArgumentException(phi+","+lam+" is not in any elements in the cell at "+phiC+","+lamC+" "+elements);
	}
}
