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

import java.util.Arrays;
import java.util.HashSet;


/**
 * A simple interface to make these generics easier with which to deal.
 * Represents all of the mesh vertices from the same point on the globe.
 * 
 * @author Justin Kunimune
 */
public class VertexSet extends HashSet<Vertex> {
	
	private static final long serialVersionUID = -8067299072982490603L;
	
	public static final int NORTHEAST = 0, NORTHWEST = 1, SOUTHWEST = 2, SOUTHEAST = 3;
	
	private final Vertex[] attatchedTo = new Vertex[4]; //which vertex is attatched to each sector (there must be exactly one)
	
	
	public VertexSet(Vertex vertex) {
		super(Arrays.asList(vertex));
		for (int i = 0; i < 4; i ++)
			attatchedTo[i] = vertex; //only one vertex
	}
	
	
	public VertexSet(Vertex ne, Vertex nw, Vertex sw, Vertex se) {
		super(Arrays.asList(ne, nw, sw, se));
		attatchedTo[NORTHEAST] = ne;
		attatchedTo[NORTHWEST] = nw;
		attatchedTo[SOUTHWEST] = sw;
		attatchedTo[SOUTHEAST] = se;
	}
}
