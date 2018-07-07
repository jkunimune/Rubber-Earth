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
package view;

import java.util.HashMap;
import java.util.Map;

import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.shape.Line;
import model.Mesh;
import model.Vertex;

/**
 * A class devoted to making the mesh visible and pretty.
 * 
 * @author Justin Kunimune
 */
public class Renderer {
	
	private final Group entities;
	private final Map<Vertex, Line[]> lines;
	
	private final int size;
	private final double scale, offset, frameTime;
	
	private Mesh mesh;
	
	private boolean rendering = false;
	private long lastFrame = 0;
	
	
	
	public Renderer(int size, double frameRate, Mesh mesh) {
		this.entities = new Group();
		this.lines = new HashMap<Vertex, Line[]>();
		this.mesh = mesh;
		this.size = size;
		this.scale = size/(2*Math.PI);
		this.offset = size/2.;
		this.frameTime = 1000/frameRate;
	}
	
	
	
	/**
	 * Get the canvas on which we will draw stuff
	 * @return canvas
	 */
	public Scene getScene() {
		Scene scene = new Scene(entities, size, size, true);
		return scene;
	}
	
	
	/**
	 * Draw the current thing to canvas
	 */
	public void render() {
		long now = System.currentTimeMillis();
		if (now-lastFrame < frameTime)
			return; //don't render more quickly than is necessary
		else
			lastFrame = now;
		
		if (rendering) //If we are already rendering,
			return; //just give up at the expense of the frame rate
		rendering = true;
		
		for (Vertex v0: mesh) { // for every vertex
			if (!lines.containsKey(v0))
				lines.put(v0, new Line[4]); // make sure it's in the map
			
			for (int i = Vertex.NNE; i < Vertex.SSW; i ++) { // for half the directions
				if (v0.getNeighbor(i) != null) { // assuming there is a connection that way,
					Vertex v1 = v0.getNeighbor(i); // take the neighbor
					if (lines.get(v0)[i-Vertex.NNE] == null) {
						lines.get(v0)[i-Vertex.NNE] = new Line();
						entities.getChildren().add(lines.get(v0)[i-Vertex.NNE]);
					}
					Line l = lines.get(v0)[i-Vertex.NNE]; // and the line between them
					
					l.setStartX(offset+scale*v0.getX()); //update the line to its current position
					l.setStartY(offset-scale*v0.getY());
					l.setEndX(  offset+scale*v1.getX());
					l.setEndY(  offset-scale*v1.getY());
				}
				else if (lines.get(v0)[i-Vertex.NNE] != null) { // if there is no connection, but we have a line
					throw new RuntimeException("We have an extra line!");
				}
			}
		}
		rendering = false;
	}
	
	
	/**
	 * Save an image of the current thing to disk
	 */
	public void saveFrame() {
		// TODO: Implement this
		
	}
}
