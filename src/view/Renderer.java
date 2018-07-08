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
import model.Cell;
import model.Mesh;
import model.Vertex;

/**
 * A class devoted to making the mesh visible and pretty.
 * 
 * @author Justin Kunimune
 */
public class Renderer {
	
	private final Group entities;
	private final Map<Cell, Line[]> lines;
	
	private final int size;
	private final double scale, offset;
	
	private Mesh mesh;
	
	
	public Renderer(int size, Mesh mesh) {
		this.entities = new Group();
		this.lines = new HashMap<Cell, Line[]>();
		this.mesh = mesh;
		this.size = size;
		this.scale = size/(2*Math.PI);
		this.offset = size/2.;
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
		for (Cell c: mesh.getCellsUnmodifiable()) { // for every vertex
			if (!lines.containsKey(c))
				lines.put(c, new Line[4]); // make sure it's in the map
			
			for (int i = 0; i < 4; i ++) { // for each edge
				if (lines.get(c)[i] == null) {
					lines.get(c)[i] = new Line(); // make sure the line exists for it
					entities.getChildren().add(lines.get(c)[i]);
				}
				Vertex v0 = c.getCorner(i); // get the endpoints
				Vertex v1 = c.getCorner((i+1)%4);
				
				lines.get(c)[i].setStartX(offset+scale*v0.getX()); //update the line to its current position
				lines.get(c)[i].setStartY(offset-scale*v0.getY());
				lines.get(c)[i].setEndX(  offset+scale*v1.getX());
				lines.get(c)[i].setEndY(  offset-scale*v1.getY());
			}
		}
	}
	
	
	/**
	 * Save an image of the current thing to disk
	 */
	public void saveFrame() {
		// TODO: Implement this
		
	}
}
