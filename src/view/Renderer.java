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

import javafx.geometry.VPos;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.shape.Line;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
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
	private final Text readout;
	private final Map<Cell, Line[]> lines;
	
	private final int size;
	private final double scale, offset;
	private final double decayTime; // in milliseconds
	private final boolean saveImages; // save frames to disk?
	
	private Mesh mesh;
	private long lastRender;
	
	
	public Renderer(int size, Mesh mesh, double decayTime, boolean saveImages) {
		this.lines = new HashMap<Cell, Line[]>();
		this.mesh = mesh;
		this.size = size;
		this.scale = size/(2*Math.PI);
		this.offset = size/2.;
		this.decayTime = decayTime;
		this.saveImages = saveImages;
		
		this.readout = new Text(10, 0, "");
		this.readout.setTextOrigin(VPos.TOP);
		this.readout.setFont(Font.font(20));
		this.entities = new Group(this.readout);
		
		this.lastRender = System.currentTimeMillis();
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
		for (Cell c: mesh.getCellsUnmodifiable()) { // for every vertex
			if (!lines.containsKey(c))
				lines.put(c, new Line[4]); // make sure it's in the map
			
			for (int i = 0; i < 4; i ++) { // for each edge
				if (lines.get(c)[i] == null) {
					lines.get(c)[i] = new Line(offset,offset,offset,offset); // make sure the line exists for it
					entities.getChildren().add(lines.get(c)[i]);
				}
				Vertex v0 = c.getCorner(i); // get the endpoints
				Vertex v1 = c.getCorner((i+1)%4);
				Line l = lines.get(c)[i];
				
				double c1 = Math.exp((lastRender-now)/decayTime);
				double c2 = 1 - c1;
				
				l.setStartX(c1*l.getStartX() + c2*(offset+scale*v0.getX())); // update the line to its current position
				l.setStartY(c1*l.getStartY() + c2*(offset-scale*v0.getY())); // but slowly
				l.setEndX(c1*l.getEndX() + c2*(offset+scale*v1.getX()));
				l.setEndY(c1*l.getEndY() + c2*(offset-scale*v1.getY()));
			}
		}
		
		this.readout.setText(String.format("%.2fJ", mesh.getElasticEnergy()));
		this.lastRender = now;
		
		if (saveImages) {
			// TODO: rendering
		}
	}
	
	
	/**
	 * Save an image of the current thing to disk
	 */
	public void saveFrame() {
		// TODO: Implement this
		
	}
}
