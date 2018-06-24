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

import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import model.Mesh;
import model.Vertex;

/**
 * A class devoted to making the mesh visible and pretty.
 * 
 * @author Justin Kunimune
 */
public class Renderer {
	
	private Canvas canvas;
	private GraphicsContext g;
	private Mesh mesh;
	private final double scale, offset, frameTime;
	
	private long lastFrame = 0;
	
	
	
	public Renderer(int size, double frameRate, Mesh mesh) {
		this.canvas = new Canvas(size, size);
		this.g = canvas.getGraphicsContext2D();
		this.mesh = mesh;
		this.scale = size/(2*Math.PI);
		this.offset = size/2.;
		this.frameTime = 1000/frameRate;
	}
	
	
	
	/**
	 * Get the canvas on which we will draw stuff
	 * @return canvas
	 */
	public Canvas getCanvas() {
		return this.canvas;
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
		
		g.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
		for (Vertex v0: mesh) {
			for (int i = 0; i < 2; i ++) {
				Vertex v1 = v0.getNeighbor(i);
				if (v1 == null) 	continue;
				g.moveTo(v0.getX()*scale+offset, v0.getY()*scale+offset);
				g.lineTo(v1.getX()*scale+offset, v1.getY()*scale+offset);
			}
		}
		g.stroke();
	}
	
	
	/**
	 * Save an image of the current thing to disk
	 */
	public void saveFrame() {
		// TODO: Implement this
		
	}
}
