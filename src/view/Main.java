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

import javafx.application.Application;
import javafx.concurrent.Task;
import javafx.stage.Stage;
import model.Mesh;
import model.Mesh.InitialConfiguration;


/**
 * The main class for running and stuff.
 * 
 * @author Justin Kunimune
 */
public final class Main extends Application {
	
	public static final int MESH_RESOLUTION = 3;
	public static final int VIEW_SIZE = 600;
	public static final double MAX_FRAME_RATE = 30;
	public static final double TIME_STEP = 0.001;
	
	private final Mesh mesh;
	private final Renderer renderer;
	private Task<Void> task;
	
	
	public Main() {
		mesh = new Mesh(MESH_RESOLUTION, InitialConfiguration.SINUSOIDAL);
		renderer = new Renderer(VIEW_SIZE, MAX_FRAME_RATE, mesh);
	}
	
	
	@Override
	public void start(Stage root) throws Exception {
		root.setTitle("Creating the perfect map…");
		root.setScene(renderer.getScene());
		renderer.render();
		root.show();
		
		task = new Task<Void>() {
			protected Void call() throws Exception {
				int i = 0;
				while (!isCancelled() && !mesh.isDone()) {
					mesh.update(TIME_STEP);
					renderer.render();
					if (i%10 == 0)
						renderer.saveFrame();
					i ++;
				}
				return null;
			}
			
			protected void succeeded() {
				super.succeeded();
				System.out.println("Done!");
				root.setTitle("Introducing the Danseiji I projection!");
			}
			
			protected void cancelled() {
				super.cancelled();
				System.out.println("Cancelled!");
			}
			
			protected void failed() {
				super.failed();
				task.getException().printStackTrace(System.err);
				System.out.println("Failed!");
			}
		};
		
		new Thread(task).start();
	}
	
	
	@Override
	public void stop() {
		task.cancel();
	}
	
	
	public static void main(String[] args) {
		launch(args);
	}
}
