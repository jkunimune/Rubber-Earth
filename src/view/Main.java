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
import javafx.application.Platform;
import javafx.concurrent.ScheduledService;
import javafx.concurrent.Task;
import javafx.stage.Stage;
import javafx.util.Duration;
import model.Mesh;
import model.Mesh.InitialConfig;


/**
 * The main class for running and stuff.
 * 
 * @author Justin Kunimune
 */
public final class Main extends Application {
	
	public static final double LAMBDA = 10., MU = 1.; // material properties
	public static final int MESH_RESOLUTION = 12; // the number of nodes from the equator to the pole NOTE: takes about 60 seconds to visibly converge at res 12
	public static final double PRECISION = 1e-4; // if the mean squared speed does not exceed this, we're done
	public static final double TEAR_LENGTH = 1.5; // the total allowable amount of tearing
	public static final int VIEW_SIZE = 600; // size of the viewing window
	public static final double MAX_FRAME_RATE = 30; // don't render more frames than this per second
	public static final double DECAY_TIME = 500; // the number of milliseconds that it smoothes
	public static final boolean SAVE_IMAGES = false; // save renderings as images for later processing
	
	public static final double INITIAL_DAMP_FACTOR = .9; // I've found this to work well experimentally
	
	private final Mesh mesh;
	private final Renderer renderer;
	private Task<Void> modelWorker;
	private ScheduledService<Void> viewWorker;
	
	
	public Main() {
		mesh = new Mesh(MESH_RESOLUTION, InitialConfig.SINUSOIDAL, LAMBDA, MU, PRECISION);
		renderer = new Renderer(VIEW_SIZE, mesh, DECAY_TIME, SAVE_IMAGES);
	}
	
	
	@Override
	public void start(Stage root) throws Exception {
		root.setTitle("Creating the perfect map…");
		root.setScene(renderer.getScene());
		
		modelWorker = new Task<Void>() {
			protected Void call() throws Exception {
				long start = System.currentTimeMillis();
				while (!isCancelled()){
					while (!isCancelled() && mesh.update(INITIAL_DAMP_FACTOR)) {} // descend until you can descend no more
					while (!isCancelled() && mesh.update(0)) {} // make sure you didn't miss anything
					if (mesh.getTotalTearLength() >= TEAR_LENGTH || !mesh.rupture()) // then tear
						break;
				}
				long end = System.currentTimeMillis();
				System.out.println("Final convergence of "+Math.round((mesh.getElasticEnergy())));
				System.out.println("Finished in "+Math.round((end-start)/100.)/10.+"s.");
				return null;
			}
			
			protected void succeeded() {
				super.succeeded();
				root.setTitle("Introducing the Danseiji IV projection!");
			}
			
			protected void failed() {
				super.failed();
				this.getException().printStackTrace(System.err);
			}
		};
		
		viewWorker = new ScheduledService<Void>() {
			protected Task<Void> createTask() {
				return new Task<Void>() {
					protected Void call() throws Exception {
						renderer.render();
						return null;
					}
				};
			}
			
			protected void failed() {
				super.failed();
				this.getException().printStackTrace(System.err);
			}
		};
		viewWorker.setPeriod(Duration.seconds(1./MAX_FRAME_RATE));;
		viewWorker.setExecutor(Platform::runLater);
		
		new Thread(modelWorker).start();
		viewWorker.start();
		root.show();
	}
	
	
	@Override
	public void stop() {
		modelWorker.cancel();
	}
	
	
	public static void main(String[] args) {
		launch(args);
	}
}
