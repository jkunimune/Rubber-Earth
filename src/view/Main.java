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

import java.util.Timer;
import java.util.TimerTask;

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
	
	public static final double LAMBDA = 1e0, MU = 1.; // material properties
	public static final int MESH_RESOLUTION = 12; // the number of nodes from the equator to the pole NOTE: takes about 60 seconds to visibly converge at res 12
	public static final double PRECISION = 1e-5; // if the energy changes by less than this in one step, we're done
	public static final double TEAR_LENGTH = 2*Math.PI; // the total allowable amount of tearing
	public static final int VIEW_SIZE = 600; // size of the viewing window
	public static final double MAX_FRAME_RATE = 30; // don't render more frames than this per second
	public static final double DECAY_TIME = 1000; // the number of milliseconds that it smoothes
	public static final boolean SAVE_IMAGES = false; // save renderings as images for later processing
	public static final String[] GEO_DATA_SOURCES = {
			"ne_110m_admin_0_countries", "ne_110m_graticules_15"};
	
	public static final double INITIAL_DAMP_FACTOR = .9; // I've found this to work well experimentally
	
	private final Mesh mesh;
	private final Renderer renderer;
	private Task<Void> modelWorker;
	private ScheduledService<Void> viewWorker;
	
	
	public Main() {
		mesh = new Mesh(MESH_RESOLUTION, InitialConfig.SINUSOIDAL, LAMBDA, MU, PRECISION);
		renderer = new Renderer(VIEW_SIZE, mesh, DECAY_TIME, SAVE_IMAGES, GEO_DATA_SOURCES);
	}
	
	
	@Override
	public void start(Stage root) throws Exception {
		root.setTitle("Creating the perfect map̤…");
		root.setScene(renderer.getScene());
		
		modelWorker = new Task<Void>() {
			protected Void call() throws Exception {
				long start = System.currentTimeMillis();
				while (!isCancelled()){
					while (!isCancelled() && mesh.update()) {} // make sure you didn't miss anything
					if (mesh.getTotalTearLength() >= TEAR_LENGTH || !mesh.rupture()) // then tear
						break;
				}
				if (!isCancelled()) {
					long end = System.currentTimeMillis();
					System.out.println(String.format("It finished in %.1fs.", (end-start)/1000.));
					System.out.println(String.format("The final convergence is %.3fJ.", mesh.getTotEnergy()));
				}
				return null;
			}
			
			protected void succeeded() {
				super.succeeded();
				root.setTitle("Introducing the Danseiji IV projection!");
				new Timer().schedule(new TimerTask() {
					public void run() {
						Platform.runLater(viewWorker::cancel); // tell the viewer to stop updating after giving it a moment to settle
					}
				}, (long)(2*DECAY_TIME));
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
