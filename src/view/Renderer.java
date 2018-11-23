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

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.geometry.jts.Geometries;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.filter.Filter;

import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.VPos;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.image.Image;
import javafx.scene.paint.Color;
import javafx.scene.shape.LineTo;
import javafx.scene.shape.MoveTo;
import javafx.scene.shape.Path;
import javafx.scene.shape.Polygon;
import javafx.scene.shape.Rectangle;
import javafx.scene.shape.Shape;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import model.Mesh;
import model.Vertex;

/**
 * A class devoted to making the mesh visible and pretty.
 * 
 * @author Justin Kunimune
 */
public class Renderer {
	
	private final Group entities;
	private final Polygon border;
	private final Rectangle background;
	private final Text readout;
	private final Map<Geometry, Path> shapes;
	
	private final int size;
	private final double decayTime; // in milliseconds
	private final boolean saveImages; // save frames to disk?
	
	private double scale, offset;
	private Mesh mesh;
	private long lastRender;
	private int frameNum = 0;
	
	
	public Renderer(int size, Mesh mesh, double decayTime, boolean saveImages, String[] shpFiles) {
		this.mesh = mesh;
		this.size = size;
		this.scale = 2*Math.PI/size;
		this.offset = size/2.;
		this.decayTime = decayTime;
		this.saveImages = saveImages;
		this.entities = new Group();
		
		List<Geometry> geoData = tryLoadShapefile(shpFiles);
		this.shapes = createShapes(geoData);
		for (Geometry geom: geoData) // make sure to go from bottom to top here
			this.entities.getChildren().add(shapes.get(geom));
		
		this.border = new Polygon(); // this will later get subtracted from background
		
		this.background = new Rectangle(size, size, Color.WHITE);
		this.entities.getChildren().add(background);
		
		this.readout = new Text(10, 0, "");
		this.readout.setTextOrigin(VPos.TOP);
		this.readout.setFont(Font.font(20));
		this.entities.getChildren().add(readout);
		
		this.lastRender = System.currentTimeMillis();
		this.render();
	}
	
	
	/**
	 * Attempt to load the shapefile located at the filename and return the DataStore therein.
	 * @param pathname - The full pathname of the shapefile. Pretty self-explanatory.
	 * @return the DataStore of the shapefile or an empty DataStore if you couldn't load it.
	 */
	private List<Geometry> tryLoadShapefile(String[] filenames) {
		List<Geometry> output = new LinkedList<Geometry>();
		for (String filename: filenames) {
			DataStore dataStore;
			String[] typeNames;
			try {
				Map<String, Object> params = Collections.singletonMap(
						"url", new File("data\\"+filename+".shp").toURI().toURL());
				dataStore = DataStoreFinder.getDataStore(params);
				typeNames = dataStore.getTypeNames();
			} catch (IOException e) {
				System.err.println("Could not read from data\\"+filename+".shp: "+e.getMessage());
				continue;
			}
			
			for (String typeName: typeNames) {
				try {
					SimpleFeatureCollection features =
							dataStore.getFeatureSource(typeName).getFeatures(Filter.INCLUDE);
					try (SimpleFeatureIterator iterator = features.features()) {
						while (iterator.hasNext()) {
							SimpleFeature f = iterator.next();
							Geometry topGeom = (Geometry)f.getDefaultGeometry();
							for (int i = 0; i < topGeom.getNumGeometries(); i ++) { // why is this iteration not built in
								topGeom.getGeometryN(i).setUserData(f.getID());
								output.add(topGeom.getGeometryN(i));
							}
						}
					}
				} catch (IOException e) {
					System.err.println("Could not read "+typeNames+"' from data\\"+filename+".shp: "+e.getMessage());
					continue;
				} catch (RuntimeException e) {
					System.err.println("Could not read "+typeNames+"' from data\\"+filename+".shp: "+e.getMessage());
				}
			}
		}
		return output;
	}
	
	
	private Map<Geometry, Path> createShapes(Collection<Geometry> geometries) {
		Map<Geometry, Path> shapes = new HashMap<Geometry, Path>();
		for (Geometry geom: geometries) {
			double[] points = new double[2*geom.getNumPoints()];
			for (int i = 0; i < points.length; i ++)
				points[i] = offset;
			
			Path pgon = new Path();
			pgon.getElements().add(new MoveTo(0, 0));
			for (int i = 1; i < geom.getNumPoints(); i ++)
				pgon.getElements().add(new LineTo(0, 0)); // TODO: do I need a closepath?
			
			if (Geometries.get(geom) == Geometries.POLYGON) { // formatting depends on whether its a polygon
				pgon.setStroke(Color.BLACK);
				pgon.setStrokeWidth(.1);
				pgon.setFill(randomColor((String)geom.getUserData()));
			}
			else if (Geometries.get(geom) == Geometries.LINESTRING) { // or polyline
				pgon.setStroke(Color.gray(.2));
				pgon.setStrokeWidth(.5);
				pgon.setFill(null);
			}
//			case POINT:
//				continue;
			else {
				System.err.println("I have not accounted for "+Geometries.get(geom)+"s");
			}
			
			shapes.put(geom, pgon);
		}
		return shapes;
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
		
		double mapSize = 0;
		for (Vertex v: mesh.getVerticesUnmodifiable()) // adjust the scale based on the current map size
			mapSize = Math.max(mapSize, 2*Math.max(
					Math.abs(v.getX()), Math.abs(v.getY())));
		double c1 = Math.exp((lastRender-now)/decayTime); // the time scaling coefficient
		this.scale = Math.min(c1*this.scale + (1-c1)*mapSize*1.01/size, this.scale*1.1);
		
		for (Geometry geom: shapes.keySet()) {
			Path shape = shapes.get(geom);
			List<double[]> cartCoords = Arrays.stream(geom.getCoordinates())
					.map(this::mapToMesh).collect(Collectors.toList()); // map all of the geometries to the mesh
			
			for (int i = 0; i < cartCoords.size(); i ++) {
				if (shape.getElements().get(i) instanceof MoveTo) {
					((MoveTo)shape.getElements().get(i)).setX(cartCoords.get(i)[0]);
					((MoveTo)shape.getElements().get(i)).setY(cartCoords.get(i)[1]);
				}
				else if (shape.getElements().get(i) instanceof LineTo) {
					((LineTo)shape.getElements().get(i)).setX(cartCoords.get(i)[0]);
					((LineTo)shape.getElements().get(i)).setY(cartCoords.get(i)[1]);
				}
			}
		}
		
		int i = 0;
		for (Vertex v: mesh.getEdge()) { // plot the edge
			if (border.getPoints().size() < 2*i+2)
				border.getPoints().addAll(0., 0.); // make sure the polygon is big enough
			double[] xy = transform(v.getX(), v.getY());
			this.border.getPoints().set(2*i+0, xy[0]);
			this.border.getPoints().set(2*i+1, xy[1]);
			i ++;
		}
		
		Shape maskedRect = Path.subtract(background, border); // the background masks stuff outside the map
		maskedRect.setFill(Color.WHITE);
		maskedRect.setStroke(Color.BLACK);
		maskedRect.setStrokeWidth(2.);
		this.entities.getChildren().set(this.entities.getChildren().size()-2, maskedRect);
		
		this.readout.setText(String.format("%.3fJ", mesh.getTotEnergy()));
		this.lastRender = now;
		
		if (saveImages) {
			try {
				try {
					saveFrame();
				} catch (FileNotFoundException e) { // how
					System.err.println("Could not save frame: "+e.getMessage()); // can
				} // it
			} catch (IOException e) { // possibly
				System.err.println("Could not save frame: "+e.getMessage()); // be throwing
				e.printStackTrace(); // a FileNotFoundException
			} // how does it abort
		} // without saving an exception
	}
	
	
	/**
	 * Save an image of the current thing to disk
	 * @param filepath - The name of the file to which to save the image
	 * @throws IOException if there is a problem writing to disk
	 */
	public void saveImage(String filepath) throws IOException {
		Image frame = entities.snapshot(null, null);
		BufferedImage bimg = SwingFXUtils.fromFXImage(frame, null);
		bimg = bimg.getSubimage((bimg.getWidth()-size)/2, (bimg.getHeight()-size)/2, size, size); // crop it to size
		ImageIO.write(bimg, "png", new File(filepath));
	}
	
	
	/**
	 * Save an image of the thing at this instant in time and add it to the frame collection
	 * @throws IOException if there is a problem writing to disk
	 */
	public void saveFrame() throws IOException {
		if (entities.getScene() == null)
			return; // wait for the scene if it doesn't exist yet
		saveImage(String.format("frames/frame%04d.png", frameNum));
		frameNum ++;
	}
	
	
	public double[] mapToMesh(Coordinate coords) {
		double[] cartesian = this.mesh.map(
				Math.toRadians(coords.y),
				Math.toRadians(coords.x));
		return transform(cartesian[0], cartesian[1]);
	}
	
	
	/**
	 * Convert math coordinates to screen coordinates
	 * @param mathX - The x-coordinate where right is positive and order unity is 1
	 * @param mathY - The y-coordinate where up is positive and order unity is 1
	 * @return double[2] { screenX, screenY };
	 */
	public double[] transform(double mathX, double mathY) {
		return new double[] {
				offset + mathX/scale,
				offset - mathY/scale };
	}
	
	
	/**
	 * How many frames have we saved?
	 * @return the number of frames
	 */
	public int getNumFrames() {
		return this.frameNum;
	}
	
	
	/**
	 * Generates a random colour in a pseudorandom fashion.
	 * @param seed - The same seed will always return the same colour
	 * @return
	 */
	private Color randomColor(String seed) {
		char s0 = seed.charAt(seed.length()-2), s1 = seed.charAt(seed.length()-1);
		double r = Math.pow(2, 16)*s0 + s1;
		double hue = 9000.*Math.tan(Math.pow(r*2.7, 7.2))%1;
		double sat = Math.sin(r);
		double brt = Math.sin(r*2.7);
		return Color.hsb(180+hue*180, sat/20+.5, brt/20.+.75);
	}
}
