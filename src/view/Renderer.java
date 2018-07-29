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
import org.jcodec.api.awt.AWTSequenceEncoder;
import org.jcodec.common.io.NIOUtils;
import org.jcodec.common.io.SeekableByteChannel;
import org.jcodec.common.model.Rational;
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
import javafx.scene.shape.Circle;
import javafx.scene.shape.Path;
import javafx.scene.shape.Polygon;
import javafx.scene.shape.Polyline;
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
	private final Shape background;
	private final Text readout;
	private final Map<Geometry, Shape> shapes;
	
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
		
		this.border = new Polygon();
//		this.border.setStroke(Color.BLACK);
//		this.border.setStrokeWidth(2);
//		this.border.setFill(null);
//		this.entities.getChildren().add(border);
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
	
	
	private Map<Geometry, Shape> createShapes(Collection<Geometry> geometries) {
		Map<Geometry, Shape> shapes = new HashMap<Geometry, Shape>();
		for (Geometry geom: geometries) {
			double[] points = new double[2*geom.getNumPoints()];
			for (int i = 0; i < points.length; i ++)
				points[i] = offset;
			
			switch (Geometries.get(geom)) {
			case POLYGON:
				Polygon pgon = new Polygon(points);
				pgon.setStroke(Color.BLACK);
				pgon.setStrokeWidth(.1);
				pgon.setFill(randomColor((String)geom.getUserData()));
				shapes.put(geom, pgon);
				break;
			case LINESTRING:
				Polyline pline = new Polyline(points);
				pline.setStroke(Color.gray(.2));
				pline.setStrokeWidth(.5);
				pline.setFill(null);
				shapes.put(geom, pline);
				break;
			case POINT:
				Circle circ = new Circle(geom.getCoordinate().x, geom.getCoordinate().y, 2);
				circ.setStroke(null);
				circ.setFill(Color.BLACK);
				shapes.put(geom, circ);
				break;
			default:
				System.err.println("I have not accounted for "+Geometries.get(geom)+"s");
				break;
			}
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
	public void render() { //TODO: render for arbitrary amounts of time (decouple the time component)
		long now = System.currentTimeMillis();
		
		double mapSize = 0;
		for (Vertex v: mesh.getVerticesUnmodifiable()) // adjust the scale based on the current map size
			mapSize = Math.max(mapSize, 2*Math.max(
					Math.abs(v.getX()), Math.abs(v.getY())));
		double c1 = Math.exp((lastRender-now)/decayTime); // the time scaling coefficient
		this.scale = c1*this.scale + (1-c1)*mapSize*1.01/size;
		
		for (Geometry geom: shapes.keySet()) {
			Shape shape = shapes.get(geom);
			List<double[]> cartCoords = Arrays.stream(geom.getCoordinates())
					.map(this::mapToMesh).collect(Collectors.toList()); // map all of the geometries to the mesh
			
			if (shape instanceof Polygon || shape instanceof Polyline) {
				List<Double> points = (shape instanceof Polygon) ?
						((Polygon)shape).getPoints() : ((Polyline)shape).getPoints();
				for (int i = 0; i < cartCoords.size(); i ++) {
					points.set(2*i+0, cartCoords.get(i)[0]);
					points.set(2*i+1, cartCoords.get(i)[1]);
				}
			}
			else if (shape instanceof Circle) {
				Circle circ = (Circle) shape;
				circ.setCenterX(cartCoords.get(0)[0]);
				circ.setCenterY(cartCoords.get(0)[1]);
			}
			else {
				assert false : "What is this";
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
	 * @throws IOException if there is a problem writing to disk
	 */
	public void saveFrame() throws IOException {
		if (entities.getScene() == null)
			return; // wait for the scene if it doesn't exist yet
		Image frame = entities.snapshot(null, null);
		BufferedImage bimg = SwingFXUtils.fromFXImage(frame, null);
		ImageIO.write(bimg, "png", new File(String.format("frames/frame%04d.png", frameNum)));

		frameNum ++;
	}
	
	
	/**
	 * Compile all frames saved this run into an AVI and save that to disk
	 * @throws IOException if there is a problem reading or writing the disk
	 */
	public void compileFrames() throws IOException {
		if (!saveImages) 	return; // don't bother if there are no frames to save
		System.out.println("Compiling frames...");
		int imgSize = this.size - 4;
		SeekableByteChannel out = null;
		try {
			out = NIOUtils.writableFileChannel(String.format("frames/convergence %s.mp4", "I"));
			AWTSequenceEncoder encoder = new AWTSequenceEncoder(out, Rational.R(25, 1));
			for (int i = 0; i < frameNum; i ++) {
				System.out.println(i);
				BufferedImage image = ImageIO.read(
						new File(String.format("frames/frame%04d.png", i)));
				image = image.getSubimage((image.getWidth()-imgSize)/2,
						(image.getHeight()-imgSize)/2, imgSize, imgSize); // crop it to size
				encoder.encodeImage(image);
			}
			encoder.finish();
			System.out.printf("Successfully compiled frames into frames/convergence %s.mp4\n", "I");
		} finally {
			NIOUtils.closeQuietly(out);
		}
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
