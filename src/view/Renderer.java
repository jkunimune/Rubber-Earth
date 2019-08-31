/**
 * This is free and unencumbered software released into the public domain.
 * 
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 * 
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 * 
 * For more information, please refer to <http://unlicense.org>
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
import java.util.LinkedHashMap;
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
import javafx.geometry.Rectangle2D;
import javafx.geometry.VPos;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.image.Image;
import javafx.scene.paint.Color;
import javafx.scene.shape.ClosePath;
import javafx.scene.shape.LineTo;
import javafx.scene.shape.MoveTo;
import javafx.scene.shape.Path;
import javafx.scene.shape.PathElement;
import javafx.scene.shape.Polygon;
import javafx.scene.shape.Rectangle;
import javafx.scene.shape.Shape;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import model.Element;
import model.Mesh;
import model.Vertex;

/**
 * A class devoted to making the mesh visible and pretty.
 * 
 * @author Justin Kunimune
 */
public class Renderer {
	
	private static final double MAX_SEGMENT_LENGTH = 20;
	
	private final Group entities;
	private final Rectangle background;
	private final Group mask;
	private final Text leftText, rightText;
	private final Map<Element, Polygon> meshShapes;
	private final Map<Geometry, Path> geoShapes;
	
	private final int size, margin; // the window size and margin width
	private final double decayTime; // in seconds
	private final boolean saveImages; // save frames to disk?
	
	private double viewX, viewY, viewTh, viewW; // the dimensions of the viewbox
	private Mesh mesh;
	private long lastRender;
	private int frameNum = 0;
	
	
	public Renderer(int size, int margin, Mesh mesh, double decayTime, double initViewSize,
			boolean drawMesh, boolean saveImages, String[] shpFiles,
			double lambda, double mu, double maxTear) {
		this.mesh = mesh;
		this.size = size;
		this.margin = margin;
		this.decayTime = decayTime;
		this.saveImages = saveImages;
		
		this.viewX = 0;
		this.viewY = 0;
		this.viewTh = 0;
		this.viewW = initViewSize;
		
		this.entities = new Group();// make sure to go from bottom to top here
		
		List<Geometry> geoData = tryLoadShapefile(shpFiles);
		this.geoShapes = createGeoShapes(geoData);
		this.entities.getChildren().addAll(geoShapes.values());
		
		this.background = new Rectangle(size+2*margin, size, Color.WHITE);
		this.mask = new Group();
		this.entities.getChildren().add(mask);
		
//		Rectangle left = new Rectangle(0, 0, margin, size);
//		left.setFill(Color.BLACK);
//		Rectangle right = new Rectangle(margin+size, 0, margin, size);
//		right.setFill(Color.BLACK);
//		this.entities.getChildren().addAll(left, right);
		
		if (drawMesh) {
			this.meshShapes = createMeshShapes(mesh.getElementsUnmodifiable());
			this.entities.getChildren().addAll(meshShapes.values());
		}
		else {
			this.meshShapes = Collections.emptyMap();
		}

		this.leftText = new Text(10, 0, "");
		this.leftText.setTextOrigin(VPos.TOP);
		this.leftText.setFont(Font.font(24));
		this.leftText.setText(String.format("R = %.3f m\nλ = %.3f Pa\nμ = %.3f Pa\nL* = %.3f m",
				1., lambda, mu, maxTear));
		this.rightText = new Text(margin+size+10, 0, "");
		this.rightText.setTextOrigin(VPos.TOP);
		this.rightText.setFont(Font.font(24));
		this.entities.getChildren().addAll(leftText, rightText);
		
		this.lastRender = System.currentTimeMillis();
		this.render();
	}
	
	
	/**
	 * Attempt to load the shapefile located at the filename and return the DataStore therein.
	 * @param pathname - The full pathname of the shapefile. Pretty self-explanatory.
	 * @return the DataStore of the shapefile or an empty DataStore if you couldn't load it.
	 */
	private static List<Geometry> tryLoadShapefile(String[] filenames) {
		List<Geometry> output = new LinkedList<Geometry>();
		for (int i = filenames.length-1; i >= 0; i --) { // for each shapefile to read (reversed because of the end of this method)
			String filename = filenames[i];
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
			
			for (String typeName: typeNames) { // for each type of feature in it
				try {
					SimpleFeatureCollection features =
							dataStore.getFeatureSource(typeName).getFeatures(Filter.INCLUDE);
					try (SimpleFeatureIterator iterator = features.features()) { // iterate over all the features
						while (iterator.hasNext()) {
							SimpleFeature f = iterator.next();
							Geometry topGeom = (Geometry)f.getDefaultGeometry();
							for (int j = 0; j < topGeom.getNumGeometries(); j ++) { // for ecah entiti in that entity
								topGeom.getGeometryN(j).setUserData(f.getID()); // (why is this iteration not built in)
								output.add(topGeom.getGeometryN(j));
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
		Collections.reverse(output); // do this just because Antarctica starts with A, so this puts it on the bottom
		return output;
	}
	
	
	private static Map<Geometry, Path> createGeoShapes(Collection<Geometry> geometries) {
		Map<Geometry, Path> shapes = new LinkedHashMap<Geometry, Path>();
		for (Geometry geom: geometries) {
			double[] points = new double[2*geom.getNumPoints()];
			for (int i = 0; i < points.length; i ++)
				points[i] = 0;
			
			Path pgon = new Path();
			pgon.getElements().add(new MoveTo(0, 0));
			for (int i = 1; i < geom.getNumPoints(); i ++)
				pgon.getElements().add(new LineTo(0, 0));
			
			if (Geometries.get(geom) == Geometries.POLYGON) { // formatting depends on whether its a polygon
				pgon.getElements().add(new ClosePath());
				pgon.setStroke(Color.BLACK);
				pgon.setStrokeWidth(.1);
				pgon.setFill(randomColor((String)geom.getUserData()));
			}
			else if (Geometries.get(geom) == Geometries.LINESTRING) { // or polyline
				pgon.setStroke(Color.gray(.2));
				pgon.setStrokeWidth(.5);
				pgon.setFill(null);
			}
			else {
				System.err.println("I have not accounted for "+Geometries.get(geom)+"s");
			}
			
			shapes.put(geom, pgon);
		}
		return shapes;
	}
	
	
	private static Map<Element, Polygon> createMeshShapes(Iterable<Element> elements) {
		Map<Element, Polygon> shapes = new HashMap<Element, Polygon>();
		for (Element elem: elements) {
			Polygon pgon = new Polygon(0,0, 0,0, 0,0);
			pgon.setFill(null);
			pgon.setStroke(Color.BLACK);
			pgon.setStrokeWidth(0.5);
			shapes.put(elem, pgon);
		}
		return shapes;
	}
	
	
	
	/**
	 * Get the canvas on which we will draw stuff
	 * @return canvas
	 */
	public Scene getScene() {
		Scene scene = new Scene(entities, size+2*margin, size, true); // TODO: expand window
		return scene;
	}
	
	
	/**
	 * Draw the current thing to canvas
	 */
	public void render() {
		long now = System.currentTimeMillis();
		double c1 = 1 - Math.exp((lastRender-now)/(decayTime*1000)); // compute the time scaling coefficient
		double[] meshBox = mesh.getBoundingBox(!mesh.isActive());

		this.viewX = (1-c1)*viewX + c1*meshBox[0]; // figure out the current camera coordinates
		this.viewY = (1-c1)*viewY + c1*meshBox[1];
		this.viewTh = (1-c1)*viewTh + c1*meshBox[2];
		this.viewW = (1-c1)*viewW + c1*Math.max(meshBox[3], meshBox[4])*1.01;
		
		for (Geometry geom: geoShapes.keySet()) {
			Path shape = geoShapes.get(geom);
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
			
			if (shape.getElements().get(shape.getElements().size()-1) instanceof LineTo) { // and if it is an open Path
				for (int i = 1; i < cartCoords.size(); i ++) {
					PathElement ele0 = shape.getElements().get(i-1), ele1 = shape.getElements().get(i);
					double x0 = (ele0 instanceof LineTo) ? ((LineTo)ele0).getX() : ((MoveTo)ele0).getX();
					double y0 = (ele0 instanceof LineTo) ? ((LineTo)ele0).getY() : ((MoveTo)ele0).getY();
					double x1 = (ele1 instanceof LineTo) ? ((LineTo)ele1).getX() : ((MoveTo)ele1).getX();
					double y1 = (ele1 instanceof LineTo) ? ((LineTo)ele1).getY() : ((MoveTo)ele1).getY();
					if (ele1 instanceof LineTo &&
							Math.hypot(x1 - x0, y1 - y0) > MAX_SEGMENT_LENGTH) {
						shape.getElements().set(i, new MoveTo(x1, y1)); // cut any components that are too long
					}
				}
			}
		}
		
		for (Element elem: meshShapes.keySet()) { // display the mesh directly, if desired
			Polygon shape = meshShapes.get(elem);
			for (int i = 0; i < 3; i ++) {
				double[] coords = transform(elem.getVertex(i).getX(), elem.getVertex(i).getY());
				shape.getPoints().set(2*i+0, coords[0]);
				shape.getPoints().set(2*i+1, coords[1]);
			}
		}
		
		Polygon border = new Polygon();
		for (Vertex v: mesh.getEdge()) { // plot the edge
			double[] xy = transform(v.getX(), v.getY());
			border.getPoints().addAll(xy[0], xy[1]);
		}
		
		Shape maskedRect = Path.subtract(background, border); // the background masks stuff outside the map
		maskedRect.setFill(Color.WHITE);
		maskedRect.setStroke(Color.BLACK);
		maskedRect.setStrokeWidth(2.);
		this.mask.getChildren().setAll(maskedRect);
		
		this.rightText.setText(String.format("U = %.3f J\nL = %.3f m",
				mesh.getTotEnergy(), mesh.getTotalTearLength()));
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
	 * @param includeText - Should the margins be included?
	 * @throws IOException if there is a problem writing to disk
	 */
	public void saveImage(String filepath, boolean includeText) throws IOException {
		SnapshotParameters params = new SnapshotParameters();
		if (includeText)
			params.setViewport(new Rectangle2D(0, 0, size+2*margin, size));
		else
			params.setViewport(new Rectangle2D(margin, 0, size, size));
		Image frame = entities.snapshot(params, null);
		BufferedImage bimg = SwingFXUtils.fromFXImage(frame, null);
		ImageIO.write(bimg, "png", new File(filepath));
	}
	
	
	/**
	 * Save an image of the thing at this instant in time and add it to the frame collection
	 * @throws IOException if there is a problem writing to disk
	 */
	public void saveFrame() throws IOException {
		if (entities.getScene() == null)
			return; // wait for the scene if it doesn't exist yet
		saveImage(String.format("frames/frame%04d.png", frameNum), true);
		frameNum ++;
	}
	
	
	public double[] mapToMesh(Coordinate coords) {
		double[] cartesian = this.mesh.map(
				Math.toRadians(coords.y),
				Math.toRadians(coords.x));
		return transform(cartesian[0], cartesian[1]);
	}
	
	
	/**
	 * Convert math coordinates to screen coordinates, going through the mesh's linear transform
	 * @param mathX - The x-coordinate where right is positive and order unity is 1
	 * @param mathY - The y-coordinate where up is positive and order unity is 1
	 * @return double[2] { screenX, screenY };
	 */
	public double[] transform(double mathX, double mathY) {
		return new double[] {
				margin+size/2. + size/viewW*( (mathX-viewX)*Math.cos(viewTh) + (mathY-viewY)*Math.sin(viewTh)),
				       size/2. - size/viewW*(-(mathX-viewX)*Math.sin(viewTh) + (mathY-viewY)*Math.cos(viewTh)) };
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
	private static Color randomColor(String seed) {
		char s0 = seed.charAt(seed.length()-2), s1 = seed.charAt(seed.length()-1);
		double r = Math.pow(2, 16)*s0 + s1;
		double hue = 9000.*Math.tan(Math.pow(r*2.7, 7.2))%1;
		double sat = Math.sin(r);
		double brt = Math.sin(r*2.7);
		return Color.hsb(180+hue*180, sat/20+.5, brt/20.+.75);
	}
}
