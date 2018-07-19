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

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.geometry.jts.Geometries;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.Feature;
import org.opengis.filter.Filter;

import javafx.geometry.VPos;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.LineTo;
import javafx.scene.shape.MoveTo;
import javafx.scene.shape.Path;
import javafx.scene.shape.Shape;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import model.Mesh;

/**
 * A class devoted to making the mesh visible and pretty.
 * 
 * @author Justin Kunimune
 */
public class Renderer {
	
	private final Group entities;
	private final Text readout;
	private final Map<Feature, Shape> shapes;
	
	private final int size;
	private final double scale, offset;
	private final double decayTime; // in milliseconds
	private final boolean saveImages; // save frames to disk?
	
	private Mesh mesh;
	private long lastRender;
	
	
	public Renderer(int size, Mesh mesh, double decayTime, boolean saveImages, String[] shpFiles) {
		this.mesh = mesh;
		this.size = size;
		this.scale = size/(2*Math.PI);
		this.offset = size/2.;
		this.decayTime = decayTime;
		this.saveImages = saveImages;
		
		this.readout = new Text(10, 0, "");
		this.readout.setTextOrigin(VPos.TOP);
		this.readout.setFont(Font.font(20));
		this.entities = new Group();
		List<Feature> geoData = tryLoadShapefile(shpFiles);
		this.shapes = createShapes(geoData);
		for (Feature f: geoData)
			this.entities.getChildren().add(shapes.get(f));
		this.entities.getChildren().add(readout);
		
		this.lastRender = System.currentTimeMillis();
		this.render();
	}
	
	
	/**
	 * Attempt to load the shapefile located at the filename and return the DataStore therein.
	 * @param pathname - The full pathname of the shapefile. Pretty self-explanatory.
	 * @return the DataStore of the shapefile or an empty DataStore if you couldn't load it.
	 */
	private List<Feature> tryLoadShapefile(String[] filenames) {
		List<Feature> output = new LinkedList<Feature>();
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
							output.add(iterator.next());
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
	
	
	private Map<Feature, Shape> createShapes(Collection<Feature> features) {
		Map<Feature, Shape> shapes = new HashMap<Feature, Shape>();
		for (Feature f: features) {
			Geometry geom = (Geometry)f.getDefaultGeometryProperty().getValue();
			Geometries geomType = Geometries.get(geom);
			
			switch (geomType) {
			case LINESTRING:
			case MULTILINESTRING:
			case POLYGON:
			case MULTIPOLYGON:
				Path display = new Path();
				for (int i = 0; i < geom.getNumGeometries(); i ++) { // why is this iteration not built in?
					for (int j = 0; j < geom.getGeometryN(i).getNumPoints(); j ++) {
						Coordinate c = geom.getGeometryN(i).getCoordinates()[j];
						if (j == 0)
							display.getElements().add(new MoveTo(c.x, c.y));
						else
							display.getElements().add(new LineTo(c.x, c.y));
					}
				}
				if (geomType == Geometries.LINESTRING || geomType == Geometries.MULTILINESTRING) {
					display.setStroke(Color.gray(.2));
					display.setFill(null);
				}
				else {
					display.setStroke(null);
					display.setFill(Color.hsb(Math.random()*360, .45+Math.random()*.1, .7+Math.random()*.1));
				}
				shapes.put(f, display);
				break;
				
			case POINT:
			case MULTIPOINT:
				for (Coordinate c: geom.getCoordinates())
					shapes.put(f, new Circle(c.x, c.y, 2));
				break;
				
			default:
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
		double c1 = Math.exp((lastRender-now)/decayTime);
		double c2 = 1. - c1;
		
//		for (Cell c: mesh.getCellsUnmodifiable()) { // for every vertex
//			if (!lines.containsKey(c))
//				lines.put(c, new Line[4]); // make sure it's in the map
//			
//			for (int i = 0; i < 4; i ++) { // for each edge
//				if (lines.get(c)[i] == null) {
//					lines.get(c)[i] = new Line(offset,offset,offset,offset); // make sure the line exists for it
//					entities.getChildren().add(lines.get(c)[i]);
//				}
//				Vertex v0 = c.getCorner(i); // get the endpoints
//				Vertex v1 = c.getCorner((i+1)%4);
//				Line l = lines.get(c)[i];
//				
//				l.setStartX(c1*l.getStartX() + c2*(offset+scale*v0.getX())); // update the line to its current position
//				l.setStartY(c1*l.getStartY() + c2*(offset-scale*v0.getY())); // but slowly
//				l.setEndX(c1*l.getEndX() + c2*(offset+scale*v1.getX()));
//				l.setEndY(c1*l.getEndY() + c2*(offset-scale*v1.getY()));
//			}
//		}
		
		for (Feature f: this.shapes.keySet()) {
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
