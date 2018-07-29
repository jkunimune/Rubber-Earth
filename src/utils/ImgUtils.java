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
package utils;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import org.apache.commons.imaging.ImageReadException;
import org.apache.commons.imaging.formats.tiff.TiffImageParser;

/**
 * A set of functions for dealing with TIFF input and output
 * 
 * @author Justin Kunimune
 */
public class ImgUtils {
	/**
	 * Load the grayscale TIFF file into an array.
	 * @param filename - The name of the TIFF file, which is in ../data/
	 * @param resolution - The resolution of the mesh for which this shall be used
	 * @param logBase - The base of the logarithm used to store this data, or 0 if it is linear
	 * @param minVal - Values will be rescaled to go [minVal, 1]
	 * @return the array of values read from the blue channel.
	 */
	public static double[][] loadTiffData(String filename, int resolution, double logBase,
			double minVal) {
		BufferedImage bimg = null;
		if (filename != null) {
			try {
				bimg = new TiffImageParser().getBufferedImage(
						new File("data/"+filename+".tif"), null);
			} catch (ImageReadException e) {
				System.err.println("Warning: could not load data/"+filename+".tif; "+e.getMessage());
			} catch (IOException e) {
				System.err.println("Warning: could not load data/"+filename+".tif; "+e.getMessage());
			}
		}
		double[][] data = new double[2*resolution][4*resolution];
		if (bimg == null) {
			for (int i = 0; i < data.length; i ++)
				for (int j = 0; j < data[i].length; j ++)
					data[i][j] = 1; // default to all ones
		}
		else {
			for (int i = 0; i < data.length; i ++) {
				int ii = bimg.getHeight()*i/data.length;
				int n = bimg.getHeight()*(i+1)/data.length - ii;
				for (int j = 0; j < data[i].length; j ++) {
					int ji = bimg.getWidth()*j/data[i].length;
					int m = bimg.getWidth()*(j+1)/data[i].length - ji;
					for (int di = 0; di < n; di ++) {
						for (int dj = 0; dj < m; dj ++) {
							int argb = bimg.getRGB(ji+dj, ii+di);
							double val = (argb&0x0000FF)/255.;
							if (logBase > 0) 	val = Math.pow(logBase, val)/logBase;
							data[i][j] += (logBase != 0) ? Math.pow(logBase, val)/logBase : val;
						}
					}
					data[i][j] = minVal + (1-minVal)*data[i][j]/(n*m);
				}
			}
		}
		return data;
	}
	
	
	/**
	 * Create a uniform array of ones
	 * @param size - Half the height and a quarter the width
	 * @return the ones
	 */
	public static double[][] uniform(int size) {
		double[][] output = new double[2*size][4*size];
		for (int i = 0; i < output.length; i ++)
			for (int j = 0; j < output[i].length; j ++)
				output[i][j] = 1;
		return output;
	}
}
