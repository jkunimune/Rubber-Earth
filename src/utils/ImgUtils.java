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
package utils;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import javax.imageio.ImageIO;

import org.apache.commons.imaging.ImageFormats;
import org.apache.commons.imaging.ImageReadException;
import org.apache.commons.imaging.ImageWriteException;
import org.apache.commons.imaging.Imaging;
import org.jcodec.api.awt.AWTSequenceEncoder;
import org.jcodec.common.io.NIOUtils;
import org.jcodec.common.io.SeekableByteChannel;
import org.jcodec.common.model.Rational;

/**
 * A set of functions for dealing with TIFF input and output
 * 
 * @author Justin Kunimune
 */
public class ImgUtils {
	/**
	 * Load the grayscale TIFF file into an array.
	 * @param filename - The name of the .tif file, which is in ../data/
	 * @param resolution - The resolution of the mesh for which this shall be used
	 * @param logBase - The base of the logarithm used to store this data, or 0 if it is linear
	 * @param minVal - Values will be rescaled to go [minVal, 1]
	 * @return the array of values read from the blue channel.
	 * @throws IOException if there is a problem getting the data from disk
	 * @throws ImageReadException if there is a problem with the Tiff image
	 */
	public static double[][] loadTiffData(String filename, int resolution, double logBase,
			double maxVal, double minVal) throws ImageReadException, IOException {
		BufferedImage bimg = Imaging.getBufferedImage(
				new File(String.format("data/%s.tif", filename)));
		if (resolution == 0)
			resolution = bimg.getHeight()/2;
		double[][] data = new double[2*resolution][4*resolution];
		
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
						if (logBase > 0 && val > 0) 	val = Math.pow(logBase, val)/logBase;
						data[i][j] += (logBase != 0) ? Math.pow(logBase, val)/logBase : val;
					}
				}
				data[i][j] /= n*m;
			}
		}
		
		data = normalised(data);
		for (int i = 0; i < data.length; i ++) { // do a linear interpolation type thing to account for the minval
			for (int j = 0; j < data[i].length; j ++) {
				data[i][j] = minVal + (maxVal-minVal)*data[i][j];
			}
		}
		return data;
	}
	
	
	/**
	 * Save the double array as a grayscale Tiff.
	 * @param data - The array of doubles to be saved as grayscale
	 * @param filename - The name of the .tif file, which will be in ../data/
	 * @param logBase - The base of the logarithm used to store this data, or 0 if it is linear
	 * @return the array of values read from the blue channel.
	 * @throws IOException if there is a problem getting the data to disk
	 * @throws ImageWriteException if there is a problem saving the Tiff image
	 */
	public static void saveTiffData(double[][] data, String filename, double logBase) throws ImageWriteException, IOException {
		BufferedImage bimg = new BufferedImage(data[0].length, data.length, BufferedImage.TYPE_INT_ARGB);
		for (int i = 0; i < data.length; i ++) {
			for (int j = 0; j < data[i].length; j ++) {
				double dVal = (logBase > 0 && data[i][j] > 0) ?
						Math.log(data[i][j]*logBase)/Math.log(logBase) : data[i][j];
				int val = 0x0000FF&(int)Math.round(Math.min(Math.max(0, dVal)*255, 255));
				bimg.setRGB(j, i, (val<<16)|(val<<8)|(val<<0));
			}
		}
		Imaging.writeImage(bimg, new File(String.format("data/%s.tif", filename)),
				ImageFormats.TIFF, new HashMap<String, Object>());
	}
	
	
	/**
	 * Compile all frames saved this run into an MP4 and save that to disk
	 * @throws IOException if there is a problem reading or writing the disk
	 */
	public static void compileFrames(String dir, String filename, int numFrames) throws IOException {
		System.out.println("Compiling frames...");
		SeekableByteChannel out = null;
		try {
			out = NIOUtils.writableFileChannel(String.format("%s/%s.mp4", dir, filename));
			AWTSequenceEncoder encoder = new AWTSequenceEncoder(out, Rational.R(25, 1));
			for (int i = 0; i < numFrames; i ++) {
				System.out.println(i);
				BufferedImage image = ImageIO.read(
						new File(String.format("%s/frame%04d.png", dir, i)));
				encoder.encodeImage(image);
			}
			encoder.finish();
			System.out.printf("Successfully compiled frames into %s/%s.mp4\n", dir, filename);
		} finally {
			NIOUtils.closeQuietly(out);
		}
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
	
	
	/**
	 * Scale an array to a new size, maintaining its general scales and properties
	 * @param raw - The image to be scaled
	 * @param width - The desired new width; must be less than or equal to raw[0].length
	 * @param height - The desired new height; must be less than or equal to raw.length
	 * @return the scaled array
	 */
	public static double[][] resize(double[][] raw, int width, int height) {
		double[][] out = new double[height][width];
		for (int is = 0; is < height; is ++) {
			for (int js = 0; js < width; js ++) {
				int irmin = is*raw.length/height;
				int irmax = (is+1)*raw.length/height;
				int jrmin = js*raw[0].length/width;
				int jrmax = (js+1)*raw[0].length/width;
				for (int ir = irmin; ir < irmax; ir ++)
					for (int jr = jrmin; jr < jrmax; jr ++)
						out[is][js] += raw[ir][jr]/((irmax-irmin)*(jrmax-jrmin));
			}
		}
		return out;
	}
	
	
	/**
	 * Elementwise maximum.
	 * @param arrays - the arrays of the same size over which to find the max
	 * @return an array the same size as the inputs containing in each element the maximum of
	 * the corresponding elements of the inputs.
	 */
	public static double[][] max(double[][]... arrays) {
		for (int k = 1; k < arrays.length; k ++)
			if (arrays[k].length != arrays[k-1].length || arrays[k][0].length != arrays[k-1][0].length)
				throw new IllegalArgumentException("Array dimensions do not match!");
		double[][] out = new double[arrays[0].length][arrays[0][0].length];
		for (int i = 0; i < arrays[0].length; i ++)
			for (int j = 0; j < arrays[0][0].length; j ++)
				for (int k = 0; k < arrays.length; k ++)
					if (arrays[k][i][j] > out[i][j])
						out[i][j] = arrays[k][i][j];
		return out;
	}
	
	
	/**
	 * Matrix maximum.
	 * @param array - array of which to find the max.
	 * @return the maximum element of the array.
	 */
	public static double max(double[][] array) {
		double max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < array.length; i ++)
			for (int j = 0; j < array.length; j ++)
				if (array[i][j] > max)
					max = array[i][j];
		return max;
	}
	
	
	/**
	 * Apply a gaussian blur to the double[][], taking the curvature of the sphere into account.
	 * @param raw - The focused array of values to be blurred.
	 * @param sigma - The wavelength of the blur, in radians of course
	 * @return the blurred array
	 */
	public double[][] gaussianBlur(double[][] raw, double sigma) {
		return gaussianBlur(raw, sigma, 1.);
	}
	
	
	/**
	 * Apply a gaussian blur to the double[][], taking the curvature of the sphere into account.
	 * @param raw - The focused array of values to be blurred.
	 * @param sigma - The wavelength of the blur, in radians of course
	 * @param kappa - A multiplier to be applied to the blurred matrix. Note that values will be
	 * 		capped at 1.0.
	 * @return the blurred array
	 */
	public static double[][] gaussianBlur(double[][] raw, double sigma, double kappa) {
		int idxSig = (int)Math.ceil(sigma*raw.length/Math.PI); // sigma in index units
		double sig2 = sigma*sigma; // the variance is more useful than the std
		double scaleFac = 1 - Math.exp(-2/sig2);
		double dp = Math.PI/raw.length; // for the area term
		double dl = 2*Math.PI/raw[0].length;
		double[][] blurred = new double[raw.length][raw[0].length];
		
		for (int ib = 0; ib < blurred.length; ib ++) { // index height for the blurred one
			double pb = Math.PI/2*(1 - 2*(ib + 0.5)/blurred.length); // latitude for the blurred one
			for (int jb = 0; jb < blurred[ib].length; jb ++) { // index distance for the blurred one
				double lb = Math.PI*(2*(jb + 0.5)/blurred[ib].length - 1); // longitude for the blurred one
				System.out.printf("%03d / %03d,  %03d / %03d\n", ib, blurred.length, jb, blurred[ib].length);
				
				int irMin = Math.max(0, ib-4*idxSig); // set bounds, since we don't need to do the whole raw image
				int irMax = Math.min(ib+4*idxSig, raw.length);
				for (int ir = irMin; ir < irMax; ir ++) { // index height for the raw one
					double pr = Math.PI/2*(1 - 2*(ir + 0.5)/raw.length); // latitude for the raw one
					for (int jr = 0; jr < raw[ir].length; jr ++) { // index distance for the raw one
						double lr = Math.PI*(2*(jr + 0.5)/raw[ir].length - 1); // longitude for the raw one
						
						double mdv = Math.sin(pb)*Math.sin(pr) + // mean dot v
								Math.cos(pb)*Math.cos(pr)*Math.cos(lb-lr); // or the cosine of the distance
						double Nbr = Math.exp((mdv-1)/sig2)/(2*Math.PI*sig2)/scaleFac; // 2D Gaussian distribution
						double dA = Math.cos(pr)*dp*dl; // area of <ir,jr>
						blurred[ib][jb] += kappa*Nbr*raw[ir][jr]*dA; // Riemann sum
					}
				}
				blurred[ib][jb] = Math.min(blurred[ib][jb], 1);
			}
		}
		return blurred;
	}
	
	
	/**
	 * I couldn't find a name for this type of blur, so I gave it one. It's a Gaussian blur over a set of
	 * wavelengths from zero to a fixed value, where the result is the maximum of all the different
	 * wavelengths. That makes sure the coasts don't get washed out while also preventing the huge spikes
	 * I would get from small islands if I used a nearest distance-type-thing. There's probably a more
	 * efficient way to do this.
	 * @param raw - The focused array of values to be blurred.
	 * @param maxSig - The maximum wavelength of the blur, in radians of course
	 * @param kappa - A multiplier to be applied to the blurred matrix. Note that values will be
	 * 		capped at 1.0.
	 * @return the blurred array
	 */
	public static double[][] kunimuneanBlur(double[][] raw, double maxSig, double kappa) {
		int idxSig = (int)Math.ceil(maxSig*raw.length/Math.PI); // sigma in index units
		double dp = Math.PI/raw.length; // for the area term
		double dl = 2*Math.PI/raw[0].length;
		double[][] blurred = new double[raw.length][raw[0].length];
		
		for (int ib = 0; ib < blurred.length; ib ++) { // height index for the blurred one
			double pb = Math.PI/2*(1 - 2*(ib + 0.5)/blurred.length); // latitude for the blurred one
			for (int jb = 0; jb < blurred[ib].length; jb ++) { // longitude index for the blurred one
				double lb = Math.PI*(2*(jb + 0.5)/blurred[ib].length - 1); // longitude for the blurred one
				System.out.printf("%03d / %03d,  %03d / %03d\n", ib, blurred.length, jb, blurred[ib].length);
				
				blurred[ib][jb] = Math.min(1, kappa*raw[ib][jb]); // in Kunimunean blur, we only care about the max value over several wavelengths, so we'll update this as we find higher ones
//				for (double sig = maxSig; sig >= Math.PI/raw.length; sig -= Math.PI/raw.length) { // otherwise, we have to nest a loop
				for (double sig = Math.PI/raw.length; sig <= maxSig; sig += Math.PI/raw.length) {
					double sig2 = sig*sig; // the variance
					double scaleFac = 1 - Math.exp(-2/sig2); // over many variances
					double valAtThisSig = 0;
					
					int irMin = Math.max(0, ib-4*idxSig); // set bounds, since we don't need to do the whole raw image
					int irMax = Math.min(ib+4*idxSig, raw.length);
					for (int ir = irMin; ir < irMax; ir ++) { // index height for the raw one
						double pr = Math.PI/2*(1 - 2*(ir + 0.5)/raw.length); // latitude for the raw one
						for (int jr = 0; jr < raw[ir].length; jr ++) { // index distance for the raw one
							double lr = Math.PI*(2*(jr + 0.5)/raw[ir].length - 1); // longitude for the raw one
							
							double mdv = Math.sin(pb)*Math.sin(pr) + // mean dot v
									Math.cos(pb)*Math.cos(pr)*Math.cos(lb-lr); // or the cosine of the distance
							double Nbr = Math.exp((mdv-1)/sig2)/(2*Math.PI*sig2)/scaleFac; // 2D Gaussian distribution
							double dA = Math.cos(pr)*dp*dl; // area of <ir,jr>
							valAtThisSig += Nbr*raw[ir][jr]*dA; // Riemann sum
						}
					}
					
					if (kappa*valAtThisSig > blurred[ib][jb])
						blurred[ib][jb] = Math.min(1, kappa*valAtThisSig);
				}
			}
		}
		return blurred;
	}
	
	
	/**
	 * Scale this array by the given factor
	 * @param matrix - The double array to normalise
	 * @param scaled - The double value by which to scale it
	 * @return a copy whose maximum value is 1.0
	 */
	public static double[][] scaled(double[][] matrix, double factor) {
		double[][] scaled = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; i ++)
			for (int j = 0; j < matrix[i].length; j ++)
				scaled[i][j] = matrix[i][j]*factor;
		return scaled;
	}
	
	
	/**
	 * Scale this array such that its maximum value is unity
	 * @param unnormalised - The double array to normalise
	 * @return a copy whose maximum value is 1.0
	 */
	public static double[][] normalised(double[][] unnormalised) {
		return scaled(unnormalised, 1/max(unnormalised));
	}
	
	
	/**
	 * Scale this array such that its mean is one
	 * @param unnormalised - The double array to normalise
	 * @return a copy whose maximum value is 1.0
	 */
	public static double[][] standardised(double[][] unnormalised) {
		double mean = 0;
		for (int i = 0; i < unnormalised.length; i ++)
			for (int j = 0; j < unnormalised[i].length; j ++)
				mean += unnormalised[i][j]/(unnormalised.length*unnormalised[i].length);
		double[][] normalised = new double[unnormalised.length][unnormalised[0].length];
		for (int i = 0; i < unnormalised.length; i ++)
			for (int j = 0; j < unnormalised[i].length; j ++)
				normalised[i][j] = unnormalised[i][j]/mean;
		return normalised;
	}
	
	
	public static final void main(String[] args) throws IOException, ImageReadException, ImageWriteException {
		String filename = "SRTM_RAMP2_TOPO_2000-02-11_gs_3600x1800";
		System.out.println("loading...");
		double[][] raw = loadTiffData(filename, 0, 0, 1, 0);
		System.out.println("resizing...");
		double[][] small = resize(raw, 360, 180);
		System.out.println("blurring...");
//		double[][] blurred = max(small, gaussianBlur(small, .125, 2));
//		double[][] blurred = gaussianBlur(small, .09375, 1);
		double[][] blurred = kunimuneanBlur(small, .125, 2);
		System.out.println("normalising...");
		double[][] normed = normalised(blurred);
		System.out.println("saving...");
		saveTiffData(normed, filename+"_blur", 0);
		System.out.println("done!");
	}
}
