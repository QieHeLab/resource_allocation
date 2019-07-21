/******************************************************************************
 *  Compilation:  javac ReadFile.java
 *  Execution:    java ReadFile
 *  
 *  A class to preprocess the data 
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.*;

public abstract class ReadFile {

	public static class Data_format {
		int class_index = 0;
		int sample_index = 0;
		int rap_nc_index = 0;
		double[] features = null;
		int[] stats = null;

		public Data_format(int i, int j, double[] f, int rap) {
			this.class_index = i;
			this.sample_index = j;
			this.features = f;
			this.rap_nc_index = rap;
		}

		public void print_as_string() {
			System.out.print(class_index);
			System.out.print(" , ");
			System.out.print(sample_index);
			System.out.print(" , ");
			System.out.print(rap_nc_index);
			System.out.print(" , ");
			System.out.println(Arrays.toString(features));
		}
	}

	public static void main(String[] args) {
		return;
	}

	public abstract List<Data_format> readdata (int size);

	public static void standardlize(double[][] data_output) {
		//standardlize the data
		for (int i = 0; i < data_output[0].length - 1; i++) {
			double mean = 0;
			double variance = 0;
			for (int j = 0; j < data_output.length; j++) {
				mean += data_output[j][i];
				variance += data_output[j][i] * data_output[j][i];
			}
			mean /= data_output.length;

			variance /= data_output.length; 
			variance -= mean * mean;
			/*
			for (int j = 0; j < data_output.length; j++) {
				variance += (data_output[j][i] - mean) * (data_output[j][i] - mean);
			}
			variance /= data_output.length;
			*/

			double sigma =  Math.sqrt(variance);
			for (int j = 0; j < data_output.length; j++) {
				data_output[j][i] = (data_output[j][i] - mean);
				if (sigma > 1e-5) {
					data_output[j][i] /= sigma; 
				}
			}
		}
	}

}
