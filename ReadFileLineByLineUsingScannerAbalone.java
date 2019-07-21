/******************************************************************************
 *  Compilation:  javac ReadFileLineByLineUsingScannerAbalone.java
 *  Execution:    java ReadFileLineByLineUsingScannerAbalone
 *  
 *  A class to preprocess the data 
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.*;

public class ReadFileLineByLineUsingScannerAbalone extends ReadFile {
	/*
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
	*/

	public static void main(String[] args) {
		double[][] data_output = new double[4177][11];

		try {
			
			Scanner scanner = new Scanner(new File("./data/abalone.data"));
			int i = 0;
			while (scanner.hasNextLine() && i < 4177) {
				String[] data_sample = scanner.nextLine().split(",");
				//System.out.println(scanner.nextLine());
				if (data_sample[0].equals("M")) {
					data_output[i][0] = 1;
				} else if (data_sample[0].equals("F")) {
					data_output[i][1] = 1;
				} else {
					data_output[i][2] = 1;
				}

				for (int j = 3; j < 11; j++) {
					data_output[i][j] = Double.parseDouble(data_sample[j - 2]);
				}

				//System.out.println(data_output[i][10]);
				//System.out.println(data_sample[0]);
				//System.out.println(data_sample[0].length());
				i++;
			}
			scanner.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}


		//sort the data based on their output
		java.util.Arrays.sort(data_output, (a, b) -> Double.compare(a[10], b[10]));


		//count the sample in each class
		long[] data_stats = new long[29];
		int count = 0;
		int class_index = 0;
		for (int i = 0; i < 4177; i++) {
			if (Math.abs(data_output[i][10] - class_index - 1) <= 0.01) {
				count++;
			} else {
				data_stats[class_index] = count;
				count = 1;
				while (Math.abs(data_output[i][10] - class_index - 1) > 0.01) {
					class_index++;
				}
			}
		}
		data_stats[class_index] = count;

		
		//Test code to verify the results
		/*
		System.out.println(java.util.Arrays.toString(data_stats));
		System.out.println(data_output[4174][10]);
		System.out.println(data_output[4175][10]);
		System.out.println(data_output[4176][10]);
		*/

		//Create index for the data
		List<Data_format> data_index = new ArrayList<>();
		count = 0;
		for (int i = 0; i < data_stats.length; i++) {
			for (int j = 0; j < data_stats[i]; j++) {
				data_index.add(new Data_format(i + 1, j + 1, data_output[count], count));
				count++;
				data_index.get(data_index.size()-1).print_as_string();
			}
		}
		//System.out.println(count);


	}

	public List<Data_format> readdata (int size) {
		//double[][] data_output = new double[4177][11];
		//int size = 1000;
		double[][] data_output = new double[size][11];

		try {
			//Scanner scanner = new Scanner(new File("/Users/Genius/Dropbox/java2/RAP-NC/data/abalone.data"));
			Scanner scanner = new Scanner(new File("/Users/Zeyang/Dropbox/java2/RAP-NC/data/abalone.data"));
			//Scanner scanner = new Scanner(new File("./data/abalone.data"));
			int i = 0;
			while (scanner.hasNextLine() && i < size) {
				String[] data_sample = scanner.nextLine().split(",");
				//System.out.println(scanner.nextLine());

				if (data_sample[0].equals("M")) {
					data_output[i][0] = 1;
					//System.out.println(data_sample[0]);
					//System.out.println(Arrays.toString(data_output[i]));
				} else if (data_sample[0].equals("F")) {
					data_output[i][1] = 1;
				} else {
					data_output[i][2] = 1;
				}

				for (int j = 3; j < 11; j++) {
					data_output[i][j] = Double.parseDouble(data_sample[j - 2]);
				}

				//System.out.println(data_output[i][10]);
				//System.out.println(data_sample[0]);
				//System.out.println(data_sample[0].length());
				i++;
			}
			scanner.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		//standardlize the data
		for (int i = 0; i < data_output[0].length - 1; i++) {
			double mean = 0;
			for (int j = 0; j < data_output.length; j++) {
				mean += data_output[j][i];
			}
			mean /= data_output.length;

			double variance = 0;
			for (int j = 0; j < data_output.length; j++) {
				variance += (data_output[j][i] - mean) * (data_output[j][i] - mean);
			}
			variance /= data_output.length;

			double sigma =  Math.sqrt(variance);
			for (int j = 0; j < data_output.length; j++) {
				data_output[j][i] = (data_output[j][i] - mean);
				if (sigma > 1e-5) {
					data_output[j][i] /= sigma; 
				}
			}
		}

		//sort the data based on their output
		java.util.Arrays.sort(data_output, (a, b) -> Double.compare(a[10], b[10]));

		//discretized the labels into 5 classes
		for (int i = 0; i < data_output.length; i++) {
			if (data_output[i][10] <= 7) {
				data_output[i][10] = 1;
			} else if (data_output[i][10] <= 9) {
				data_output[i][10] = 2;
			} else if (data_output[i][10] <= 11) {
				data_output[i][10] = 3;
			} else if (data_output[i][10] <= 14) {
				data_output[i][10] = 4;
			} else {
				data_output[i][10] = 5;
			}
		}

		//count the sample in each class
		int[] data_stats = new int[5];
		int count = 0;
		int class_index = 0;
		for (int i = 0; i < size; i++) {
			if (Math.abs(data_output[i][10] - class_index - 1) <= 0.01) {
				count++;
			} else {
				data_stats[class_index] = count;
				count = 1;
				while (Math.abs(data_output[i][10] - class_index - 1) > 0.01) {
					class_index++;
				}
			}
		}
		data_stats[class_index] = count;

		
		//Test code to verify the results
		/*
		System.out.println(java.util.Arrays.toString(data_stats));
		System.out.println(data_output[4174][10]);
		System.out.println(data_output[4175][10]);
		System.out.println(data_output[4176][10]);
		*/

		
		//Create index for the data
		List<Data_format> data_index = new ArrayList<>();
		count = 0;
		for (int i = 0; i < data_stats.length; i++) {
			for (int j = 0; j < data_stats[i]; j++) {
				data_index.add(new Data_format(i + 1, j + 1, data_output[count], count));
				count++;
				//data_index.get(data_index.size()-1).print_as_string();
			}
		}
		
		Data_format stats = new Data_format(0, 0, null, 0);
		stats.stats = data_stats;
		data_index.add(stats);

		return data_index;
	}

}
