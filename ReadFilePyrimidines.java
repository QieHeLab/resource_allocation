/******************************************************************************
 *  Compilation:  javac ReadFilePyrimidines.java
 *  Execution:    java ReadFilePyrimidines
 *  
 *  A class to preprocess the data 
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.*;

public class ReadFilePyrimidines extends ReadFile {
	

	public static void main(String[] args) {
		int size = 10;
		System.out.println("10");

		ReadFile test = new ReadFilePyrimidines();
		test.readdata(50);
	}

	public List<ReadFile.Data_format> readdata (int size) {
		int dimension = 27;
		double[][] data_output = new double[size][dimension + 1];

		try {
			
			Scanner scanner = new Scanner(new File("./data/pyrim.data"));
			int i = 0;
			while (scanner.hasNextLine() && i < size) {
				String[] data_sample = scanner.nextLine().split(",");
				//System.out.println(scanner.nextLine());

				for (int j = 0; j < dimension + 1; j++) {
					data_output[i][j] = Double.parseDouble(data_sample[j]);
				}

				i++;
			}
			scanner.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		System.out.println(Arrays.toString(data_output[0]));
		System.out.println("Computer activities");

		//sort the data based on their output
		java.util.Arrays.sort(data_output, (a, b) -> Double.compare(a[dimension], b[dimension]));

		int[] data_stats = new int[5];
		double mean = 0;
		for (int i = 0; i < size; i++) {
 			mean += data_output[i][dimension];
		}
		mean /= size;
		System.out.println(mean);
		for (int i = 0; i < size; i++) {
			if (i == data_output.length) {
				break;
			}
			/*
			data_output[i][8] -= mean;
			data_output[i][8] /= 10000;
			if (data_output[i][8] < -9) {
				data_output[i][8] = 1;
				data_stats[0]++;
			} else if (data_output[i][8] < -5) {
				data_output[i][8] = 2;
				data_stats[1]++;
			} else if (data_output[i][8] < 0) {
				data_output[i][8] = 3;
				data_stats[2]++;
			} else if (data_output[i][8] < 6) {
				data_output[i][8] = 4;
				data_stats[3]++;
			} else {
				data_output[i][8] = 5;
				data_stats[4]++;
			}
			*/
			
			if (i < size / 5) {
				data_output[i][dimension] = 1;
				data_stats[0]++;
			} else if (i < size * 2 / 5) {
				data_output[i][dimension] = 2;
				data_stats[1]++;
			} else if (i < size * 3 / 5) {
				data_output[i][dimension] = 3;
				data_stats[2]++;
			} else if (i < size * 4 / 5) {
				data_output[i][dimension] = 4;
				data_stats[3]++;
			} else {
				data_output[i][dimension] = 5;
				data_stats[4]++;
			}
			
		}

		System.out.println(Arrays.toString(data_stats));

		standardlize(data_output);

		//sort the data based on their output
		//java.util.Arrays.sort(data_output, (a, b) -> Double.compare(a[8], b[8]));
		/*
		double[][] test = new double[][]{{1,20},{3,4}};
		java.util.Arrays.sort(test, (a, b) -> Double.compare(a[1], b[1]));
		System.out.println(Arrays.toString(test[0]));
		System.out.println(Arrays.toString(test[1]));
		*/

		//Create index for the data
		List<ReadFile.Data_format> data_index = new ArrayList<>();
		int count = 0;
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
		
		System.out.println(Arrays.toString(data_output[size / 5 - 1]));
		System.out.println(Arrays.toString(data_output[size / 5]));
		System.out.println(Arrays.toString(data_output[size - 1]));
		
		return data_index;
	}
}