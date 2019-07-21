/******************************************************************************
 *  Compilation:  javac ReadFileBostonHousing.java
 *  Execution:    java ReadFileBostonHousing
 *  
 *  A class to preprocess the data 
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.*;

public class ReadFileBostonHousing extends ReadFile {
	

	public static void main(String[] args) {
		int size = 10;
		System.out.println("10");

		ReadFileBostonHousing test = new ReadFileBostonHousing();
		test.readdata(300);
	}

	public List<ReadFile.Data_format> readdata (int size) {
		
		double[][] data_output = new double[size][14];

		try {
			
			Scanner scanner = new Scanner(new File("./data/housing.data"));
			int i = 0;
			while (scanner.hasNextLine() && i < size) {
				String[] data_sample = scanner.nextLine().split(",");
				//System.out.println(scanner.nextLine());

				for (int j = 0; j < 14; j++) {
					data_output[i][j] = Double.parseDouble(data_sample[j]);
				}

				i++;
			}
			scanner.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		System.out.println(Arrays.toString(data_output[0]));
		System.out.println("Boston housing");

		int[] data_stats = new int[5];
		for (int i = 0; i < size; i++) {
			if (i == data_output.length) {
				break;
			}

			if (data_output[i][13] < 18) {
				data_output[i][13] = 1;
				data_stats[0]++;
			} else if (data_output[i][13] < 22) {
				data_output[i][13] = 2;
				data_stats[1]++;
			} else if (data_output[i][13] < 25) {
				data_output[i][13] = 3;
				data_stats[2]++;
			} else if (data_output[i][13] < 35) {
				data_output[i][13] = 4;
				data_stats[3]++;
			} else {
				data_output[i][13] = 5;
				data_stats[4]++;
			}
		}

		System.out.println(Arrays.toString(data_stats));

		standardlize(data_output);

		//sort the data based on their output
		java.util.Arrays.sort(data_output, (a, b) -> Double.compare(a[13], b[13]));

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

		return data_index;
	}
}