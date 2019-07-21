/******************************************************************************
 *  Compilation:  javac TestInseparable.java
 *  Execution:    java TestInseparable
 *  
 *  A class to run the main test
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/
import java.util.*;
import java.io.*;

public class TestInseparable {
	public static void main(String[] args) throws FileNotFoundException {
		TestInseparable test = new TestInseparable();

		//ReadFile get_data_1 = new ReadFileLineByLineUsingScannerAbalone();
		//test.main2(get_data_1, 1000, 100, "Abalone");
		
		ReadFile get_data_2 = new ReadFileBank32nh();
		test.main2(get_data_2, 3000, 100, "Bank");
		
		//ReadFile get_data_3 = new ReadFileBostonHousing();
		//test.main2(get_data_3, 300, 100, "BostonHousing");
		
		//ReadFile get_data_4 = new ReadFileCaliforniaHousing();
		//test.main2(get_data_4, 5000, 30, "CaliforniaHousing");

		//ReadFile get_data_5 = new ReadFileCensus();
		//test.main2(get_data_5, 6000, 100, "Census");

		//ReadFile get_data_6 = new ReadFileComputer();
		//test.main2(get_data_6, 4000, 100, "Computer");

		//ReadFile get_data_7 = new ReadFileMachineCPU();
		//test.main2(get_data_7, 150, 100, "MachineCPU");

		//ReadFile get_data_8 = new ReadFilePyrimidines();
		//test.main2(get_data_8, 50, 100, "Pyrimidines");
	}

	public static void run_a_test(TestSVM_Continuous2 test_one, ReadFile get_data) {
		//read data
		long start_time = System.currentTimeMillis();
		test_one.read_data_general(test_one.max_input_size, get_data);
		long time_readdata = System.currentTimeMillis() - start_time; 
		

		//
		//initilization of upper bound and lowerbound 
		test_one.initilization_bounds();

		//initilization of variables and precompute the nested sum 
		test_one.initilization_variables();

		//initialization of the discriminant function
		test_one.initilization_discriminant();

		//initialization of stopping criteria
		double initial_tau = test_one.stop_criterion();
		System.out.println("initial_obj: " + test_one.compute_obj_val());
		System.out.println("initial_tau: " + initial_tau);

		test_one.projected_gradient_ascent();

		test_one.display_time();
		
		long time = System.currentTimeMillis() - start_time;

		test_one.time_total = time; 

		System.out.println(
					"Total time: "
					+ String.format("%10s", ((double) time) / 1000));
		
	}
	public static void main2(ReadFile get_data, int max_input_size, double C, String name) throws FileNotFoundException {
		long currentTime = System.currentTimeMillis();
		PrintStream o = new PrintStream(new File("Non-separable Dataset name " + name + "  " + "DCA--MDA" + Long.toString(currentTime) + ".txt"));
		System.setOut(o);

		test_MDA_DCA(get_data, max_input_size, C);
	}

	public static void test_MDA_DCA(ReadFile get_data, int max_input_size, double C) {

		int[] size_set = new int[]{2,4,6,8,10};

		int[] total_iteration_DCA = new int[5];
		int[] total_iteration_MDA = new int[5];
		long[] total_time_DCA = new long[5];
		long[] total_time_MDA = new long[5];
		long[] total_time_RAP_NC_DCA = new long[5];
		long[] total_time_RAP_NC_MDA = new long[5];
		int[] total_number_RAP_MDA = new int[5];
		int[] total_number_RAP_DCA = new int[5];

		//int max_input_size = 5000;

		for (int i = 0; i < 5; i++) {
			TestSVM_Continuous2 test_dca = new TestSVM_Continuous2();
			test_dca.isMDA = false;
			test_dca.display_obj_flag = false;
			test_dca.max_input_size = max_input_size;
			test_dca.working_set_size = size_set[i];
			test_dca.constant_C = C;

			run_a_test(test_dca, get_data);

			total_iteration_DCA[i] = test_dca.curt_iter + 1;

			total_time_DCA[i] = test_dca.time_total;

			total_time_RAP_NC_DCA[i] = test_dca.time_solving_RAPNC;

			total_number_RAP_DCA[i] = test_dca.number_of_RAP;
			
			TestSVM_Continuous2 test_mda = new TestSVM_Continuous2();
			test_mda.isMDA = true;
			test_mda.display_obj_flag = false;
			test_mda.max_input_size =max_input_size;
			test_mda.working_set_size = size_set[i];
			test_mda.constant_C = C;

			run_a_test(test_mda, get_data);

			total_iteration_MDA[i] = test_mda.curt_iter + 1;

			total_time_MDA[i] = test_mda.time_total;

			total_time_RAP_NC_MDA[i] = test_mda.time_solving_RAPNC;

			total_number_RAP_MDA[i] = test_mda.number_of_RAP;

		} 
	
		for (int i = 0; i < 5; i++) {
			System.out.println(
				"&  &  &" 
				+ size_set[i] 
				+ " & " 
				+ total_iteration_MDA[i] 
				+ " & "
				+ String.format("%10s", ((double) total_time_MDA[i]) / 1000)
				+ " & "
				+ String.format("%10s", ((double) total_time_RAP_NC_MDA[i]) / 1000)
				+ " & "
				+ total_number_RAP_MDA[i]
				+ " & "
				+ total_iteration_DCA[i]
				+ " & "
				+ String.format("%10s", ((double) total_time_DCA[i]) / 1000)
				+ " & "
				+ String.format("%10s", ((double) total_time_RAP_NC_DCA[i]) / 1000)
				+ " & "
				+ total_number_RAP_DCA[i]
				+ "\\\\");
		}		
	}
}