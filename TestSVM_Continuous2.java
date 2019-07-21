/******************************************************************************
 *  Compilation:  javac -cp ./ TestSVM_Continuous2.java
 *  Execution:    java TestSVM_Continuous2
 *  
 *  A test class to evaluate the numerical performance of DCA, MDA when they are used as a subroutine in the projected gradient descent of COSVM
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/
import java.util.*;
import java.io.*;

public class TestSVM_Continuous2 {
	/**
     * An internal class to define the Objective functions
     * It is a function oracle for quadratic function
     * @param a, b, etc parameter of the function
     */
	private static class QuadraticFunction extends Function {
		double a;
		double b;
		public QuadraticFunction(double a, double b) {
			this.a = a;
			this.b = b;
		}
		public double getValue(double x) {
			return a * x * x + b * x;
		}	
	}

	boolean debugging = false;
	boolean debugging_stop_criteria = false;
	boolean isMDA = true;
	boolean isStochasticBatchGradient = false; //do not converge at some times
	boolean display_obj_flag = true;
	boolean find_most_violated_samples_KKT_not_each_hyperplane = false;
	

	//data related
	List<ReadFile.Data_format> data_raw;
	int[] data_stats;
	int[] data_stats_cumulative;
	int data_size;
	double[][] coefficients;

	//stopping criterion
	double[] discriminant_function;

	//index related
	int[] alpha_index_global; 
	int[] beta_index_global; // alpha_star + beta = constant_C 
	int[] class_global;
	int[] index_in_class;
	int[] variable_original_index;

	//subproblem related
	double[] lbNested;
	double[] ubNested;

	//the solution is organized in the following way: beta_i^1, alpha_i^1, beta_i^2, alpha_i^2, \cdots
	double[] alpha; 

	//accumulated sum of variables on, indexed on the end of beta_i^j's
	double[] sum_beta_alpha; // sum of beta and alpha: used in simplified the subproblem
	double[] sum_alpha; // sum of alpha_star and alpha: used in computing the stopping criterion

	//gradient descent related
	double[] old_alpha; // store the variables to update the discriminant function incrementally
	double[] old_alpha_for_incre; //store the variables to update the discriminant function incrementally for the variables in the working set. 
	double constant_C = 30;
	double kernel_kappa = 0.2;
	int max_input_size = 5000;
	int max_iteration = 500000;
	double step_size = 0.2;
	int batch_size = 1; // deprecated for now
	int n_grad = 20;
	int working_set_size = 2;
	int[] working_set;
	int batch_curt_index = 0;
	double[] gradient_container; //Since we are using the projected gradient ascent, we cannot update the gradient in an online fashion (before the projection) as it may not be feasible.
	int[] most_vio;
	Set<Integer> most_vio_test = null;
	int curt_iter;
	double print_frequency = 0.001;
	double obj_val = 0;
	double accuracy_in_stopping_criteria = 1e-7;
	double accuracy_to_stop = 0.001;

	//use to perform the gradient descent and construct the subproblem
	int[] working_set_stats;
	int[] working_set_stats_cumulative;
	int[] working_set_variable_index;
	int[] working_set_variable_isDummy;
	int[] working_set_variable_isBeta; //different in the gradient (1 - discriminant_f(x) or -1 - discriminant_f(x));

	//projection subproblem
	RAPNC_Continuous projection_subproblem;

	//time related
	long time_gradient_update = 0;
	long time_simplfying = 0;
	long time_solving_RAPNC = 0;
	long time_test_convergence = 0;
	long time_total = 0;
	int number_of_RAP = 0;

	/**
     * read_data_test
     * With this operation we input the data Abalone
     * We compute the commonly used index and statistics. 
     * <p>
     * Time-Complexity: O(size)
     * @param size: the maximum number of sample
     */
	public void read_data_test(int size) {
		//read data
		//ReadFileLineByLineUsingScannerAbalone get_data = new ReadFileLineByLineUsingScannerAbalone();
		//ReadFileBank32nh get_data = new ReadFileBank32nh();
		//ReadFileBostonHousing get_data = new ReadFileBostonHousing();
		ReadFile get_data = new ReadFileCaliforniaHousing();
		//ReadFile get_data = new ReadFileComputer();
		//ReadFile get_data = new ReadFilePyrimidines();
		//ReadFile get_data = new ReadFileCensus();
		
		//
		data_raw = get_data.readdata(size);
		data_stats = data_raw.get(data_raw.size() - 1).stats;
		
		//cumulative 
		data_stats_cumulative = new int[data_stats.length];
		int cumulative = 0;
		for (int i =0; i < data_stats.length; i++) {
			cumulative += data_stats[i];
			data_stats_cumulative[i] = cumulative;
		}

		System.out.println("\n Data Statistics:");
		System.out.println(Arrays.toString(data_stats));
		data_size =  data_raw.size() - 1;

		//initialized the index 
		//precompute the alpha and alpha star index for each sample i
		alpha_index_global = new int[data_size];
		beta_index_global = new int[data_size];
		class_global = new int[data_size];
		index_in_class = new int[data_size];
		variable_original_index = new int[data_size*2];

		for (int i = 0; i < data_size; i++) {
			int class_i = (int) data_raw.get(i).class_index;
			class_global[i] = class_i;
			index_in_class[i] = (int) data_raw.get(i).sample_index; 

			if (class_i > 1) {
				beta_index_global[i] = (int) (data_stats_cumulative[class_i - 2] + i);
			} else {
				beta_index_global[i] = i;
			}
			alpha_index_global[i] = (int) (beta_index_global[i] + data_stats[class_i - 1]);

			variable_original_index[beta_index_global[i]] = i;
			variable_original_index[alpha_index_global[i]] = i;
		}

		//precompute the inner product
		coefficients = new double[data_size][data_size];
		for (int i = 0; i < data_size; i++) {
			ReadFile.Data_format data_i = data_raw.get(i);
			for (int j = 0; j < i; j++) {
				//System.out.println("Check here: computing coefficients: " + j);
				ReadFile.Data_format data_j = data_raw.get(j);
				double sum = 0;
				for (int k = 0; k < data_j.features.length - 1; k++) {
					//System.out.println("Check here: computing coefficients: " + sum);
					sum += - 0.5 * (data_j.features[k] - data_i.features[k]) * (data_j.features[k] - data_i.features[k]);
				}
				//System.out.println("Check here: computing coefficients: " + sum);
				//System.out.println(Arrays.toString(data_i.features));
				//System.out.println(Arrays.toString(data_j.features));
				coefficients[i][j] = Math.exp(kernel_kappa * sum);
			}
		} 

		for (int i = 0; i < data_size; i++) {
			coefficients[i][i] = 1;
			for (int j = i + 1; j < data_size; j++) {
				coefficients[i][j] = coefficients[j][i];
			}
		} 



		/*For debug*/
		if (debugging) {
			System.out.println("\n beta index:");
			System.out.println(Arrays.toString(beta_index_global));
			System.out.println("\n alpha index:");
			System.out.println(Arrays.toString(alpha_index_global));
			System.out.println("\n Class wrt sample index");
			System.out.println(Arrays.toString(class_global));
			System.out.println("\n Index with in a class wrt sample index");
			System.out.println(Arrays.toString(index_in_class));
			System.out.println("\n variable_original_index:");
			System.out.println(Arrays.toString(variable_original_index));

			//System.out.println("\n coefficients");
			for (int i = 0; i < 0; i++) {
				System.out.println(Arrays.toString(coefficients[i]));
			}
		}
		
	}

	/**
     * initilization_bounds
     * With this operation we initialize the nested bounds. The bounds will not change over time. 
     * <p>
     * @param no param
     */
	public void initilization_bounds() {
		ubNested = new double[data_stats.length];
		lbNested = new double[data_stats.length];
		for (int i = 0; i < data_stats_cumulative.length; i++) {
			lbNested[i] = data_stats_cumulative[i] * constant_C;
			ubNested[i] = data_size * constant_C;
		}

		if (debugging) {
			/*For debug*/
			System.out.println("\n Nested bounds on variables:");
			System.out.println(Arrays.toString(lbNested));
			System.out.println(Arrays.toString(ubNested));
		}
	}
	
	/**
     * initilization_bounds
     * With this operation we initialize the decision variables.
     * <p>
     * @param no param
     */
	public void initilization_variables() {
		alpha = new double[2 * data_size];
		sum_beta_alpha = new double[data_stats.length];
		sum_alpha = new double[data_stats.length];
		
		for (int i = 0; i < data_stats[0]; i++) {
			// initialize beta_i^1 = constant_C or (alpha_star_1_i = 0)
			alpha[i] = constant_C;
			sum_beta_alpha[0] += constant_C;
		}

		for (int j = 1; j < data_stats.length; j++) {
			// initialize beta_i^j = constant_C or (alpha_star_1_j = 0)
			int start_index = 2*data_stats_cumulative[j-1];
			int end_index = start_index + data_stats[j];
			
			for (int i = start_index; i < end_index; i++) {
				alpha[i] = constant_C;
				sum_beta_alpha[j] += constant_C;
			}

			sum_beta_alpha[j] += sum_beta_alpha[j-1];
		}

		for (int j = 0; j < data_stats.length; j++) {
			sum_alpha[j] = sum_beta_alpha[j] - data_stats_cumulative[j] * constant_C;
			if (Math.abs(sum_alpha[j]) > 1e-9) {
				//Make sure everything is correct
				System.out.println(sum_alpha[j]);
				System.out.println("The array sum_alpha (mu^j) is not initialized correctly!");
			}
		}

		if (debugging) {
			/*For debug*/
			System.out.println("\n Initial variables:");
			System.out.println(Arrays.toString(alpha));
			System.out.println("\n Initial nested num of variables (the transformed variables):");
			System.out.println(Arrays.toString(sum_beta_alpha));
			System.out.println("\n Initial nested num of variables (the initial variables):");
			System.out.println(Arrays.toString(sum_alpha));
		}

	}

	/**
     * initilization_discriminant()
     * The operation initialize the discriminant function for each sample
     * Then we can update the gradient incrementally
     * <p>
     * @param 
     */
	public void initilization_discriminant() {
		//preprocess the initial discriminant function
		discriminant_function = new double[data_size];
		for (int i = 0; i < data_size; i++) {
			for (int j = 0; j < data_size; j++) {
				discriminant_function[i] += (constant_C - alpha[beta_index_global[j]] - alpha[alpha_index_global[j]]) * coefficients[i][j]; 
			}
		}

		if (debugging) {
			/*For debug*/
			System.out.println("\n discriminant_function:");
			System.out.println(Arrays.toString(discriminant_function));
		}
	}


	/**
     * stop_criterion()
     * With this operation we can compute tau, (the stopping criteria)
     * At the same time, it will collected the most violated sample pair for each hyperplane 
     * <p>
     * @param discriminant
     * @param alpha
     * @param sum_alpha
     * The detail can be found in Support Vector Ordinal Regression (Wei Chu and S. Sathiya Keerthi 2007)
     */
	public double stop_criterion() {
		double tau = -1e10;
		double accuracy = accuracy_in_stopping_criteria;

		double[] b_up = new double[data_stats.length+1];
		double[] b_low = new double[data_stats.length+1];
		int[] b_low_index = new int[data_stats.length+1];
		int[] b_up_index = new int[data_stats.length+1];

		for (int i =0; i < b_low.length; i++) {
			b_low[i] = -Double.MAX_VALUE/2;
			b_up[i] = Double.MAX_VALUE/2;
			b_low_index[i] = -1;
			b_up_index[i] = -1;
		}

		for (int i = 0; i < data_stats_cumulative[0]; i++) {
			if (b_up[0] > discriminant_function[i] - 1) {
				b_up[0] = discriminant_function[i] - 1;
				b_up_index[0] = i;
			}
		}

		for (int i = 0; i < data_size; i++) {
			int class_i = class_global[i];
			double val = alpha[beta_index_global[i]];

			if (Math.abs(val - constant_C) < accuracy) {				
				// if alpha_star_i^j is close to zero	
				if (b_up[class_i - 1] > discriminant_function[i] - 1) {
					b_up[class_i - 1] = discriminant_function[i] - 1;
					b_up_index[class_i - 1] = i;
				}
			} else if (Math.abs(val) < accuracy) {
				// if alpha_star_i^j is close to constant_C	
				if (b_low[class_i - 1] < discriminant_function[i] - 1) {
					b_low[class_i - 1] = discriminant_function[i] - 1;
					b_low_index[class_i - 1] = i;
				}
			} else {
				// if alpha_star_i^j is not null	
				if (b_up[class_i - 1] > discriminant_function[i] - 1) {
					b_up[class_i - 1] = discriminant_function[i] - 1;
					b_up_index[class_i - 1] = i;
				}
				if (b_low[class_i - 1] < discriminant_function[i] - 1) {
					b_low[class_i - 1] = discriminant_function[i] - 1;
					b_low_index[class_i - 1] = i;
				}
			}

			val = alpha[alpha_index_global[i]];

			if (Math.abs(val) < accuracy) {				

				// if alpha_i^j is close to zero
				if (b_low[class_i] < discriminant_function[i] + 1) {
					b_low[class_i] = discriminant_function[i] + 1;
					b_low_index[class_i] = i;
				}				
			} else if (Math.abs(val - constant_C) < accuracy) {
				// if alpha_i^j is close to constant_C	
				if (b_up[class_i] > discriminant_function[i] + 1) {
					b_up[class_i] = discriminant_function[i] + 1;
					b_up_index[class_i] = i;
				}
			} else {
				// if alpha_i^j is not null	
				if (b_low[class_i] < discriminant_function[i] + 1) {
					b_low[class_i] = discriminant_function[i] + 1;
					b_low_index[class_i] = i;
				}
				if (b_up[class_i] > discriminant_function[i] + 1) {
					b_up[class_i] = discriminant_function[i] + 1;
					b_up_index[class_i] = i;
				}
			}
		} 

		if (debugging_stop_criteria) {
			/*For debug*/
			System.out.println("\n b_up:");
			System.out.println(Arrays.toString(b_up));
			System.out.println("\n b_up_index:");
			System.out.println(Arrays.toString(b_up_index));
			System.out.println("\n b_low:");
			System.out.println(Arrays.toString(b_low));
			System.out.println("\n b_low_index:");
			System.out.println(Arrays.toString(b_low_index));
		}

		double[] B_tilta_low = new double[data_stats.length + 1];
		double[] B_tilta_up = new double[data_stats.length + 1];

		int[] B_tilta_low_index = new int[data_stats.length + 1];
		int[] B_tilta_up_index = new int[data_stats.length + 1];


		for (int i = 0; i < B_tilta_up.length; i++) {
			B_tilta_low[i] = -Double.MAX_VALUE/2;
			B_tilta_up[i] = Double.MAX_VALUE/2;
		}
		
		for (int j = 1; j < data_stats.length + 1; j++) {
			for (int k = 1; k <= j; k++) {
				if (B_tilta_low[k] < b_low[k]) {
					B_tilta_low[k] = b_low[k];
					B_tilta_low_index[k] = b_low_index[k];
				}
			}
		}

		for (int j = 0; j < data_stats.length; j++) {
			for (int k = j; k < data_stats.length; k++) {
				if (B_tilta_up[k] > b_up[k]) {
					B_tilta_up[k] = b_up[k];
					B_tilta_up_index[k] = b_up_index[k];
				}
			}
		}

		if (debugging_stop_criteria) {
			/*For debug*/
			System.out.println("\n B_tilta_up:");
			System.out.println(Arrays.toString(B_tilta_up));
			System.out.println("\n B_tilta_up_index:");
			System.out.println(Arrays.toString(B_tilta_up_index));
			System.out.println("\n B_tilta_low:");
			System.out.println(Arrays.toString(B_tilta_low));
			System.out.println("\n B_tilta_low_index:");
			System.out.println(Arrays.toString(B_tilta_low_index));
		}

		double[] B_low = new double[data_stats.length];
		double[] B_up = new double[data_stats.length];
		int[] B_low_index = new int[data_stats.length];
		int[] B_up_index = new int[data_stats.length];

		for (int j = 1; j < data_stats.length; j++) {
			if (sum_alpha[j] > accuracy) {
				B_low[j] = B_tilta_low[j + 1];
				B_low_index[j] = B_tilta_low_index[j + 1];
			} else {
				B_low[j] = B_tilta_low[j];
				B_low_index[j] = B_tilta_low_index[j];
			}

			if (sum_alpha[j - 1] > accuracy) {
				B_up[j] = B_tilta_up[j - 1];
				B_up_index[j] = B_tilta_up_index[j - 1];
			} else {
				B_up[j] = B_tilta_up[j];
				B_up_index[j] = B_tilta_up_index[j];
			}
		}

		if (debugging_stop_criteria) {
			/*For debug*/
			System.out.println("\n B_up:");
			System.out.println(Arrays.toString(B_up));
			System.out.println("\n B_low:");
			System.out.println(Arrays.toString(B_low));
		}
			

		most_vio = new int[2];
		Set<Integer> most_vio_test_temp = new HashSet<Integer>();
		for (int j = 1; j < data_stats.length; j++) {
			if (tau < B_low[j] - B_up[j]) {
				tau = B_low[j] - B_up[j];
				most_vio[0] = B_low_index[j];
				most_vio[1] = B_up_index[j];
			} 

			if (B_low[j] - B_up[j] > 0) {
				most_vio_test_temp.add(B_low_index[j]);
				most_vio_test_temp.add(B_up_index[j]);
			}
		}

		most_vio_test = new LinkedHashSet<Integer>();
		most_vio_test.add(most_vio[0]);
		most_vio_test.add(most_vio[1]);
		for (Integer k : most_vio_test_temp) {
			most_vio_test.add(k);
		}


		if (debugging_stop_criteria) {
			/*For debug*/
			System.out.println("\n tau:");
			System.out.println(tau);
			System.out.println("\n working set:");
			System.out.println(Arrays.toString(most_vio));
		}
		return tau;
	}

	/**
     * projected_gradient_ascent
     * With this operation we can perform a projected gradient descent step
     */
	public void projected_gradient_ascent() {

		for (int iter = 0; iter < max_iteration; iter++) {
			curt_iter = iter;
			//initialize the working set
			working_set = new int[working_set_size];
			for (int i = 0; i < working_set_size; i++) {
				working_set[i] = (int) (data_size * Math.random());
			}

			
			if (most_vio != null) {
				working_set[0] = most_vio[0];
				working_set[1] = most_vio[1];
			}
			
			
			if (most_vio_test != null) {
				int i = 0;
				for (Integer val : most_vio_test) {
					if (i >= working_set.length) {
						break;
					}
					working_set[i] = val;
					i++;
				}
			}
			
			//System.out.println(Arrays.toString(most_vio));
			//System.out.println(Arrays.toString(working_set));


			Arrays.sort(working_set);
			int length = removeDuplicates(working_set);
			int[] temp_work_set = new int[length];
			for (int i = 0; i < length; i++) {
				temp_work_set[i] = working_set[i];
			}
			working_set = temp_work_set;
			Arrays.sort(working_set);

			//compute necessary information for creating the subproblem
			this.preprocess_working_set();

			old_alpha = alpha.clone();
			//gradient descent n_grads steps
			for (int inner_iter = 0; inner_iter < n_grad; inner_iter++) {

				
				//Compute the gradient
				long time_gradient_update_start = System.currentTimeMillis();
				
				this.graient_ascent_single_step();

				//move on to the other batch
				batch_curt_index += batch_size;
				if (batch_curt_index>= data_size) {
					batch_curt_index -= data_size;
				}
				
				time_gradient_update += System.currentTimeMillis() - time_gradient_update_start;

				/*Create the projection subproblem
				*Require:
				* working_set_variable_index
				* gradient_container
				*
				*/
				long time_simplfying_start = System.currentTimeMillis();
				
				this.create_projection_subproblem();

				time_simplfying += System.currentTimeMillis() - time_simplfying_start;


				//solve the projection subproblem
				long time_solving_RAPNC_start = System.currentTimeMillis();
				
				//
				ResultTypeRAPNC_Continuous res_simplified = null;

				if (!isMDA) {
					res_simplified = projection_subproblem.solveContinuousDCA();
				} else {
					/*MDA*/
					ResultTypeMDA_Continuous res_MDA = projection_subproblem.FastMDA();
					res_simplified = new ResultTypeRAPNC_Continuous(true, res_MDA.aa);
				}

				number_of_RAP += projection_subproblem.number_subproblem;
				
				if (!res_simplified.feasible) {
					System.out.println(Arrays.toString(projection_subproblem.lbVar));
					System.out.println(Arrays.toString(projection_subproblem.ubVar));
					System.out.println(Arrays.toString(projection_subproblem.lbNested));
					System.out.println(Arrays.toString(projection_subproblem.ubNested));
					double lower = 0;
					double upper = 0;
					System.out.println("Wrong at ");
					for (int i = 0; i < working_set.length * 2; i++) {
						lower += projection_subproblem.lbVar[i];
						upper += projection_subproblem.ubVar[i];

						if (lower > projection_subproblem.ubNested[i] || upper < projection_subproblem.lbNested[i]) {
							System.out.println("Wrong at " + i);
						}
					}
					System.out.println("One subproblem is infeasible");
				}

				if (false || debugging) {
					for (int i = 0; i < working_set_variable_index.length; i++) {
						int index = working_set_variable_index[i];
						System.out.println("\nIndex " + index + " : ");
						System.out.println(alpha[index] + gradient_container[i]);
						System.out.println(res_simplified.sol[i]);
					}
				}
				
				
				//ResultTypeRAPNC_Continuous res_simplified = projection_subproblem.force_feasibility(res_simplified);
				//enforce a feasible solution
				{
					double sum_vari = 0;
					double min_var = 1e9;
					double max_var = -1e9;
					int min_var_index = -1;
					int max_var_index = -1;

					for (int i = 0; i < res_simplified.sol.length; i++) {

						sum_vari += res_simplified.sol[i];

						if (Math.abs(projection_subproblem.lbVar[i] - projection_subproblem.ubVar[i]) < 1e-9){
							continue;
						} 

						if (res_simplified.sol[i] > max_var) {
							max_var_index = i;
							max_var = res_simplified.sol[i];
						} 

						if (res_simplified.sol[i] < min_var) {
							min_var_index = i;
							min_var = res_simplified.sol[i];
						} 
						
					}

					double B = projection_subproblem.ubNested[projection_subproblem.ubNested.length - 1];
					//System.out.println("Here Here!!:" + sum_vari);
					//System.out.println(B);
					//System.out.println(sum_vari == B);
					//System.out.println(Math.abs(sum_vari - B) < 1e-9);
					if (sum_vari - B > 0) {
						res_simplified.sol[max_var_index] -= sum_vari - B;
					} else {
						res_simplified.sol[min_var_index] -= sum_vari - B;
					}
				}

				if (false && debugging) {
					for (int i = 0; i < working_set_variable_index.length; i++) {
						int index = working_set_variable_index[i];
						System.out.println(alpha[index] + gradient_container[i]);
						System.out.println(res_simplified.sol[i]);
					}
				}

				time_solving_RAPNC += System.currentTimeMillis() - time_solving_RAPNC_start;

				//update the alpha's beta's sum's
				long time_gradient_update_start2 = System.currentTimeMillis();

				if (true) {
					old_alpha_for_incre = new double[data_size * 2];
					for (int i = 0; i < working_set_variable_index.length; i++) {
						int index = working_set_variable_index[i];
						old_alpha_for_incre[index] = alpha[index];
					}
				}	

				this.update_variable_with_projected_solution(res_simplified);
				
				
				if (true) {
					this.update_discriminant_function_working_set(); // update the gradient incrementally
				}

				time_gradient_update += System.currentTimeMillis() - time_gradient_update_start2;		
			}

			long time_test_convergence_start = System.currentTimeMillis();
			
			//update the discriminant function incrementally
			this.update_discriminant_function();

			//compute the stopping criterion, test the convergence
			double tau = 1e5;
			if (find_most_violated_samples_KKT_not_each_hyperplane) {
				tau = this.stop_criterion_select_most_violated();
			} else {
				tau = this.stop_criterion();
			}

			if (tau < 0.001) {
				System.out.println("Terminate by the stopping criteria.");
				System.out.println("tau at iteration " + iter + " : " + tau);
				System.out.println("optimal obj at stopped: " + this.compute_obj_val());
				return;
			}
			
			if (display_obj_flag && Math.random() < print_frequency) {
				System.out.println("tau at iteration " + iter + " : " + tau);
				double obj_val_new = this.compute_obj_val();
				if (Math.abs(obj_val - obj_val_new) < 1e-5) {
					//System.out.println("optimal obj at stopped: " + this.compute_obj_val());
					//return;
				}
				obj_val = obj_val_new;
				System.out.println("obj: " + obj_val);
			}
			
			time_test_convergence += System.currentTimeMillis() - time_test_convergence_start;

		}	

		double tau = this.stop_criterion();
		System.out.println("Terminate by the the max_iteration.");
		System.out.println("tau at iteration " + max_iteration + " : " + tau);
		System.out.println("obj at stopped: " + this.compute_obj_val());
	}

	/**
     * preprocess_working_set()
     * With this operation we compute the necessary information.
     */
	public void preprocess_working_set() {
		working_set_stats = new int[data_stats.length];
		working_set_stats_cumulative = new int[data_stats.length];
		working_set_variable_index = new int[2 * working_set.length];
		working_set_variable_isDummy = new int[2 * working_set.length];
		working_set_variable_isBeta = new int[2 * working_set.length];


		//the number of samples by classes
		for (int i = 0; i < working_set.length; i++) {
			int sample_index = working_set[i];
			working_set_stats[class_global[sample_index] - 1]++;
		}

		//accumulate
		working_set_stats_cumulative[0] = working_set_stats[0];
		for (int i = 1; i < data_stats.length; i++) {
			working_set_stats_cumulative[i] += working_set_stats[i];
			working_set_stats_cumulative[i] += working_set_stats_cumulative[i - 1];
		}

		//the variable indice
		for (int i = 0; i < working_set_stats[0]; i++) {
			int sample_index = working_set[i];
			working_set_variable_index[i] = beta_index_global[sample_index];
			working_set_variable_index[i + working_set_stats[0]] = alpha_index_global[sample_index];

			working_set_variable_isBeta[i] = 1;
			working_set_variable_isDummy[i] = 1;
		}
		for (int j = 1; j < working_set_stats.length; j++) {
			
			int start_index = working_set_stats_cumulative[j - 1];
			int end_index = working_set_stats_cumulative[j];
			
			for (int i = start_index; i < end_index; i++) {
				int temp_index = working_set[i];
				working_set_variable_index[i + working_set_stats_cumulative[j - 1]] = beta_index_global[temp_index];
				working_set_variable_index[i + working_set_stats_cumulative[j]] = alpha_index_global[temp_index];

				working_set_variable_isBeta[i + working_set_stats_cumulative[j - 1]] = 1;
			}
		}
		for (int i = 0; i < working_set_stats[working_set_stats.length - 1]; i++) {
			working_set_variable_isDummy[2 * working_set.length - 1 - i] = 1;
		}

		if (debugging) {
			/*For debug*/
			System.out.println("\n Working set");
			System.out.println(Arrays.toString(working_set));
			System.out.println("\n Working set stats");
			System.out.println(Arrays.toString(working_set_stats));
			System.out.println("\n Working set stats cumulative");
			System.out.println(Arrays.toString(working_set_stats_cumulative));
			System.out.println("\n Working set variable index");
			System.out.println(Arrays.toString(working_set_variable_index));
			System.out.println("\n Working set variable isDummy");
			System.out.println(Arrays.toString(working_set_variable_isDummy));
			System.out.println("\n Working set variable isBeta");
			System.out.println(Arrays.toString(working_set_variable_isBeta));
		}
		
	}

	/**
     * Compute the gradient for the variables in the working sets.
     */
	public void graient_ascent_single_step() {

		gradient_container = new double[working_set_variable_index.length];

 		int[] batch_gradient_index = new int[batch_size];
		for (int i = 0; i < batch_size; i++) {
			if (batch_curt_index + i < data_size) {
				batch_gradient_index[i] = batch_curt_index + i;
			} else {
				batch_gradient_index[i] = batch_curt_index + i - data_size;
			}
		}

		if (debugging) {
			/*For Debug*/
			//System.out.println(Arrays.toString(batch_gradient_index));
		}

		if (isStochasticBatchGradient) {
			batch_gradient_index = new int[batch_size];
			for (int i = 0; i < batch_size; i++) {
				batch_gradient_index[i] = (int) (data_size * Math.random());
			}
		}
		
		
		for (int i = 0; i < working_set_variable_index.length; i++) {

			int index = working_set_variable_index[i];

			//skip that variable if it is a dummy variables
			if (working_set_variable_isDummy[i] == 1) {
				continue;
			}

			double gradient = 0;
			if (working_set_variable_isBeta[i] == 1) {
				gradient = -1; 
			} else {
				gradient = 1;
			}

			/* compute the gradient naively 
			int last_batch_index = batch_curt_index + batch_size;
			//for (int k : batch_gradient_index) {
			for (int k = 0; k < data_size; k++) {	
				if (debugging) {
					//System.out.println("sample in the batch k = :" + k);
				}
				
				//Please notice that in this step, you can not use the online alpha's as they may be infeasible!.
				double temp_diff_alpha = alpha[beta_index_global[k]] + alpha[alpha_index_global[k]] - constant_C;
				
				if (false && debugging) {
					System.out.println("difference in beta alpha" + temp_diff_alpha);
					System.out.println("coefficients" + coefficients[variable_original_index[index]][k]);
				}
				

				temp_diff_alpha *= coefficients[variable_original_index[index]][k];

				gradient -= temp_diff_alpha;

			}
			*/


			if (debugging) {
				/*For Debug
				System.out.println("i = :" + i);
				System.out.println("For sample of index: " + variable_original_index[index]);
				System.out.println("variable index:" + index);
				System.out.println("It is a beta variable: " + working_set_variable_isBeta[i]);
				System.out.println("Gradient in a single iteration: " + gradient);
				*/
			}		
			
			//compute the gradient with the discriminant function
			gradient += discriminant_function[variable_original_index[index]];

			//alpha[index] += gradient;
			//gradient ascent
			gradient_container[i] = step_size * gradient;

			/*For debug
			System.out.println("i = :" + i);
			System.out.println("Gradient: " + gradient);
			System.out.println("alpha:" + alpha[index]);
			*/
		}
	}

	/**
	 * create_projection_subproblem() 
     * After the gradient step, create the projection subproblem.
     */
	public void create_projection_subproblem() {
		List<Function> obj = new ArrayList<>(); 
		double[] lbVar_simplified = new double[working_set_variable_index.length];
		double[] ubVar_simplified = new double[working_set_variable_index.length];

		//obj function L_2-norm; lower and upper bound for a single variable
		for (int i = 0; i < working_set_variable_index.length; i++) {
			int index = working_set_variable_index[i];
			double a = 1000;
			double b = - 2000 * (alpha[index] + gradient_container[i]);

			obj.add(new QuadraticFunction(a, b));

			int sample_index = variable_original_index[index];
			lbVar_simplified[i] = 0;
			ubVar_simplified[i] = constant_C;
			if (class_global[sample_index] == 1 && working_set_variable_isDummy[i] == 1) {
				lbVar_simplified[i] = constant_C;
				ubVar_simplified[i] = constant_C;
			} else if (class_global[sample_index] == working_set_stats.length && working_set_variable_isDummy[i] == 1) {
				lbVar_simplified[i] = 0;
				ubVar_simplified[i] = 0;
			}
		}

		if (debugging) {
			System.out.println("\nSimplifing the projection subproblem:");
			System.out.println("working_set_variable_index: " + Arrays.toString(working_set_variable_index));
			System.out.println("lbVar_simplified" + Arrays.toString(lbVar_simplified));
			System.out.println("ubVar_simplified" + Arrays.toString(ubVar_simplified));
		}


		//nested bound
		double[] lbNested_simplified = new double[working_set_variable_index.length];
		double[] ubNested_simplified = new double[working_set_variable_index.length];
		for (int i = 0; i < lbNested_simplified.length; i++) {
			lbNested_simplified[i] = 0;
			ubNested_simplified[i] = 2 * constant_C * working_set.length;
		}
		double adjust_fixed_sum = 0; //sum of the current values of all non-fixed variables
		
		int start_index = 0;		
		int end_index = working_set_stats_cumulative[0];

		for (int j = 0; j < working_set_stats.length; j++) {

			if (j > 0) {
				start_index = end_index;
				end_index = working_set_stats_cumulative[j] + working_set_stats_cumulative[j-1];
			} 
			
			if (debugging) {
				System.out.println("\nStart_index: " + start_index);
				System.out.println("End_index: " + end_index);
			}
			
			
			for (int i = start_index; i < end_index; i++) {
				int index = working_set_variable_index[i];
				adjust_fixed_sum += alpha[index];
			} 

			if (end_index > 0) {
				//System.out.println(adjust_fixed_sum);
				double fixed_var_sum = sum_beta_alpha[j] - adjust_fixed_sum;

				if (debugging) {
					System.out.println("\nClass: " + j);
					System.out.println("sum: " + fixed_var_sum);
					System.out.println("lb_or: " + lbNested[j]);
					System.out.println("ub_or: " + ubNested[j]);
					System.out.println("lb: " + round(lbNested[j] - fixed_var_sum));
					System.out.println("ub: " + round(ubNested[j] - fixed_var_sum));
				}
				
				lbNested_simplified[end_index - 1] = Math.max(lbNested_simplified[end_index - 1], round(lbNested[j] - fixed_var_sum));
				ubNested_simplified[end_index - 1] = Math.min(ubNested_simplified[end_index - 1], round(ubNested[j] - fixed_var_sum));
			}
		}	

		{
			start_index = end_index;
			end_index = 2 * working_set_stats_cumulative[working_set_stats.length - 1];
			//System.out.println("Dummy variables" + end_index + "??" + start_index);
			for (int i = start_index; i < end_index; i++) {
				int index = working_set_variable_index[i];
				adjust_fixed_sum += alpha[index];

				//System.out.println("Dummy variables" + (alpha[index] == 0));
			} 
			
			lbNested_simplified[end_index - 1] = Math.max(lbNested_simplified[end_index - 1], round(adjust_fixed_sum));
			ubNested_simplified[end_index - 1] = Math.min(ubNested_simplified[end_index - 1], round(adjust_fixed_sum));
		}

		if (debugging) {
			System.out.println("\nSimplifing the projection subproblem:");
			System.out.println("working_set_stats_cumulative: " + Arrays.toString(working_set_stats_cumulative));
			System.out.println("lbNested_simplified" + Arrays.toString(lbNested_simplified));
			System.out.println("ubNested_simplified" + Arrays.toString(ubNested_simplified));
		}
		

		// final nested bound
		for (int i = 1; i < 2*working_set.length; i++) {
			if (lbNested_simplified[i] < lbNested_simplified[i-1]) {
				lbNested_simplified[i] = lbNested_simplified[i-1];
			}
		}
		for (int i = 2*working_set.length-2; i > -1; i--) {
			if (ubNested_simplified[i] > ubNested_simplified[i+1]) {
				ubNested_simplified[i] = ubNested_simplified[i+1];
			}
		}

		//force feasibility
		for (int i = 0; i < 2*working_set.length; i++) {
			for (int j = i; j < 2*working_set.length; j++) {
				if (lbNested_simplified[i] > ubNested_simplified[j]) {
					lbNested_simplified[i] = ubNested_simplified[j];
				}
			} 

			if (lbNested_simplified[i] > constant_C * (i + 1)) {
				lbNested_simplified[i] = constant_C * (i + 1);
			}
		}

		if (debugging) {
			System.out.println("\nFill redundant nested bounds:");
			System.out.println("lbNested_simplified" + Arrays.toString(lbNested_simplified));
			System.out.println("ubNested_simplified" + Arrays.toString(ubNested_simplified));
		}

		projection_subproblem = new RAPNC_Continuous(obj, lbVar_simplified, ubVar_simplified, lbNested_simplified, ubNested_simplified);

	}

	/**
     * With the gradient, create the projection subproblem.
     */
	public void update_variable_with_projected_solution(ResultTypeRAPNC_Continuous res_simplified) {
		double[] sol = res_simplified.sol;
		/*
		for (int j = 0; j < working_set_stats.length; j++) {

			int start_index = 0; 
			int end_index = 2 * working_set_stats_cumulative[j];
			if (j > 0) {
				start_index = 2 * working_set_stats_cumulative[j - 1];
			}

			for (int i = start_index; i < end_index; i++) {
				int index = working_set_variable_index[i];
				alpha[index] = sol[i];
			}
		}
		*/

		double adjust_fixed_sum = 0; //sum of the current values of all non-fixed variables
		double adjust_fixed_sum_alpha = 0; //sum of the orginal values of all non-fixed variables
		
		int start_index = 0;		
		int end_index = working_set_stats_cumulative[0];

		for (int j = 0; j < working_set_stats.length; j++) {

			if (j > 0) {
				start_index = end_index;
				end_index = working_set_stats_cumulative[j] + working_set_stats_cumulative[j-1];
			} 
			
			if (debugging) {
				System.out.println("\nAdjusting nested sum of variables: ");
				System.out.println("\nStart_index: " + start_index);
				System.out.println("End_index: " + end_index);
			}
			
			for (int i = start_index; i < end_index; i++) {
				int index = working_set_variable_index[i];

				adjust_fixed_sum += sol[i];
				adjust_fixed_sum_alpha += alpha[index];

				alpha[index] = sol[i];
			} 

			if (end_index > 0) {
				//System.out.println(adjust_fixed_sum);
				sum_beta_alpha[j] += adjust_fixed_sum - adjust_fixed_sum_alpha;
			}
		}	

		{
			start_index = end_index;
			end_index = 2 * working_set_stats_cumulative[working_set_stats.length - 1];
			for (int i = start_index; i < end_index; i++) {
				int index = working_set_variable_index[i];
				adjust_fixed_sum += alpha[index];
				adjust_fixed_sum_alpha += alpha[index];
			} 
			//System.out.println("adjust_fixed_sum: " + adjust_fixed_sum);
			//System.out.println("adjust_fixed_sum_alpha: " + adjust_fixed_sum_alpha);
			//System.out.println("adjust_fixed_sum_alpha: " + sum_beta_alpha[data_stats.length - 1]);
			sum_beta_alpha[data_stats.length - 1] += adjust_fixed_sum - adjust_fixed_sum_alpha;
			//System.out.println("adjust_fixed_sum_alpha: " + sum_beta_alpha[data_stats.length - 1]);   
			//double fixed_var_sum = sum_beta_alpha[working_set_stats.length - 1] - adjust_fixed_sum;
		}

		for (int j = 0; j < data_stats.length; j++) {
			if (sum_beta_alpha[j] - data_stats_cumulative[j] * constant_C > 1e-9) {
				sum_alpha[j] = sum_beta_alpha[j] - data_stats_cumulative[j] * constant_C;
			} else {
				sum_alpha[j] = 0;
			}
		}

		if (debugging) {
			System.out.println(Arrays.toString(sum_alpha));
		}

	}


	/**
	 * update the discriminant_function()
	 * Update the discriminant function for each sample incrementally 
	 * For the samples in the working set, there is no need to update as their values has already been updated
     */
	public void update_discriminant_function() {
		int indicator_working_set = 0;
		indicator_working_set = 0;
		for (int i = 0; i < data_size; i++) {
			if (indicator_working_set < working_set.length && i == working_set[indicator_working_set]) {
				indicator_working_set++;
				continue;
			}
			double change = 0;
			for (int j = 0; j < working_set.length; j++) {
				int index = working_set[j];
				double temp = - alpha[beta_index_global[index]] - alpha[alpha_index_global[index]];
				temp -= - old_alpha[beta_index_global[index]] - old_alpha[alpha_index_global[index]];

				temp *= coefficients[i][index];
				change += temp;
			}
			//System.out.println(change);
			discriminant_function[i] += change;
		}
	}

	/**
	 * update_discriminant_function_working_set()
     * Update the discriminant_function for variables in the working set. This routine will be called in every inner iteration
     */
	public void update_discriminant_function_working_set() {
		for (int i = 0; i < working_set.length; i++) {
			int curt_index = working_set[i];
			double change = 0;
			for (int j = 0; j < working_set.length; j++) {
				int index = working_set[j];
				double temp = - alpha[beta_index_global[index]] - alpha[alpha_index_global[index]];
				temp -= - old_alpha_for_incre[beta_index_global[index]] - old_alpha_for_incre[alpha_index_global[index]];

				temp *= coefficients[curt_index][index];
				change += temp;
			}
			//System.out.println(change);
			discriminant_function[curt_index] += change;
		}
	}

	//main function
	public static void main(String[] args) {

		TestSVM_Continuous2 test_one = new TestSVM_Continuous2();

		//read data
		long start_time = System.currentTimeMillis();
		test_one.read_data_test(test_one.max_input_size);
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
		System.out.println(
					"Total time: "
					+ String.format("%10s", ((double) time) / 1000));
		
	}

	/**
	 * display_time()
     * Print all the information about this run 
     */
	public void display_time() {
		if (!isMDA) {
			System.out.println(
					"Solution time while using DCA"
			);
		} else {
			System.out.println(
					"Solution time while using MDA"
			);
		}
		if (find_most_violated_samples_KKT_not_each_hyperplane) {
			System.out.println(
					"The working_set contains the most violated sample pairs for KKT conditions"
			);
		} else {
			System.out.println(
					"The working_set contains the most violated sample pairs for KKT conditions for each hyperplane."
			);
		}
		System.out.println(
					"Number of samples: "
					+ data_size
		);
		System.out.println(
					"Maximum Iteration: "
					+ max_iteration
		);
		System.out.println(
					"Regularizer coefficient: C = "
					+ constant_C
		);
		System.out.println(
					"working_set_size: "
					+ working_set_size
		);
		System.out.println(
					"Iteration when the program stopped: "
					+ (curt_iter + 1)
		);
		System.out.println(
					"Total time in solving RAPNC: "
					+ String.format("%10s", ((double) time_solving_RAPNC) / 1000) 
		);
		System.out.println(
					"Total time in gradient update: "
					+ String.format("%10s", ((double) time_gradient_update) / 1000)
		);
		System.out.println(
					"Total time in simplifying RAP: "
					+ String.format("%10s", ((double) time_simplfying) / 1000)
		);
		System.out.println(
					"Total time in checking convergence: "
					+ String.format("%10s", ((double) time_test_convergence) / 1000)
		);
	}

	/**
	 * compute_obj_val()
     * This method compute the current obj value to make sure the gradient ascent algorithm is working
     * It is an O(n^2) operation and should be turn off in the numerical experiment 
     */
	public double compute_obj_val () {
		double sum = 0;

		for (int i = 0; i < data_size; i++) {
			sum += alpha[alpha_index_global[i]] + constant_C - alpha[beta_index_global[i]];
			for (int k = 0; k < data_size; k++) {
				sum -= 0.5 * (constant_C - alpha[alpha_index_global[i]] - alpha[beta_index_global[i]]) * (constant_C - alpha[alpha_index_global[k]] - alpha[beta_index_global[k]]) * coefficients[i][k];
			}
		}

		return sum;
	}
	
	/**
	 * removeDuplicates()
     * A simple helper function to move the duplicate elements to the end of the array 
     */
	public static int removeDuplicates(int[] A) {
		if (A.length < 2)
			return A.length;
	 
		int j = 0;
		int i = 1;
	 
		while (i < A.length) {
			if (A[i] != A[j]) {
				j++;
				A[j] = A[i];
			}
	 
	                i++;
		}
	 
		return j + 1;
	}

	/**
	 * round()
     * A simple helper function. It rounds any val to 1e-9 precision to avoid floating point error
     */
	public static double round(double val) {
		val *= 1e9;
		val = Math.round(val);
		val /= 1e9;

		return val;
	}

	/**
	 * An helper class for recording the most violated sample pairs for the KKT conditions
	 */
	public class ValIndexPair{
		int index;
		double val;
		public ValIndexPair(double val, int index) {
			this.index = index;
			this.val = val;
		}
	}
	/**
     * stop_criterion
     * With this operation we can compute tau, the stopping criteria
     * It will also record the set of most violated samples pairs for the KKT conditions
     * Notice that it is different from the stop_criterion(). 
     * <p>
     * @param discriminant
     * @param alpha
     * @param sum_alpha
     * The detail can be found in Support Vector Ordinal Regression (Wei Chu and S. Sathiya Keerthi 2007)
     */
	
	public double stop_criterion_select_most_violated() {
		double tau = -1e10;
		double accuracy = accuracy_in_stopping_criteria;

		double[] b_up = new double[data_stats.length+1];
		double[] b_low = new double[data_stats.length+1];
		int[] b_low_index= new int[data_stats.length+1];
		int[] b_up_index = new int[data_stats.length+1];
		List<TreeSet<ValIndexPair>> b_low_queue= new ArrayList<>();
		List<TreeSet<ValIndexPair>> b_up_queue = new ArrayList<>();

		for (int i = 0; i < data_stats.length + 1; i++) { 
            TreeSet<ValIndexPair> b_low_queue_temp = new TreeSet<ValIndexPair>(new Comparator<ValIndexPair>() {
    			@Override
    			public int compare(ValIndexPair a, ValIndexPair b) {
        			if (a.val - b.val > 0) {
        				return 1;
        			} else if (a.val - b.val < 0) {
        				return -1;
        			} else {
        				return 0;
        			}
    			}
			}); 
            b_low_queue.add(b_low_queue_temp);
            TreeSet<ValIndexPair> b_up_queue_temp = new TreeSet<ValIndexPair>(new Comparator<ValIndexPair>() {
    			@Override
    			public int compare(ValIndexPair a, ValIndexPair b) {
        			if (a.val - b.val > 0) {
        				return -1;
        			} else if (a.val - b.val < 0) {
        				return 1;
        			} else {
        				return 0;
        			}
    			}
			}); 
			b_up_queue.add(b_up_queue_temp);
        }

		for (int i =0; i < b_low.length; i++) {
			b_low[i] = -Double.MAX_VALUE/2;
			b_up[i] = Double.MAX_VALUE/2;
			b_low_index[i] = -1;
			b_up_index[i] = -1;
		}

		for (int i = 0; i < data_stats_cumulative[0]; i++) {
			if (b_up[0] > discriminant_function[i] - 1) {
				b_up[0] = discriminant_function[i] - 1;
				b_up_index[0] = i;
			}

			this.add_index_for_up(b_up_queue.get(0), discriminant_function[i] - 1, i);
		}

		/*for debug*/
		if (b_up_queue.get(0).last().index != b_up_index[0]) {
			System.out.println("Wrong in the TreeSet!");
			System.out.println(b_up_index[0]);
			System.out.println(b_up_queue.get(0).last().index);
		}
		/*for debug*/
		

		for (int i = 0; i < data_size; i++) {
			int class_i = class_global[i];
			double val = alpha[beta_index_global[i]];

			if (Math.abs(val - constant_C) < accuracy) {				
				// if alpha_star_i^j is close to zero	
				if (b_up[class_i - 1] > discriminant_function[i] - 1) {
					b_up[class_i - 1] = discriminant_function[i] - 1;
					b_up_index[class_i - 1] = i;
				}

				this.add_index_for_up(b_up_queue.get(class_i - 1), discriminant_function[i] - 1, i);

			} else if (Math.abs(val) < accuracy) {
				// if alpha_star_i^j is close to constant_C	
				if (b_low[class_i - 1] < discriminant_function[i] - 1) {
					b_low[class_i - 1] = discriminant_function[i] - 1;
					b_low_index[class_i - 1] = i;
				}

				this.add_index_for_low(b_low_queue.get(class_i - 1), discriminant_function[i] - 1, i);
			} else {
				// if alpha_star_i^j is not null	
				if (b_up[class_i - 1] > discriminant_function[i] - 1) {
					b_up[class_i - 1] = discriminant_function[i] - 1;
					b_up_index[class_i - 1] = i;
				}
				
				this.add_index_for_up(b_up_queue.get(class_i - 1), discriminant_function[i] - 1, i);

				if (b_low[class_i - 1] < discriminant_function[i] - 1) {
					b_low[class_i - 1] = discriminant_function[i] - 1;
					b_low_index[class_i - 1] = i;
				}

				this.add_index_for_low(b_low_queue.get(class_i - 1), discriminant_function[i] - 1, i);
			}

			val = alpha[alpha_index_global[i]];

			if (Math.abs(val) < accuracy) {				

				// if alpha_i^j is close to zero
				if (b_low[class_i] < discriminant_function[i] + 1) {
					b_low[class_i] = discriminant_function[i] + 1;
					b_low_index[class_i] = i;
				}				

				this.add_index_for_low(b_low_queue.get(class_i), discriminant_function[i] + 1, i);

			} else if (Math.abs(val - constant_C) < accuracy) {
				// if alpha_i^j is close to constant_C	
				if (b_up[class_i] > discriminant_function[i] + 1) {
					b_up[class_i] = discriminant_function[i] + 1;
					b_up_index[class_i] = i;
				}

				this.add_index_for_up(b_up_queue.get(class_i), discriminant_function[i] + 1, i);
			} else {
				// if alpha_i^j is not null	
				if (b_low[class_i] < discriminant_function[i] + 1) {
					b_low[class_i] = discriminant_function[i] + 1;
					b_low_index[class_i] = i;
				}
				if (b_up[class_i] > discriminant_function[i] + 1) {
					b_up[class_i] = discriminant_function[i] + 1;
					b_up_index[class_i] = i;
				}

				this.add_index_for_up(b_up_queue.get(class_i), discriminant_function[i] + 1, i);

				this.add_index_for_low(b_low_queue.get(class_i), discriminant_function[i] + 1, i);

			}
		}

		/* For debug */
		for (int i = 0; i < b_up_index.length; i++) {
			//System.out.println("Set index: " + i);
			//System.out.println(b_up_index[i]);
			if (b_up_queue.get(i).isEmpty() && b_up_index[i] == -1) {
				continue;
			}
			if (b_up_queue.get(i).last().index != b_up_index[i]) {
				System.out.println("Wrong in the TreeSet!");
				System.out.println(b_up_index[i]);
				System.out.println(b_up_queue.get(i).last().index);
			}
		}

		for (int i = 0; i < b_up_index.length; i++) {
			if (b_low_queue.get(i).isEmpty() && b_low_index[i] == -1) {
				continue;
			}
			if (b_low_queue.get(i).last().index != b_low_index[i]) {
				System.out.println("Wrong in the TreeSet!2");
				System.out.println(b_low_index[i]);
				System.out.println(b_low_queue.get(i).last().index);
			}
		}

		tau = this.find_most_violated_samples(b_up_queue, b_low_queue);

		return tau;
	}

	public void add_index_for_up(TreeSet<ValIndexPair> set, double val, int index) {
		if (set.size() < working_set_size) {
			set.add(new ValIndexPair(val, index));
		} else if (set.first().val > val) {
			set.pollFirst();
			set.add(new ValIndexPair(val, index));
		}
	}

	public void add_index_for_low(TreeSet<ValIndexPair> set, double val, int index) {
		if (set.size() < working_set_size) {
			set.add(new ValIndexPair(val, index));
		} else if (set.first().val < val) {
			set.pollFirst();
			set.add(new ValIndexPair(val, index));
		}
	}
	
	public double find_most_violated_samples(List<TreeSet<ValIndexPair>> b_up_queue, List<TreeSet<ValIndexPair>> b_low_queue) {
		Set<Integer> res = new LinkedHashSet<Integer>();
		double tau_val = -1e10;
		boolean first_iteration = true;


		while (res.size() < working_set_size) {
			double tau = -1e10;
			double accuracy = accuracy_in_stopping_criteria;
			double[] b_up = new double[data_stats.length+1];
			double[] b_low = new double[data_stats.length+1];
			int[] b_low_index= new int[data_stats.length+1];
			int[] b_up_index = new int[data_stats.length+1];

			for (int i = 0; i < b_up_queue.size(); i++) {
				//System.out.println("Set index: " + i);
				//System.out.println(b_up_index[i]);
				if (b_up_queue.get(i).isEmpty()) {
					b_up_index[i] = -1;
					b_up[i] = Double.MAX_VALUE/2;
				} else {
					b_up_index[i] = b_up_queue.get(i).last().index;
					b_up[i] = b_up_queue.get(i).last().val;
				}

				if (b_low_queue.get(i).isEmpty()) {
					b_low_index[i] = -1;
					b_low[i] = -Double.MAX_VALUE/2;
				} else {
					b_low_index[i] = b_low_queue.get(i).last().index;
					b_low[i] = b_low_queue.get(i).last().val;
				}
			}

			double[] B_tilta_low = new double[data_stats.length + 1];
			double[] B_tilta_up = new double[data_stats.length + 1];

			int[] B_tilta_low_index = new int[data_stats.length + 1];
			int[] B_tilta_up_index = new int[data_stats.length + 1];


			for (int i = 0; i < B_tilta_up.length; i++) {
				B_tilta_low[i] = -Double.MAX_VALUE/2;
				B_tilta_up[i] = Double.MAX_VALUE/2;
			}
			
			for (int j = 1; j < data_stats.length + 1; j++) {
				for (int k = 1; k <= j; k++) {
					if (B_tilta_low[k] < b_low[k]) {
						B_tilta_low[k] = b_low[k];
						B_tilta_low_index[k] = b_low_index[k];
					}
				}
			}

			for (int j = 0; j < data_stats.length; j++) {
				for (int k = j; k < data_stats.length; k++) {
					if (B_tilta_up[k] > b_up[k]) {
						B_tilta_up[k] = b_up[k];
						B_tilta_up_index[k] = b_up_index[k];
					}
				}
			}

			if (debugging_stop_criteria) {
				/*For debug*/
				System.out.println("\n B_tilta_up:");
				System.out.println(Arrays.toString(B_tilta_up));
				System.out.println("\n B_tilta_up_index:");
				System.out.println(Arrays.toString(B_tilta_up_index));
				System.out.println("\n B_tilta_low:");
				System.out.println(Arrays.toString(B_tilta_low));
				System.out.println("\n B_tilta_low_index:");
				System.out.println(Arrays.toString(B_tilta_low_index));
			}

			double[] B_low = new double[data_stats.length];
			double[] B_up = new double[data_stats.length];
			int[] B_low_index = new int[data_stats.length];
			int[] B_up_index = new int[data_stats.length];

			for (int j = 1; j < data_stats.length; j++) {
				if (sum_alpha[j] > accuracy) {
					B_low[j] = B_tilta_low[j + 1];
					B_low_index[j] = B_tilta_low_index[j + 1];
				} else {
					B_low[j] = B_tilta_low[j];
					B_low_index[j] = B_tilta_low_index[j];
				}

				if (sum_alpha[j - 1] > accuracy) {
					B_up[j] = B_tilta_up[j - 1];
					B_up_index[j] = B_tilta_up_index[j - 1];
				} else {
					B_up[j] = B_tilta_up[j];
					B_up_index[j] = B_tilta_up_index[j];
				}
			}

			if (debugging_stop_criteria) {
				/*For debug*/
				System.out.println("\n B_up:");
				System.out.println(Arrays.toString(B_up));
				System.out.println("\n B_low:");
				System.out.println(Arrays.toString(B_low));
			}
				

			most_vio = new int[2];
			most_vio_test = new HashSet<Integer>();
			for (int j = 1; j < data_stats.length; j++) {
				if (tau < B_low[j] - B_up[j]) {
					tau = B_low[j] - B_up[j];
					most_vio[0] = B_low_index[j];
					most_vio[1] = B_up_index[j];
				}
			}
			
			if (debugging_stop_criteria) {
				/*For debug*/
				System.out.println("\n tau:");
				System.out.println(tau);
				System.out.println("\n working set:");
				System.out.println(Arrays.toString(most_vio));
			}
			
			if (tau > 0) {
				res.add(most_vio[0]);
				res.add(most_vio[1]);
			} else {
				break;
			}
			

			for (int i = 0; i < b_up_queue.size(); i++) {
				//System.out.println("Set index: " + i);
				//System.out.println(b_up_index[i]);
				if (!b_up_queue.get(i).isEmpty() 
					&& res.contains(b_up_queue.get(i).last().index)
					){ 
					b_up_queue.get(i).pollLast();
				}

				if (!b_low_queue.get(i).isEmpty() 
					&& res.contains(b_low_queue.get(i).last().index)
					){ 
					b_low_queue.get(i).pollLast();
				}
			}

			if (first_iteration) {
				tau_val = tau;
			}
			first_iteration = false;

			if (debugging) {
				System.out.println("\nMost violated samples");
				System.out.println(Arrays.toString(most_vio));
			}
		}

		most_vio_test = res;
		
		return tau_val;

	}
	

	public void read_data_general(int size, ReadFile get_data) {
		//read data
		//ReadFileLineByLineUsingScannerAbalone get_data = new ReadFileLineByLineUsingScannerAbalone();
		//ReadFileBank32nh get_data = new ReadFileBank32nh();
		//ReadFileBostonHousing get_data = new ReadFileBostonHousing();
		//ReadFileCaliforniaHousing get_data = new ReadFileCaliforniaHousing();
		
		//
		data_raw = get_data.readdata(size);
		data_stats = data_raw.get(data_raw.size() - 1).stats;
		
		//cumulative 
		data_stats_cumulative = new int[data_stats.length];
		int cumulative = 0;
		for (int i =0; i < data_stats.length; i++) {
			cumulative += data_stats[i];
			data_stats_cumulative[i] = cumulative;
		}

		System.out.println("\n Data Statistics:");
		System.out.println(Arrays.toString(data_stats));
		data_size =  data_raw.size() - 1;

		//initialized the index 
		//precompute the alpha and alpha star index for each sample i
		alpha_index_global = new int[data_size];
		beta_index_global = new int[data_size];
		class_global = new int[data_size];
		index_in_class = new int[data_size];
		variable_original_index = new int[data_size*2];

		for (int i = 0; i < data_size; i++) {
			int class_i = (int) data_raw.get(i).class_index;
			class_global[i] = class_i;
			index_in_class[i] = (int) data_raw.get(i).sample_index; 

			if (class_i > 1) {
				beta_index_global[i] = (int) (data_stats_cumulative[class_i - 2] + i);
			} else {
				beta_index_global[i] = i;
			}
			alpha_index_global[i] = (int) (beta_index_global[i] + data_stats[class_i - 1]);

			variable_original_index[beta_index_global[i]] = i;
			variable_original_index[alpha_index_global[i]] = i;
		}

		//precompute the inner product
		coefficients = new double[data_size][data_size];
		for (int i = 0; i < data_size; i++) {
			ReadFile.Data_format data_i = data_raw.get(i);
			for (int j = 0; j < i; j++) {
				//System.out.println("Check here: computing coefficients: " + j);
				ReadFile.Data_format data_j = data_raw.get(j);
				double sum = 0;
				for (int k = 0; k < data_j.features.length - 1; k++) {
					//System.out.println("Check here: computing coefficients: " + sum);
					sum += - 0.5 * (data_j.features[k] - data_i.features[k]) * (data_j.features[k] - data_i.features[k]);
				}
				//System.out.println("Check here: computing coefficients: " + sum);
				//System.out.println(Arrays.toString(data_i.features));
				//System.out.println(Arrays.toString(data_j.features));
				coefficients[i][j] = Math.exp(sum);
			}
		} 

		for (int i = 0; i < data_size; i++) {
			coefficients[i][i] = 1;
			for (int j = i + 1; j < data_size; j++) {
				coefficients[i][j] = coefficients[j][i];
			}
		} 

		/*For debug*/
		if (debugging) {
			System.out.println("\n beta index:");
			System.out.println(Arrays.toString(beta_index_global));
			System.out.println("\n alpha index:");
			System.out.println(Arrays.toString(alpha_index_global));
			System.out.println("\n Class wrt sample index");
			System.out.println(Arrays.toString(class_global));
			System.out.println("\n Index with in a class wrt sample index");
			System.out.println(Arrays.toString(index_in_class));
			System.out.println("\n variable_original_index:");
			System.out.println(Arrays.toString(variable_original_index));

			//System.out.println("\n coefficients");
			for (int i = 0; i < 0; i++) {
				System.out.println(Arrays.toString(coefficients[i]));
			}
		}
		
	}
}