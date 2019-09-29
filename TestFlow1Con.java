/******************************************************************************
 *  Compilation:  javac TestFlow1Con.java
 *  Execution:    java TestFlow1Con
 *  
 *  A test class to evaluate the numerical performance of DCA, MDA, and FlowS
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/
import java.util.*;
import java.io.*;

public class TestFlow1Con {

	/**
     * An internal class to define the Objective functions
     * The name is QuadraticFunction, you can change it to any (convex) function 
     * by overwriting the abstract method getValue(double x)  
     * @param a, b, etc. parameter of the function
     */
	private static class QuadraticFunction extends Function {
		double a;
		double b;
		public QuadraticFunction(double a, double b) {
			this.a = a;
			this.b = b;
			this.b = b * 10; //Crash function
		}
		public double getValue(double x) {
			return b * x;
			//return 10*b + a / x; //Crash function
			//return a * b * b / x / x / x; // Fuel function
			//return x * x * x * x / 4 +  b * x ; // F function
			//return b * x;
		}
	}


	/**
     * An internal class to define the penalty objective functions defined in MDA
     * The name is QuadraticFunction, you can change it to any (convex) function 
     * by overwriting the abstract method getValue(double x)  
     * @param a, b, etc parameter of the function
     * @param c, d, boundary of the variables
     * @param L, Lipschitz constant
     */
	private static class QuadraticLOnePenalty extends Function {
		double a;
		double b;
		double c;
		double d;
		double L = 2;

		public QuadraticLOnePenalty(double a, double b, double c, double d) {
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
		}
		public double getValue(double x) {
			if (x < c) {
				return a * c + L*(c - x);
			} else if (x > d) {
				return a * d + L*(x - d);
			}
			return a * x ;
		}
	}

	/**
     * main method
     * initial the test and record the test result in a text file
     * You can modify the name, the header of the text file as well as the test cases
     */
	public static void main(String[] args) throws FileNotFoundException{
		//Set output environment
		long currentTime = System.currentTimeMillis();
		PrintStream o = new PrintStream(new File("DCA----FlowS----MDA----GreedyS" + " Convex_LINEAR_SUB " + Long.toString(currentTime) + "_Vb_100" + ".txt"));
		//PrintStream o = new PrintStream(new File("DCA----FlowS----MDA" + " Con " + "TTTTTTESSSTTTTT100" + ".txt"));
		System.setOut(o);

		
		//test settings (dimension)
		int[] testSize = new int[]{800, 1600, 3200, 6400, 12800, 25600, 51200, 51200<<1};
		testSize = new int[]{800, 1600, 3200, 6400, 12800, 25600, 51200, 51200<<1, 51200<<2, 51200<<3, 51200<<4,  51200<<5, 51200<<6, 51200<<7};
		//testSize = new int[]{51200<<7, 51200<<8};
		//testSize = new int[]{52100<<2, 51200<<3, 51200<<4,  51200<<5, 51200<<6, 51200<<7, 51200<<8};
		//int x = 51200;
		//testSize = new int[]{5};
		int rep = 10;
		int varBound = 100;

		System.out.println("Test: Variable bound " + varBound + ":General Convex Function");
		System.out.println("Dimension 		DCA 	FlowS 	FlowSWOS 	MDA 	FastMDA 	number_subproblem_DCA  		number_subproblem_MDA		GreedyS");

		for (int i = 0; i < testSize.length; i++) {
			for (int j = 0; j < rep; j++) {
				test(testSize[i], varBound);
			}
			System.out.println(" ");
		}
		
		return;
	}
    
    /**
     * test method
     * This is a static method that performs the tests.
     * You can modify the parameters to perform any tests 
     * The generation procedure ensures that every test instance is feasible
     * Default:
     * @param size dimension of the problem
     * @param varBound related to the upper bound of the box constraints: d_i < 1.3 * varBound
     */
    public static void test (int size, int varBound) {
    	List<Function> obj = new ArrayList<Function>();
		List<Function> objMDA = new ArrayList<Function>();
		int dimension = size; 
		long[] lbVar = new long[dimension];
		long[] ubVar = new long[dimension];
		long[] lbNested = new long[dimension];
		long[] ubNested = new long[dimension];
		long B = (long) (0.75 * dimension * dimension);
		String[] varNames = new String[dimension];


		// Create box constraints on variables
      	for (int i = 0; i < dimension; i++) {
			long b = (long) ((Math.random() + 0.3)* varBound);
			long c = (long) (Math.random() * varBound);
			if (b < c) {
				lbVar[i] = b;
				ubVar[i] = c;
			} else {
				lbVar[i] = c;
				ubVar[i] = b;
			}				
			lbVar[i] = 0;
		}		
			
      		
      	// Set objective      		
      	for (int i = 0; i < dimension; i++) {
      		double a = Math.random();
      		double b = Math.random();
      		if (Math.random() < 0.5) {
      			b = b;
      		}
      		QuadraticFunction fi = new QuadraticFunction(a, b);
      		QuadraticLOnePenalty fiMDA = new QuadraticLOnePenalty(a, b, lbVar[i], ubVar[i]);
			obj.add(fi);
			objMDA.add(fiMDA);
      	}
            

      	//Nested constraints
      	double a = 0;
      	double b = 0;
      	//reset from 1 to dimension
      	for (int i = 0; i < dimension - 1; i++) {
      		a += lbVar[i] + Math.random() * (ubVar[i] - lbVar[i]);
      		b += lbVar[i] + Math.random() * (ubVar[i] - lbVar[i]);
      		if (a > b) {
      			lbNested[i] = (long) b;
      			ubNested[i] = (long) a;
      		} else {
      			lbNested[i] = (long) a;
      			ubNested[i] = (long) b;
      		}
      	}
            
        // Add constraint:\sum x_i = B
      	a += lbVar[dimension - 1] + Math.random() * (ubVar[dimension - 1] - lbVar[dimension - 1]);
      	b += lbVar[dimension - 1] + Math.random() * (ubVar[dimension - 1] - lbVar[dimension - 1]);
      	B = (long) Math.max(a, b);
      	// Greedy Algorithm
      	lbNested[dimension - 1] = B;
      	ubNested[dimension - 1] = B;
      	
      	/*
      	System.out.println("ubNested");        
        System.out.println(Arrays.toString(ubNested));
        System.out.println("lbNested");
        System.out.println(Arrays.toString(lbNested));	
      	*/
      	
        /********************************************************************************
        **
        **DCA: solve the problem by DCA and record the time
        **
        ********************************************************************************/
      	//
      	RAPNC test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
      	long startTime = System.currentTimeMillis();
		ResultTypeRAPNC res = test.solveIntegerDCA();
		//System.out.println("DCA Done!");
		//ResultTypeRAPNC res = test.solveIntegerLinearDCA();
		//solve the problem by DCA 100 times
		for (int i = 0; i < 0; i++) {
			test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
			//res = test.solveIntegerDCALinear();
			res = test.solveIntegerLinearDCA();
		}	

		//time 
		long endTime = System.currentTimeMillis();
		long timeDCA = endTime - startTime;
		//number of subproblems
		long number_subproblem_DCA = test.number_subproblem;
		/*
    	System.out.println("Our algorithm took " + dca + " milliseconds");

		if (res.feasible) {
			System.out.println("yes, the problem is feasible");
			//System.out.println(Arrays.toString(res.sol));
		} else {
			System.out.println("no, the problem is infeasible");
		}
		*/
		
		/********************************************************************************
        **
        **FlowS: solve the problem by SFA and record the time
        **
        ********************************************************************************/
      	//
		
		startTime = System.currentTimeMillis();
		ResultTypeRAPNC resFlow = new ResultTypeRAPNC(true, null);
		//for general convex objective functions
		test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
		resFlow.sol = test.createLotSizing().solve();
		//System.out.println("FLOWS Done!");
		//resFlow.sol = test.createLotSizing().solveLinear();
		//solve the problem by  100 times
		//resFlow2.sol = test.createLotSizingWOSegment().solveLinear();
		for (int i = 0; i < 0; i++) {
			test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
			//resFlow.sol = test.createLotSizing().solve();
			//resFlow.sol = test.createLotSizing().solveLinear();
		}	
		
		//time 
		endTime = System.currentTimeMillis();
		//System.out.println(Arrays.toString(resFlow.sol));

		long timeFlowS = (endTime - startTime);

		/********************************************************************************
        **
        **FlowS: solve the problem by FlowS and record the time Without Segment tree
        **
        ********************************************************************************/
      	//
		
		startTime = System.currentTimeMillis();
		test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
		ResultTypeRAPNC resFlow2 = new ResultTypeRAPNC(true, null);
		//for general convex objective functions
		resFlow2.sol = test.createLotSizingWOSegment().solve();
		//resFlow2.sol = test.createLotSizingWOSegment().solveLinear();
		//solve the problem by 100 times
		for (int i = 0; i < 0; i++) {
			test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
			//resFlow2.sol = test.createLotSizingWOSegment().solve();
			//resFlow2.sol = test.createLotSizingWOSegment().solveLinear();
		}	
		
		//time 
		endTime = System.currentTimeMillis();
		//System.out.println(Arrays.toString(resFlow.sol));

		long timeFlowSWOSegment= (endTime - startTime);

		
		/********************************************************************************
        **
        **MDA: solve the problem by MDA and record the time
        **
        ********************************************************************************/
      	//
		//RAPNC testMDA = new RAPNC(objMDA, lbVar, ubVar, lbNested, ubNested);
		startTime = System.currentTimeMillis();
		//ResultTypeMDA resMDA = testMDA.MDA();
		ResultTypeMDA resMDA = new ResultTypeMDA(resFlow.sol, resFlow.sol, resFlow.sol, resFlow.sol);

		//time 
		endTime = System.currentTimeMillis();
		long timeMDA = endTime - startTime;
		

		/********************************************************************************
        **
        **MDA: solve the problem by FastMDA and record the time
        **
        ********************************************************************************/
      	//
		RAPNC testFastMDA = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
		startTime = System.currentTimeMillis();
		ResultTypeMDA resFastMDA = testFastMDA.FastMDA();
		//System.out.println("MDA Done!");
		//ResultTypeMDA resFastMDA = new ResultTypeMDA(resFlow.sol, resFlow.sol, resFlow.sol, resFlow.sol);

		//time 
		endTime = System.currentTimeMillis();
		long timeFastMDA = endTime - startTime;

		//number of subproblems
		long number_subproblem_MDA = testFastMDA.number_subproblem;


		/********************************************************************************
        **
        **FlowS: solve the problem with GreedyS and record the time
        **
        ********************************************************************************/
      	//
		
		startTime = System.currentTimeMillis();
		ResultTypeRAPNC resGreedyS = new ResultTypeRAPNC(true, null);
		//for general convex objective functions
		test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
		resGreedyS.sol = test.createLotSizing().greedy_s_solve();
		
		//time 
		endTime = System.currentTimeMillis();
		//System.out.println(Arrays.toString(resFlow.sol));

		long timeGreedyS = (endTime - startTime);

    	//System.out.println("MDA took " + mda + " milliseconds");
    	//System.out.println("Flow: " + Arrays.toString(resFlow.sol));
    	//System.out.println("Our:" + Arrays.toString(res.sol));
    	//check whether the three algorithm produces the same solution
    	double sumOur = 0;
    	double sumFlow = 0;
    	long sumDemandMda = 0;
    	long sumDemandOur = 0;
    	for (int i = 0; i < dimension; i++) {
    		sumDemandOur += res.sol[i];
    		sumDemandMda += resFlow.sol[i];    		
			if (Math.abs(resFlow.sol[i] - res.sol[i]) >= 0.01 
				|| Math.abs(resFlow.sol[i] - resMDA.aa[i]) >= 0.01 
				|| Math.abs(resFlow.sol[i] - resFastMDA.aa[i]) >= 0.01 
				|| Math.abs(resFlow.sol[i] - resGreedyS.sol[i]) >= 0.01
			) {
				System.out.println("No, the solution x[" + i + "] is different.");
				System.out.println(resFlow.sol[i] + " " + res.sol[i] + " " + resMDA.aa[i] + " " + resFastMDA.aa[i] + " " + resGreedyS.sol[i]);
				System.out.println("Coefficient: " + obj.get(i).getValue(3) + obj.get(i).getValue(2));
			}
			//sumOur += obj.get(i).getValue(res.sol[i]);
			//sumFlow += obj.get(i).getValue(resFlow.sol[i]);
		}
	
		System.out.println(
				String.format("%10s", dimension) 
				+ String.format("%10s", ((double) timeDCA) / 1000) 
				+ String.format("%10s", ((double) timeFlowS) / 1000) 
				+ String.format("%10s", ((double) timeFlowSWOSegment) / 1000) 
				+ String.format("%10s", ((double) timeMDA) / 1000) 
				+ String.format("%10s", ((double) timeFastMDA) / 1000) 
				+ String.format("%20s", (number_subproblem_DCA)) 
				+ String.format("%20s", (number_subproblem_MDA)) 
				+ String.format("%10s", ((double) timeGreedyS) / 1000) 
		);
    }
}
