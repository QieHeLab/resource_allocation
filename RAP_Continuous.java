/******************************************************************************
 *  Compilation:  javac RAP_Continuous.java
 *  Execution:    java RAP_Continuous
 *  
 *  A class for simple resource allocation problem(RAP) written by Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/

/******************************************************************************
* This is a class for simple resource allocation with separable convex quadratic objectives
*      
* Method: solveRAP_Continous()
*   
* Bisection is used to find the optimal dual variable. 
* 
* @author Zeyang Wu
*/

import java.util.*;
public class RAP_Continuous {
	//list of function oracles: The user should write their own function classes which extends the abstract class Function .
	List<Function> obj;
	double B;
	double[] lbVar;
	double[] ubVar;
	int dimension;

	public RAP_Continuous(List<Function> obj, double B, double[] lbVar, double[] ubVar) {
		this.obj = obj;
		this.B = B;
		this.lbVar = lbVar;
		this.ubVar = ubVar;
		this.dimension = lbVar.length;	
	}


	/**
     * SolveRAP_Continuous operations.
     * With this operation you can solve the simple resource allocation problem with quadratic objectives
     * It use bisection to find the optimal dual variable
     * It can be generalized to general separable convex objectives if the derivative is provided.
     * <p>
     * 
     * @param no param
     */
	ResultTypeRAP_Continuous solveRAP_Continuous() {
		ResultTypeRAP_Continuous res = null;
		double[] a = new double[dimension];
		double[] b = new double[dimension];

		double left = Double.MAX_VALUE/2;
		double right = -Double.MAX_VALUE/2;
		for (int i = 0; i < dimension; i++) {
			a[i] = obj.get(i).getValue(1) + obj.get(i).getValue(-1);
			b[i] = obj.get(i).getValue(1) - obj.get(i).getValue(-1);
			b[i] /= 2;

			right = Math.max(right, - a[i] * lbVar[i] -b[i]);
			left = Math.min(left, - a[i] * ubVar[i] - b[i]);
		}
		
		double sum = 0;
		double accuracy = 1e-9;

		while (Math.abs(right - left) > accuracy || Math.abs(sum - B) > accuracy)  {
			double mid = left + (right - left) / 2;
			sum = 0;
			for (int i = 0; i < dimension; i++) {
				if (mid > - a[i] * lbVar[i] - b[i]) {
					sum += lbVar[i];
				} else if (mid < - a[i] * ubVar[i] - b[i]) {
					sum += ubVar[i];
				} else {
					sum += - (mid + b[i]) / a[i];
				}
			}

			if (sum >= B) {
				left = mid;
			} else {
				right = mid;
			}
			//System.out.print("Current mid: ");
			//System.out.println(mid);
			//System.out.print("Current sum: ");
			//System.out.println(sum);
		}

		if (Math.abs(right - left) > accuracy || Math.abs(sum - B) > accuracy) {
			System.out.println("Wrong at RAP");
			return new ResultTypeRAP_Continuous(false, null);
		}

		double[] sol = new double[dimension];
		for (int i = 0; i < dimension; i++) {
			if (left > - a[i] * lbVar[i] -b[i]) {
				sol[i] = lbVar[i];
			} else if (left < - a[i] * ubVar[i] - b[i]) {
				sol[i] = ubVar[i];
			} else {
				sol[i] = - (left + b[i]) / a[i];
			}
		}

		return new ResultTypeRAP_Continuous(true, sol);
	}
}