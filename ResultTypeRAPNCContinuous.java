/******************************************************************************
* This a class for storing the solution to RAP-NC with continuous variables
****************************************************************************/
public class ResultTypeRAPNCContinuous {
	boolean feasible;
	double[] sol;
	public ResultTypeRAPNCContinuous(boolean feasible, double[] sol) {
		this.feasible = feasible;
		this.sol = sol;
	}
}
