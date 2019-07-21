/******************************************************************************
* This is a class for storing the solution to an instance continuous RAP
****************************************************************************/
public class ResultTypeRAP_Continuous {
	boolean feasible;
	double[] sol;
	public ResultTypeRAP_Continuous(boolean feasible, double[] sol) {
		this.feasible = feasible;
		this.sol = sol;
	}
}