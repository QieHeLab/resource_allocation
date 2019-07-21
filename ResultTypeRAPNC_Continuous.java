/******************************************************************************
* This is a class for storing the solution to an instance continuous RAP-NC  
****************************************************************************/
public class ResultTypeRAPNC_Continuous {
	boolean feasible;
	double[] sol;
	public ResultTypeRAPNC_Continuous(boolean feasible, double[] sol) {
		this.feasible = feasible;
		this.sol = sol;
	}
}