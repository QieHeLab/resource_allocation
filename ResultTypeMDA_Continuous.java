/******************************************************************************
* This is a class for storing the solution to the subproblems of RAP-NC in the MDA method
****************************************************************************/

public class ResultTypeMDA_Continuous {
	double[] aa;
	double[] ab;
	double[] ba;
	double[] bb;

	public ResultTypeMDA_Continuous(double[] aa, double[] ab, double[] ba, double[] bb) {
		this.aa = aa;
		this.ab = ab;
		this.ba = ba;
		this.bb = bb;
	}
} 