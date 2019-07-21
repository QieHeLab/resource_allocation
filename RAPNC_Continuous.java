/*******************************************************************************
 *  Compilation:  javac RAPNC_continuous.java
 *  Execution:    java RAPNC_continuous
 *  
 *  A class for resource allocation problem with lower and upper nested constraints (RAPNC) 
 *  @author Zeyang Wu @ University of Minnesota
 *
 * Now this class can only solve RAP-NC with quadratic objective as the method 
 * RAP_Continuous.solve_RAP_Continuous()
 * can only handle the quadratic case 
 *
 * It can be easily generalized to general convex function as long as the derivatives can be query.  
 *
 *******************************************************************************/


import java.util.*;

public class RAPNC_Continuous {
	//problem dimension
	int dimension;

    //list of function oracles 
	List<Function> obj;

	//list of parameters
	double[] ubVar;
	double[] lbVar;
	double[] ubNested;
	double[] lbNested;
    double[] lbCopyMDA; //reserve bounds for MDA
    double[] ubCopyMDA;
	private double scaleFactor; //deprecated
    int number_subproblem; // number of RAP subproblem solved

    //normal constructor
    public RAPNC_Continuous(int K) {
    	//normal constructor
    	this.obj = new ArrayList<Function>();
		this.lbVar = new double[K];
		this.ubVar =  new double[K];
		this.lbNested =  new double[K];
		this.ubNested =  new double[K];
		this.dimension = lbVar.length;
		this.scaleFactor = 1;
        this.number_subproblem = 0;
    }

    //normal constructor
	public RAPNC_Continuous(List<Function> obj, double[] lbVar, double[] ubVar, double[] lbNested, double[] ubNested) {
		this.obj = obj;
		this.lbVar = lbVar;
		this.ubVar = ubVar;
		this.lbNested = lbNested;
		this.ubNested = ubNested;
		this.dimension = lbVar.length;
		this.scaleFactor = 1;		
        this.number_subproblem = 0;
	}

	//This is used in solving continuous RAP-NC    
    public void setScaleFactor(double scaleFactor) {
    	this.scaleFactor = scaleFactor;
    }

     /**
     * createRAP()
     * By this method we can create an instance of RAP by relaxing the nested constraints.
     * Time-Complexity: O(n)     *
     * @param no argument
     * @return RAP relaxation of RAP-NC
     */
    public RAP_Continuous createRAP() {
    	//The resource bound is the bound of the last nested constraint.
    	double rapB = ubNested[dimension - 1];
    	double[] raplb = lbVar.clone();
    	double[] rapub = ubVar.clone();
    	RAP_Continuous res = new RAP_Continuous(obj, rapB, raplb, rapub);
    	return res;
    }

    /**
     * createRAPNC
     * By this method we can create two RAPNC subproblems by a index K.
     * RAP-NC(s, K) and RAP-NC(K + 1, e) 
     * Time-Complexity: O(n) 
     * 
     * @return List<RAPNC> containing the two subproblems
     */
    public List<RAPNC_Continuous> createRAPNC(int K) {
    	//The first subproblem contains K + 1 variables and the second subproblem contains (n - K - 1) variables
    	//return a list containing the two problems

    	//copy array
    	RAPNC_Continuous left = new RAPNC_Continuous(K + 1);
    	RAPNC_Continuous right = new RAPNC_Continuous(dimension - K - 1);
    	//setup left obj
    	System.arraycopy(this.ubNested, 0, left.ubNested, 0, K + 1);
    	System.arraycopy(this.lbNested, 0, left.lbNested, 0, K + 1);
    	System.arraycopy(this.ubVar, 0, left.ubVar, 0, K + 1);
    	System.arraycopy(this.lbVar, 0, left.lbVar, 0, K + 1);
    	left.dimension = K + 1;
    	left.obj = this.obj;
    	left.scaleFactor = this.scaleFactor;


    	//setup right
    	right.dimension = dimension - K - 1;
    	System.arraycopy(this.ubVar, K + 1, right.ubVar, 0, dimension - K - 1);
    	System.arraycopy(this.lbVar, K + 1, right.lbVar, 0, dimension - K - 1);
    	right.obj = new ArrayList<Function>();
    	for (int i = K + 1; i < dimension; i++) {
    		right.ubNested[i - K - 1] = ubNested[i] - ubNested[K];
    		right.lbNested[i - K - 1] = lbNested[i] - ubNested[K];
    		right.obj.add(this.obj.get(i));
    	}
    	right.scaleFactor = this.scaleFactor;

    	List<RAPNC_Continuous> res = new ArrayList<RAPNC_Continuous>();
    	res.add(left);
    	res.add(right);
    	return res;
    }


    /**
     * solveContinuousDCA() 
     * This is the implementation of DCA method
     * Time-Complexity: O(n^2 log(B)) 
     * 
     * @return ResultTypeRAPNC_Continuous containing the solution and feasibility
     */
    public ResultTypeRAPNC_Continuous solveContinuousDCA() {
    	//Solve the RAPNC with integer variables
    	//1. solve the relaxation problem
    	//2. find the maximum violation and then divide the problem into two subproblems
    	//3. solve the two subproblem recursively and conquer the results

        this.number_subproblem++; //record the subproblems
        //double accuracy = 1e-5;
    	//Trivial case
    	if (dimension == 1) {
    		//check feasibility
    		if (ubNested[0] >= lbVar[0] - 1e-5 && ubNested[0] <= ubVar[0] + 1e-5) {
    			double[] sol = new double[1];
    			sol[0] = ubNested[0];
    			return new ResultTypeRAPNC_Continuous(true, sol);
    		} else {
                System.out.println("Wrong at RAP-NC");
    			return new ResultTypeRAPNC_Continuous(false, null);
    		}
    	}


    	//call Greedy Algorithm to solve the relaxed instance with greedy algorithm.
    	RAP_Continuous re = createRAP();
    	ResultTypeRAP_Continuous solRAP = re.solveRAP_Continuous(); 

    	if (!solRAP.feasible) {
    		return new ResultTypeRAPNC_Continuous(false, null);
    	}
    	//obtain a solution to RAP
    	double[] solRe = solRAP.sol.clone();
    	/*Debug
        
		System.out.println("yes, the problem is feasible");
		System.out.println(Arrays.toString(solRe));
		*/

    	//find the maximum violation
    	double sum = 0;
    	int maxIndex = -1;
    	double maxVio = 1e-4;
    	int maxFlag = 0;//excess 1, shortage 0
    	double[] violation = new double[dimension];
    	
    	for (int i = 0; i < dimension - 1; i++) {
    		sum += solRe[i];
    		int flag = 1;
    		if (sum > ubNested[i]) {
    			violation[i] = sum - ubNested[i];
    		} else if (sum < lbNested[i]) {
    			violation[i] = lbNested[i] - sum;
    			flag = 0;
    		}
    		if (violation[i] > maxVio) {
    			maxIndex = i;
    			maxVio = violation[i];
    			maxFlag = flag;
    		}
    	}

        //If the solution to RAP satisfies all nested constraints. 
    	if (maxIndex == -1) {
    		//.out.println(Arrays.toString(solRe));
    		return new ResultTypeRAPNC_Continuous(true, solRe);
    	} 
    	/*Debug
    	System.out.print("Split at:");
        System.out.print(maxIndex);
        System.out.print(" with ");
        System.out.println(maxFlag);
		*/

    	//else divide the problem into two small problems. 
    	if (maxFlag == 0) {
    		ubNested[maxIndex] = lbNested[maxIndex];
    	} else {
    		lbNested[maxIndex] = ubNested[maxIndex];
    	}
    	
    	/*Debug
    	System.out.print("TightenNested:");
        System.out.print(ubNested[maxIndex]);
        System.out.print(" with ");
        System.out.println(ubNested[maxIndex]);
    	*/  	
    	
    	List<RAPNC_Continuous> divide = createRAPNC(maxIndex);
    	//System.out.println(maxIndex);
    	if (maxIndex == dimension - 1) {
    		System.out.println(Arrays.toString(violation));
    		int[] ary = new int[9];
    		ary[9]++;
    	}
    	ResultTypeRAPNC_Continuous left = divide.get(0).solveContinuousDCA(); 
    	ResultTypeRAPNC_Continuous right = divide.get(1).solveContinuousDCA();

        this.number_subproblem += divide.get(0).number_subproblem;
        this.number_subproblem += divide.get(1).number_subproblem;

    	if (!(left.feasible && right.feasible)) {
    		return new ResultTypeRAPNC_Continuous(false, null);
    	} else {
    		System.arraycopy(left.sol, 0, solRe, 0, left.sol.length);
    		System.arraycopy(right.sol, 0, solRe, maxIndex + 1, right.sol.length);
    		return new ResultTypeRAPNC_Continuous(true, solRe);
    	}
    } 



    // New implementation of MDA 
    /**
     * FastMDA() 
     * This is an implementation of MDA method，it initiate the Main recursion FastMDA(u,v)
     * We discuss with vidal at informs annual meeting 2018 to have a better understanding of the detailed implementation.
     * Time-Complexity: O(n log(n) log(B)) 
     * @param obj the function oracles should be transformed to corresponding obj with penalties
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */
    public ResultTypeMDA_Continuous FastMDA() {
        //lbCopyMDA = Arrays.copyOfRange(lbVar, 0, dimension);
        //ubCopyMDA = Arrays.copyOfRange(ubVar, 0, dimension);
        lbCopyMDA = new double[dimension];
        ubCopyMDA = new double[dimension];
        //Arrays.fill(lbCopyMDA, - 2 * ubNested[dimension - 1]);
        //Arrays.fill(ubCopyMDA, 2 * ubNested[dimension - 1]);
        return FastMDA(0, dimension - 1);
    }

    /**
     * FastMDA() 
     * This is the implementation of MDA method
     * Time-Complexity: O(n log(n) log(B)) 
     * @param obj 
     * @param v start index
     * @param u end index
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */

    public ResultTypeMDA_Continuous FastMDA(int v, int w) {
        
        this.number_subproblem += 4; // record the number of subproblems
        //initialize return of the f
        double[] aa = new double[w - v + 1];
        double[] ab = new double[w - v + 1];
        double[] ba = new double[w - v + 1];
        double[] bb = new double[w - v + 1];

        //simple cases
        if (v == w) {
            if (v == 0) {
                aa[0] = lbNested[v];
                ab[0] = ubNested[v];
                ba[0] = lbNested[v];
                bb[0] = ubNested[v];
            } else {
                aa[0] = lbNested[v] - lbNested[v - 1];
                ab[0] = ubNested[v] - lbNested[v - 1];
                ba[0] = lbNested[v] - ubNested[v - 1];
                bb[0] = ubNested[v] - ubNested[v - 1];
            }
            return new ResultTypeMDA_Continuous(aa, ab, ba, bb);            
        }

        //divide 
        int u = v + (w - v) / 2;
        ResultTypeMDA_Continuous left = FastMDA(v, u);
        ResultTypeMDA_Continuous right = FastMDA(u + 1, w);


        //conquer
        //update the bounds
        //aa
        //check the conditions
        boolean flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        boolean flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            aa = subproblemRAPSolveFastMDA(v, w, lbNested[w]);
        } else {
            aa = subproblemRAPSolveFastMDA(v, w, lbNested[w] - lbNested[v - 1]);
        }

        //avoid the last computation
        if (v == 0 && w == dimension - 1) {
            //System.out.println("yes");
            return new ResultTypeMDA_Continuous(aa, aa, aa, aa);
        }

        //ab
        flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            ab = subproblemRAPSolveFastMDA(v, w, ubNested[w]);
        } else {
            ab = subproblemRAPSolveFastMDA(v, w, ubNested[w] - lbNested[v - 1]);
        }

        //ba
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            ba = subproblemRAPSolveFastMDA(v, w, lbNested[w]);
        } else {
            ba = subproblemRAPSolveFastMDA(v, w, lbNested[w] - ubNested[v - 1]);
        }

        //bb
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            bb = subproblemRAPSolveFastMDA(v, w, ubNested[w]);
        } else {
            bb = subproblemRAPSolveFastMDA(v, w, ubNested[w] - ubNested[v - 1]);
        }


        //return 
        return new ResultTypeMDA_Continuous(aa, ab, ba, bb); 
    }



    /**
     * subproblemRAPSolveFastMDA
     * This is an method used to solve the special RAP subproblems in MDA method
     * Time-Complexity: O(n log n) 
     * @param v start index
     * @param u end index
     * @param LR the total amount of resource in the subproblems
     * @return solution to the subproblem. Here we simply assume that the problem is always feasible
     */
    public double[] subproblemRAPSolveFastMDA(int v, int w, double LR) { 
        //set up checkers
        double sum_c_prime = 0;
        double sum_d_prime = 0;
        double sum_c_bar = 0;
        double sum_d_bar = 0;

        double rapB = LR;
        double[] raplb = new double[w + 1 - v];
        double[] rapub = new double[w + 1 - v];
        List<Function> rapObj = obj.subList(v, w + 1);

        for (int i = v; i < w + 1; i++) {

            raplb[i - v] = Math.max(lbCopyMDA[i], Math.min(lbVar[i], ubCopyMDA[i]));
            rapub[i - v] = Math.min(ubCopyMDA[i], Math.max(ubVar[i], lbCopyMDA[i]));
            
            sum_c_bar += lbCopyMDA[i];
            sum_d_bar += ubCopyMDA[i];
            sum_c_prime += raplb[i - v];      
            sum_d_prime += rapub[i - v];
        }

        //case one: the problem is infeasible w.r.t. sum_c
        // This can be viewed as a special implementation of the greedy algorithms for RAP. (Because the penalty cost are all the same and linear.)
        if (sum_c_prime > LR) {
            // return some solution
            double[] construct_sol = Arrays.copyOfRange(lbCopyMDA, v, w + 1);

            for (int i = v; i < w + 1; i++) {
                
                if (LR - sum_c_bar == 0) {
                    break;
                }

                double delta = Math.max(Math.min(raplb[i - v] - lbCopyMDA[i], LR - sum_c_bar), 0);
                construct_sol[i - v] += delta;
                sum_c_bar += delta;
            }

            return construct_sol;
        }

        //case two:
        if (sum_d_prime < LR) {
            // return some solution
            double[] construct_sol = Arrays.copyOfRange(ubCopyMDA, v, w + 1);
            //System.out.println("Case 2 Checked");
            for (int i = v; i < w + 1; i++) {
                
                if (sum_d_bar - LR == 0) {
                    break;
                }
                
                double delta = Math.max(Math.min(ubCopyMDA[i] - rapub[i - v], sum_d_bar - LR), 0);
                construct_sol[i - v] -= delta;
                sum_d_bar -= delta;
            }

            return construct_sol;
        }


        //case three: the RAP problem is feasible under the original bound
        RAP_Continuous rap = new RAP_Continuous(rapObj, rapB, raplb, rapub);
        double sum_c = 0;
        double sum_d = 0;
        for (int i = v; i < w + 1; i++) {
            sum_c += raplb[i - v];
            sum_d += rapub[i - v];
            if (raplb[i - v] > rapub[i - v]) {
                System.out.println("Something may be wrong");
            }
        }
        if (sum_c > rapB) {
            System.out.println("lower bound is wrong");
        } else if (sum_d < rapB) {
            System.out.println("upper bound is wrong");
        }
        //System.out.println("Checked");
        //rap.scaleFactor = this.scaleFactor;
        ResultTypeRAP_Continuous res = rap.solveRAP_Continuous();
        //ResultTypeRAP res = rap.solveRAPLinear();
        if (!res.feasible) {
            System.out.println("Subproblem" + v + " " + w + "Infeasible");
        }
        return res.sol;
    }

    /**
     * check method 
     * This is an method check if the component-wise vector inequality nums1 <= nums2 is violated 
     * Time-Complexity: O(n) 
     * @param nums1 
     * @param nums2 
     * @param v start index
     * @param u end index
     * @return boolean variable if nums1 <= nums2 then FALSE
     */
    public boolean check(double[] nums1, double[] nums2, int v, int u) {
        //check if nums1 <= nums2
        boolean res = false;
        for (int i = v; i <= u; i++) {
            if (nums1[i - v] > nums2[i - v]) {
                res = true;
            }
        }
        return res;
    }

    /**
     * adjust method 
     * This is an method adjust the two vector such that nums1 <= nums2 is satisfied in a way that
     * the modified arrays are still optimal to the corresponding subproblems
     * Time-Complexity: O(n) 
     * @param nums1 
     * @param nums2 
     * @param v start index
     * @param u end index
     * @return no return 
     */
    public void adjust(double[] nums1, double[] nums2, int v, int u) {
        double delta = 0;
        for (int i = v; i <= u; i++) {
             if (nums1[i - v] > nums2[i - v]) {
                delta += nums1[i - v] - nums2[i - v];
                nums1[i - v] = nums2[i - v];
             }
        }

        for (int i = v; i <= u; i++) {
             if (nums1[i - v] < nums2[i - v]) {
                double ins = Math.min(nums2[i - v] - nums1[i - v], delta);
                delta -= ins;
                nums1[i - v] += ins;
             }
        }
    }

    /**
     * update method 
     * This is an method used to update the variable bounds in the MDA method
     * Time-Complexity: O(n) 
     * @param nums1 
     * @param nums2 
     * @param v start index
     * @param u end index
     * @return no
     */
    public void updateBounds(double[] nums1, double[] nums2, int v, int u) {
        for (int i = v; i <= u; i++) {
             lbCopyMDA[i] = nums1[i - v];
             ubCopyMDA[i] = nums2[i - v];
        }
    }

}