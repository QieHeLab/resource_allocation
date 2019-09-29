/******************************************************************************
 *  Compilation:  javac CapacitatedLotSizing.java
 *  Execution:    java CapacitatedLotSizing
 *  
 *  A class for CapacitatedLotSizing problem written by Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/

/******************************************************************************
* This a class for capacitated lot-sizing problem with convex production cost and no holding cost. 
* Backorder is not allowed.      
* Method: solve()
*         ssp(int, long[])
* A scaling algorithm inspired by Hochbaum 2008 paper is implemented. The worse cases running time is 
* O(nlognlonB) where n is the number of periods and B is total amount of demands.  
* @author Zeyang Wu
*/

import java.util.*;

public class CapacitatedLotSizing {
	
	//instance variables
	//number of periods: n
	int numPeriod;

	//demands
	long[] demands;
    
    //production cost for period 0, 1, ..., n - 1, the length is n 
    List<Function> obj;

    //no holding cost

    //capacity
    //production capacity
    long[] prodCap;

    //inventory capacity: 0 -> 1 -> ... -> n - 1, the length is n - 1
    long[] inventCap;


    //Constructor
	public CapacitatedLotSizing() {
		this.numPeriod = 0; 
		this.demands = null;
		this.obj = null;
		this.prodCap = null;
		this.inventCap = null;
	}	
    
    //The DataType class represents a path with its id and its cost/capacity
    private class DataType{
    	int id;
    	double value;
    	public DataType(int id, double value) {
    		this.id = id;
    		this.value = value;
    	}

    	/*
    	@Override
    	public int compareTo(Object o) {
    		DataType a = (DataType) o;
    		return new Integer(this.id).compareTo(a.id);
    	}

    	@Override
    	public boolean equals(Object obj) {
    		if (obj == this) {
    			return true;
    		}

    		if (obj instanceof DataType) {
    			DataType a = (DataType) obj;
    			if (this.id == a.id) {
    				return true;
    			}
    		}

    		return false;
    	}
    	*/
    }

    private class DataTypeComparator implements Comparator<DataType> {
    	//compare the cost
    	public int compare(DataType a, DataType b) {    		
    		if (a.id == b.id) {
    			return 0;
    		}    		
    		if (a.value > b.value || (a.value == b.value && a.id < b.id)) {
    			return -1;
    		}
    		return 1;
    	}
    }

    /**
     * Successive Shortest Path Algorithm with step size s.
     * The is a special implementation for solving minimum network flow problem with convex arc cost.
     * Balanced Binary Search tree is used to store the possible paths and supports O(log n) add/remove/min operations. 
     * Segment tree is used to maintain the capacities and supports O(log n) RangeMinQuery and RangeAdd operations.
     * <p>
     * Time-Complexity: O(n + m log(n)) where m is the total number of increments made during this process 
     *
     * @param s  step size
     * @param pardFlow a initial pseudo flow, can be 0 or any component-wise lower of the optimal solution  
     */


	public void ssp(long s, long[] prodFlow) {

		//Make copies of the parameters 
		double[] cost = new double[numPeriod];
		long[] u0 = Arrays.copyOfRange(prodCap, 0, numPeriod);
		long[] u = Arrays.copyOfRange(inventCap, 0, numPeriod - 1);
		long[] de = Arrays.copyOfRange(demands, 0, numPeriod);
		
		//preprocess the prodFlow to update the capacities O(n)
		long  sum = 0; 
		for (int i = 0; i < numPeriod; i++) {
			cost[i] = obj.get(i).getValue(prodFlow[i] + 1) - obj.get(i).getValue(prodFlow[i]);
			u0[i] -= prodFlow[i];
			if (sum + prodFlow[i] <= de[i]) {				
				de[i] -= (prodFlow[i] + sum);
				sum = 0;
			} else {
				sum += (prodFlow[i] - de[i]);
				de[i] = 0;				
			}

			if (i < numPeriod - 1 && sum > 0) {
				u[i] -= sum;
			}
		}
		

		//main iteration
		//A forward tree store the cost of possible paths
		List<DataType> hashCost = new ArrayList<>();

		for (int i = 0; i < numPeriod; i++) {
			hashCost.add(new DataType(i, cost[i]));
		}

		//curt iteration
		int curt = 0;
		//smallest index of production period
		int a = 0;
		TreeSet<DataType> forwardTree = new TreeSet<>(new DataTypeComparator());
		forwardTree.add(hashCost.get(0));

		//segment tree data structure used to store the capacities of horizontal arcs
		SegmentTree1 bo = new SegmentTree1(u);
		

		while(curt < numPeriod) {
			//satisfies the demand at node curt sequentially
			while (de[curt] > 0) {		

				//For debugging use. If the problem is feasible, this will not pop up
				if (forwardTree.isEmpty()) {					
					System.out.println(forwardTree.isEmpty());
					System.out.println(Arrays.toString(prodFlow));
					System.out.println("Remaining Demand ");
					System.out.println(Arrays.toString(de));
					System.out.println(curt);
					break;
				}
				
				//find the minimum unit production cost path 
				int minIndex = forwardTree.pollLast().id;

				//range minimum query to find the bottleneck
				long bottleneck = Integer.MAX_VALUE;
				int bottleneckIndex = -1;	
				
				if (minIndex < curt) {
					bottleneck = bo.rMinQ(minIndex, curt - 1).value;
					bottleneckIndex = bo.rMinQ(minIndex, curt - 1).id;
				}
				
				//step size
				long delta = Math.min(u0[minIndex], s);
				delta = Math.min(delta, bottleneck);
				delta = Math.min(delta, de[curt]);

				//update the capacities
				de[curt] -= delta;
				u0[minIndex] -= delta;
				prodFlow[minIndex] += delta;
				
				//update the cost		
				hashCost.get(minIndex).value = obj.get(minIndex).getValue(prodFlow[minIndex] + 1) - obj.get(minIndex).getValue(prodFlow[minIndex]);			
				forwardTree.add(hashCost.get(minIndex));
				
				//rangeupdate the capacities of the horizontal arcs
				if (minIndex < curt) {
					bo.rangeAdd(minIndex, curt - 1, -delta);
				}				
				
				//remove path if necessary
				if (u0[minIndex] == 0) {
					//System.out.print("Remove arc by prodcap ");
					//System.out.print(minIndex);
					forwardTree.remove(hashCost.get(minIndex));
				}

				if (delta == bottleneck) {
					for (int i = a; i <= bottleneckIndex; i++) {
						forwardTree.remove(hashCost.get(i));
					}
					a = bottleneckIndex + 1;
				}
			}

			curt++;
			if (curt == numPeriod) {
				break;
			}
			
			//Add new path
			forwardTree.add(hashCost.get(curt));
		}
	}

	/**
     * Solve operations.
     * With this operation you can solve the lot sizing problem with convex production cost
     * The is a scaled version of Hochbaum 2008 Solving Linear Cost Dynamic Lot-Sizing Problems in O(n log n) Time
     * <p>
     * Time-Complexity: O(n log(n) log(B)) where B is the total demand
     *
     * @param no param
     */
	public long[] solve() {
		long[] prodFlow = new long[numPeriod];

		//ssp(1, prodFlow);
		
		long B = 0;
		for (int i = 0; i < numPeriod; i++) {
			B += demands[i];
		}
		
		long s = (long) Math.ceil(((double) B) / numPeriod / 2);
		while (s > 1) {
			//greedy step with size s
			ssp(s, prodFlow);			
			
			
 			//System.out.println("Sum of flows: " + java.util.stream.LongStream.of(prodFlow).sum());
			//undo the last step of greedy(s)
			for (int i = 0; i < numPeriod; i++) {
				//avoid infeasibility				
				prodFlow[i] = Math.max(prodFlow[i] - s, 0);
			}
            //System.out.println("Pseudo flow after adjustment with step size " + s + ": " + Arrays.toString(prodFlow));

			if (s % 2 == 0) {
				s = s / 2;
			} else {
				s = s / 2 + 1;
			}
		}
		
		ssp(1, prodFlow);
		
		return prodFlow;
	}

	/**
     * Solve operations for LSP with linear production cost.
     * With this operation you can solve the lot sizing problem with linear production cost
     * The is an algorithm derived in Hochbaum 2008 Solving Linear Cost Dynamic Lot-Sizing Problems in O(n log n) Time
     * <p>
     * Time-Complexity: O(n log(n))
     *
     * @param no param
     */
	public long[] solveLinear() {
		long[] prodFlow = new long[numPeriod];

		//ssp(1, prodFlow);
		
		long B = 0;
		for (int i = 0; i < numPeriod; i++) {
			B += demands[i];
		}
		
		long s = B;
		
		
		ssp(s, prodFlow);
		
		return prodFlow;
	}


	/**
     * A greedy algorithm with step size s.
     * The is a special implementation of Hochbaum's algorithm in her 1994 paper for DRAP with Network constraints
     * <p>
     *
     * @param s  step size
     * @param pardFlow a initial pseudo flow, can be 0 or any component-wise lower of the optimal solution  
     */

	public void greedy_s(long s, long[] prodFlow) {

		//Make copies of the parameters 
		double[] cost = new double[numPeriod];
		long[] u0 = Arrays.copyOfRange(prodCap, 0, numPeriod);
		long[] u = Arrays.copyOfRange(inventCap, 0, numPeriod - 1);
		long[] de = Arrays.copyOfRange(demands, 0, numPeriod);
		
		//preprocess the prodFlow to update the capacities O(n)
		long  sum = 0; 
		long total_demand = 0;
		for (int i = 0; i < numPeriod; i++) {
			cost[i] = obj.get(i).getValue(prodFlow[i] + 1) - obj.get(i).getValue(prodFlow[i]);
			u0[i] -= prodFlow[i];
			if (sum + prodFlow[i] <= de[i]) {				
				de[i] -= (prodFlow[i] + sum);
				sum = 0;
			} else {
				sum += (prodFlow[i] - de[i]);
				de[i] = 0;				
			}

			if (i < numPeriod - 1 && sum > 0) {
				u[i] -= sum;
			}
			total_demand += de[i];
		}
		

		//main iteration
		//A forward tree store the cost of possible paths
		TreeSet<DataType> forwardTree = new TreeSet<>(new DataTypeComparator());
		List<DataType> hashCost = new ArrayList<>();

		for (int i = 0; i < numPeriod; i++) {
			hashCost.add(new DataType(i, cost[i]));
			forwardTree.add(hashCost.get(i));
		}

		while (total_demand > 0) {
			int minIndex = forwardTree.pollLast().id;

			//check how much we could increase in that minIndex
			long maxIncrease = de[numPeriod - 1];
			for (int i = numPeriod - 2; i >= minIndex; i--) {
				maxIncrease = Math.min(maxIncrease, u[i]);
				maxIncrease += de[i];
			}

			maxIncrease = Math.min(maxIncrease, u0[minIndex]);

			if (maxIncrease >= s) {
				maxIncrease = s;
			} 

			
			if (maxIncrease >= 1) {	
				//update the prodFlow
				prodFlow[minIndex] += maxIncrease;
				total_demand -= maxIncrease;
				u0[minIndex] -= maxIncrease;

				//update the cost if it is feasible
				hashCost.get(minIndex).value = obj.get(minIndex).getValue(prodFlow[minIndex] + 1) - obj.get(minIndex).getValue(prodFlow[minIndex]);			
				forwardTree.add(hashCost.get(minIndex));

				//update the capacities
				long temp_sum = maxIncrease;
				for (int i = minIndex; i < numPeriod; i++) {
					if (de[i] > maxIncrease) {
						de[i] -= maxIncrease;
						break;
					} else {
						maxIncrease -= de[i];
						de[i] = 0;
						if (i < numPeriod - 1) {
							u[i] -= maxIncrease;
						}
					}
				}
			}
		}
	}

	public long[] greedy_s_solve() {
		long[] prodFlow = new long[numPeriod];

		//ssp(1, prodFlow);
		
		long B = 0;
		for (int i = 0; i < numPeriod; i++) {
			B += demands[i];
		}
		
		long s = (long) Math.ceil(((double) B) / numPeriod / 2);
		while (s > 1) {
			//greedy step with size s
			greedy_s(s, prodFlow);			
			
			
 			//System.out.println("Sum of flows: " + java.util.stream.LongStream.of(prodFlow).sum());
			//undo the last step of greedy(s)
			for (int i = 0; i < numPeriod; i++) {
				//avoid infeasibility				
				prodFlow[i] = Math.max(prodFlow[i] - s, 0);
			}
            //System.out.println("Pseudo flow after adjustment with step size " + s + ": " + Arrays.toString(prodFlow));

			if (s % 2 == 0) {
				s = s / 2;
			} else {
				s = s / 2 + 1;
			}
		}
		
		greedy_s(1, prodFlow);
		
		return prodFlow;
	}
}