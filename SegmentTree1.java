/******************************************************************************
 *  Compilation:  javac SegmentTree1.java
 *  Execution:    java SegmentTree1
 *  
 *  A segment tree data structure written by Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/

import java.util.Arrays;
import java.util.*;

/**
 * The {@code SegmentTree} class is an structure for efficient search of cumulative data.
 * It performs  Range Minimum Query O(log(n)) time.
 * <p>
 * Also it has been develop with  {@code LazyPropagation} for range updates, which means
 * when you perform update operations over a range, the update process affects the least nodes as possible
 * so that the bigger the range you want to update the less time it consumes to update it. Eventually those changes will be propagated
 * to the children and the whole array will be up to date.
 * <p>
 * Example:
 * <p>
 * SegmentTreeHeap st = new SegmentTreeHeap(new Integer[]{1,3,4,2,1, -2, 4});
 * st.update(0,3, 1)
 * In the above case only the node that represents the range [0,3] will be updated (and not their children) so in this case
 * the update task will be less than n*log(n)
 *
 * Memory usage:  O(n)
 *
 * @author Zeyang Wu
 */


public class SegmentTree1 {
	public class DataType {
		//An internal class with two instance variables
		long value;
		int id;

		public DataType(int id, long value) {
			this.id = id;
			this.value = value;
		}
	}

	private class Node {
		long min;
		int minIndex;

		Long pendingVal = null;

		int from;
		int to;

		int size() {
			return to - from + 1;
		}
	}


	private Node[] heap;
	private DataType[] array;
    private int size;

    /**
     * Time-Complexity:  O(n*log(n))
     *
     * @param array the Initialization array
     */
    public SegmentTree1(long[] array) {
        //this.array = Arrays.copyOf(array, array.length);
        this.array = new DataType[array.length];
        for (int i = 0; i < array.length; i++) {
            this.array[i] = new DataType(i, array[i]);
        }
        
        //The max size of this array is about 2 * 2 ^ log2(n) + 1
        size = (int) (2 * Math.pow(2.0, Math.floor((Math.log((long) array.length) / Math.log(2.0)) + 1)));
        heap = new Node[size];
        build(1, 0, array.length);
    }

    public int size() {
        return array.length;
    }

    //Initialize the Nodes of the Segment tree
    private void build(int v, int from, int size) {
        heap[v] = new Node();
        heap[v].from = from;
        heap[v].to = from + size - 1;
        
        if (size == 1) {
            heap[v].min = array[from].value;
            heap[v].minIndex = array[from].id;
        } else {
            //Build childs
            build(2 * v, from, size / 2);
            build(2 * v + 1, from + size / 2, size - size / 2);

            //min = min of the children
            //heap[v].min = Math.min(heap[2 * v].min, heap[2 * v + 1].min);
            if (heap[2 * v].min > heap[2 * v + 1].min) {
                heap[v].min = heap[2 * v + 1].min;
                heap[v].minIndex = heap[2 * v + 1].minIndex;
            } else {
                heap[v].min = heap[2 * v].min;
                heap[v].minIndex = heap[2 * v].minIndex;
            }
        }
    }

	/**
     * Range Min Query
     * 
     * Time-Complexity: O(log(n))
     *
     * @param  from from index
     * @param  to to index
     * @return min
     */
    public DataType rMinQ(int from, int to) {
        return rMinQ(1, from, to);
    }

    private DataType rMinQ(int v, int from, int to) {
        Node n = heap[v];        

        if (contains(from, to, n.from, n.to)) {
            return new DataType(heap[v].minIndex, heap[v].min);
        }

        //If you did a range update that contained this node, you can infer the Min value without going down the tree
        /*
        if (n.pendingVal != null && contains(n.from, n.to, from, to)) {
            return n.pendingVal;
        }
        */

        if (intersects(from, to, n.from, n.to)) {
            propagate(v);
            DataType leftMin = rMinQ(2 * v, from, to);
            DataType rightMin = rMinQ(2 * v + 1, from, to);

            //return Math.min(leftMin, rightMin);
            if (leftMin.value > rightMin.value) {
                return rightMin;
            } else {
                return leftMin;
            }
        }

        return new DataType(Integer.MAX_VALUE, Integer.MAX_VALUE);
    }

    //Test if the range1 contains range2
    private boolean contains(int from1, int to1, int from2, int to2) {
        return from2 >= from1 && to2 <= to1;
    }

    //check inclusive intersection, test if range1[from1, to1] intersects range2[from2, to2]
    private boolean intersects(int from1, int to1, int from2, int to2) {
        return from1 <= from2 && to1 >= from2   //  (.[..)..] or (.[...]..)
                || from1 >= from2 && from1 <= to2; // [.(..]..) or [..(..)..]
    }


     /**
     * Range Add Operation.
     * With this operation you can add either one position or a range of positions with a given number.
     * The rangeAdd operations will add the less it can to add the whole range (Lazy Propagation).
     * The values will be propagated lazily from top to bottom of the segment tree.
     * This behavior is really useful for adds on portions of the array
     * <p>
     * Time-Complexity: O(log(n))
     *
     * @param from  from index
     * @param to    to index
     * @param value value
     */
    public void rangeAdd(int from, int to, long value) {
        rangeAdd(1, from, to, value);
    }

    private void rangeAdd(int v, int from, int to, long value) {

        //The Node of the heap tree represents a range of the array with bounds: [n.from, n.to]
        Node n = heap[v];

        /**
         * If the updating-range contains the portion of the current Node  We lazily update it.
         * This means We do NOT update each position of the vector, but update only some temporal
         * values into the Node; such values into the Node will be propagated down to its children only when they need to.
         */
        if (contains(from, to, n.from, n.to)) {
            change(n, value);
            return;
        }

        if (n.size() == 1) return;

        if (intersects(from, to, n.from, n.to)) {
            /**
             * Before keeping going down to the tree We need to propagate the
             * the values that have been temporally/lazily saved into this Node to its children
             * So that when We visit them the values  are properly updated
             */
            propagate(v);

            rangeAdd(2 * v, from, to, value);
            rangeAdd(2 * v + 1, from, to, value);

            //heap[v].min = Math.min(heap[2 * v].min, heap[2 * v + 1].min);
            if (heap[2 * v].min > heap[2 * v + 1].min) {
                n.min = heap[2 * v + 1].min;
                n.minIndex = heap[2 * v + 1].minIndex;
            } else {
                n.min = heap[2 * v].min;
                n.minIndex = heap[2 * v].minIndex;
            }
        }
    }

    //Propagate temporal values to children
    private void propagate(int v) {
        Node n = heap[v];

        if (n.pendingVal != null) {
            change(heap[2 * v], n.pendingVal);
            change(heap[2 * v + 1], n.pendingVal);
            n.pendingVal = null; //unset the pending propagation value
        }
    }

    //Save the temporal values that will be propagated lazily
    private void change(Node n, long value) {
    	if (n.pendingVal == null) {
    		n.pendingVal = Long.valueOf(0);
    	}
        n.pendingVal += value;
        n.min += value;
        //array[n.from].value += value;
    }

    
    public static void main(String[] args) {
    	SegmentTree1 st = new SegmentTree1(new long[]{1,3,4,2,1, -2, 4,8});

    	//test whether the tree is built correctly
    	
    	for (Node i : st.heap) {
    		if (i != null) { 
    			System.out.println("range: " + i.from + " to " + i.to + " with min " + i.min + " at position " + i.minIndex);
    	 	}
    	}
    	   	

    	//test whether the minQuery is correct
    	//System.out.println("rangeQuery: 1 to 4 with min " + st.rMinQ(1, 4).value + " at position " + st.rMinQ(1, 4).id);
    	//System.out.println("rangeQuery: 13 to 9 with min " + st.rMinQ(13, 9).value + " at position " + st.rMinQ(13, 9).id);
    	//for (DataType i : st.array) {
    		//System.out.println(i.id + " " + i.value);
    	//}

    	//test rangeAdd

    	System.out.println("rangeQuery: 1 to 4 with min " + st.rMinQ(1, 4).value + " at position " + st.rMinQ(1, 4).id);
    	
    	/*
    	for (Node i : st.heap) {
    		if (i != null && i.size() == 1) { 
    			System.out.println("range: " + i.from + " to " + i.to + " with min " + i.min + " at position " + i.minIndex);
    	 	}
    	}
    	*/
    	st.rangeAdd(0, 3, 2);
    	//st.rangeAdd(0, 0, 4);
    	//st.rangeAdd(4, 8, -3);
    	//st.rangeAdd(5, 5, 22);
    	//st.rangeAdd(7, 7, -78);
    	System.out.println("After rangeAdd, rangeQuery: 0 to 7 with min " + st.rMinQ(0, 7).value + " at position " + st.rMinQ(0, 7).id);
    	System.out.println("After rangeAdd, rangeQuery: 4 to 6 with min " + st.rMinQ(4, 6).value + " at position " + st.rMinQ(4, 6).id);
    	System.out.println("After rangeAdd, rangeQuery: 2 to 2 with min " + st.rMinQ(2, 2).value + " at position " + st.rMinQ(2, 2).id);
    	/*
    	st.rangeAdd(1, 4, 3);
    	System.out.println("After rangeAdd, rangeQuery: 1 to 4 with min " + st.rMinQ(1, 4).value + " at position " + st.rMinQ(1, 4).id);
    	st.rangeAdd(4, 4, 67);
    	System.out.println("After rangeAdd, rangeQuery: 1 to 4 with min " + st.rMinQ(1, 4).value + " at position " + st.rMinQ(1, 4).id);
    	st.rangeAdd(4, 5, -67);
    	System.out.println("After rangeAdd, rangeQuery: 1 to 4 with min " + st.rMinQ(1, 4).value + " at position " + st.rMinQ(1, 4).id);
    	System.out.println("After rangeAdd, rangeQuery: 1 to 4 with min " + st.rMinQ(1, 5).value + " at position " + st.rMinQ(1, 5).id);
    	*/
    	
    }
    
}