import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 * This class implements six different comparison sorts (and an optional seventh
 * sort for extra credit):
 * <ul>
 * <li>selection sort</li>
 * <li>insertion sort</li>
 * <li>merge sort</li>
 * <li>quick sort</li>
 * <li>heap sort</li>
 * <li>selection2 sort</li>
 * <li>(extra credit) insertion2 sort</li>
 * </ul>
 * It also has a method that runs all the sorts on the same input array and
 * prints out statistics.
 */

public class ComparisonSort {

	private static int moves = 0;

	/**
	 * Sorts the given array using the selection sort algorithm. You may use
	 * either the algorithm discussed in the on-line reading or the algorithm
	 * discussed in lecture (which does fewer data moves than the one from the
	 * on-line reading). Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void selectionSort(E[] A) {

		int j, k, minIndex;
		E min;
		int N = A.length;

		for (k = 0; k < N; k++) {
			min = A[k];
			moves++;
			minIndex = k;
			for (j = k + 1; j < N; j++) {
				if (A[j].compareTo(min) < 0) {
					min = A[j];
					minIndex = j;
					moves++;
				}
			}
			A[minIndex] = A[k];
			moves++;
			A[k] = min;
			moves++;
		}

	}

	/**
	 * Sorts the given array using the insertion sort algorithm. Note: after
	 * this method finishes the array is in sorted order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void insertionSort(E[] A) {
		int k, j;
		E tmp;
		int N = A.length;

		for (k = 1; k < N; k++) {
			tmp = A[k];
			moves++;
			j = k - 1;
			while ((j >= 0) && (A[j].compareTo(tmp) > 0)) {
				A[j + 1] = A[j];
				moves++;
				j--;
			}
			A[j + 1] = tmp;
			moves++;
		}
	}

	/**
	 * Sorts the given array using the merge sort algorithm. Note: after this
	 * method finishes the array is in sorted order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */

	public static <E extends Comparable<E>> void mergeSort(E[] A) {
		mergeAux(A, 0, A.length - 1); // call the aux. function to do all the
										// work
	}

	private static <E extends Comparable<E>> void mergeAux(E[] A, int low,
			int high) {
		// base case
		if (low == high)
			return;

		// recursive case

		// Step 1: Find the middle of the array (conceptually, divide it in
		// half)
		int mid = (low + high) / 2;

		// Steps 2 and 3: Sort the 2 halves of A
		mergeAux(A, low, mid);
		mergeAux(A, mid + 1, high);

		// Step 4: Merge sorted halves into an auxiliary array
		E[] tmp = (E[]) (new Comparable[high - low + 1]);
		moves = moves + high - low + 1;
		int left = low; // index into left half
		int right = mid + 1; // index into right half
		int pos = 0; // index into tmp

		while ((left <= mid) && (right <= high)) {
			// choose the smaller of the two values "pointed to" by left, right
			// copy that value into tmp[pos]
			// increment either left or right as appropriate
			// increment pos
			if (A[left].compareTo(A[right]) <= 0) {
				tmp[pos] = A[left];
				left++;
				moves++;
			} else {
				tmp[pos] = A[right];
				right++;
				moves++;
			}
			pos++;
		}

		// when one of the two sorted halves has "run out" of values, but
		// there are still some in the other half, copy all the remaining
		// values to tmp
		// Note: only 1 of the next 2 loops will actually execute
		while (left <= mid) {
			tmp[pos] = A[left];
			left++;
			pos++;
			moves++;
		}
		while (right <= high) {
			tmp[pos] = A[right];
			right++;
			pos++;
			moves++;
		}

		// all values are in tmp; copy them back into A
		System.arraycopy(tmp, 0, A, low, tmp.length);
	}

	/**
	 * Sorts the given array using the quick sort algorithm, using the median of
	 * the first, last, and middle values in each segment of the array as the
	 * pivot value. Note: after this method finishes the array is in sorted
	 * order.
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * 
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void quickSort(E[] A) {
		quickAux(A, 0, A.length - 1);
	}

	private static <E extends Comparable<E>> void quickAux(E[] A, int low,
			int high) {
		if (high - low < 2) {
			if (high - low == 1) {
				if (A[low].compareTo(A[high]) > 0)
					swap(A, low, high);
			}

		} else {
			int right = partition(A, low, high);
			quickAux(A, low, right);
			quickAux(A, right + 2, high);
		}
	}

	private static <E extends Comparable<E>> int partition(E[] A, int low,
			int high) {
		// precondition: A.length > 3

		E pivot = medianOfThree(A, low, high); // this does step 1

		int left = low + 1;
		int right = high - 2;
		while (left <= right) {
			while (A[left].compareTo(pivot) < 0)

				left++;
			while (A[right].compareTo(pivot) > 0)
				right--;
			if (left <= right) {
				swap(A, left, right);
				left++;
				right--;
			}
		}

		swap(A, right + 1, high - 1); // step 4
		return right;
	}

	private static <E extends Comparable<E>> void swap(E[] A, int low, int high) {

		E temp = A[low];
		A[low] = A[high];
		A[high] = temp;
		moves = moves + 3;
	}

	private static <E extends Comparable<E>> E medianOfThree(E[] A, int low,
			int high) {
		int mid = (low + high) / 2;
		if (A[low].compareTo(A[mid]) > 0)
			swap(A, low, mid);
		if (A[low].compareTo(A[high]) > 0)
			swap(A, low, high);
		if (A[mid].compareTo(A[high]) > 0)
			swap(A, mid, high);
		swap(A, mid, high - 1);

		return A[high - 1];
	}

	/**
	 * Sorts the given array using the heap sort algorithm outlined below. Note:
	 * after this method finishes the array is in sorted order.
	 * <p>
	 * The heap sort algorithm is:
	 * </p>
	 * 
	 * <pre>
	 * for each i from 1 to the end of the array
	 *     insert A[i] into the heap (contained in A[0]...A[i-1])
	 *     
	 * for each i from the end of the array up to 1
	 *     remove the max element from the heap and put it in A[i]
	 * </pre>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void heapSort(E[] A) {
		if (A == null || A.length <= 1) {
			return;
		}

		buildMaxHeap(A);

		for (int i = A.length - 1; i >= 1; i--) {
			swap(A, 0, i);

			maxHeap(A, i, 0);
		}
	}

	private static <E extends Comparable<E>> void buildMaxHeap(E[] A) {
		if (A == null || A.length <= 1) {
			return;
		}

		int half = A.length / 2;
		for (int i = half; i >= 0; i--) {
			maxHeap(A, A.length, i);
		}
	}

	private static <E extends Comparable<E>> void maxHeap(E[] A, int heapSize,
			int index) {
		int left = index * 2 + 1;
		int right = index * 2 + 2;

		int largest = index;
		if (left < heapSize && A[left].compareTo(A[index]) > 0) {
			largest = left;
		}

		if (right < heapSize && A[right].compareTo(A[largest]) > 0) {
			largest = right;
		}

		if (index != largest) {
			swap(A, index, largest);

			maxHeap(A, heapSize, largest);
		}
	}


	/**
	 * Sorts the given array using the selection2 sort algorithm outlined below.
	 * Note: after this method finishes the array is in sorted order.
	 * <p>
	 * The selection2 sort is a bi-directional selection sort that sorts the
	 * array from the two ends towards the center. The selection2 sort algorithm
	 * is:
	 * </p>
	 * 
	 * <pre>
	 * begin = 0, end = A.length-1
	 * 
	 * // At the beginning of every iteration of this loop, we know that the 
	 * // elements in A are in their final sorted positions from A[0] to A[begin-1]
	 * // and from A[end+1] to the end of A.  That means that A[begin] to A[end] are
	 * // still to be sorted.
	 * do
	 *     use the MinMax algorithm (described below) to find the minimum and maximum 
	 *     values between A[begin] and A[end]
	 *     
	 *     swap the maximum value and A[end]
	 *     swap the minimum value and A[begin]
	 *     
	 *     ++begin, --end
	 * until the middle of the array is reached
	 * </pre>
	 * <p>
	 * The MinMax algorithm allows you to find the minimum and maximum of N
	 * elements in 3N/2 comparisons (instead of 2N comparisons). The way to do
	 * this is to keep the current min and max; then
	 * </p>
	 * <ul>
	 * <li>take two more elements and compare them against each other</li>
	 * <li>compare the current max and the larger of the two elements</li>
	 * <li>compare the current min and the smaller of the two elements</li>
	 * </ul>
	 * 
	 * @param <E>
	 *            the type of values to be sorted
	 * @param A
	 *            the array to sort
	 */
	public static <E extends Comparable<E>> void selection2Sort(E[] A) {
		int begin = 0;
		int end = A.length - 1;

		while (begin <= end) {
			int min = begin;
			int max = end;

			int left = begin;
			int right = end;
			while (left <= right) {
				if (A[left].compareTo(A[right]) < 0) {
					if (A[left].compareTo(A[min]) < 0) {
						min = left;
					}
					if (A[right].compareTo(A[max]) > 0) {
						max = right;
					}
				} else {
					if (A[left].compareTo(A[max]) > 0) {
						max = left;
					}
					if (A[right].compareTo(A[min]) < 0) {
						min = right;
					}
				}
				left++;
				right--;
			}
			
			if ((begin == max) && (end == min)) {
				swap(A, begin, end);
			} else {
				if (min != begin) {
					swap(A, min, begin);
				}
				if (max != end) {
					swap(A, max, end);
				}
			}

			begin++;
			end--;
		}

	}


	/**
	 * <b>Extra Credit:</b> Sorts the given array using the insertion2 sort
	 * algorithm outlined below. Note: after this method finishes the array is
	 * in sorted order.
	 * <p>
	 * The insertion2 sort is a bi-directional insertion sort that sorts the
	 * array from the center out towards the ends. The insertion2 sort algorithm
	 * is:
	 * </p>
	 * 
	 * <pre>
	 * precondition: A has an even length
	 * left = element immediately to the left of the center of A
	 * right = element immediately to the right of the center of A
	 * if A[left] > A[right]
	 *     swap A[left] and A[right]
	 * left--, right++ 
	 *  
	 * // At the beginning of every iteration of this loop, we know that the elements
	 * // in A from A[left+1] to A[right-1] are in relative sorted order.
	 * do
	 *     if (A[left] > A[right])
	 *         swap A[left] and A[right]
	 *  
	 *     starting with with A[right] and moving to the left, use insertion sort 
	 *     algorithm to insert the element at A[right] into the correct location 
	 *     between A[left+1] and A[right-1]
	 *     
	 *     starting with A[left] and moving to the right, use the insertion sort 
	 *     algorithm to insert the element at A[left] into the correct location 
	 *     between A[left+1] and A[right-1]
	 *  
	 *     left--, right++
	 * until left has gone off the left edge of A and right has gone off the right 
	 *       edge of A
	 * </pre>
	 * <p>
	 * This sorting algorithm described above only works on arrays of even
	 * length. If the array passed in as a parameter is not even, the method
	 * throws an IllegalArgumentException
	 * </p>
	 * 
	 * @param A
	 *            the array to sort
	 * @throws IllegalArgumentException
	 *             if the length or A is not even
	 */
	public static <E extends Comparable<E>> void insertion2Sort(E[] A) {
		if (A.length % 2 != 0) {
			throw new IllegalArgumentException();
		}
		int N = A.length;
		int left = N / 2 - 1;
		int right = N / 2;
		E tmp;

		if (A[left].compareTo(A[right]) > 0) {
			swap(A, left, right);
			moves++;
		}
		left--;
		right++;
		while (right < N && left >= 0) {
			if (A[left].compareTo(A[right]) > 0) {
				swap(A, left, right);
			}
			moves++;
			tmp = A[right];
			int ind = right - 1;
			while ((ind >= left + 1) && (A[ind].compareTo(tmp) > 0)) {
				moves++;
				A[ind + 1] = A[ind];
				ind--;
			}
			moves++;
			A[ind + 1] = tmp;
			moves++;
			tmp = A[left];
			int ind2 = left + 1;
			while ((ind2 < right) && (A[ind2].compareTo(tmp) < 0)) {
				moves++;
				A[ind2 - 1] = A[ind2];
				ind2++;
			}
			moves++;
			A[ind2 - 1] = tmp;
			left--;
			right++;
		}
	}

	private static <E extends Comparable<E>> int getmoves() {
		return moves;
	}

	private static <E extends Comparable<E>> void resetMoves() {
		moves = 0;
	}

	/**
	 * Internal helper for printing rows of the output table.
	 * 
	 * @param sort
	 *            name of the sorting algorithm
	 * @param compares
	 *            number of comparisons performed during sort
	 * @param moves
	 *            number of data moves performed during sort
	 * @param milliseconds
	 *            time taken to sort, in milliseconds
	 */
	private static void printStatistics(String sort, int compares, int moves,
			long milliseconds) {
		System.out.format("%-23s%,15d%,15d%,15d\n", sort, compares, moves,
				milliseconds);
	}

	/**
	 * Sorts the given array using the six (seven with the extra credit)
	 * different sorting algorithms and prints out statistics. The sorts
	 * performed are:
	 * <ul>
	 * <li>selection sort</li>
	 * <li>insertion sort</li>
	 * <li>merge sort</li>
	 * <li>quick sort</li>
	 * <li>heap sort</li>
	 * <li>selection2 sort</li>
	 * <li>(extra credit) insertion2 sort</li>
	 * </ul>
	 * <p>
	 * The statistics displayed for each sort are: number of comparisons, number
	 * of data moves, and time (in milliseconds).
	 * </p>
	 * <p>
	 * Note: each sort is given the same array (i.e., in the original order) and
	 * the input array A is not changed by this method.
	 * </p>
	 * 
	 * @param A
	 *            the array to sort
	 */
	static public void runAllSorts(SortObject[] A) {
		System.out.format("%-23s%15s%15s%15s\n", "algorithm", "data compares",
				"data moves", "milliseconds");
		System.out.format("%-23s%15s%15s%15s\n", "---------", "-------------",
				"----------", "------------");

		// TODO: run each sort and print statistics about what it did

		// Selection Sort Output
		SortObject[] arr1 = new SortObject[A.length];
		System.arraycopy(A, 0, arr1, 0, A.length);

		long startSelection = System.currentTimeMillis();
		selectionSort(arr1);
		long endSelection = System.currentTimeMillis();
		printStatistics("selection", SortObject.getCompares(), getmoves(),
				endSelection - startSelection);
		resetMoves();
		SortObject.resetCompares();

		// Insertion Sort output
		SortObject[] arr2 = new SortObject[A.length];
		System.arraycopy(A, 0, arr2, 0, A.length);

		long startInsertion = System.currentTimeMillis();
		insertionSort(arr2);
		long endInsertion = System.currentTimeMillis();
		printStatistics("insertion", SortObject.getCompares(), getmoves(),
				endInsertion - startInsertion);
		resetMoves();
		SortObject.resetCompares();

		// Merge Sort output
		SortObject[] arr3 = new SortObject[A.length];
		System.arraycopy(A, 0, arr3, 0, A.length);

		long startMerge = System.currentTimeMillis();
		mergeSort(arr3);
		long endMerge = System.currentTimeMillis();
		printStatistics("merge", SortObject.getCompares(), getmoves(), endMerge
				- startMerge);
		resetMoves();
		SortObject.resetCompares();

		// Quick Sort output
		SortObject[] arr4 = new SortObject[A.length];
		System.arraycopy(A, 0, arr4, 0, A.length);

		long startQuick = System.currentTimeMillis();
		quickSort(arr4);
		long endQuick = System.currentTimeMillis();
		printStatistics("quick", SortObject.getCompares(), getmoves(), endQuick
				- startQuick);
		resetMoves();
		SortObject.resetCompares();

		// Heap Sort output
		SortObject[] arr5 = new SortObject[A.length];
		System.arraycopy(A, 0, arr5, 0, A.length);

		long startHeap = System.currentTimeMillis();
		heapSort(arr5);
		long endHeap = System.currentTimeMillis();
		printStatistics("heap", SortObject.getCompares(), getmoves(), endHeap
				- startHeap);
		resetMoves();
		SortObject.resetCompares();

		// Selection2Sort Output
		SortObject[] arr6 = new SortObject[A.length];
		System.arraycopy(A, 0, arr6, 0, A.length);

		long startSelection2 = System.currentTimeMillis();
		selection2Sort(arr6);
		long endSelection2 = System.currentTimeMillis();
		printStatistics("selection2", SortObject.getCompares(), getmoves(),
				endSelection2 - startSelection2);
		resetMoves();
		SortObject.resetCompares();

		// Insertion2Sort Output
		SortObject[] arr7 = new SortObject[A.length];
		System.arraycopy(A, 0, arr7, 0, A.length);

		long startInsertion2 = System.currentTimeMillis();
		insertion2Sort(arr7);
		long endInsertion2 = System.currentTimeMillis();
		printStatistics("insertion2", SortObject.getCompares(), getmoves(),
				endInsertion2 - startInsertion2);
		resetMoves();
		SortObject.resetCompares();
	}
}