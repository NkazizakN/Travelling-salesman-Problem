//To use a different file to Read data in 
//Please update the public variable <fileName> and include the extension ".txt"

package cst3170.cw1;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
public class TravellingSalesProblem {
	
	//UPDATE ME
	public static String fileName ="test1tsp.txt";
	
	public static int NumberOfCities;
	public static Point[] cities = new Point[500];
	public static double[][] adjacencyMatrix = new double[500][500];	
	public static  int start;
	public static List<Integer> tour = new ArrayList<>();
	public static double minTourCost = Double.MAX_VALUE;
	public static boolean ranSolver = false;

	public static void main(String[] args) {
		
		//starting the clock
		long timeStart = System.nanoTime();
        String dataFile = System.clearProperty("user.dir") + File.separator + "CityData" + File.separator + fileName;

        File file = new File(dataFile);
        Scanner sc;

        sc = null;
        try {
            sc = new Scanner(file);
        } catch (FileNotFoundException f) {
            System.out.println("File Not Found ...");
            System.exit(0);
        }
        int count = 0;
//reading City data from File
        while (sc.hasNextLine()) {

        	String line = sc.nextLine().trim();
        	if(line == "")
        	{
        		continue;
        	}
        	//try to read file with spaces
        	try {
        	String temp[] = line.split(" ");
        	int position = Integer.parseInt(temp[0]);
        	Point city = new Point();
        	city.setX(Double.parseDouble(temp[1]));
        	city.setY(Double.parseDouble(temp[2]));
        	cities[position] = city;
        	count++;
        	}
        	//catch to read file with tabs
        	catch(NumberFormatException ex)
        	{
        		String temp[] = line.split("	");
            	int position = Integer.parseInt(temp[0]);
            	Point city = new Point();
            	city.setX(Double.parseDouble(temp[1]));
            	city.setY(Double.parseDouble(temp[2]));
            	cities[position] = city;
            	count++;
        	}
        }

        NumberOfCities = count;
        //populating the entire array with 0;
        for (double[] row : adjacencyMatrix) java.util.Arrays.fill(row, 0);
        createAdjacencyMatrix();
        
        start = 1;
        List<Integer> finalPath = getTour();   
        System.out.print("Tour :");
        for(int i : finalPath) 
        {
        	System.out.print(" -> "+ (i+1));
        }
        System.out.println("\n\nTour cost: " + getTourCost());
        long timeEnd = System.nanoTime();
        System.out.println("\n\nTotal time taken (in Nano second):" + (timeEnd - timeStart));
	}//end of Main()
	
	public static List<Integer> getTour()
	{
		if(!ranSolver)
			solve();
		return tour;
	}
	
	public static double getTourCost()
	{
		if(!ranSolver)
			solve();
		return minTourCost;
	}
	public static void solve() {

	    if (ranSolver) return;

	    final int END_STATE = (1 << NumberOfCities) - 1;
	    Double[][] memo = new Double[NumberOfCities][1 << NumberOfCities];

	    // Add all outgoing edges from the starting node to memo table.
	    for (int end = 0; end < NumberOfCities; end++) 
	    {
	      if (end == start) continue;
	      memo[end][(1 << start) | (1 << end)] = adjacencyMatrix[start][end];
	    }

	    for (int r = 3; r <= NumberOfCities; r++) 
	    {
	      for (int subset : combinations(r, NumberOfCities)) 
	      {
	        if (notIn(start, subset)) continue;
	        for (int next = 0; next < NumberOfCities; next++) 
	        {
	          if (next == start || notIn(next, subset)) 
	        	  {
	        	  	continue;
	        	  }
	          int subsetWithoutNext = subset ^ (1 << next);
	          double minDist = Double.MAX_VALUE;
	          for (int end = 0; end < NumberOfCities; end++) 
	          {
	            if (end == start || end == next || notIn(end, subset)) continue;
	            double newDistance = memo[end][subsetWithoutNext] + adjacencyMatrix[end][next];
	            if (newDistance < minDist) 
	            {
	              minDist = newDistance;
	            }
	          }
	          memo[next][subset] = minDist;
	        }
	      }
	    }

	    // Connect tour back to starting node and minimize cost.
	    for (int i = 0; i < NumberOfCities; i++) 
	    {
	      if (i == start) continue;
	      double tourCost = memo[i][END_STATE] + adjacencyMatrix[i][start];
	      if (tourCost < minTourCost) 
	      {
	        minTourCost = tourCost;
	      }
	    }

	    int lastIndex = start;
	    int state = END_STATE;
	    tour.add(start);

	    // Reconstruct TSP path from memo table.
	    for (int i = 1; i < NumberOfCities; i++) 
	    {
	      int bestIndex = -1;
	      double bestDist = Double.MAX_VALUE;
	      for (int j = 0; j < NumberOfCities; j++) 
	      {
	        if (j == start || notIn(j, state)) continue;
	        double newDist = memo[j][state] + adjacencyMatrix[j][lastIndex];
	        if (newDist < bestDist) 
	        {
	          bestIndex = j;
	          bestDist = newDist;
	        }
	      }
	      tour.add(bestIndex);
	      state = state ^ (1 << bestIndex);
	      lastIndex = bestIndex;
	    }

	    tour.add(start);
	    Collections.reverse(tour);

	    ranSolver = true;
	  }
	  private static boolean notIn(int elem, int subset) 
	  {
		    return ((1 << elem) & subset) == 0;
	  }
	  // This method generates all bit sets of size n where r bits
	  // are set to one. The result is returned as a list of integer masks.
	  public static List<Integer> combinations(int r, int n) 
	  {
	    List<Integer> subsets = new ArrayList<>();
	    combinations(0, 0, r, n, subsets);
	    return subsets;
	  }
	  
	  
	  // To find all the combinations of size r we need to recurse until we have
	  // selected r elements (r = 0), otherwise if r != 0 then we still need to select an element which is found after the position of our last selected element
 
	  private static void combinations(int set, int at, int r, int n, List<Integer> subsets) 
	  {

	    // Return early if there are more elements left to select than what is available.
	    int elementsLeftToPick = n - at;
	    if (elementsLeftToPick < r) return;

	    // We selected 'r' elements so we found a valid subset!
	    if (r == 0) 
	    {
	      subsets.add(set);
	    } else {
	      for (int i = at; i < n; i++) 
	      {
	        // Try including this element
	        set ^= (1 << i);

	        combinations(set, i + 1, r - 1, n, subsets);

	        // Backtrack and try the instance where we did not include this element
	        set ^= (1 << i);
	      }
	    }
	  }
	
	//a function to calculate the distance between any two point	
	public static double getDistance(double x1, double y1, double x2, double y2)
	{
			
		double dist = Math.sqrt( Math.abs((x2-x1)*(x2 - x1)) + Math.abs(( y2 - y1) *(y2-y1)));
		return dist;
	}//end of getDistance()
	
	
	public static void createAdjacencyMatrix() 
	{
		for(int row = 0; row < NumberOfCities; row++)
		{
			for(int column = 0; column < NumberOfCities; column++)
			{
				//distance between same city is 0 thus we are skipping this iteration
				if(row == column)
				{
					continue;
				}
				if(adjacencyMatrix[row][column] == 0)
				{
					//if the distance between city 1 and city 2 is 100 then the distance between city 2 and city 1 should be the same
					adjacencyMatrix[row][column] = getDistance(cities[row+1].getX(), cities[row+1].getY(),cities[column+1].getX(),cities[column+1].getY());
					adjacencyMatrix[column][row] = adjacencyMatrix[row][column];
				}
				else 
				{
					continue;
				}
						
			}
		}//end of loop for populating adjacency Matrix

	}

}
