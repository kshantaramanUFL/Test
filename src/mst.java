import java.io.*;
import java.util.*;

/**
 * @author = Karthik Shantaraman
 * 
 * The program implements Prim's algorithm 
 * on the adjacency list representation of a graph using Simple scheme and Fibonacci Heap
 * 
 */

/** Class Name : mst
 * Function : Starting point of the program
 * Data Members :
 *  g: Graph -> Object of class Graph to perform all the operations on the Graph
 */
public class mst {

	/**
	 * Function : main
	 * Parameter : args
	 * args[0]:
	 * -r : Random Mode
	 * -s : Simple Scheme
	 * -f : Fibonacci Heap Scheme
	 * args[1]
	 * n - Number of Vertices for Random Mode
	 * FileName - For User Input Mode
	 * args[2]
	 * d - Edge density of the graph
	 * 
	 */

	static Graph g;

	public static void main(String args[]) throws IOException
	{

		int numOfVertices=0,numOfEdges=0;
		if(args[0].equals("-r"))
		{
			int edgeCounter=0,density=0;
			mstEdge[] result;
			// Converting from String to Int : args[1] to numOfVertices
			for(int i=0;i<args[1].length();i++)
			{
				numOfVertices+=(args[1].charAt(i)-48)*Math.pow(10,args[1].length()-i-1);
			}
			result=new mstEdge[numOfVertices-1];
			for(int i=0;i<result.length;i++)
				result[i]=null;

			/* Allocating memory to the Graph object
				Since this is an adjacency list implementation, Memory allocation is done 
							based on number of vertices */
			g=new Graph(numOfVertices+1);

			// Converting from String to Int : args[2] to density
			for(int i=0;i<args[2].length();i++)
			{
				density+=(args[2].charAt(i)-48)*Math.pow(10,args[2].length()-i-1);
			}
			// Initializing number of edges based on edge density
			numOfEdges=(numOfVertices*(numOfVertices-1)/2)*density/100;
			Random rand=new Random();

			// Building the graph in random mode

			while(g.edgeCounter<=numOfEdges)
			{
				// Generating Random values for Edge properties (source,destination) and Cost
				int source=rand.nextInt(numOfVertices);
				int destination=rand.nextInt(numOfVertices);
				int cost=rand.nextInt(1000)+1;
				if(g.edgeCounter==numOfEdges)
				{
					break;
				}
				//Calling function to add Edge to Graph
				if(source!=destination)
					g.buildEdge(source, destination, cost);
			}
			System.out.println("Graph Constructed");
			g.checkConnected();
			long sTart=System.currentTimeMillis();
			//Calling function to compute MST using f-heap scheme
			g.prims_fibonacci();
			long eNd=System.currentTimeMillis();
			System.out.println(" Fibonacci : " + (eNd-sTart));
			sTart=System.currentTimeMillis();
			//Calling function to compute MST using Simple scheme
			g.prims();
			eNd=System.currentTimeMillis();
			System.out.println(" Simple : " + (eNd-sTart));
		}
		else
		{
			// User input mode reading the file
			// Input File should be present in same folder as the source code

			InputStream is = mst.class.getClassLoader().getResourceAsStream(args[1]);
			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(is));
			ArrayList<String> fileline=new ArrayList<String>();
			String lines=" ";
			int count=1;
			while ((lines= bufferedReader.readLine()) != null) {
				fileline.add(lines);
				count++;
			}
			bufferedReader.close();
			String[] startLine = fileline.get(0).split("\\s+");

			//Obtaining number of Vertices

			for(int i=0;i<startLine[0].length();i++)
				numOfVertices+=(startLine[0].charAt(i)-48)*Math.pow(10,startLine[0].length()-i-1);

			//Obtaining number of Edges

			for(int i=0;i<startLine[1].length();i++)
				numOfEdges+=(startLine[1].charAt(i)-48)*Math.pow(10,startLine[1].length()-i-1);
			mstEdge[] result;
			result=new mstEdge[numOfVertices+1];
			for(int i=0;i<result.length;i++)
				result[i]=null;
			g=new Graph(numOfVertices+1);

			for(int i=1;i<=numOfEdges;i++)
			{
				String[] eachline=fileline.get(i).split("\\s+");
				int x=0,y=0,cost=0;

				//Source Vertice
				for(int j=0;j<eachline[0].length();j++)
					x+=(eachline[0].charAt(j)-48)*Math.pow(10,eachline[0].length()-j-1);

				//Destination Vertice
				for(int j=0;j<eachline[1].length();j++)
					y+=(eachline[1].charAt(j)-48)*Math.pow(10,eachline[1].length()-j-1);

				//Cost of (Source,Destination)
				for(int j=0;j<eachline[2].length();j++)
					cost+=(eachline[2].charAt(j)-48)*Math.pow(10,eachline[2].length()-j-1);

				//Calling function to add edge to Graph
				g.buildEdge(x, y, cost);
			}

			// Simple scheme
			if(args[0].equals("-s"))
			{


				// Calling function that computes and returns the mst
				result = g.prims();

				//Denotes the MinimumSpanningCost
				int minCost=0; 
				System.out.println("Prims Tree is : ");

				//Printing out tree
				for(int i=0;i<numOfVertices-1;i++)
				{
					minCost+=result[i].getCost();
					System.out.println("("+result[i].getMstSource()+","+result[i].getMstDest()+")-->"+result[i].getCost());
				}
				System.out.println("Cost is : "+minCost);
			}
			// f-heap scheme
			else
			{
				int mincost=0;
				//g.printList();

				// Calling function that computes and returns the mst
				result=g.prims_fibonacci();
				for(int i=0;i<numOfVertices-1;i++)
				{
					mincost+=result[i].getCost();
					System.out.println("("+result[i].getMstSource()+","+result[i].getMstDest()+")-->"+result[i].getCost());
				}
				System.out.println(mincost);
			}
		}
	}
}

/**
 * class Name : Vertex
 * Used in : Representation of an Element in Vertex Adjacency List
 * Data Members :
 * 	destination : Destination Vertex
 * 	Weight : Cost to destination from Source
 * 
 */

class Vertex
{
	public int destination;
	public int weight;

	// Constructor that assigns value to the data members when a Vertex Object is created
	public Vertex(int destination, int weight) {
		this.destination = destination;
		this.weight = weight;
	}

	// Default Constructor
	public Vertex() {
		// TODO Auto-generated constructor stub
	}
}

/*
 * class Name : Graph
 * Denotes : Adjacency List  representation of Graph
 * Data Members : 
 * 		vertices : ArrayList<Vertex> -> Represents ArrayList of Vertex Objects
 * 		edgeCounter : int -> Keeps a count of the number of edges added to the graph
 *  
 * Member Functions:
 * 		buildEdge(int,int,int) : void -> Adds the Edge to Graph 
 * 		checkConnected() : Boolean -> Checks if Graph is connected or not
 * 		prims() : mstEdge[] -> Computes MST using Prim's using Simple Data Structures 
 * 		prims_fibonacci() : mstEdge[] -> Computes MST using Prim's using f-heap
 *
 */

class Graph
{	
	ArrayList<Vertex>[] vertices;
	int edgeCounter;

	/** 
	 * Parameterized Constructor of class Graph
	 * @ param : numofVertices - int
	 * Constructor used to initialize data members of the class 
	 */

	Graph(int numOfVertices)
	{
		vertices=new ArrayList[numOfVertices];
		edgeCounter=0;
		for(int i=0;i<numOfVertices;i++)
		{
			vertices[i]=new ArrayList<Vertex>();
		}
	}
	/** 
	 * Method: buildEdge(int,int,int) : void
	 * 
	 * @param: source - int
	 * @param: destination - int
	 * @param: cost - int
	 * 	Constructs an Edge from source->destination with Cost=cost
	 * 
	 * Return: Nothing 
	 */
	public void buildEdge(int source,int destination,int cost)
	{
		//Flag to represent whether Edge has to be added or not
		Boolean addFlag=true; 

		// If the adjacency list of the corresponding vertex is not empty
		if(!vertices[source].isEmpty())
			for (Iterator<Vertex> iterator = vertices[source].iterator(); iterator.hasNext();) 
			{
				Vertex v= (Vertex) iterator.next();
				//Destination Vertex already present in Graph
				if(v.destination==destination)
				{
					addFlag=false;
					break;
				}

				// The (destination,source) edge is not present
				else if(v.destination==source)
				{
					addFlag=false;
					break;
				}
			}
		// Add the edge
		if(addFlag==true)
		{
			// Since graph is undirected
			// adding from source -> destination
			// adding from destination -> source
			vertices[source].add(new Vertex(destination, cost));
			vertices[destination].add(new Vertex(source, cost));
			edgeCounter++;
		}
	}

	/** 
	 * Method: checkConnected() : Boolean
	 * 		Checks whether Graph is connected or not
	 * Return: 
	 * 			True : Graph is connected
	 * 			False : Graph not connected 
	 */
	public Boolean checkConnected()
	{
		for(int i=0;i<vertices.length;i++)
		{
			// If any one of the vertices list does not have an adjacent verTex 
			//Graph is not connected

			if(vertices[i].isEmpty())
			{
				return false;
			}
		}
		return true;
	}

	/** 
	 * Method: prims() : mstEdge[]
	 * Finds the minimum spanning tree of the Graph using Simple Data Structure
	 * Return: Array of spanning edges
	 * 
	 * Working:
	 * 
	 * 1. Set source Vertex
	 * 2. Set visited Flag of starting vertex is true
	 * 3. Loop till # of vertices visited = Total # of vertices-1
	 * 4. Add  the source vertex to visited Vertices List
	 * 5. Loop till size of the visited vertices list
	 * 6. Iterate over the list and Find minimum cost edge
	 * 7. Set the destination as next source 
	 * 8. Update visited flag of that vertex and add it to visited vertices list
	 * 9. Continue looping
	 *  
	 */  
	public mstEdge[] prims()
	{
		ArrayList<Vertex> starting_vertex=vertices[0]; // List of starting vertices
		ArrayList<Integer> visitedVertices= new ArrayList<Integer>(); //list of visited vertices
		Vertex minKeyNode= null;
		boolean visitFlag[]= new boolean[vertices.length-1];
		int source=0,visitCounter=1;
		for(int i=0;i<visitFlag.length;i++)
			visitFlag[i]=false;
		visitedVertices.add(source); //Starting from Source Vertex 0
		visitFlag[source]=true;
		mstEdge[] mstedge = new mstEdge[vertices.length-1];
		for(int i=0;i<vertices.length-1;i++)
			mstedge[i]=null;
		//Loop till number of vertices visited
		while(visitCounter<vertices.length)  
			//equals numberofVertices in Graph
		{
			int mincost=Integer.MAX_VALUE;          
			for(int i : visitedVertices) //Loop to find Minimum Cost Vertex among all Visited Vertices
			{
				starting_vertex=vertices[i];	
				for (Iterator<Vertex> iterator = starting_vertex.iterator(); iterator.hasNext();)
				{
					Vertex v=(Vertex) iterator.next();

					//If the destination vertex has not yet been visited from this source
					if(!visitFlag[v.destination])
					{
						if(v.weight<mincost)
						{
							minKeyNode=v;
							mincost=v.weight;
							source=i;
						}
					}
				}
			}
			// Add the visited edge to the set of Spanning edges
			mstedge[visitCounter-1]=new mstEdge(source, minKeyNode.destination, minKeyNode.weight);
			visitCounter++;

			//Updating VisitedVertices Array
			visitedVertices.add(minKeyNode.destination);

			//Updating flag
			visitFlag[minKeyNode.destination]=true;
		}

		//Return the set of spanning edges
		return mstedge;
	}

	/** 
	 * Method: prims_fibonacci() : mstEdge[]
	 * Finds the minimum spanning tree of the Graph using f-heap
	 * Return: Array of spanning edges
	 * 
	 * Working:
	 * 
	 * 1. Set source Vertex
	 * 2. Set visited Flag of starting vertex is true
	 * 3. Loop till # of visited vertices = # of vertices- -1
	 * 4. Add adjacent vertices of source vertex to F-heap. 
	 * 		Before adding check whether the adjacency list of that vertex has already been added 
	 * 5. Extract min from the fibonacci heap
	 * 6. Add the source,destination and cost of min vertex to Spanning tree
	 * 7. Update destination as next source and set the visited flag to true
	 * 
	 */

	public mstEdge[] prims_fibonacci() {
		// TODO Auto-generated method stub

		//
		int visitCount=1,minCostDestination=0,mincost,minCostSource=0;
		FibonacciHeap<Integer> fHeap=new FibonacciHeap<Integer>();
		int source=0;
		//ArrayList<Integer> visitedVertices= new ArrayList<Integer>();
		ArrayList<Vertex> starting_vertex=vertices[source];
		boolean visitFlag[]=new boolean[vertices.length-1];
		for(int i=0;i<visitFlag.length;i++)
			visitFlag[i]=false;
		FibonacciHeapNode<Integer> fHeapNode;
		FibonacciHeapNode<Integer> fminKeyNode=new FibonacciHeapNode<Integer>(0,source,0);
		visitFlag[source]=true;
		mstEdge[] mstedge=new mstEdge[vertices.length-1];
		while(visitCount<vertices.length-1)
		{
			starting_vertex=vertices[source];
			for (Iterator<Vertex> iterator = starting_vertex.iterator(); iterator.hasNext();)
			{
				Vertex v=(Vertex) iterator.next();
				if(visitCount==1)
				{
					fHeapNode=new FibonacciHeapNode<Integer>(v.destination,source,v.weight);
					fHeap.insert(fHeapNode, v.weight);
				}				
				else
				{
					//Boolean addFlag=true;
					fHeapNode=new FibonacciHeapNode<Integer>(v.destination,source,v.weight);
					if(!visitFlag[v.destination])
						fHeap.insert(fHeapNode, v.weight);
				}

			}
			//System.out.println(fHeap.toString());
			minCostDestination=(int) fHeap.getMin().value;
			mincost=(int)fHeap.getMin().key;
			minCostSource = (int) fHeap.getMin().source;
			//fMinNode=fHeap.min();
			fminKeyNode.value=minCostDestination;
			fminKeyNode.key=mincost;
			
			//If the MinCostVertex has already been added before to Spanning tree, do not add it again
			if(!visitFlag[minCostDestination])
			{
				mstedge[visitCount-1]=new mstEdge(minCostSource, minCostDestination, mincost);
				//System.out.println("(" +minCostSource+","+minCostDestination+")->"+mincost);
				visitCount++;
				visitFlag[minCostDestination]=true;
				source=minCostDestination;
			}
			fHeap.removeMin();
		}
		return mstedge;
	}
}

/**
 * class Name : mstEdge
 * Denotes : Wrapper class for Minimum Spanning Edge
 * Data Members : 
 * 		mstSource: int -> Represents source Vertex of Spanning Edge
 *		mstDestination: int -> Represents source Vertex of Spanning Edge
 *		cost: int -> Weight of mstEdge 
 * Member Functions :
 * 		getMstSource() : int -> Returns the source Edge
 *  	getMstDest() : int -> Returns Destination Edge
 *  	getCost()    : int - Returns Cost
 *  
 */

class mstEdge
{
	int mstSource,mstDest,cost;

	/** Constructor mstEdge(int,int,int)
	 * 
	 * @param : mstSource : int - Source vertex of minimum Spanning edge 
	 * @param : mstDest : int - Destination vertex of minimum Spanning edge	
	 * @param : cost : int - Cost minimum Spanning edge
	 * 
	 * Sets the value of the mstEdge
	 * 		
	 */
	
	public mstEdge(int mstSource, int mstDest, int cost) {
		this.mstSource = mstSource;
		this.mstDest = mstDest;
		this.cost = cost;
	}
	
	/** 
	 * Method: getMstSource() : int
	 * 	Retrieves the value of Source Vertex of the Edge
	 * Return: index of source Edge
	 * 
	 */
	
	public int getMstSource() {
		return mstSource;
	}
	/** 
	 * Method: getMstDest() : int
	 * 	Retrieves the value of Destination Vertex of the Edge
	 * Return: index of Destination Edge
	 * 
	 */
	public int getMstDest() {
		return mstDest;
	}
	/** 
	 * Method: getCost() : int
	 * Retrieves the cost from Source to Destination
	 * Return: Cost of Edge
	 * 
	 */
	public int getCost() {
		return cost;
	}
}

/**
 * class Name : FibonacciHeap<Integer>
 * Denotes : Represents an Fibonacci Heap of Integers
 * Data Members : 
 * 		HeapSize : Double -> Represents Size of Heap
 * 		minKeyNode : FibonacciHeapNode -> Represents the node with the least Key in the Heap
 * 		numNodes : int -> Represents the number of nodes in the heap
 *  
 * Member Functions:
 * 		insert(FibonacciHeapNode<Integer> fHeapNode,int key) : void - Inserts a node into the heap with given Key
 * 		combine(FibonacciHeapNode<Integer> fHeapNode1,FibonacciHeapNode<Integer> fHeapNode2) : void -> Links two nodes by making one a child of other
 * 		getMin() : FibonacciHeapNode<Integer> -> Returns the node with the least key from the heap 
 * 		consolidate() : void -> Consolidates the tree when a remove min happens
 * 		removeMin() : void -> Removes the node with the least key from the heap
 *
 */

class FibonacciHeap<Integer>
{
	private static final double HeapSize =
			1.0 / Math.log((1.0 + Math.sqrt(5.0)) / 2.0);

	//Fibonacci Heap Node with lowest key
	private FibonacciHeapNode<Integer> minKeyNode;

	//Number of nodes in the heap.
	private int numNodes;

	/** 
	 * Method: insert(FibonacciHeapNode<Integer> node, int key) : void
	 * 
	 * Inserts a new node into the heap with Key Value key
	 * 
	 * Return: Nothing 
	 * 
	 */
	
	public void insert(FibonacciHeapNode<Integer> node, int key)
	{
		node.key = key;

		// Add the node into the Minimum Heap list
		if (minKeyNode != null) {
			node.leftptr = minKeyNode;
			node.rightptr = minKeyNode.rightptr;
			minKeyNode.rightptr = node;
			node.rightptr.leftptr = node;

			if (key < minKeyNode.key) {
				minKeyNode = node;
			}
		} else {
			minKeyNode = node;
		}

		numNodes++;
	}

	/** 
	 * Method: getMin() : FibonacciHeapNode<Integer>
	 * 
	 * Returns the element with the minimum key
	 * 
	 * Return: FibonacciHeapNode<Integer> - Minimum Node 
	 * 
	 */
	
	public FibonacciHeapNode<Integer> getMin()
	{
		return minKeyNode;
	}

	/** 
	 * Method: removeMin() : FibonacciHeapNode<Integer>
	 * 
	 * Removes the smallest element from the heap.
	 *
	 *  Return : Nothing
	 */
	
	public void removeMin()
	{
		FibonacciHeapNode<Integer> tempNode = minKeyNode;

		if (tempNode != null) {
			int numKids = tempNode.degree;
			FibonacciHeapNode<Integer> node1 = tempNode.child;
			FibonacciHeapNode<Integer> tempRight;

			while (numKids > 0) {
				tempRight = node1.rightptr;

				// Reassigning pointers
				node1.leftptr.rightptr = node1.rightptr;
				node1.rightptr.leftptr = node1.leftptr;
				node1.leftptr = minKeyNode;
				node1.rightptr = minKeyNode.rightptr;
				minKeyNode.rightptr = node1;
				node1.rightptr.leftptr = node1;

				// Point parent of the root to null
				node1.parent = null;
				node1 = tempRight;
				numKids--;
			}

			// remove tempNode from root list of heap
			tempNode.leftptr.rightptr = tempNode.rightptr;
			tempNode.rightptr.leftptr = tempNode.leftptr;

			if (tempNode == tempNode.rightptr) {
				minKeyNode = null;
			} else {
				minKeyNode = tempNode.rightptr;
				consolidate();
			}

			// decrementing heap size as a node has been removed
			numNodes--;
		}
	}

	/** 
	 * Method: consolidate() : void
	 * 
	 * Consolidate the heap after a removeMin has been performed
	 * 
	 * Return: Nothing 
	 **/
	
	protected void consolidate()
	{
		int arraySize =
				((int) Math.floor(Math.log(numNodes) * HeapSize)) + 1;

		List<FibonacciHeapNode<Integer>> array =
				new ArrayList<FibonacciHeapNode<Integer>>(arraySize);

		for (int i = 0; i < arraySize; i++) {
			array.add(null);
		}

		int numRoots = 0;
		FibonacciHeapNode<Integer> node1 = minKeyNode;

		if (node1 != null) {
			numRoots++;
			node1 = node1.rightptr;

			while (node1 != minKeyNode) {
				numRoots++;
				node1 = node1.rightptr;
			}
		}

		
		while (numRoots > 0) {
			
			int d = node1.degree;
			FibonacciHeapNode<Integer> next = node1.rightptr;

			for (;;) {
				FibonacciHeapNode<Integer> node2 = array.get(d);
				if (node2 == null)
					break;
				
				if (node1.key > node2.key) {
					FibonacciHeapNode<Integer> temp = node2;
					node2 = node1;
					node1 = temp;
				}

				combine(node2, node1);

				array.set(d, null);
				d++;
			}

			array.set(d, node1);

			// Iterate through list.
			node1 = next;
			numRoots--;
		}

		// Point minKeyNode to null and reconstruct the root list
		minKeyNode = null;

		for (int i = 0; i < arraySize; i++) {
			FibonacciHeapNode<Integer> node2 = array.get(i);
			if (node2 == null) {
				continue;
			}

			// Adding minKeyNode to Heap.
			if (minKeyNode != null) {
				node2.leftptr.rightptr = node2.rightptr;
				node2.rightptr.leftptr = node2.leftptr;
				node2.leftptr = minKeyNode;
				node2.rightptr = minKeyNode.rightptr;
				minKeyNode.rightptr = node2;
				node2.rightptr.leftptr = node2;

				if (node2.key < minKeyNode.key) 
					minKeyNode = node2;
			} else
				minKeyNode = node2;
		}
	}

	/** 
	 * Method: combine(FibonacciHeapNode<Integer> node2, FibonacciHeapNode<Integer> node1) : void
	 * 	Combines node node1 and node2 by making node2 a child of node1.
	 * 
	 * @param node2 : FibonacciHeapNode<Integer> - Node that becomes a child
	 * @param node1 : FibonacciHeapNode<Integer> - Node that becomes the parent
	 * 
	 * Return: Nothing
	 *  
	 **/
	
	protected void combine(FibonacciHeapNode<Integer> node2, FibonacciHeapNode<Integer> node1)
	{
		// remove node2 from root list of heap
		node2.leftptr.rightptr = node2.rightptr;
		node2.rightptr.leftptr = node2.leftptr;

		// make node2 a child of node1
		node2.parent = node1;

		if (node1.child == null) {
			node1.child = node2;
			node2.rightptr = node2;
			node2.leftptr = node2;
		} else {
			node2.leftptr = node1.child;
			node2.rightptr = node1.child.rightptr;
			node1.child.rightptr = node2;
			node2.rightptr.leftptr = node2;
		}

		node1.degree++;
		node2.mark = false;
	}
}

/**
 * 
 * class Name : FibonacciHeapNode<Integer>
 * Denotes : Represents a node of the Fibonacci Heap of Integers
 * Data Members : 
 * 		data : Integer -> Value of the Node. Denotes Destination vertice of an Edge from source
 * 		source : Integer -> Variable to keep track of source vertice of an edge
 * 		child : FibonacciHeapNode -> first Child Node
 * 		leftptr : FibonacciHeapNode -> left sibling Node
 * 		parent : FibonacciHeapNode -> parent Node
 * 		rightptr : FibonacciHeapNode -> right sibling Node
 * 		mark : Boolean -> Fag to keep track whether the node has lost a child
 * 		key : int -> Represents the cost of an Edge from source->data. Key value of F-Heap Node
 * 		degree : int -> Number of children of this node
 *  
 * Member Functions:
 * 		getKey() : int -> Obtains the key value of this Node
 * 		getdata() : int -> Obtains the data of this Node
 *
 */

class FibonacciHeapNode<Integer>
{
	//Node value.
	int value;
	
	//Node source to obtain source in Prim's
	int source;
	
	//first child 
	FibonacciHeapNode<Integer> child;
	
	//left sibling 
	FibonacciHeapNode<Integer> leftptr;
	
	//parent 
	FibonacciHeapNode<Integer> parent;
	
	//right sibling
	FibonacciHeapNode<Integer> rightptr;
	
	// Flag to track cascading cut operations
	boolean mark;
	
	//key value for this node. Represents the costs of an Edge
	int key;
	
	//number of children of this node 
	int degree;
	
	/**
	 *  
	 * Parameterized Constructor of class FibonacciHeapNode.
	 * Initializes the pointers of the Node, and makes this a circularly doubly-linked list.
	 * 
	 * @ param : value - int : value for this node. Denotes Destination vertice
	 * @ param : data1 - int : Sets source vertice. Useful for Prims
	 * @ param : key- int : initial key for this node. Denotes Cost of Edge between (data1,value)
	 * Constructor used to initialize data members of the class
	 *  
	 */
	
	public FibonacciHeapNode(int value, int data1,int key)
	{
		rightptr = this;
		leftptr = this;
		this.value = value;
		this.key = key;
		this.source =data1;
	}

	/**
	 *  
	 * Method: getKey() : int
	 * 
	 * Returns the cost of the edge for Prim's
	 * Returns key of the node for combine and removeMin and min operations of FibonacciHeap class
	 * 
	 * Return: key 
	 * 
	 */
	
	public int getKey()
	{
		return key;
	}

	/** 
	 * Method: getvalue() : int
	 * 
	 * Returns the destination of the edge from Prim's
	 * Returns value of the node for combine,removeMin and min operations of FibonacciHeap class
	 * 
	 * Return: key 
	 * 
	 */
	
	public int getvalue()
	{
		return value;
	}
}