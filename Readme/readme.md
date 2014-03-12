MINIMUM SPANNING TREE USING PRIM’S ALGORITHM

I.	COMPILATION INSTRUCTIONS

COMPILER USED:  Eclipse Compiler for Java / OpenJDK for Ubuntu/ Javac for Windows
HOW TO COMPILE:
The program compiles normally in all the cases. However, in random mode, there is a specific when the number of vertices is 5000 and the edge density is >50%. The program has to be compiled in Eclipse IDE with the following VM arguments
-Xmx2048M
This increases the heap space of the Eclipse Compiler.

Normally, the program can be compiled from command prompt using
javac mst.java
During execution time, the input file and the program have to present in the same folder.








II.	Program Structure
The file mst.java consists of the following 6 classes
1)	public class mst
2)	class Graph
3)	class Vertex
4)	class mstEdge
5)	class FibonacciHeap<Integer>
6)	class FibonacciHeapNode<Integer>

1)	public class mst
The primary class of the program that contains the main method which is the entry point of the program.
Data Members: 
	g: Graph -> Object of class Graph to perform all the operations on the Graph
Member functions:
public static void main(String args[])
This function is the starting point of the program. It receives the command line arguments parses them and decided whether to perform Random Mode operations of User Input Mode operations. In the User Input mode, the function switches between Simple Scheme and F-Heap Scheme based on the command line argument. 
-r for Random mode
-s for Simple scheme in User Input Mode
-f for F-Heap scheme in User Input Mode






2) class Graph
This class is the Adjacency List implementation of the Graph
Data Members
	vertices: ArrayList<Vertex> -> Represents an ArrayList of objects of Class Vertex
	edgeCounter: int -> Keeps a count of the number of edges added to the graph
Member functions
	Graph(int numOfVertices)
Parameterized Constructor of class Graph that is used to determine the size of the ArrayList vertices based of the number of vertices and initializes counter value to Zero.
	public void buildEdge(int source,int destination,int cost)

This functions adds an edge between source and destination with cost. Before adding an edge the function checks the following

	The reverse edge (destination,source) is not present already
	The edge (source,destination) is not already present

If both the above checks pass, then the edge is added to the graph.

	public Boolean checkConnected()
This functions checks whether the Graph is connected. If any of the individual Array List of a Vertex is empty it means Graph is not connected. 
Functions returns False in such a case, else function returns true
	public mstEdge[] prims()
This is the function that computes the Minimum Spanning Tree using Simple Data structure. This is how the function works
1)	Loop runs until number of Vertices Visited = total number of vertices-1
2)	Picks up a starting Vertex. Set visitedFlag of that vertice to true and adds it to the list of Visited Vertices
3)	Iterate through the adjacency list of the visited Vertices and pick up minimum cost element
4)	Set this destination as next source. Set visitedFlag as true. Add vertex to set of visited Vertices. 
5)	Add the minimum Cost Vertex to MST by creating object of class mstEdge.
Function returns list of spanning edges of type mstEdge[] – Array of objects of type mstEdge


	public mstEdge[] prims_fibonacci()

This function computes the Minimum Spanning Tree using F-Heap Data Structure. This is how the function works.

1)	Set source Vertex
2)	Set visited Flag of starting vertex is true
3)	Loop till # of visited vertices = # of vertices- -1
4)	Add adjacent vertices of source vertex to F-heap. Before adding check whether the adjacency list of that vertex has already been added
5)	Extract minimum from the fibonacci heap.
6)	Add the source,destination and cost of min vertex to Spanning tree
7)	Update destination as next source and set the visited flag to true
8)	If the Minimum Cost Vertex has already been added before to Spanning tree, do not add it again
9)	Remove the minimum element from the heap.

Function returns list of spanning edges of type mstEdge[] – Array of objects of type mstEdge
3) class Vertex

This class denotes the representation of an element of the Adjacency List.
Data Members
	destination: int -> Destination Vertex
	Weight: int -> Cost to destination from Source
Member Functions
	public Vertex(int destination, int weight)
Parameterized Constructor of class Vertex that is used to set the values of the destination and weight of an Edge for a vertex in the graph adjacency list.


4) class mstEdge
This class is the representation of the Minimum Spanning Tree obtained using Prims.
Data Members
	mstSource : int -> Source Vertex
	mstDest : int -> Destination Vertex
	cost : int -> Cost of Minimum Spanning Edge
Member Functions
	public mstEdge(int mstSource, int mstDest, int cost)
Parameterized Constructor of class mstEdge that is used to set the values of the source, destination and cost of the mstEdge.
	public int getMstSource()

This functions returns the source Vertex of an edge of the Minimum Spanning Tree.

	public int getMstDest()

This functions returns the destination Vertex of an edge of the Minimum Spanning Tree

	public int getCost()

This functions returns the cost of an edge of the Minimum Spanning Tree















5)	class FibonacciHeapNode<Integer>
Represents a node in the Fibonacci Heap of Integers. The <Integer> indicates that, the value stored in the Node is an Integer.
Data Members
	value: int -> Value of the Node. Denotes Destination vertice of an Edge from source
	source: int -> Variable to keep track of source vertice of an edge
	child: FibonacciHeapNode -> first Child Node
	leftptr: FibonacciHeapNode -> left sibling Node
	parent: FibonacciHeapNode -> parent Node
	rightptr: FibonacciHeapNode -> right sibling Node
	mark: Boolean -> Fag to keep track whether the node has lost a child
	key: int -> Represents the cost of an Edge from source->value. Key value of F-Heap Node
	degree: int -> Number of children of this node
Member Functions
	public FibonacciHeapNode(int value, int data1,int key)

Parameterized Constructor of class FibonacciHeapNode. Initializes the right and left pointers, making this a circular doubly-linked list.

value indicates destination of Edge
data1 indicates source of Edge
key indicates cost between value and data1

	  public int getKey()
Returns the cost of the edge for Prim's and key of the node for link and removeMin and min operations of the Fibonacci Heap.
	public int getValue()
Returns the destination of the edge from Prim's and value of the node for link, removeMin and min operations of the Fibonacci class.




6)	class FibonacciHeap<Integer>
Represents an Fibonacci Heap of Integers.
Data Members
	HeapSize : Double -> Represents Size of Heap
	minNode : FibonacciHeapNode -> Represents the node with the least Key in the Heap
	numNodes : int -> Represents the number of nodes in the heap
Member Functions
	public void insert(FibonacciHeapNode<Integer> node, int key)
Inserts a new Node element into the heap with Key Value key
	protected void combine(FibonacciHeapNode<Integer> node2, FibonacciHeapNode<Integer> node1)
Combines nodes node1 and node2 by making node2 a child of node1
	protected void consolidate()

Consolidate the heap after a removeMin has been performed.

	public FibonacciHeapNode<Integer> getMin()
Returns the element with the minimum key. This is used in Prim’s to obtain the edge with minimum cost.
	public void removeMin()

Removes the element with the minimum weight from heap. This is used in Prim’s to remove the extract the Minimum Spanning Edge







III.	Comparison of Results

The program runs in two modes, Simple Scheme and F-Heap scheme. The simple scheme has been implemented using an ArrayList.
Expectation
For sparse graphs, the performance offered by both these schemes is similar.
For Dense Graph, Graphs having more than a thousand vertices and 1000 edges, the Fibonacci heap scheme is expected to run faster than the simple scheme.
For a graph G= (V,E), consisting of N vertices, the complexity of Simple Scheme is O(N2) whereas for Fibonacci Scheme, complexity is of the order O(N log N + E), because the amortized complexity of removeMin, insertion are O(1) and O(log N) respectively.
Results of Execution
From the table below it can be seen that the values generate corroborate with our expectation.
The values have been averaged over 5 trials for each value of Number of Vertices and Edge Density.
				1000 Vertices		3000 Vertices	5000 Vertices
Edge Densities	Simple	Fibonacci	Simple	Fibonacci	Simple	Fibonacci
10				533		109			12815	1263		57359	2951
20				893		170			24215	1553		110843	7102
30				1302	161			36004	3396		163814	10864
40				1609	193			46053	3868		187700	11193
50				2055	409			57977	5073		262154	16306
60				2479	1097		67639	6110		319117	16970
70				2973	1087		77809	6635		510265	23749
80				3295	1103		88502	11337		599707	27398
90				3770	1185		104837	12044		607238	31974
100				4070	1200		126226	12188		759466	34382

