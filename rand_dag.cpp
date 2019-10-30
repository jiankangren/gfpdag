//Generate a random connected directed-acyclic-graph
//David Ferry
//Sept 23, 2011
// Minor modification by Son Dinh for G-FP analysis (Oct 2018).


#include <iostream>
#include <list>
#include <cstdlib>
#include "rand_dag.h"
#include "common.h"

using namespace std;

//Suppose a graph has no cycles in it, and we want to check if adding a
//particular edge introduces a cycle. Since there are no previous cycles
//in the graph this new edge must be a part of any introduced cycle.
//Thus, if we want to add a new edge originating in "start" and 
//terminating in "end" then we need only do a DFS from "end" to see
//if "start" is reachable. 

//If this method returns "false" then a cycle exists, if it returns
//"true" then no cycle exists
bool acyclic_check(vector<vector<unsigned int> >& am, 
				   unsigned int size, 
				   unsigned int start, 
				   unsigned int end){
	
	if( start == end){
		return false;
	}

	vector<bool> visited (size,false);
	
	unsigned int current;

	//DFS code from here to end of routine
	list<unsigned int> stack;
	stack.push_front(end);
	
	while(!stack.empty()){
		current = stack.front();
		stack.pop_front();
		visited[current] = true;
	
		//If the start node is reachable from any node on the stack then there
		//is a cycle in the graph
		if(am[current][start] == 1){
			return false;
		}

		//If not we push all unvisited children onto the stack
		for(unsigned int i=0; i<size; i++){
			if(am[current][i] == 1 && visited[i] == false){
				stack.push_front(i);
			}
		}
	} 

	//If we've searched the whole graph and not gotten back to start then
	//there is no cycle in the graph. 
	return true;
}

//A digraph is weakly connected if a DFS along the undirected analog can reach
//every node in the graph
bool weak_conn_test(vector<vector<unsigned int> >& am){
	
	unsigned int size = am.size();

	vector<bool> visited (size,false);
	
	list<unsigned int> stack;
	stack.push_front(0);	

	unsigned int current;

	//DFS traversal with undirected edges
	while(!stack.empty()){
    current = stack.front();
    stack.pop_front();
    visited[current] = true;

    //We push all unvisited children onto the stack
    for(unsigned int i=0; i<size; i++){
      //Forward edges
			if(am[current][i] == 1 && visited[i] == false){
        stack.push_front(i);
      }
			//Backward edges
			if(am[i][current] == 1 && visited[i] == false){
				stack.push_front(i);
			}
			
    }
  }

	//Test to see if every node was visited or not	
	bool connected = true;
	for(vector<bool>::iterator it = visited.begin();
			 it != visited.end();
			 it++){

		if( *it == false ){
			connected = false;
		}
	}

return connected;
}

//Generates an adjacency matrix representation of a DAG
void gen_adj_matrix(vector<vector<unsigned int> >& am, unsigned int size ){

	//Get size of DAG to generate
	unsigned int n = size;

	//Declare an adjacency matrix and initialize it to zeroes
	am.resize(n);
	for(unsigned int i = 0; i < n; i++){
		am[i].resize(n);
	}

	for(unsigned int i = 0; i < n; i++){
		for(unsigned int j = 0; j < n; j++){
			am[i][j] = 0;
		}
	}

	#ifdef VERBOSE
	cout << "Generating directed graph of size " << n << endl;
	#endif

	unsigned int node;
	unsigned int target;
	bool graph_is_connected = false;

	//We continue to add edges (i.e. loop) until we have a DAG that weakly
	//connects all vertices
	while(!graph_is_connected){
		//		node = rand()%n;
		//		target = rand()%n;

		// Son -- begin
		node = Common::uniform_int_gen(0, n-1);
		target = Common::uniform_int_gen(0, n-1);
		// Son -- end
		
		//we've already looked at this edge
		if(am[node][target] != 0) {
			continue;
		}
		
		switch (am[node][target]){
			//Case: we haven't looked at this edge yet
		case 0:
			//Check to see if adding this edge adds a cycle
			if(acyclic_check(am,n,node,target)){
				//If no cycle is introduced, add the edge to the graph
				am[node][target] = 1;
				
				//If the graph is now weakly connected then we terminate
				if(weak_conn_test(am)){
					graph_is_connected = true;
				}	
				
			} else {
				//If a cycle is introduced, mark it as such
				am[node][target] = 2;
			}	
			break;
	
			//Case: this edge already exists
		case 1:
			break;
			
			//Case: we already looked at this edge and discarded it
		case 2:
			break;
		default:
			cout << "Found bad value in adjacency matrix" << endl;
			abort();
		}
		
	
	}
	
	/** For debugging **
		cout << "Resultant graph:" << endl;
		for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
		cout << am[i][j] << " ";
		}
		cout << endl;
		}
	*/
	
}

//A digraph is weakly connected if a DFS along the undirected analog can reach
//every node in the graph
int weak_conn_add(vector<vector<unsigned int> >& am){
	
	unsigned int size = am.size();
	vector<bool> visited (size,false);
	list<unsigned int> stack;
	stack.push_front(0);	
	unsigned int current;

	//DFS traversal with undirected edges
	while(!stack.empty()){
		current = stack.front();
		stack.pop_front();
		visited[current] = true;

		//We push all unvisited children onto the stack
		for(unsigned int i=0; i<size; i++){
			//Forward edges
			if(am[current][i] == 1 && visited[i] == false){
				stack.push_front(i);
			}
			//Backward edges
			if(am[i][current] == 1 && visited[i] == false){
				stack.push_front(i);
			}

		}
	}

	//Test to see if every node was visited or not	
	int connected = 0;
	for(unsigned int it = 0; it < size; it++){
		if( visited[it] == false ){
			connected = (int) it;
		}
	}
	return connected;
}

//Generates an adjacency matrix representation of a DAG
void gen_erdos_matrix(vector<vector<unsigned int> >& am, unsigned int size, double p ){

	//Get size of DAG to generate
	unsigned int n = size;
	
	//Declare an adjacency matrix and initialize it to zeroes
	am.resize(n);
	for(unsigned int i = 0; i < n; i++){
		am[i].resize(n);
	}

	for(unsigned int i = 0; i < n; i++){
		for(unsigned int j = 0; j < n; j++){
			am[i][j] = 0;
		}
	}

	unsigned int node;
	unsigned int target;
	int graph_is_connected;

	for(unsigned int i = 0; i < n; i++){
		for(unsigned int j = i + 1; j < n; j++){
			//			double tmp = rand()/(RAND_MAX+1.0);
			// Son -- begin
			double tmp = Common::uniform_real_gen(0, 1);
			// Son -- end

			if(tmp <= p)
			{
				//Check to see if adding this edge adds a cycle
				//if(!acyclic_check(am,n,i,j)){
				//	cout << "cycle occurs" << endl; abort();}
				am[i][j] = 1;
			}
		}
	}
	graph_is_connected = weak_conn_add(am);	
	//	if(graph_is_connected != 0){cout<<".";}
	//We continue to add edges (i.e. loop) until we have a DAG that weakly
	//connects all vertices
	while(graph_is_connected != 0){
		target = graph_is_connected;
		//		node = rand()%target;
		// Son -- begin
		node = Common::uniform_int_gen(0, target-1);
		// Son -- end

		//Check to see if adding this edge adds a cycle
		//if(acyclic_check(am,n,node,target)){
			//If no cycle is introduced, add the edge to the graph
			am[node][target] = 1;
			//If the graph is now weakly connected then we terminate
			graph_is_connected = weak_conn_add(am);
		//} else {
			//If a cycle is introduced, mark it as such
		//	cout << "cycle occurs" << endl; abort();
		//}	
	}

	/** For debugging **
		cout << "Resultant graph:" << endl;
		for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
		cout << am[i][j] << " ";
		}
		cout << endl;
		}
	*/
	
}

//Generates a set of random computational work values to go with 
//generated adjacency matrix to form a complete task DAG
//Returns the total computational load of all nodes
unsigned long gen_node_lengths(vector<unsigned long>& pa, 
							   unsigned int size, 
							   unsigned int min, 
							   unsigned int max){
	unsigned int temp;
	unsigned int C = 0;
	
	if( min == 0 ){
		cout << "Minimum period length must be positive" << endl;
		abort();
	}

	pa.resize(size);

	for(unsigned int i = 0; i < size; i++){
		//		temp = min + rand()%(max-min+1);
		// Son -- begin
		temp = Common::uniform_int_gen(min, max);
		// Son -- end

		pa[i] = temp;
		C += temp;
	}
	
	return C;
}

unsigned long gen_periods_NP(vector<unsigned long>& pa, 
							 unsigned int size, 
							 unsigned int min, 
							 unsigned int max){
	unsigned int temp;
	unsigned int C = 0;
	
	if( min == 0 ){
		cout << "Minimum period length must be positive" << endl;
		abort();
	}

	pa.resize(size);

	int NP = (int)(max/min);
	for(unsigned int i = 0; i < size; i++){
		temp = min *(1 + rand()%NP);
		pa[i] = temp;
		C += temp;
	}

	return C;
}

//Calculates the span of a DAG
unsigned long calc_span(vector<vector<unsigned int> >& am, vector<unsigned long>& pa){

	unsigned int size = pa.size();
	
	list<unsigned int> stack;
	unsigned int current;
	list<unsigned int> source_nodes;
	vector<unsigned int> dist (size, 0); 
	unsigned int span = 0;	
	
	//Find all source nodes in the graph
	//If a column in the adjacency matrix is zeroes, then that node is a source
	for(unsigned int j=0; j<size; j++){
		bool source = true;
		for(unsigned int i=0; i<size; i++){
			if(am[i][j] == 1){
				source = false;
			}
		}
		if(source){
			source_nodes.push_back(j);
		}
	}
	
	//Initialize the starting nodes to be source nodes
	for(list<unsigned int>::iterator iter = source_nodes.begin();
		 iter != source_nodes.end();
		 iter++){
		
		int a_source = *iter;
		stack.push_back(a_source);
		dist[a_source] = pa[a_source];
	}
	
	unsigned int temp;
	
	while(!stack.empty()){
		current = stack.front();
		stack.pop_front();
		
		for(unsigned int i = 0; i < size; i++){
			if(am[current][i] == 1 ){
				
				temp = dist[current] + pa[i];      
				
				//Push back the starting times of children nodes
				//This is not optimal because we process multiple
				//nodes more than once
				if(temp > dist[i]){
					dist[i] = temp;
					stack.push_back(i);
					if(span < temp){
						span = temp;
					}
				}
			}
		}
	}
	
	//	for(unsigned int i = 0; i < size; i++){
	//		cout << "Node " << i << " is dist " << dist[i] << endl;
	//	}	
	return span;
}


//Generates the DAG used as an example in Preemptive and Non-Preemptive 
//Scheduling for Parallel Tasks of DAG Model by Abu and Kunal
//Note that nodes here are numbered one less than the paper
void test_DAG(vector<vector<unsigned int> >& am, 
			  vector<unsigned long>& pa){

	am.resize(10);
	for (unsigned int i = 0; i < 10; i++){
		am[i].resize(10);
		for(unsigned int j = 0; j < 10; j++){
			am[i][j] = 0;
		}
	}	
	
	pa.resize(10);
	
	//Specify outgoing edges
	am[0][3] = 1;
	am[0][6] = 1;
	am[0][7] = 1;
	am[0][9] = 1;
	
	am[1][3] = 1;
	am[1][6] = 1;
	am[1][4] = 1;
	
	am[2][5] = 1;
	am[2][8] = 1;

	am[3][7] = 1;

	am[4][6] = 1;
	am[4][5] = 1;

	am[5][8] = 1;

	am[6][7] = 1;

	am[7][8] = 1;
	am[7][9] = 1;

	//Specify periods
	pa[0] = 4;	
	pa[1] = 2;
	pa[2] = 4;
	pa[3] = 5;
	pa[4] = 3;
	pa[5] = 4;
	pa[6] = 2;
	pa[7] = 2;
	pa[8] = 3;
	pa[9] = 3;
}

