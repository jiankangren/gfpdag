#ifndef __RAND_DAG
#define __RAND_DAG

#include <vector>

using namespace std;

// Function prototypes
bool acyclic_check(vector<vector<unsigned int> > &am, unsigned int size,
				   unsigned int start, unsigned int end);

bool weak_conn_test(vector<vector<unsigned int> > &am);

void gen_adj_matrix(vector<vector<unsigned int> > &am, unsigned int size);

int weak_conn_add(vector<vector<unsigned int> > &am);

void gen_erdos_matrix(vector<vector<unsigned int> > &am, unsigned int size, 
					  double p);

unsigned long gen_node_lengths(vector<unsigned long> &pa, unsigned int size,
							   unsigned int min, unsigned int max);

unsigned long gen_periods_NP(vector<unsigned long> &pa, unsigned int size, 
							 unsigned int min, unsigned int max);

unsigned long calc_span(vector<vector<unsigned int> > &am, 
						vector<unsigned long> &pa);

void test_DAG(vector<vector<unsigned int> > &am, vector<unsigned long> &pa);

#endif
