#include <iostream>
#include <cstdlib>
#include <cassert>
#include <stack>
#include "task.h"
#include "rand_dag.h"
#include "common.h"

// A constant to control the utilization of the task. 
// The utilization is generated uniformly in [beta, C/L]
//const double beta = 0.1;

// The minimum and maximum number of vertices in a task
const int MinNodes = 10;
const int MaxNodes = 20;
const int DefaultNumNodes = 10; // Default is 10 nodes/DAG

Task::Task(unsigned int id, unsigned long wmin, unsigned long wmax, double p, unsigned int m, double beta) {
	this->id = id;
	this->size = gen_num_nodes();
	this->wmin = wmin;
	this->wmax = wmax;
	this->p = p;
	this->m = m;
	this->beta = beta;

	// Generate DAG task using Erdos-Renyi method
	gen_erdos_matrix(am, size, p);
	this->work = gen_node_lengths(node_works, size, wmin, wmax);
	this->span = calc_span(am, node_works);

	// Generate period
	this->period = gen_period();

	// Generate relative deadline & calculate utilization
	this->deadline = gen_deadline();
	this->util = (work*1.0)/period;
	
	// Set initial response time and slack value
	this->rtime = 0;
	this->slack = 0;

	// Init the caches for carry-in and carry-out workloads
	ci_workloads.resize(span, 0);
	co_workloads.resize(span, 0);
	// The first elements of the two vectors are 0 since they are 
	// for carry-in and carry-out window of length 0
	ci_workloads[0] = 0;
	co_workloads[0] = 0;

	// Compute the adjacency matrix of the reverse DAG
	calc_reverse_am();

	// Compute the adjacency lists representation of the two DAGs
	calc_adjacency_list(am, al);
	calc_adjacency_list(reverse_am, reverse_al);

	// Compute lists of paths for the two DAGs.
	// Must pass the adjacency lists.
	calc_paths(al, reverse_al, paths);
	calc_paths(reverse_al, al, reverse_paths);

	// Calculate the earliest start times.
	start_times.resize(size);
	calc_start_times();
}

// Period is generated uniformly in [L, C/beta]
unsigned long Task::gen_period() const {
	int min = span;
	int max = (int) (work/beta);
	
	int period = Common::uniform_int_gen(min, max);
	return period;
}

// Generate relative deadline
unsigned long Task::gen_deadline() const {
	int min = span;
	int max = period;

	/*
	double mean = (max+min)/2;
	double stddev = (max-min)/4;
	int deadline;

	do {
		deadline = (int)Common::normal_gen(mean, stddev);
	} while (deadline < min || deadline > max);
	*/

	// Deadline is generated using normal distribution with mean = period, 
	// standard deviation = (period - span)/2. 
	double mean = max;
	double stddev = (max - min)/2;
	int deadline;

	do {
		deadline = (int)Common::normal_gen(mean, stddev);
	} while (deadline < min || deadline > (2*max-min));

	if (deadline > max) {
		deadline = max - (deadline - max);
	}

	return deadline;
}

// Generate the number of nodes uniformly in range [MinNodes, MaxNodes]
unsigned int Task::gen_num_nodes() const {
	return DefaultNumNodes;
	//return Common::uniform_int_gen(MinNodes, MaxNodes);
}

unsigned int Task::get_id() const {
	return id;
}

unsigned long Task::get_work() const {
	return work;
}

unsigned long Task::get_span() const {
	return span;
}

unsigned long Task::get_period() const {
	return period;
}

unsigned long Task::get_deadline() const {
	return deadline;
}

unsigned int Task::get_size() const {
	return size;
}

unsigned int Task::get_procs() const {
	return m;
}

vector<unsigned int> Task::get_sources() const {	
	return sources;
}

double Task::get_util() const {
	return util;
}

// Access & modify response time
unsigned long Task::get_rtime() const {
	return rtime;
}

void Task::set_rtime(unsigned long rtime) {
	this->rtime = rtime;
}

// Access & modify slack value
unsigned long Task::get_slack() const {
	return slack;
}

void Task::set_slack(unsigned long slack) {
	this->slack = slack;
}

// Change task's period. Also update its utilization. 
// Must re-generate deadline after updating period.
void Task::set_period(unsigned long period) {
	this->period = period;
	this->util = (work*1.0)/period;
	this->deadline = gen_deadline();
}

// Change task's util. Also update its period.
// Must re-generate deadline after updating util.
void Task::set_util(double util) {
	this->util = util;
	this->period = (int) (work/util);
	this->deadline = gen_deadline();
}

// Calculate the earliest start times of the nodes.
void Task::calc_start_times() {
	for (int i = 0; i < size; i++) {
		if (paths.find(i) == paths.end()) {
			// Node with no paths is a source node and can start immediately.
			start_times[i] = 0;
		} else {
			// Otherwise, the earliest start time is equal to the longest 
			// distance to the node from one of the sources.
			vector<vector<unsigned int> >& my_paths = paths[i];
			unsigned long stime = 0;
			for (vector<unsigned int> path : my_paths) {
				unsigned long dist = 0;
				for (unsigned int d : path) {
					dist += d;
				}
				if (stime < dist) stime = dist;
			}
			start_times[i] = stime;
		}
	}
}

// Calculate the lists of paths to vertices, each for a vertex.
// @al: the adjacency list of the input DAG
// @reverse_al: the adjacency list of the reverse DAG
// @paths: the output lists of paths for the DAG @al
void Task::calc_paths(const vector<vector<unsigned int> > &al, 
					  const vector<vector<unsigned int> > &reverse_al, 
					  map<unsigned int, vector<vector<unsigned int> > > &paths) {
	// There is a list of paths for each target vertex
	//	paths.resize(size);
	
	// The source vertices of a DAG are the ones with 
	// zero out-degree in its reverse DAG. 
	//vector<unsigned int> sources;
	for (unsigned int i = 0; i < size; i++) {
		if (reverse_al[i].empty()) {
			sources.push_back(i);
		}
	}

	// Make sure there are at least 1 source vertices
	assert(sources.size() > 0);
	
	// Traverse the DAG from each source vertex.
	// Note that for source vertices, their list of paths are empty since 
	// there is no vertex preceding them.
	for (int i = 0; i < sources.size(); i++) {
		unsigned int source = sources[i];
		const vector<unsigned int> &neighbors = al[source];
		
		// Save the paths as we traverse the graph.
		// Its content is the same as the stack below.
		vector<unsigned int> path;
		path.push_back(source);
	
		for (int j = 0; j < neighbors.size(); j++) {
			visit(neighbors[j], path, al, paths);
		}
	}
}

// Helper method to visit vertices and add paths
void Task::visit(unsigned int target, vector<unsigned int> &path, 
				 const vector<vector<unsigned int> > &al, 
				 map<unsigned int, vector<vector<unsigned int> > > &paths) {
	paths[target].push_back(path);
	path.push_back(target);
	
	for (int i = 0; i < al[target].size(); i++) {
		visit(al[target][i], path, al, paths);
	}

	path.pop_back();
}

// Compute the adjacency matrix for the reverse DAG
void Task::calc_reverse_am() {
	reverse_am.resize(size);

	for (int i = 0; i < size; i++) {
		reverse_am[i].resize(size);
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			reverse_am[i][j] = am[j][i];
		}
	}
}

// Convert adjacency matrix to adjacency list
// @am: the complete adjacency matrix of the DAG
// @al: the desired adjacency list of the DAG
void Task::calc_adjacency_list(const vector<vector<unsigned int> > &amatrix, 
							   vector<vector<unsigned int> > &alist) {
	int n = size; // Number of vertices of the DAG
	alist.resize(n);
	
	// Each internal vector contains the list of adjacent 
	// vertices of the corresponding vertex
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (amatrix[i][j] == 1) {
				alist[i].push_back(j);
			}
		}
	}
}

// Get the WCETs of the nodes in the DAG
vector<unsigned long> Task::get_node_works() const {
	return node_works;
}

// Get the earliest start times of the nodes in relative 
// to the whole task's release time.
vector<unsigned long> Task::get_start_times() const {
	return start_times;
}

// Get the lists of paths, each for a vertex
map<unsigned int, vector<vector<unsigned int> > > Task::get_paths() const {
	return paths;
}

// The same lists for the reverse DAG
map<unsigned int, vector<vector<unsigned int> > > Task::get_reverse_paths() const {
	return reverse_paths;
}

// Caller must make sure the input @len is in [1, span-1].
// Return the corresponding carry-in workload.
unsigned long Task::get_ci_workload(unsigned long len) const {
	return ci_workloads[len];
}

void Task::set_ci_workload(unsigned long len, unsigned long workload) {
	if (len < 0 || len > (span-1)) {
		cout << "Warning: incorrect carry-in window length: " << len << endl;
		return;
	}

	ci_workloads[len] = workload;
}

// Caller must make sure the input @len is in [1, span-1].
// Return the corresponding carry-out workload.
unsigned long Task::get_co_workload(unsigned long len) const {
	return co_workloads[len];
}

void Task::set_co_workload(unsigned long len, unsigned long workload) {
	if (len < 0 || len > (span-1)) {
		cout << "Warning: incorrect carry-out window length: " << len << endl;
		return;
	}

	co_workloads[len] = workload;
}

map<unsigned long, unsigned long> Task::get_sum_cico_workload() const {
	return sum_cico_workloads;
}

void Task::set_sum_cico_workload(unsigned long len, unsigned long workload) {
	sum_cico_workloads[len] = workload;
}
