#ifndef _TASK_H_
#define _TASK_H_

#include <vector>
#include <map>

using namespace std;

class Task {
 public:
	// Inputs: 
	// @id: ID of the task
	// @weight_min: minimum weight of each subtask
	// @weight_max: maximum weight of each subtask
	// @p: probability of adding an edge in its DAG
	// @m: number of processors of the system
	// @beta: minimum task utilization
	// @size: number of subtasks. If not provided, the task will generate the size itself.
	Task(unsigned int id, unsigned long weight_min, unsigned long weight_max, double p, unsigned int m, double beta);

	// Get task's basic parameters
	unsigned int get_id() const;
	unsigned long get_work() const;
	unsigned long get_span() const;
	unsigned long get_period() const;
	unsigned long get_deadline() const;
	unsigned int get_size() const; // Number of subtasks
	double get_util() const;
	unsigned int get_procs() const;

	// Get the source nodes of the graph
	vector<unsigned int> get_sources() const;
	
	// Set task's period or utilization
	void set_period(unsigned long period);
	void set_util(double util);

	// Access task's response-time
	unsigned long get_rtime() const;
	void set_rtime(unsigned long rtime);

	// Access task's slack value
	unsigned long get_slack() const;
	void set_slack(unsigned long slack);

	// Get the task's DAG information
	vector<unsigned long> get_node_works() const;
	map<unsigned int, vector<vector<unsigned int> > > get_paths() const;
	map<unsigned int, vector<vector<unsigned int> > > get_reverse_paths() const;
	vector<unsigned long> get_start_times() const;

	// Get the cached value for carry-in or carry-out workload.
	// Return -1 if the value has not been in the cache yet.
	unsigned long get_ci_workload(unsigned long len) const;
	void set_ci_workload(unsigned long len, unsigned long workload);
	unsigned long get_co_workload(unsigned long len) const;
	void set_co_workload(unsigned long len, unsigned long workload);

	// Store the sum of carry-in and carry-out workloads
	map<unsigned long, unsigned long> get_sum_cico_workload() const;
	void set_sum_cico_workload(unsigned long len, unsigned long workload);
	
 private:
	unsigned long gen_period() const;
	unsigned long gen_deadline() const;
	unsigned int gen_num_nodes() const;

	void calc_start_times();

	// Compute the adjacency matrix for the reverse DAG
	void calc_reverse_am();

	// Convert adjacency matrix representation to adjacency list
	void calc_adjacency_list(const vector<vector<unsigned int> > &am,
							 vector<vector<unsigned int> > &al);

	// Compute the lists of paths for a given DAG
	void calc_paths(const vector<vector<unsigned int> > &al, 
					const vector<vector<unsigned int> > &reverse_al, 
					map<unsigned int, vector<vector<unsigned int> > > &paths);
	void visit(unsigned int target, vector<unsigned int> &path, 
			   const vector<vector<unsigned int> > &al, 
			   map<unsigned int, vector<vector<unsigned int> > > &paths);

 private:
	unsigned int id;
	unsigned int size;
	unsigned long period;
	unsigned long deadline; // relative deadline
	unsigned long work;
	unsigned long span;
	unsigned long wmin;
	unsigned long wmax;
	unsigned int m; // Number of processors in the system
	double beta; // minimum task utilization

	unsigned long rtime; // response-time bound
	unsigned long slack; // lower bound on the slack of this task

	double p; // probability for edge in Erdos method
	double util;

	vector<unsigned int> sources; // list of ids of source nodes.

	vector<vector<unsigned int> > am; // adjacency matrix of the DAG
	vector<vector<unsigned int> > al; // adjacency list ot the DAG

	vector<unsigned long> node_works; // works of the nodes
	vector<unsigned long> start_times; // earliest start times of the nodes
	vector<vector<unsigned int> > reverse_am; // adjacency matrix of the reverse DAG
	vector<vector<unsigned int> > reverse_al; // adjacency list of the reverse DAG

	// Paths from sources to specific vertices
	map<unsigned int, vector<vector<unsigned int> > > paths; // list of paths from sources to vertices
	map<unsigned int, vector<vector<unsigned int> > > reverse_paths; // list of paths in the reverse DAG task

	// Cache the values for carry-in and carry-out workloads to be reused.
	// Use the length of the carry-in or carry-out window as index.
	// Each stores values for length in range [0, span-1].
	vector<unsigned long> ci_workloads;
	vector<unsigned long> co_workloads;
	map<unsigned long, unsigned long> sum_cico_workloads;
};

#endif
