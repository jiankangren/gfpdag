#include <iostream>
#include <algorithm>
#include <map>
#include <sstream>
#include <unordered_map>
#include <omp.h>
#include "ilcplex/ilocplex.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "taskset.h"

ILOSTLBEGIN

// Construct a task set by adding tasks until it reaches the 
// desired total utilization.
TaskSet::TaskSet(double desired_util, unsigned int num_procs, 
				 unsigned long weight_min, unsigned long weight_max, double p, double beta) {
	m = num_procs;
	double sum_util = 0;
	unsigned int id = 0; // Task ID starts from 0.

	while (sum_util < desired_util) {
		Task task(id, weight_min, weight_max, p, num_procs, beta);
		double util = task.get_util();
		double remain_util = desired_util - sum_util;

		if ((util > remain_util) || (remain_util-util < 0.1)) {
			// Change the period of the task to match the remain util.
			task.set_util(remain_util);
			util = remain_util;
		}

		sum_util += util;
		id += 1;
		tasks.push_back(task);
	}
	
	// Save the number of tasks.
	n = id;

	print_taskset();
}

// Construct a task set with a given set of tasks' utilizations.
TaskSet::TaskSet(vector<double> utils, unsigned int num_procs, 
				 unsigned long weight_min, unsigned long weight_max, double p) {
	m = num_procs;
	n = utils.size();

	for (int i=0; i<n; i++) {
		Task task(i, weight_min, weight_max, p, num_procs, 0.1);
		task.set_util(utils[i]);
		tasks.push_back(task);
	}
}

// Compare two task in ascending deadline.
bool TaskSet::comp_dm(Task i, Task j) {
	return (i.get_deadline() < j.get_deadline());
}

// Assign priorities in Deadline monotonic and 
// sort the list of tasks in decreasing order of priorities.
void TaskSet::sort_dm() {
	sort(tasks.begin(), tasks.end(), comp_dm);
	//	print_taskset();
}

bool TaskSet::comp_ws(Task i, Task j) {
	unsigned long work_i = i.get_work();
	unsigned long deadline_i = i.get_deadline();
	unsigned long work_j = j.get_work();
	unsigned long deadline_j = j.get_deadline();
	unsigned int procs = i.get_procs();
	
	// Make sure the number of processors is correct.
	assert(i.get_procs() == j.get_procs());

	return ((deadline_i - work_i/procs) < (deadline_j - work_j/procs));
}

// Assign priorties using workload slack and sort tasks in decreasing priority order.
void TaskSet::sort_ws() {
	sort(tasks.begin(), tasks.end(), comp_ws);
}

bool TaskSet::comp_ls(Task i, Task j) {
	unsigned long span_i = i.get_span();
	unsigned long deadline_i = i.get_deadline();
	unsigned long span_j = j.get_span();
	unsigned long deadline_j = j.get_deadline();

	return ((deadline_i - span_i) < (deadline_j - span_j));
}

// Assign priorities using critical-path length slack and 
// sort tasks in decreasing priority order.
void TaskSet::sort_ls() {
	sort(tasks.begin(), tasks.end(), comp_ls);
}

// Perform response-time analysis with Deadline Monotonic priority order.
bool TaskSet::rta_dm() {
	sort_dm();
	return rta(OUR_RTA);
}

// Perform RTA with Workload Slack priority order.
bool TaskSet::rta_ws() {
	sort_ws();
	return rta(OUR_RTA);
}

// Perform RTA with Critical-Path Length Slack priority order.
bool TaskSet::rta_ls() {
	sort_ls();
	return rta(OUR_RTA);
}

// Default sched test is RTA with DM priority order.
bool TaskSet::is_schedulable() {
	return rta_dm();
}

// The response-time analysis by Melani et al., "Response-time analysis conditional DAG tasks in 
// multiprocessor systems", ECRTS 2015. Apply the test by Melani et al., for different priority assignments.
bool TaskSet::melani_dm() {
	sort_dm();
	return rta(MELANI_RTA);
}

bool TaskSet::melani_ws() {
	sort_ws();
	return rta(MELANI_RTA);
}

bool TaskSet::melani_ls() {
	sort_ls();
	return rta(MELANI_RTA);
}

// Response-time analysis. Assuming tasks are already sorted in decreasing priority order.
// (First task has the highest priority.)
// @type: either our method or Melani et al's method.
bool TaskSet::rta(RtaType type) {
	unsigned int size = tasks.size();
	for (int i = 0; i < size; i++) {
		unsigned int id = tasks[i].get_id();
		unsigned long work = tasks[i].get_work();
		unsigned long span = tasks[i].get_span();
		unsigned long temp = span + (work - span)/m;
		tasks[i].set_rtime(temp);
	}

	// If the response time of any task is larger than its 
	// initial deadline, then the task set is unschedulable.
	for (int i = 0; i < size; i++) {
		if (tasks[i].get_rtime() > tasks[i].get_deadline()) {
			return false;
		}
	}

	// Calculate the response times for the other tasks.
	for (int i = 1; i < size; i++) {
		unsigned long rtime = calc_rt_bound(i, type);
		if (rtime > tasks[i].get_deadline()) {
			return false;
		}
		
		// Store the response time bound because lower-priority
		// tasks need this to calculate the interfering workload 
		// caused by this task.
		tasks[i].set_rtime(rtime);
	}

	return true;
}


// Calculate the response time bound for a task with a 
// specified index in the list of tasks sorted in decreasing 
// priority order (i.e., fixed-point iteration calculation).
// @i: index of the considered task.
// @type: OUR_RTA for our method, MELANI_RTA for Melani et al's method.
unsigned long TaskSet::calc_rt_bound(unsigned int i, RtaType type) {
	Task &task = tasks[i];
	unsigned long work = task.get_work();
	unsigned long span = task.get_span();
	unsigned long deadline = task.get_deadline();

	unsigned long temp = span + (work - span)/m;
	integral_t interfering_load = 0;

	unsigned long rtime_curr = 0;
	unsigned long rtime_next = temp;

	while (rtime_curr < rtime_next) {
		// No need to calculate further if the current value of 
		// response time is larger than deadline.
		if (rtime_next > deadline) {
			return rtime_next;
		}

		rtime_curr = rtime_next;
		interfering_load = 0;
		
		// Accumulating interfering workload from higher-priority tasks.
		for (int k = 0; k < i; k++) {
			interfering_load += calc_workload(k, rtime_curr, type);
		}
		interfering_load /= m;
		rtime_next = temp + interfering_load.get_ui();
		//rtime_next = temp + interfering_load/m;
	}
	
	return rtime_curr;
}

// A general method to compute workload of a task at index k.
// @type: 1 if using our method, 2 if using Melani et al's method.
unsigned long TaskSet::calc_workload(unsigned int k, unsigned long len, RtaType type) {
	if (type == OUR_RTA) {
		return calc_workload_ours(k, len);
	} else if (type == MELANI_RTA) {
		return calc_workload_melani(k, len);
	} else {
		cout << "Incorrect type for the RTA method used !!!" << endl;
		exit(-1);
	}
}


// Calculate bound for workload generated by task with index k 
// in an interval of a specific length using our method.
unsigned long TaskSet::calc_workload_ours(unsigned int k, unsigned long len) {
	cout << "======= Calculating intefering workload from task: " << k << " with length: " << len << endl;
	Task &task = tasks[k];
	unsigned long span = task.get_span();
	unsigned long period = task.get_period();
	unsigned long deadline = task.get_deadline();
	unsigned long rtime = task.get_rtime();
	unsigned long work = task.get_work();

	// Body workload
	unsigned long body = calc_body_workload(k, len);
	if (body == 0) {
		// Must account this case separately since it means the window 
		// is too small and thus contains no body jobs. 
		// The distance between the carry-in and carry-out jobs.
		if (len <= period - rtime) {
			// The window cannot cross both carry-in and carry-out jobs.
			if (len <= deadline - rtime) {
				// The window cannot have carry-in job, only carry-out job.
				return calc_carryout_workload(k, len);
			} else {
				unsigned long co = calc_carryout_workload(k, len);
				unsigned long ci = calc_carryin_workload(k, len - (deadline-rtime));
				return max(ci, co);
			}
		} else {
			// The window can cross both carry-in and carry-out jobs.
			// So we must consider 3 cases: (i) there is only carry-in job, 
			// (ii) there is only carry-out job, and (iii) there are both.
			// Return the maximum workload among these 3 cases.
			unsigned long case1 = calc_carryin_workload(k, len - (deadline-rtime));
			unsigned long case2 = calc_carryout_workload(k, len);
			
			// Case 3
			unsigned long len_sum = len - (period - rtime);
			unsigned long len_ci = len_sum;
			unsigned long ci = 0, co = 0;
			unsigned long case3 = 0;
			while (len_ci >= 0) {
				if (len_ci > len_sum) {
					cout << "xxxxxxxxxxxx WRONG: Carry-in window length too big: " << len_ci << endl;
				}
				ci = calc_carryin_workload(k, len_ci);
				co = calc_carryout_workload(k, len_sum - len_ci);

				if (ci + co > case3) {
					case3 = ci + co;
				}
				
				if (len_ci > 0) {
					// Only decrease when length is greater than 0, otherwise it will overflow.
					len_ci--;
				} else if (len_ci == 0) {
					break;
				}
			}
			
			return max(max(case1, case2), case3);
		}
	}


	// If body workload is larger than 0.
	unsigned long carryin = 0, carryout = 0;
	unsigned long sum_ci_co = 0; // Sum of carry-in and carry-out workloads.

	// Calculate the sum of carry-in and carry-out windows.
	unsigned long len_sum = max(len - ((len-span+rtime)/period)*period + rtime, (unsigned long)0);

	if (len_sum == 0) {
		carryin = 0;
		carryout = 0;
		sum_ci_co = 0;
	} else if (len_sum >= 2*span) {
		carryin = work;
		carryout = work;
		sum_ci_co = min(2*work, len_sum*m);
	} else {
		// Compute maximum sum of carry-in and carry-out workloads.
		if (task.get_sum_cico_workload().find(len_sum) != 
			task.get_sum_cico_workload().end()) {
			// If this sum is already computed, just reuse it.
			sum_ci_co = task.get_sum_cico_workload()[len_sum];
		} else {
			// Otherwise, need to compute it and store for later reuses.
			unsigned long max_len_ci = min(len_sum, span);
			unsigned long min_len_ci = len_sum - max_len_ci;
			unsigned long len_ci = max_len_ci;
			
			while (len_ci >= min_len_ci) {
				if (len_ci > len_sum) {
					cout << "xxxxxxxxxxxx WRONG: Carry-in window length too big: " << len_ci << endl;
				}
				carryin = calc_carryin_workload(k, len_ci);
				carryout = calc_carryout_workload(k, len_sum - len_ci);
				
				if (carryin + carryout > sum_ci_co) {
					sum_ci_co = carryin + carryout;
				}
				
				// Update carry-in (and thus carry-out window) length.
				if (len_ci > 0) {
					// Only decrease when it is greater than 0, otherwise it will overflow.
					len_ci--;
				} else if (len_ci == 0) {
					break;
				}
			}
			
			// Cache the computation for the sum of carry-in and carry-out workloads.
			task.set_sum_cico_workload(len_sum, sum_ci_co);
		}
	}

	return (body + sum_ci_co);
}

// Bound for body workload generated by an interfering task.
unsigned long TaskSet::calc_body_workload(unsigned int k, unsigned long len) {
	if (len == 0) {
		return 0;
	}

	Task &task = tasks[k];
	unsigned long work = task.get_work();
	unsigned long span = task.get_span();
	unsigned long period = task.get_period();
	unsigned long rtime = task.get_rtime();

    long body = ((len - span + rtime)/period - 1)*work;

	return max(body, (long)0);
}

// Bound for carry-in workload generated by an interfering task.
unsigned long TaskSet::calc_carryin_workload(unsigned int k, unsigned long len) {
	if (len == 0) {
		return 0;
	}

	if (len >= tasks[k].get_span()) {
		return min(tasks[k].get_work(), m*len);
	}
	
	unsigned long carryin = calc_ci_co_bound(tasks[k], len, CARRY_IN);
	
	if (carryin > tasks[k].get_work()) {
		cout << "WRONG: Carry-in workload: " << carryin << " > work: " << tasks[k].get_work() << " !" << endl;
		exit(-1);
	}

	return min(carryin, m*len);
}

// Bound for carry-out workload generated by an interfering task.
unsigned long TaskSet::calc_carryout_workload(unsigned int k, unsigned long len) {
	if (len == 0) {
		return 0;
	}

 	if (len >= tasks[k].get_span()) {
		return min(tasks[k].get_work(), m*len);
	}

	unsigned long carryout = calc_ci_co_bound(tasks[k], len, CARRY_OUT);

	if (carryout > tasks[k].get_work()) {
		cout << "WRONG: Carry-out workload: " << carryout << " > work: " << tasks[k].get_work() << " !" << endl;
		exit(-1);
	}

	return min(carryout, m*len);
}

// A general method to bound carry-in or carry-out workload.
// @type: CARRY_OUT if calculating carry-out workload, 
//        CARRY_IN if calculating carry-in workload
unsigned long TaskSet::calc_ci_co_bound(Task &task, unsigned long len, WorkloadType type) {
	// If the workload is already computed and cached, just return it.
	if (type == CARRY_OUT) {
	    unsigned long cached_value = task.get_co_workload(len);
		if (cached_value > 0) {
			//cout << "=== Cache hit for carry-out ..." << endl;
			return cached_value;
		}
	} else {
		// For carry-in workload.
		unsigned long cached_value = task.get_ci_workload(len);
		if (cached_value > 0) {
			//cout << "=== Cache hit for carry-in ..." << endl;
			return cached_value;
		}
	}

	if (type == CARRY_OUT) {
		return maximize_direct_call_scip(task, len, type);
	} else {
		return maximize_carryin_fonseca(task, len);
		/*
		unsigned long ci_bound = maximize_carryin_scip(task, len);
		if (ci_bound > task.get_work()) {
			// If somehow the solver cannot solve the improved bound, we use compute and 
			// return the bound from the general method. 
			return maximize_direct_call_scip(task, len, type);
		}
		// Otherwise, we return the improved bound.
		return ci_bound;
		*/
	}
}

// Compute the maximum carry-in workload using the method by Fonseca et al.
// We agree now that this method computes a safe upper-bound for carry-in workload, 
// thus we will use this method.
unsigned long TaskSet::maximize_carryin_fonseca(Task &task, unsigned long len) {
	vector<unsigned long> nodes_work = task.get_node_works();
	vector<unsigned long> start_times = task.get_start_times();
	unsigned int size = task.get_size();
	unsigned long span = task.get_span();

	unsigned long wrkload = 0;
	for (unsigned int i = 0; i < size; i++) {
		wrkload += (unsigned long)max((long)nodes_work[i] - max((long)span - (long)start_times[i] - (long)len, (long)0), (long)0);
	}

	// Cache this value for the next time.
	task.set_ci_workload(len, wrkload);
	
	return wrkload;
}


// The old method to compute maximum workload generated by a given task in 
// a given carry-in (or carry-out window) length.
// It creates and writes the problem to a .lp file using CPLEX interfaces 
// and then call SCIP to read and solve the problem.
// Currently we use a direct method to create and solve problem using SCIP directly.
unsigned long TaskSet::maximize_indirect_call_scip(Task &task, unsigned long len, WorkloadType type) {
	vector<unsigned long> nodes_work = task.get_node_works();
	map<unsigned int, vector<vector<unsigned int> > > paths;
	
	if (type == CARRY_OUT) { // For carry-out workload
		paths = task.get_paths();
	} else { // For carry-in workload
		paths = task.get_reverse_paths();
	}

	unsigned int size = task.get_size();
	unsigned long span = task.get_span();

	SCIP_Real sol_val;
	IloEnv env;

	try {
		IloModel model(env);
		IloRangeArray cons(env);

		// Arrays of variables for subtasks. Each element corresponds to a subtask.
		IloNumVarArray x_vars(env);
		IloNumVarArray w_vars(env);
		IloNumVarArray s_vars(env);
		IloNumVarArray m_vars(env);
		IloNumVarArray a_vars(env);
		
		// Each element in this vector corresponds to an array of variables for a subtask.
		vector<IloNumVarArray> d_vars_vec;

		// Set bounds for variables.
		for (int i = 0; i < size; i++) {
			// Each subtask has a variable of each type, except the distance variables.
			// There can be more than 1 distance variables for a subtask since there 
			// can be more than 1 path from sources to this subtask.
			x_vars.add(IloNumVar(env, 0, nodes_work[i], ILOINT));
			w_vars.add(IloNumVar(env, 0, nodes_work[i], ILOINT));
			s_vars.add(IloNumVar(env, 0, span, ILOINT));
			m_vars.add(IloNumVar(env, 0, len, ILOINT));
			a_vars.add(IloNumVar(env, 0, 1, ILOINT));
			//a_vars.add(IloNumVar(env, 0, 1, ILOFLOAT)); // Try using continuous value
			
			// Add distance variables for this subtask.
			vector<vector<unsigned int> > &my_paths = paths[i];
			IloNumVarArray d_vars(env);
			for (int j = 0; j < my_paths.size(); j++) {
				d_vars.add(IloNumVar(env, 0, span, ILOINT));
			}
			d_vars_vec.push_back(d_vars);
		}
		
		// Impose a set of constraints.
		for (int i = 0; i < size; i++) {
			// Constraint 2
			model.add(w_vars[i] - x_vars[i] <= 0);

			vector<vector<unsigned int> > &my_paths = paths[i];
			IloNumVarArray &d_vars = d_vars_vec[i];

			for (int j = 0; j < my_paths.size(); j++) {
				vector<unsigned int> &path = my_paths[j];
				IloExpr expr(env);
				expr += d_vars[j];
				for (int p = 0; p < path.size(); p++) {
					unsigned int node_id = path[p];
					expr -= x_vars[node_id];
				}

				// Constraints 3 and 4
				cons.add(expr <= 0);
				cons.add(expr >= 0);
				expr.end();
				model.add(cons);

				// Constraint 5
				model.add(s_vars[i] - d_vars[j] >= 0);
			}

			// Constraint 6
			model.add(w_vars[i] - m_vars[i] <= 0);
			
			// Constraint 8
			model.add(m_vars[i] - (IloInt)len*a_vars[i] + s_vars[i]*a_vars[i] <= 0);
		}
		
		// Objective function
		IloExpr obj(env);
		for (int i = 0; i < size; i++) {
			obj += w_vars[i];
		}

		model.add(IloMaximize(env, obj));
		obj.end();

		// Now solve the maximization problem.
		IloCplex cplex(model);

		cplex.setOut(env.getNullStream());

		/*
		if (!cplex.solve()) {
			env.error() << "Failed to maximize the workload!" << endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, w_vars);
		env.out() << "Values = " << w_vars << endl;
		*/

		int tid = omp_get_thread_num();
		stringstream ss;
		ss << "tmp/thread" << tid << ".lp";
		string file_name = ss.str();
		cplex.exportModel(file_name.c_str());
		
		SCIP *scip = NULL;
		SCIP_CALL( SCIPcreate(&scip) );

		SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
		
		SCIP_CALL( SCIPsetMessagehdlr(scip, NULL) );

		SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 10) );

		SCIP_CALL( SCIPreadProb(scip, file_name.c_str(), NULL) );
		
		SCIP_CALL( SCIPsolve(scip) );

		sol_val = SCIPgetPrimalbound(scip);
		//cout << "<== INDIRECT CALL ==> Objective value: " << sol_val << endl;
		
		SCIP_CALL( SCIPfree(&scip) );

		BMScheckEmptyMemory();

	} catch (IloException &e) {
		cerr << "Ilog concert exception: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception: " << endl;
	}

	env.end();

	// Store the computed workload to its cache.
	if (type == CARRY_OUT) {
		task.set_co_workload(len, (unsigned long)sol_val);
	} else if (type == CARRY_IN) {
		task.set_ci_workload(len, (unsigned long)sol_val);
	}

	return (unsigned long)sol_val;
}

// A general method to directly call SCIP to compute maximum carry-in (or carry-out) workload.
unsigned long TaskSet::maximize_direct_call_scip(Task &task, unsigned long len, WorkloadType type) {
	//	cout << "\t\t\tCalling SCIP for length " << len << endl;
	vector<unsigned long> nodes_work = task.get_node_works();
	map<unsigned int, vector<vector<unsigned int> > > paths;

	if (type == CARRY_OUT) {
		paths = task.get_paths();
	} else {
		// CARRY_IN
		paths = task.get_reverse_paths();
	}

	unsigned int size = task.get_size();
	unsigned long span = task.get_span();

	SCIP *scip = NULL;
	SCIP_CALL( SCIPcreate(&scip) );
	SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
	SCIP_CALL( SCIPsetMessagehdlr(scip, NULL) );
	SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 10) );

	vector<SCIP_VAR*> x_vars(size, NULL);
	vector<SCIP_VAR*> w_vars(size, NULL);
	vector<SCIP_VAR*> s_vars(size, NULL);
	vector<SCIP_VAR*> m_vars(size, NULL);
	vector<SCIP_VAR*> a_vars(size, NULL);
	vector<vector<SCIP_VAR*> > d_vars(size);

	vector<SCIP_CONS*> cons_workload(size, 0);
	vector<vector<SCIP_CONS*> > cons_distance(size);
	vector<vector<SCIP_CONS*> > cons_starttime(size);
	vector<SCIP_CONS*> cons_workload2(size, 0);
	vector<SCIP_CONS*> cons_max_upper(size, 0); // upper-bound for the intermediate variables

	string prob_name = string("Maximize ") + (type == CARRY_OUT? "Carry-Out" : "Carry-In");
	SCIP_CALL( SCIPcreateProbBasic(scip, prob_name.c_str()) );

	// Default objective sense is minimizing, so we must set it to maximizing.
	SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

	// Create variables
	for (int i = 0; i < size; i++) {
		string name = "exetime_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &x_vars[i], name.c_str(), 0.0, nodes_work[i], 0.0, /*SCIP_VARTYPE_CONTINUOUS*/ SCIP_VARTYPE_INTEGER) );
		name = "workload_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &w_vars[i], name.c_str(), 0.0, nodes_work[i], 1.0, /*SCIP_VARTYPE_CONTINUOUS*/ SCIP_VARTYPE_INTEGER) );
		name = "starttime_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &s_vars[i], name.c_str(), 0.0, span, 0.0, /*SCIP_VARTYPE_CONTINUOUS*/ SCIP_VARTYPE_INTEGER) );
		name = "max_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &m_vars[i], name.c_str(), 0.0, len, 0.0, /*SCIP_VARTYPE_CONTINUOUS*/ SCIP_VARTYPE_INTEGER) );
		name = "binary_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &a_vars[i], name.c_str(), 0.0, 1.0, 0.0, /*SCIP_VARTYPE_CONTINUOUS*/ SCIP_VARTYPE_BINARY) );

		name = "distance_" + to_string(i) + "_";
		// Number of paths for this subtask.
		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_VAR *tmp;
			string d_name = name + to_string(j);
			SCIP_CALL( SCIPcreateVarBasic(scip, &tmp, d_name.c_str(), 0.0, span, 0.0, /*SCIP_VARTYPE_CONTINUOUS*/ SCIP_VARTYPE_INTEGER) );
			d_vars[i].push_back(tmp);
		}
	}

	// Add variables to the problem.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPaddVar(scip, x_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, w_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, s_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, m_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, a_vars[i]) );

		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPaddVar(scip, d_vars[i][j]) );
		}
	}
	
	// Setup constraints.
	for (int i = 0; i < size; i++) {
		// Create set of constraints 2.
		string name = "cons_workload_" + to_string(i);
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_workload[i], name.c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload[i], w_vars[i], 1.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload[i], x_vars[i], -1.0) );

		for (int j = 0; j < paths[i].size(); j++) {
			vector<unsigned int>& path = paths[i][j];
			// Create set of constraint 3 & 4 using a single SCIP_CONS variable.
			name = "cons_distance_" + to_string(i) + "_" + to_string(j);
			SCIP_CONS *tmp;
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmp, name.c_str(), 0, NULL, NULL, 0.0, 0.0) );
			for (int p = 0; p < path.size(); p++) {
				unsigned int id = path[p];
				SCIP_CALL( SCIPaddCoefLinear(scip, tmp, x_vars[id], 1.0) );
			}
			SCIP_CALL( SCIPaddCoefLinear(scip, tmp, d_vars[i][j], -1.0) );
			cons_distance[i].push_back(tmp);

			// Create set of constraints 5.
			name = "cons_starttime_" + to_string(i) + "_" + to_string(j);
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmp, name.c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
			SCIP_CALL( SCIPaddCoefLinear(scip, tmp, d_vars[i][j], 1.0) );
			SCIP_CALL( SCIPaddCoefLinear(scip, tmp, s_vars[i], -1.0) );
			cons_starttime[i].push_back(tmp);
		}

		// Create set of constraints 6.
		name = "cons_workload2_" + to_string(i);
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_workload2[i], name.c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload2[i], w_vars[i], 1.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload2[i], m_vars[i], -1.0) );

		// Create set of constraints 8.
		name = "cons_max_" + to_string(i);
		SCIP_Real licoef = 1.0;
		SCIP_Real quadcoef = 1.0;
		SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons_max_upper[i], name.c_str(), 
												1, &m_vars[i], &licoef, 
												1, &s_vars[i], &a_vars[i], &quadcoef, -SCIPinfinity(scip), 0.0) );
		// WARNING: len is unsigned long, so cannot use (SCIP_Real)(-len) for coefficient.
		SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, cons_max_upper[i], a_vars[i], -(SCIP_Real)len) );
	}
	
	// Add all constraints to the problem.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPaddCons(scip, cons_workload[i]) );
		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPaddCons(scip, cons_distance[i][j]) );
			SCIP_CALL( SCIPaddCons(scip, cons_starttime[i][j]) );
		}
		SCIP_CALL( SCIPaddCons(scip, cons_workload2[i]) );
		SCIP_CALL( SCIPaddCons(scip, cons_max_upper[i]) );
	}

	// Release variables and constraints.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPreleaseVar(scip, &x_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &w_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &s_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &m_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &a_vars[i]) );

		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPreleaseVar(scip, &d_vars[i][j]) );
		}
	}

	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPreleaseCons(scip, &cons_workload[i]) );
		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPreleaseCons(scip, &cons_distance[i][j]) );
			SCIP_CALL( SCIPreleaseCons(scip, &cons_starttime[i][j]) );
		}
		SCIP_CALL( SCIPreleaseCons(scip, &cons_workload2[i]) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons_max_upper[i]) );
	}
		
	//SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

	//	SCIP_CALL( SCIPpresolve(scip) );
	SCIP_CALL( SCIPsolve(scip) );

	/*
	if( SCIPgetNSols(scip) > 0 ) {
		SCIPinfoMessage(scip, NULL, "\nSolution:\n");
		SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
	}
	*/
	SCIP_Real val = SCIPgetPrimalbound(scip);
	//cout << "<== DIRECT CALL ==> Objective value: " << val << endl;
	//exit(1);

	SCIP_CALL( SCIPfree(&scip) );

	BMScheckEmptyMemory();

	// Write the computed workload to its cache.
	if (type == CARRY_OUT) {
		task.set_co_workload(len, (unsigned long)val);
	} else if (type == CARRY_IN) {
		task.set_ci_workload(len, (unsigned long)val);
	}

	return (unsigned long)val;
}

// Improved bound for carry-in workload.
unsigned long TaskSet::maximize_carryin_scip(Task &task, unsigned long len) {
	vector<unsigned long> nodes_work = task.get_node_works();
	map<unsigned int, vector<vector<unsigned int> > > paths;
	paths = task.get_reverse_paths();
	unsigned int size = task.get_size();
	unsigned long work = task.get_work();
	unsigned long span = task.get_span();
	vector<unsigned int> sources = task.get_sources();
	
	SCIP *scip = NULL;
	SCIP_CALL( SCIPcreate(&scip) );
	SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
	SCIP_CALL( SCIPsetMessagehdlr(scip, NULL) );
	SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 10) );

	vector<SCIP_VAR*> x_vars(size, 0);
	vector<SCIP_VAR*> w_vars(size, 0);
	vector<SCIP_VAR*> s_vars(size, 0);
	vector<SCIP_VAR*> m_vars(size, 0);
	vector<SCIP_VAR*> a_vars(size, 0);
	vector<vector<SCIP_VAR*> > d_vars(size);
	// Map from source node IDs to their corresponding variables.
	unordered_map<unsigned int, vector<SCIP_VAR*> > b_vars;
	SCIP_VAR *lci_var; // Length of carry-in window.
	SCIP_VAR *cp_var; // Actual critical path length.
	SCIP_VAR *obj_var; // Objective function.

	vector<SCIP_CONS*> cons_workload(size, 0);
	vector<vector<SCIP_CONS*> > cons_distance(size);
	vector<vector<SCIP_CONS*> > cons_starttime(size);
	vector<SCIP_CONS*> cons_workload2(size, 0);
	vector<SCIP_CONS*> cons_max_upper(size, 0);
	SCIP_CONS *cons_lci; // A constraint for carry-in window length.
	SCIP_CONS *cons_cp; // A constraint for critical path length.
	SCIP_CONS *cons_b_vars;
	SCIP_CONS *cons_obj_1;
	SCIP_CONS *cons_obj_2;

	SCIP_CALL( SCIPcreateProbBasic(scip, "Maximize Carry-In") );
	SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
	
	// Create variables.
	for (int i = 0; i < size; i++) {
		string name = "exetime_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &x_vars[i], name.c_str(), 0.0, nodes_work[i], 0.0, SCIP_VARTYPE_INTEGER) );
		name = "workload_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &w_vars[i], name.c_str(), 0.0, nodes_work[i], 0.0, SCIP_VARTYPE_INTEGER) );
		name = "starttime_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &s_vars[i], name.c_str(), 0.0, span, 0.0, SCIP_VARTYPE_INTEGER) );
		name = "max_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &m_vars[i], name.c_str(), 0.0, len, 0.0, SCIP_VARTYPE_INTEGER) );
		name = "binary_" + to_string(i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &a_vars[i], name.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

		name = "distance_" + to_string(i) + "_";
		// Number of paths for this subtask.
		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_VAR *tmp;
			string d_name = name + to_string(j);
			SCIP_CALL( SCIPcreateVarBasic(scip, &tmp, d_name.c_str(), 0.0, span, 0.0, SCIP_VARTYPE_INTEGER) );
			d_vars[i].push_back(tmp);
		}
	}

	for (unsigned int id : sources) {
		string name = "source_" + to_string(id) + "_";
		int num_paths = paths[id].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_VAR *tmp;
			string b_name = name + to_string(j);
			SCIP_CALL( SCIPcreateVarBasic(scip, &tmp, b_name.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
			b_vars[id].push_back(tmp);
		}
	}
	SCIP_CALL( SCIPcreateVarBasic(scip, &lci_var, "l_ci", 0.0, len, 0.0, SCIP_VARTYPE_INTEGER) );
	SCIP_CALL( SCIPcreateVarBasic(scip, &cp_var, "cp_len", 0.0, span, 0.0, SCIP_VARTYPE_INTEGER) );
	SCIP_CALL( SCIPcreateVarBasic(scip, &obj_var, "objval", 0.0, work, 1.0, SCIP_VARTYPE_INTEGER) );

	// Add variables to the problem.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPaddVar(scip, x_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, w_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, s_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, m_vars[i]) );
		SCIP_CALL( SCIPaddVar(scip, a_vars[i]) );

		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPaddVar(scip, d_vars[i][j]) );
		}
	}
	
	for (unsigned int id : sources) {
		for (int j = 0; j < paths[id].size(); j++) {
			SCIP_CALL( SCIPaddVar(scip, b_vars[id][j]) );
		}
	}
	SCIP_CALL( SCIPaddVar(scip, lci_var) );
	SCIP_CALL( SCIPaddVar(scip, cp_var) );
	SCIP_CALL( SCIPaddVar(scip, obj_var) );

	// Set up constraints.
	for (int i = 0; i < size; i++) {
		// Constraints: W_{i,a} - X_{i,a} <= 0.
		string name = "cons_workload_" + to_string(i);
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_workload[i], name.c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload[i], w_vars[i], 1.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload[i], x_vars[i], -1.0) );

		for (int j = 0; j < paths[i].size(); j++) {
			vector<unsigned int>& path = paths[i][j];
			// Constraints: sum(X_{i,j}) - D^p_{i,a} = 0.
			name = "cons_distance_" + to_string(i) + "_" + to_string(j);
			SCIP_CONS *tmp;
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmp, name.c_str(), 0, NULL, NULL, 0.0, 0.0) );
			for (int p = 0; p < path.size(); p++) {
				unsigned int id = path[p];
				SCIP_CALL( SCIPaddCoefLinear(scip, tmp, x_vars[id], 1.0) );
			}
			SCIP_CALL( SCIPaddCoefLinear(scip, tmp, d_vars[i][j], -1.0) );
			cons_distance[i].push_back(tmp);

			// Constraints: D^p_{i,a} - S_{i,a} <= 0.
			name = "cons_starttime_" + to_string(i) + "_" + to_string(j);
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmp, name.c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
			SCIP_CALL( SCIPaddCoefLinear(scip, tmp, d_vars[i][j], 1.0) );
			SCIP_CALL( SCIPaddCoefLinear(scip, tmp, s_vars[i], -1.0) );
			cons_starttime[i].push_back(tmp);
		}

		// Constraints: W_{i,a} - M_{i,a} <= 0.
		name = "cons_workload2_" + to_string(i);
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_workload2[i], name.c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload2[i], w_vars[i], 1.0) );
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_workload2[i], m_vars[i], -1.0) );

		// Constraints: M_{i,a} + S_{i,a}*A_{i,a} - L^CI*A_{i,a} <= 0.
		name = "cons_max_" + to_string(i);
		SCIP_Real licoef = 1.0;
		SCIP_Real quadcoef = 1.0;
		SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons_max_upper[i], name.c_str(), 1, &m_vars[i], &licoef, 1, &s_vars[i], &a_vars[i], &quadcoef, -SCIPinfinity(scip), 0.0) );
		SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons_max_upper[i], lci_var, a_vars[i], -quadcoef) );
	}
	
	// Constraint: L^CI - sum(X_{i,a})*1/m - CP_i*(m-1)/m = len - C_i/m - L_i*(m-1)/m.
	SCIP_Real rhs = (SCIP_Real)(len*m) - (SCIP_Real)work - (SCIP_Real)(span*(m-1));
	SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_lci, "cons_lci", 0, NULL, NULL, rhs, rhs) );
	SCIP_CALL( SCIPaddCoefLinear(scip, cons_lci, lci_var, (SCIP_Real)m) );
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_lci, x_vars[i], -1.0) );
	}
	SCIP_CALL( SCIPaddCoefLinear(scip, cons_lci, cp_var, -(SCIP_Real)(m-1)) );
	
	// Constraint: CP_i - sum( (D^p_{i,a} + X_{i,a})*B^p_{i,a} ) = 0, 
	// for sum over all source nodes of \tau_i.
	SCIP_Real li_coef = 1.0;
	SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons_cp, "cons_cp", 1, &cp_var, &li_coef, 0, NULL, NULL, NULL, 0.0, 0.0) );
	for (unsigned int id : sources) {		
		for (int j = 0; j < paths[id].size(); j++) {
			SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons_cp, d_vars[id][j], b_vars[id][j], -1.0) );
			SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons_cp, x_vars[id], b_vars[id][j], -1.0) );
		}
	}
	
	// Constraint: sum(B^p_{i,a}) = 1, for sum over source nodes of \tau_i.
	SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_b_vars, "cons_b_vars", 0, NULL, NULL, 1.0, 1.0) );
	for (unsigned int id : sources) {
		for (int j = 0; j < paths[id].size(); j++) {
			SCIP_CALL( SCIPaddCoefLinear(scip, cons_b_vars, b_vars[id][j], 1.0) );
		}
	}

	// Constraint: OBJ - sum(W_{i,a}) <= 0.
	SCIP_Real obj_coef1 = 1.0;
	SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_obj_1, "cons_obj_1", 1, &obj_var, &obj_coef1, -SCIPinfinity(scip), 0.0) );
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPaddCoefLinear(scip, cons_obj_1, w_vars[i], -1.0) );
	}

	// Constraint: OBJ - m*L^CI <= 0.
	SCIP_Real obj_coef2 = 1.0;
	SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons_obj_2, "cons_obj_2", 1, &obj_var, &obj_coef2, -SCIPinfinity(scip), 0.0) );
	SCIP_CALL( SCIPaddCoefLinear(scip, cons_obj_2, lci_var, -(SCIP_Real)m) );


	// Add constraints to the problem.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPaddCons(scip, cons_workload[i]) );
		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPaddCons(scip, cons_distance[i][j]) );
			SCIP_CALL( SCIPaddCons(scip, cons_starttime[i][j]) );
		}
		SCIP_CALL( SCIPaddCons(scip, cons_workload2[i]) );
		SCIP_CALL( SCIPaddCons(scip, cons_max_upper[i]) );
	}
	SCIP_CALL( SCIPaddCons(scip, cons_lci) );
	SCIP_CALL( SCIPaddCons(scip, cons_cp) );
	SCIP_CALL( SCIPaddCons(scip, cons_b_vars) );
	SCIP_CALL( SCIPaddCons(scip, cons_obj_1) );
	SCIP_CALL( SCIPaddCons(scip, cons_obj_2) );

	// Release variables.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPreleaseVar(scip, &x_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &w_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &s_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &m_vars[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &a_vars[i]) );

		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPreleaseVar(scip, &d_vars[i][j]) );
		}
	}

	for (unsigned int id : sources) {
		for (int j = 0; j < paths[id].size(); j++) {
			SCIP_CALL( SCIPreleaseVar(scip, &b_vars[id][j]) );
		}
	}
	SCIP_CALL( SCIPreleaseVar(scip, &lci_var) );
	SCIP_CALL( SCIPreleaseVar(scip, &cp_var) );
	SCIP_CALL( SCIPreleaseVar(scip, &obj_var) );

	// Release constraints.
	for (int i = 0; i < size; i++) {
		SCIP_CALL( SCIPreleaseCons(scip, &cons_workload[i]) );
		int num_paths = paths[i].size();
		for (int j = 0; j < num_paths; j++) {
			SCIP_CALL( SCIPreleaseCons(scip, &cons_distance[i][j]) );
			SCIP_CALL( SCIPreleaseCons(scip, &cons_starttime[i][j]) );
		}
		SCIP_CALL( SCIPreleaseCons(scip, &cons_workload2[i]) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons_max_upper[i]) );
	}
	SCIP_CALL( SCIPreleaseCons(scip, &cons_lci) );
	SCIP_CALL( SCIPreleaseCons(scip, &cons_cp) );
	SCIP_CALL( SCIPreleaseCons(scip, &cons_b_vars) );
	SCIP_CALL( SCIPreleaseCons(scip, &cons_obj_1) );
	SCIP_CALL( SCIPreleaseCons(scip, &cons_obj_2) );

	// Solve the problem.
	SCIP_CALL( SCIPsolve(scip) );

	SCIP_Real val = SCIPgetPrimalbound(scip);

	SCIP_CALL( SCIPfree(&scip) );
	BMScheckEmptyMemory();

	// Cache this value for the next time.
	task.set_ci_workload(len, (unsigned long)val);

	return (unsigned long)val;
}


// Calculate workload generated by a task k in an interval 
// of a given length using Melani et al's method.
unsigned long TaskSet::calc_workload_melani(unsigned int k, unsigned long len) {
	Task &task = tasks[k];
	unsigned long work = task.get_work();
	unsigned long period = task.get_period();
	unsigned long rtime = task.get_rtime();

	unsigned long body = max(((len + rtime - work/m)/period)*work, (unsigned long)0);
	unsigned long workload = body;

	long len_co = (len + rtime - work/m) % period;
	if (len_co > 0) {
		workload += min(work, (unsigned long)m*len_co);
	}

	return workload;
}

// Print out the task set for testing.
void TaskSet::print_taskset() const {
	int size = tasks.size();
	cout << "Task set has " << size << " tasks." << endl;

	for (int i = 0; i < size; i++) {
		cout << "<Task " << i << ">:";
		cout << "\tWork: " << tasks[i].get_work() 
			 << ". Span: " << tasks[i].get_span() 
			 << ". Period: " << tasks[i].get_period() 
			 << ". Deadline: " << tasks[i].get_deadline()
			 << ". Util: " << tasks[i].get_util() << endl;
	}
}

// Slack-based iterative test with different priority assignments.
bool TaskSet::slacktest_dm() {
	sort_dm();
	return slacktest();
}

bool TaskSet::slacktest_ws() {
	sort_ws();
	return slacktest();
}

bool TaskSet::slacktest_ls() {
	sort_ls();
	return slacktest();
}

// A general slack-based iterative schedulability test.
// The tasks must be already sorted in decreasing priorities.
// This follows Figure 5 in Bertogna et al., "Schedulability 
// analysis of global scheduling algorithms on multiprocessor
// platforms", TPDS 2008.
bool TaskSet::slacktest() {
	unsigned int size = tasks.size();
	for (int i = 0; i < size; i++) {
		Task &task = tasks[i];
		unsigned long work = task.get_work();
		unsigned long span = task.get_span();
		unsigned long deadline = task.get_deadline();

		if ((span + (work-span)/m) > deadline) {
			return false;
		}
	}


	bool updated = true;
	int nround = 0;
	for (int i = 0; i < size; i++) {
		tasks[i].set_slack(0);
	}

	while (updated && nround < NroundLimit) {
		bool feasible = true;
		updated = false;

		for (int i = 0; i < size; i++) {
			long slack_new = slack_compute(i);
			if (slack_new < 0) {
				feasible = false;
			} else if (slack_new > tasks[i].get_slack()) {
				tasks[i].set_slack(slack_new);
				updated = true;
			}
		}
		nround++;
		if (feasible) {
			return true;
		}
	}

	return false;
}

// Calculate slack of a task with index i.
// The tasks are already sorted in decreasing priorities.
long TaskSet::slack_compute(unsigned int i) {
	Task &task = tasks[i];
	unsigned long work = task.get_work();
	unsigned long span = task.get_span();
	unsigned long deadline = task.get_deadline();

    long slack = deadline - span - (long)ceil((work - span)*1.0/m);
	
	// Compute interfering workloads from higher-priority tasks.
	integral_t workloads = 0;
	for (int k = 0; k < i; k++) {
		workloads += calc_workload_slack(k, deadline);
	}
	// Take ceiling of workloads/m.
	integral_t tmp;
	mpz_t ncores;
	mpz_init(ncores);
	mpz_set_ui(ncores, m);
	mpz_cdiv_q(tmp.get_mpz_t(), workloads.get_mpz_t(), ncores);
	mpz_clear(ncores);
	slack -= tmp.get_ui();

	//	slack -= (int)ceil(1.0*workloads/m);
	
	return slack;
}

// Compute workload generated by body jobs of task at index k.
unsigned long TaskSet::calc_body_workload_slack(unsigned int k, unsigned long len) {
	if (len <= 0) {
		return 0;
	}

	Task &task = tasks[k];
	unsigned long work = task.get_work();
	unsigned long span = task.get_span();
	unsigned long period = task.get_period();
	unsigned long deadline = task.get_deadline();
	unsigned long slack = task.get_slack();

	// Implicit floor.
	long body = ((len-span-slack+deadline)/period - 1)*work;

	if (body <= 0) {
		return 0;
	}

	return (unsigned long)body;
}

// Compute workload generated by a task at index k during 
// an interval of a given length, with slack accounted.
unsigned long TaskSet::calc_workload_slack(unsigned int k, unsigned long len) {
	Task &task = tasks[k];
	unsigned long span = task.get_span();
	unsigned long period = task.get_period();
	unsigned long slack = task.get_slack();
	unsigned long work = task.get_work();
	unsigned long deadline = task.get_deadline();

	// Body workload.
	unsigned long body = calc_body_workload_slack(k, len);
	unsigned long carryin = 0, carryout = 0;
	unsigned long sum_ci_co = 0; // Sum of carry-in and carry-out workloads.

	// Calculate the sum of carry-in and carry-out windows.
	long len_sum = len - slack + deadline - (int)floor(1.0*(len-span-slack+deadline)/period)*period;

	if (len_sum <= 0) {
		carryin = 0;
		carryout = 0;
		sum_ci_co = 0;
	} else if (len_sum >= 2*span) {
		carryin = work;
		carryout = work;
		sum_ci_co = 2*work;
	} else {
		// Find maximum sum of carry-in and carry-out workloads.
		unsigned long init_len_ci = min((unsigned long)len_sum, span);
		unsigned long init_len_co = len_sum - init_len_ci;
		unsigned long len_ci = init_len_ci; // Got maximum length initially.
		unsigned long len_co = init_len_co; // Got minimum length initially.

		while (len_ci >= init_len_co) {
			carryin = calc_carryin_workload(k, len_ci);
			carryout = calc_carryout_workload(k, len_co);

			if (carryin + carryout > sum_ci_co) {
				sum_ci_co = carryin + carryout;
			}

			// Update carry-in and carry-out windows.
			len_ci--;
			len_co++;
		}
	}

	return (body + sum_ci_co);	
}
