#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <omp.h>
#include <unistd.h>
#include <sched.h>
#include "taskset.h"

using namespace std;

// Min and max WCET for each subtask
const int WcetMin = 1;
const int WcetMax = 100;

// Three types of experiment
int vary_utilization(int argc, char *argv[]);
int vary_num_cores(int argc, char *argv[]);
int vary_num_tasks();
int test();

int m; // Number of processors
int NTaskSets; // Number of task sets
double Beta; // Minimum task utilization
double EdgeProb; // Probability of an edge between 2 nodes
int NThreads; // Number of OpenMP threads for this experiments
double util_frac; // Normalized total utilization

// @first_core and @last_core are the 5th and 6th arguments, respectively.
int main(int argc, char *argv[]) {
	int first_core, last_core;
	stringstream ss;
	ss << argv[5] << " " << argv[6];
	ss >> first_core; // First core on which this experiment runs
	ss >> last_core; // Last core on which this experiment runs
	//	cout << "First core: " << first_core << ", last core: " << last_core << endl;

	// Bind the experiment to the given cores
	cpu_set_t mask;
	CPU_ZERO(&mask);
	for (unsigned i = first_core; i <= last_core; i++) {
		CPU_SET(i, &mask);
	}
	
	int ret_val = sched_setaffinity(getpid(), sizeof(mask), &mask);
	if (ret_val != 0) {
		cout << "WARNING: Could not bind to all cores!" << endl;
	}

	return vary_utilization(argc, argv);
	//	return vary_num_cores(argc, argv);
}

// Experiment with fixed number of processors and varying total utilization.
int vary_utilization(int argc, char *argv[]) {
	if (argc != 7) {
		cout << "Usage: <Vary Util> M NTaskSets Beta EdgeProb FirstCore LastCore" << endl;
		exit(1);
	}

	int first_core, last_core;
	stringstream ss;
	ss << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6];
	ss >> m;
	ss >> NTaskSets;
	ss >> Beta;
	ss >> EdgeProb;
	ss >> first_core; // First core on which this experiment runs
	ss >> last_core; // Last core on which this experiment runs
	NThreads = last_core - first_core + 1;
	//	cout << "m: " << m << ", #tasksets: " << NTaskSets << ", beta: " << Beta << ", edge prob: " << EdgeProb 
	//		 << ", first core: " << first_core << ", last core: " << last_core << endl;
	//	exit(1);

	//	double utils[] = {4.0};
	double utils[] = {10.0, 11.0, 12.0, 13.0, 14.0};
	int size = sizeof(utils)/sizeof(double);

	// Each contains the numbers of schedulable task sets for each method
	vector<int> our_dm(size, 0); // Our method with DM priority
	vector<int> melani_dm(size, 0); // Melani's method with DM priority
	
	// Count number of task sets that schedulable with Melani's RTA but 
	// unschedulable with our RTA.
	vector<int> diff(size, 0);

	// Set number of OpenMP threads to the number of cores
	omp_set_dynamic(0);
	omp_set_num_threads(NThreads);

	// For each total utilization value
	for (int i = 0; i < size; i++) {
		double util = utils[i];
		
		// Count number of finished task sets
		int count_finished = 0;

		// For each generated task set
#pragma omp parallel for
		for (int j = 0; j < NTaskSets; j++) {
			TaskSet taskset(util, m, WcetMin, WcetMax, EdgeProb, Beta);
			//taskset.print_taskset();
			
			bool is_rta_dm = taskset.rta_dm();
			bool is_melani_dm = taskset.melani_dm();

			// Atomically update the shared counters
			#pragma omp parallel 
			{
				if (is_rta_dm) {our_dm[i]++;}
				if (is_melani_dm) {melani_dm[i]++;}
				if (is_melani_dm && !is_rta_dm) {diff[i]++;}

				count_finished++;
				// Write the intermediate results
				stringstream ss;
				ss << "m=" << m << "util=" << util << "vary_util.txt"; 
				ofstream ofile(ss.str().c_str());
				
				if (ofile.is_open()) {
					ofile << "Our RTA: " << our_dm[i] << "/" << count_finished 
						  << ". Melani RTA: " << melani_dm[i] << "/" << count_finished << "\n";
				}
				ofile.close();
			}
		}

		// Write the overall ratio for this utilization
		stringstream ss;
		ss << "m=" << m << "util=" << util << "vary_util.txt"; 
		ofstream ofile(ss.str().c_str());
		if (ofile.is_open()) {
			ofile << "Our RTA: " << 1.0*our_dm[i]/NTaskSets << ". Melani RTA: " 
				  << 1.0*melani_dm[i]/NTaskSets << ". Diff: " << diff[i] << "\n";
		}

		ofile.close();
	}

	//	vector<double> our_dm_ratio(size, 0);
	//	vector<double> melani_dm_ratio(size, 0);

	//	for (int i = 0; i < size; i++) {
	//		our_dm_ratio[i] = 1.0*our_dm[i]/NTaskSets;
	//		melani_dm_ratio[i] = 1.0*melani_dm[i]/NTaskSets;
	//	}
}

// For experiment with varying number of cores.
int vary_num_cores(int argc, char *argv[]) {
	if (argc != 7) {
		cout << "Usage: <Vary Cores> Util NTaskSets Beta EdgeProb FirstCore LastCore" << endl;
		exit(1);
	}

	int first_core, last_core;
	stringstream ss;
	ss << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6];
	ss >> util_frac;
	ss >> NTaskSets;
	ss >> Beta;
	ss >> EdgeProb;
	ss >> first_core; // First core on which this experiment runs
	ss >> last_core; // Last core on which this experiment runs
	NThreads = last_core - first_core + 1;
	//	cout << "Util fraction: " << util_frac << ", #tasksets: " << NTaskSets << ", beta: " << Beta << ", edge prob: " << EdgeProb 
	//	 << ", first core: " << first_core << ", last core: " << last_core << endl;
	//exit(1);

	int nprocs[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	int size = sizeof(nprocs)/sizeof(int);

	// Each contains the numbers of schedulable task sets for each method
	vector<int> our_dm(size, 0); // Our method with DM priority
	vector<int> melani_dm(size, 0); // Melani's method with DM priority
	vector<int> diff(size, 0);

	// Set number of OpenMP threads to the number of cores
	omp_set_dynamic(0);
	omp_set_num_threads(NThreads);

	// For each total utilization value
	for (int i = 0; i < size; i++) {
		int nproc = nprocs[i];
		double util = util_frac * nproc;

		int count_finished = 0;

		// For each generated task set
#pragma omp parallel for
		for (int j = 0; j < NTaskSets; j++) {
			TaskSet taskset(util, nproc, WcetMin, WcetMax, EdgeProb, Beta);
			//taskset.print_taskset();

			bool is_rta_dm = taskset.rta_dm();
			bool is_melani_dm = taskset.melani_dm();

			#pragma omp parallel 
			{
				if (is_rta_dm) {our_dm[i]++;}
				if (is_melani_dm) {melani_dm[i]++;}
				if (is_melani_dm && !is_rta_dm) {diff[i]++;}

				count_finished++;
				// Write the intermediate results
				stringstream ss;
				ss << "util=" << util << "m=" << nproc << "vary_procs.txt"; 
				ofstream ofile(ss.str().c_str());
				
				if (ofile.is_open()) {
					ofile << "Our RTA: " << our_dm[i] << "/" << count_finished 
						  << ". Melani RTA: " << melani_dm[i] << "/" << count_finished << "\n";
				}
				ofile.close();
			}
		}

		// Write the overall ratio for this utilization
		stringstream ss;
		ss << "util=" << util << "m=" << nproc << "vary_procs.txt"; 
		ofstream ofile(ss.str().c_str());
		if (ofile.is_open()) {
			ofile << "Our RTA: " << 1.0*our_dm[i]/NTaskSets << ". Melani RTA: " 
				  << 1.0*melani_dm[i]/NTaskSets << ". Diff: " << diff[i] << "\n";
		}

		ofile.close();
	}

	//	vector<double> our_dm_ratio(size, 0);
	//	vector<double> melani_dm_ratio(size, 0);

	//	for (int i = 0; i < size; i++) {
	//		our_dm_ratio[i] = 1.0*our_dm[i]/NTaskSets;
	//		melani_dm_ratio[i] = 1.0*melani_dm[i]/NTaskSets;
	//	}

	/*
	ofstream ofile("results_vary_procs.txt");
	if (ofile.is_open()) {
		ofile << "#Procs: 4\t6\t8\t10\t12\t14\t16\t18\t20\n";
		ofile << "OurDM: ";
		for (int i = 0; i < size; i++) {
			ofile << our_dm_ratio[i] << "\t";
		}
		ofile << "\n";

		ofile << "Melani: ";
		for (int i = 0; i < size; i++) {
			ofile << melani_dm_ratio[i] << "\t";
		}
		ofile << "\n";
	}
	ofile.close();
	*/
}

// Fix total utilization: 4.0, number of cores: 8.
// For each number of tasks, read individual utilizations from a file.
int vary_num_tasks(int argc, char *argv[]) {
	int ntasks[] = {4};//{2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
	double total_util = 4.0;
	int nprocs = 8;
	int size = sizeof(ntasks)/sizeof(int);

	// Each contains the numbers of schedulable task sets for each method
	vector<int> our_dm(size, 0); // Our method with DM priority
	vector<int> melani_dm(size, 0); // Melani's method with DM priority

	vector<int> diff(size, 0);

	// Set number of OpenMP threads to the number of cores
	omp_set_dynamic(0);
	omp_set_num_threads(NThreads);

	// For each total utilization value
	for (int i = 0; i < size; i++) {
		int ntask = ntasks[i];
		// Read all utilizations from a file and pass it to the task set constructor
		vector<vector<double> > utils_vec;
		stringstream tmp;
		tmp << "data/" << ntask << ".txt";
		ifstream ifile(tmp.str().c_str());
		if (ifile.is_open()) {
			while (!ifile.eof()) {
				string line;
				getline(ifile, line);
				stringstream line_ss(line);
				vector<double> utils;
				for (int i=0; i<ntask; i++) {
					double util;
					line_ss >> util;
					utils.push_back(util);
				}
				utils_vec.push_back(utils);
			}
		}
		ifile.close();
		int count_finished = 0;

		// For each generated task set
#pragma omp parallel for
		for (int j = 0; j < NTaskSets; j++) {
			for (int k=0; k<utils_vec[j].size(); k++) {
				cout << utils_vec[j][k] << " ";
			}
			cout << endl;

			TaskSet taskset(utils_vec[j], nprocs, WcetMin, WcetMax, EdgeProb);
			//taskset.print_taskset();

			bool is_rta_dm = taskset.rta_dm();
			bool is_melani_dm = taskset.melani_dm();

			#pragma omp parallel 
			{
				if (is_rta_dm) {our_dm[i]++;}
				if (is_melani_dm) {melani_dm[i]++;}
				if (is_melani_dm && !is_rta_dm) {diff[i]++;}

				count_finished++;
				// Write the intermediate results
				stringstream ss;
				ss << "n=" << ntask << "util=" << total_util << "m=" << nprocs << "vary_ntasks.txt"; 
				ofstream ofile(ss.str().c_str());
				
				if (ofile.is_open()) {
					ofile << "Our RTA: " << our_dm[i] << "/" << count_finished 
						  << ". Melani RTA: " << melani_dm[i] << "/" << count_finished << "\n";
				}
				ofile.close();
			}
		}

		// Write the overall ratio for this utilization
		stringstream ss;
		ss << "n=" << ntask << "util=" << total_util << "m=" << nprocs << "vary_ntasks.txt"; 
		ofstream ofile(ss.str().c_str());
		if (ofile.is_open()) {
			ofile << "Our RTA: " << 1.0*our_dm[i]/NTaskSets << ". Melani RTA: " 
				  << 1.0*melani_dm[i]/NTaskSets << ". Diff: " << diff[i] << "\n";
		}

		ofile.close();
	}

	vector<double> our_dm_ratio(size, 0);
	vector<double> melani_dm_ratio(size, 0);

	for (int i = 0; i < size; i++) {
		our_dm_ratio[i] = 1.0*our_dm[i]/NTaskSets;
		melani_dm_ratio[i] = 1.0*melani_dm[i]/NTaskSets;
	}

	/*
	ofstream ofile("results_vary_ntasks.txt");
	if (ofile.is_open()) {
		ofile << "#Procs: 2\t4\t6\t8\t10\t12\t14\t16\t18\t20\n";
		ofile << "OurDM: ";
		for (int i = 0; i < size; i++) {
			ofile << our_dm_ratio[i] << "\t";
		}
		ofile << "\n";

		ofile << "Melani: ";
		for (int i = 0; i < size; i++) {
			ofile << melani_dm_ratio[i] << "\t";
		}
		ofile << "\n";
	}
	ofile.close();	
	*/
}

int test() {
	TaskSet taskset(2.5, 4, 1, 100, 0.2, Beta);

	taskset.print_taskset();

	if (taskset.rta_dm()) {
		cout << "Task set is DM schedulable !!" << endl;
	} else {
		cout << "Task set is DM unschedulable !!" << endl;
	}

	if (taskset.melani_dm()) {
		cout << "Task set is MELANI_DM schedulable !!" << endl;
	} else {
		cout << "Task set is MELANI_DM unschedulable !!" << endl;
	}
	return 0;
}
