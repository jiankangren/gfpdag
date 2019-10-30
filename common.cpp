#include "common.h"

using namespace std;

// Generate an integer uniformly in range [a, b]
int Common::uniform_int_gen(int a, int b) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	uniform_int_distribution<int> distribution(a, b);
	return distribution(generator);
}

// Generate a real number uniformly in range [a, b)
double Common::uniform_real_gen(double a, double b) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	uniform_real_distribution<double> distribution(a, b);
	return distribution(generator);
}

double Common::normal_gen(double mean, double stddev) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	
	normal_distribution<double> distribution(mean, stddev);
	return distribution(generator);
}
