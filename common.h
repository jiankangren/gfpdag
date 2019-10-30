#ifndef _COMMON_H_
#define _COMMON_H_

#include <random>
#include <chrono>
#include <cstddef>
#include <gmpxx.h>

// Number of rounds for the slack-based test
const int NroundLimit = 2;

enum RtaType {
	OUR_RTA = 1,
	MELANI_RTA = 2
};

enum WorkloadType {
	CARRY_OUT = 1,
	CARRY_IN = 2
};

typedef mpz_class integral_t;
typedef mpq_class fractional_t;

class Common {
 public:
	static int uniform_int_gen(int a, int b);
	static double uniform_real_gen(double a, double b);
	static double normal_gen(double mean, double stddev);
};


#endif
