/**
 * @file util_math.c
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#include "util_math.h"

double rand1(void) {
	return (double)rand() / ((double)RAND_MAX + 1.);
}


double randn(double n) {
	return n * (double)rand() / ((double)RAND_MAX + 1.);
}


double randnk(double lower, double upper) {
	return (lower + upper) * (double)rand() / ((double)RAND_MAX + 1.) - lower;
}

long ceil_po2(double num) {
	register double temp = log(num) / M_LN2;
#ifdef __USE_ISOC99
	return (long) exp2(ceil(temp));
#endif
	//return (long) pow(2., ceil(temp));
	return (long) exp(ceil(temp) * M_LN2);
}
