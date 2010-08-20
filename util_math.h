/**
 * @file util_math.h
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <stdlib.h>
#include <time.h>
#include <math.h>

/**		Returns the square of the argument.
 * @param	num
 * @return num*num
 */
#define SQR(A) ((A)*(A))

/**  	Returns a random number between [0,1).
 * Use srand() beforhand.
 * @return the random number
 */
double rand1(void);
#define RAND1 (double)rand() / ((double)RAND_MAX + 1.)

/**		Returns a random number between [0, n).
 * Use srand() beforhand.
 * @param[in]	n	: the upper limit of the generated numbers.
 * @return the random number.
 */
double randn(double n);
#define RANDN(A) (A) * RAND1

/**		Returns a random number between [lower, upper).
 * Use srand() beforhand.
 * @param[in]	lower	: the lowest possible number
 * @param[in]	upper	: the upper limit of the generated numbers.
 * @return the random number.
 */
double randnk(double lower, double upper);
#define RANDNK(A,B) ((A) + (B)) * RAND1 - (A);

/**		Returns the smallest power of two no less than num.
 * @param[in]	num	: 
 * @return 
 */
long ceil_po2(double num);

/**		Returns the largest power of two not greater than num.
 * @param[in]	num	: 
 * @return 
 */
long floor_po2(double num);

/**		Rounds to power of two.
 * @param[in]	num	: 
 * @return 
 */
long round_po2(double num);

#endif	// UTIL_MATH_H
