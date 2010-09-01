/**
 * @file fitting.h
 * @date 2010-07-28 15.38.59
 * @author László Veréb
 * @brief Declarations of various functions to calculate the fitting of to
 * signals.
 */

#ifndef FITTING_H
#define FITTING_H

#include "datatypes.h"
#include "util_math.h"
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

extern double dt;
extern double FITTING_PI;

/**		Calculates the phase shift between the two signals.
 * Assuming that the maximal shift is not greater than \f$\pi/2\f$
 * @param[in]	signal1	: array of the two signals
 * @return	the phase shift, -1 if error accoured
 */
double calculate_Phase_Shift(Array signal[]);
double calculate_Phase_Shift1(Array signal[]);

///	The statistic
typedef struct {
	size_t size;	///< size of the statistic
	double *stat;	///< the "two" dimensional statistic
} Statistic;

///	The spin
typedef struct {
	double chi;	///< the dimensionless spin amplitude
	double phi;	///< the \f$\phi\f$ angle of the spin
	double cth;	///< the \f$\cos\theta\f$ angle of the spin
	double x;	///< the x component of the dimensionless spin
	double y;	///< the y component of the dimensionless spin
	double z;	///< the z component of the dimensionless spin
	double m;	///< mass 
} Spins;

/// Parameters
typedef struct {
	double lower;	///< the lower boundary of the spins
	double upper;	///< the upper boundary of the spins
	double step;	///< the step of the spins
	Spins spin[2];	///< the spin parameters
	double fp;		///< the plus antenna function
	double fc;		///< the cross antenna function
	double theta;	///< angle for \f$f_p\f$, \f$f_\times\f$
	double phi;		///< angle for \f$f_p\f$, \f$f_\times\f$
	double pol;		///< polariosation angle
} Params;

/**		Returns the chi statistic of the waveform with given parameters.
 * @param[in]	params	: inspiral parameters of the waveform
 * @param[in]	pparams	: ppn parameters of the waveform
 * @param[in]	par		: other parameters
 * @return	the statistic
 */
Statistic chi_Statistic(SimInspiralTable *params, PPNParamStruc *pparams, Params par);

#endif // FITTING_H
