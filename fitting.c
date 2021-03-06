/**
 * @file fitting.c
 * @date 2010-07-28 15.36.24
 * @author László Veréb
 * @brief Various methods to calculate the fitting of two signals.
 */

#include "fitting.h"
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

double FITTING_PI = M_PI;

double calculate_Phase_Shift(Array signal[]) {
	// variable definitions and initialisations
	size_t i, j;
	size_t index[][4] = {
		{0, 0, 0, 0},
		{0, 0, 0, 0}
	};

#define LOCALMAX (signal[i].data[j-1] < signal[i].data[j] && signal[i].data[j+1] < signal[i].data[j])
	// find the first two amplitudes
	for (i = 0; i < 2; i++) {
		for (j = 1; j < signal[i].length - 1; j++) {
			if (LOCALMAX) {
				if (!index[i][2]) {
					index[i][2] = j;
					index[i][3]++;
				} else {
					index[i][1] = j;
					index[i][3]++;
				}
				break;
			}
		}
	}

	// find the last two amplitudes
	while (index[i][3] == index[i][3]) {
		for (i = 0; i < 2; i++) {
			for (; j < signal[i].length - 1; j++) {
				if (LOCALMAX) {
					index[i][0] = index[i][1];
					index[i][1] = index[i][2];
					index[i][2] = j;
					index[i][3]++;
					break;
				}
			}
		}
	}

	// calculate the phase shift
#define PHASESHIFT(a,b) FITTING_PI*fabs(index[0][a]-index[1][a-b]) / ((double)(index[0][a]-index[0][a-1] + index[1][a-b]-index[1][a-b-1] )/2.)
	if (index[0][3] == index[1][3]) {
        return PHASESHIFT(2,0);
	} else if (index[0][3] < index[1][3]) {
        return PHASESHIFT(2,1);
    } else {
        return PHASESHIFT(1,-1);
	}
}

void angle_To_Component(Spins *spin) {
	short i;
	for (i = 0; i < 2; i++) {
		spin->x = spin->chi * sqrt(1. - SQR(spin->cth)) * cos(spin->phi);
		spin->y = spin->chi * sqrt(1. - SQR(spin->cth)) * sin(spin->phi);
		spin->z = spin->chi * spin->cth;
	}
}

Statistic chi_Statistic(SimInspiralTable *params, 
		PPNParamStruc *pparams,
		Params par) {
	LALStatus status;
	CoherentGW wave[2];
	memset(&status, 0, sizeof(LALStatus));
	memset(&wave[0], 0, sizeof(CoherentGW));
	memset(&wave[1], 0, sizeof(CoherentGW));
	Array signal[2];
	Statistic stat;
	stat.size = (size_t)ceil((par.upper - par.lower) / par.step);
	stat.stat = malloc(stat.size * stat.size * sizeof(double));
	double actual[2] = {par.lower, par.lower};
	/*double fp = 0.5 * (1. + par.theta * par.theta) * cos(par.phi) * cos(par.pol) -
		par.theta * sin(par.phi) * sin(par.pol);
	double fc = 0.5 * (1. + par.theta * par.theta) * cos(par.phi) * sin(par.pol) -
		par.theta * sin(par.phi) * cos(par.pol);*/

	size_t s1, s2, i, j;
	double a1, a2, phi, shift;
	for (s1 = 0; s1 < stat.size; s1++) {
		par.spin[0].chi = actual[0];
		for (s2 = 0; s2 < stat.size; s2++) {
			par.spin[1].chi = actual[1];
			for (i = 0; i < 2; i++) {
				angle_To_Component(&par.spin[i]);
				params->spin1x = par.spin[0].x;
				params->spin1y = par.spin[0].y;
				params->spin1z = par.spin[0].z;
				params->spin2x = par.spin[1].x;
				params->spin2y = par.spin[1].y;
				params->spin2z = par.spin[1].z;
				LALGenerateInspiral(&status, &wave[i], params, pparams);
				if (status.statusCode) {
					fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );
					free(stat.stat);
					stat.size = 0;
					return stat;
				}
				mallocArray(&signal[i], wave[i].f->data->length);
				for (j = 0; j < signal[i].length; j++) {
					a1  = wave[i].a->data->data[2*j];
					a2  = wave[i].a->data->data[2*j+1];
					phi     = wave[i].phi->data->data[j] - wave[i].phi->data->data[0];
					shift   = wave[i].shift->data->data[j];
					signal[i].data[j] = par.fp * (a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi)) +
										par.fc * (a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
				}
			}
			stat.stat[s1 + stat.size * s2] = calculate_Phase_Shift(signal);
			actual[1] += par.step;
		}
		actual[0] += par.step;
		actual[1] = par.lower;
	}
	return stat;
}

