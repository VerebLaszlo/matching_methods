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
#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/GenerateInspiral.h>

double dt;
double FITTING_PI = M_PI;

double calculate_Phase_Shift1(Array signal[]) {
	short found[2] = {1, 1};
	size_t ind[2] = {2, 2};
	long index[2][3] = {{0, 0, 0}, {0, 0, 0}};
	size_t i;

#define SIGNAL(i,j) signal[i].data[j]
#define LOCALMAX ( SIGNAL(i,ind[i]-2) < SIGNAL(i,ind[i]) && SIGNAL(i,ind[i]-1) < SIGNAL(i,ind[i]) &&\
				   SIGNAL(i,ind[i]+1) < SIGNAL(i,ind[i]) && SIGNAL(i,ind[i]+2) < SIGNAL(i,ind[i]))
	while (found[0] && found[1]) {
		for (i = 0; i < 2; i++) {
			for (;ind[i] < signal[i].length - 2; ind[i]++) {
				if ((found[i] = LOCALMAX)) {
					//printf("%d %d ", i, ind[i]);
					if (index[i][0]) {
						index[i][2] = index[i][1];
						index[i][1] = index[i][0];
						index[i][0] = ind[i];
					} else if (index[i][1]) {
						index[i][0] = ind[i];
					} else {
						index[i][1] = ind[i];
					}
					ind[i]++;
					break;
				}
			}
			if (!found[i]) {
				break;
			}
		}
		//puts("");
	}
#undef LOCALMAX
	if (found[0] && !found[1]) {
		index[0][0] = index[0][1];
		index[0][1] = index[0][2];
	}
	//printf("%ld %ld\n", index[0][1], index[0][0]);
	//printf("%ld %ld\n", index[1][1], index[1][0]);
    return fabs(index[0][0] - index[1][0]) * 2. * FITTING_PI / ((double)(index[0][0] - index[0][1] + index[1][0] - index[1][1]) / 2.);
}

double calculate_Phase_Shift(Array signal[]) {
	// variable definitions and initialisations
	static const size_t diff = 1;
	size_t i, j;
    /// index[i][0] - number of maximums
    /// index[i][1] - the index of last but one maximum
    /// index[i][2] - the index of last maximum
    /// index[i][3] - the latest maximum in !i
	long index[][4] = {
		{0, 0, 0, 0},
		{0, 0, 0, 0}
	};
        size_t maxlength,minlength;

        // i stores the maximum
        i = signal[0].length<signal[1].length;
        maxlength = signal[i].length;
        minlength = signal[!i].length;

#define SIGNAL(i,j) signal[i].data[j]
#define LOCALMAX ( SIGNAL(i,j-1) < SIGNAL(i,j) && SIGNAL(i,j+1) < SIGNAL(i,j) )
    for (j=1; j<minlength-1; j++) {
        for (i=0; i<2; i++) {
            if (LOCALMAX) {
                index[i][1] = index[i][2];
                index[i][2] = j;
                index[i][3] = index[!i][2];
                index[i][0]++;
#ifdef DEBUG
                printf("MAX%d %ld. %d (%ld)\n", i, index[i][0], j, index[!i][2]);
#endif
                j+=diff;
            }
        }
    }
#ifdef DEBUG
    printf("Number of maximums: %ld, %ld\n",index[0][0],index[1][0]);
    printf("Latest maximums:\n\t%ld %ld %ld\n\t%ld %ld %ld\n",
            index[0][1],index[0][2],index[0][3],
            index[1][1],index[1][2],index[1][3]);
#endif
    // if some of index is 0
    if (!(index[0][0]*index[1][0])) {
        printf("Error - no maximum!\n");
        return -1;
    }

    // if the number of maximums are non-equal
    if (abs(index[0][0]-index[1][0]) > 1) {
        // which signal has more maximum
        i = index[0][0]<index[1][0]; 
        index[i][2] = index[!i][3];
        j = index[i][2];
        while (!(LOCALMAX)) {
            j--;
        }
        index[i][1] = j;
    }

    for (i=0; i<2; i++) {
        printf("%d. %ld %ld\n", i, index[i][1], index[i][2]);
    }
    
    /// TODO
    // in previous section we guarantee that the one of the index are between
    // the other indexes.
    if (index[0][2] == index[1][2]) {
            return 0.;
    }
    if (index[0][2]<index[1][2] && index[0][2]>index[1][1]) {
        i=0;
    } else {
        i=1;
    }

    return fabs(index[i][2] - index[!i][2]) * 2. * FITTING_PI / ((double)(index[i][2] - index[i][1] + index[!i][2] - index[!i][1]) / 2.);

}

void angle_To_Component(Spins *spin) {
	short i;
	for (i = 0; i < 2; i++) {
		spin[i].x = spin[i].chi * sqrt(1. - SQR(spin[i].cth)) * cos(spin[i].phi);
		spin[i].y = spin[i].chi * sqrt(1. - SQR(spin[i].cth)) * sin(spin[i].phi);
		spin[i].z = spin[i].chi * spin[i].cth;
		printf("C: %d, x=%lg,\ty=%lg,\tz=%lg\n",i , spin[i].x, spin[i].y, spin[i].z);
	}
}

Statistic chi_Statistic(SimInspiralTable *params, PPNParamStruc *pparams, Params par) {
	LALStatus status;
	CoherentGW wave[2];
	Array signal[2];
	Statistic stat;
 	stat.size = (size_t)ceil((par.upper - par.lower) / par.step) + 1;
	stat.stat = malloc(stat.size * stat.size * sizeof(double));
	double actual[2] = {par.lower, par.lower};

	size_t s1, s2, i, j;
	double a1, a2, phi, shift;
	for (s1 = 0; s1 < stat.size; s1++) {
		par.spin[0].chi = actual[0];
		for (s2 = 0; s2 < stat.size; s2++) {
			par.spin[1].chi = actual[1];
			printf("X %lg, %lg\n", par.spin[0].chi, par.spin[1].chi);
			angle_To_Component(par.spin);
			params->spin1x = par.spin[0].x;
			params->spin1y = par.spin[0].y;
			params->spin1z = par.spin[0].z;
			params->spin2x = par.spin[1].x;
			params->spin2y = par.spin[1].y;
			params->spin2z = par.spin[1].z;
			for (i = 0; i < 2; i++) {
				memset(&status, 0, sizeof(LALStatus));
				memset(&wave[0], 0, sizeof(CoherentGW));
				memset(&wave[1], 0, sizeof(CoherentGW));
				//printf("C: %d, x=%lg,\ty=%lg,\tz=%lg\n",i , par.spin[i].x, par.spin[i].y, par.spin[i].z);
				LALGenerateInspiral(&status, &wave[i], params, pparams);
				if (status.statusCode) {
					fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );
					free(stat.stat);
					stat.size = 0;
					return stat;
				}
				mallocArray(&signal[i], wave[i].phi->data->length);
				for (j = 0; j < signal[i].length; j++) {
					a1  = wave[i].a->data->data[2*j];
					a2  = wave[i].a->data->data[2*j+1];
					phi     = wave[i].phi->data->data[j] - wave[i].phi->data->data[0];
					shift   = wave[i].shift->data->data[j];
					signal[i].data[j] = par.fp * (a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi)) +
										par.fc * (a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
					if (j > signal[i].length - 10) {
						printf("%lg ", signal[i].data[j]);
					}
				}
				XLALSQTPNDestroyCoherentGW(&wave[i]);
				puts("");
			}
			stat.stat[s1 + stat.size * s2] = calculate_Phase_Shift1(signal);
			freeArray(&signal[0]);
			freeArray(&signal[1]);
			actual[1] += par.step;
		}
		actual[0] += par.step;
		actual[1] = par.lower;
	}
	return stat;
}

