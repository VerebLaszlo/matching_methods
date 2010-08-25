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

double dt;
double FITTING_PI = M_PI;

double calculate_Phase_Shift(Array signal[]) {
	// variable definitions and initialisations
	static const size_t diff = 1;
	size_t i, j;
    /// index[i][0] - number of maximums
    /// index[i][1] - the index of last but one maximum
    /// index[i][2] - the index of last maximum
	size_t index[][4] = {
		{0, 0, 0, 0},
		{0, 0, 0, 0}
	};

#define LOCALMAX (signal[i].data[j-diff] < signal[i].data[j] && signal[i].data[j+diff] < signal[i].data[j])
    for (i=0; i<2; i++) {
        for (j=1; j<signal[i].length-1; j++) {
            if (LOCALMAX) {
                index[i][1] = index[i][2];
                index[i][2] = j,
                index[i][0]++;
#ifdef DEBUG
                printf("MAX %d\n", j);
#endif
            }
        }
#ifdef DEBUG
        printf("\n");
#endif
    }
    printf("Number of maximums: %d, %d\n",index[0][0],index[1][0]);
    // if some of index is 0
    if (!index[0][0]*index[1][0]) {
        printf("Error - no maximum!\n");
        return -1;
    }

    // if the number of maximums are non-equal
    if (index[0][0] != index[1][0]) {
        /// TODO
        /// we should decide, what we want to do in this case
    }

    // in previous section we guarantee that the one of the index are between
    // the other indexes.
    if (index[0][2]<index[1][2] && index[0][2]>index[1][1]) {
        i=0;
    } else {
        i=1;
    }
    return FITTING_PI/2 * (index[i][2]-index[!i][1]) / (index[!i][2]-index[!i][1]);

/*
	// find the first two amplitudes
	for (i = 0; i < 2; i++) {
		for (j = 1; j < signal[i].length - 1; j++) {
			if (LOCALMAX) {
				if (!index[i][1]) {
					index[i][1] = j;
					index[i][3]++;
				} else {
					index[i][2] = j;
					index[i][3]++;
					break;
				}
			}
		}
	}
	for (i = 0; i < 2; i++) {
		printf("%d: %lg, %d: %lg, %d: %lg %d\n", index[i][0], dt * index[i][0],
				index[i][1], dt * index[i][1], index[i][2], dt * index[i][2],
				index[i][3]);
	}*/
	// find the last two amplitudes
	j++;
	//size_t temp = index[0][1];
//	while (index[0][3] == index[1][3]) {
/*		for (; j < signal[i].length - 1; j++) {
			for (i = 0; i < 2; i++) {
				if (LOCALMAX) {
					index[i][0] = index[i][1];
					index[i][1] = index[i][2];
					index[i][2] = j;
					index[i][3]++; */
					/*short xxxx;
					for (xxxx = 0; xxxx < 2; xxxx++) {
						printf("%d %d:: %d: %lg, %d: %lg, %d: %lg %d\n",signal[i].length, j, index[xxxx][0], dt * index[xxxx][0],
								index[xxxx][1], dt * index[xxxx][1], index[xxxx][2], dt * index[xxxx][2],
								index[xxxx][3]);
					}*/
/*				}
			}
			if (index[0][3] != index[1][3]) {
				break;
			}
		}
//	}

	// calculate the phase shift
#define PHASESHIFT(a,b) FITTING_PI*fabs(index[0][a]-index[1][a-b]) / ((double)(index[0][a]-index[0][a-1] + index[1][a-b]-index[1][a-b-1] )/2.)
	if (index[0][3] == index[1][3]) {
        // (02-12) / (02-01+12-11)
        return PHASESHIFT(2,0);
	} else if (index[0][3] < index[1][3]) {
        // (02-11) / (02-01+11-10)
        return PHASESHIFT(2,1);
    } else {
        // (01-12) / (01-00+12-11)
        return PHASESHIFT(1,-1);
	}*/
}

void angle_To_Component(Spins *spin) {
	short i;
	for (i = 0; i < 2; i++) {
		spin[i].x = spin[i].chi * sqrt(1. - SQR(spin[i].cth)) * cos(spin[i].phi);
		spin[i].y = spin[i].chi * sqrt(1. - SQR(spin[i].cth)) * sin(spin[i].phi);
		spin[i].z = spin[i].chi * spin[i].cth;
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
	double actual[2] = {par.lower + par.step, par.lower + par.step};
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
			printf("%d %d, %lg %lg\n", s1, s2, par.spin[0].chi, par.spin[1].chi);
			angle_To_Component(par.spin);
			for (i = 0; i < 2; i++) {
				params->spin1x = par.spin[0].x = 0.077934276;
				params->spin1y = par.spin[0].y = 0.156671629 ;
				params->spin1z = par.spin[0].z = 0.967892924;
				params->spin2x = par.spin[1].x = 0.107035774;
				params->spin2y = par.spin[1].y = -0.15506005;
				params->spin2z = par.spin[1].z = 0.882409491;
				//printf("%lg %lg %lg\n", par.spin[0].x, par.spin[0].y, par.spin[0].z);
				//printf("%lg %lg %lg\n", par.spin[1].x, par.spin[1].y, par.spin[1].z);
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
				}
			}
			/*
			FILE *out = fopen("out.out", "w");
			size_t xxx;
			for (xxx=0; xxx < signal[0].length; xxx++) {
				fprintf(out, "%lg %lg\n", xxx * dt, signal[0].data[xxx]);
			}
			fclose(out);*/
			stat.stat[s1 + stat.size * s2] = calculate_Phase_Shift(signal);
			actual[1] += par.step;
		}
		actual[0] += par.step;
		actual[1] = par.lower;
	}
	return stat;
}

