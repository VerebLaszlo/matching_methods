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

void print_SimInspiralTable(SimInspiralTable *params) {
	printf("Masses: %10lg %10lg\n", params->mass1, params->mass2);
	printf("0: x=%10lg, y=%10lg, z=%10lg\t", params->spin1x, params->spin1y, params->spin1z);
	printf("1: x=%10lg, y=%10lg, z=%10lg\n", params->spin2x, params->spin2y, params->spin2z);
	printf("Inc, f: %10lg %10lg %10lg\t", params->inclination, params->f_lower, params->f_final);
	printf("L, QM:  %10lg %10lg %10lg\n", params->distance, params->qmParameter1, params->qmParameter2);
}

double calculate_Phase_Shift1(Array signal[]) {
	short found[2] = {1, 1};
	size_t ind[2] = {2, 2};
	long index[2][3] = {{0, 0, 0}, {0, 0, 0}};
	size_t i;

	size_t xxx = 0;

#define SIGNAL(i,j) signal[i].data[j]
#define LOCALMAX ( SIGNAL(i,ind[i]-2) < SIGNAL(i,ind[i]) && SIGNAL(i,ind[i]-1) < SIGNAL(i,ind[i]) &&\
				   SIGNAL(i,ind[i]+1) < SIGNAL(i,ind[i]) && SIGNAL(i,ind[i]+2) < SIGNAL(i,ind[i]))
	while (found[0] && found[1]) {
		found[0] = found[1] = 0;	// pontosan miért is kell?
		for (i = 0; i < 2; i++) {
			for (;ind[i] < signal[i].length - 2; ind[i]++) {
				if ((found[i] = LOCALMAX)) {
					xxx++;
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
//		printf("%d ", xxx / 2);
	}
//	puts("");
	if (found[0] && !found[1]) {
		index[0][0] = index[0][1];
		index[0][1] = index[0][2];
	}
    return fabs(index[0][0] - index[1][0]) * 2. * FITTING_PI / ((double)(index[0][0] - index[0][1] + index[1][0] - index[1][1]) / 2.);
}
#undef LOCALMAX

#define SIGNAL(i,j) signal[i].data[j]
#define LOCALMAX ( SIGNAL(i,j-2) < SIGNAL(i,j) && SIGNAL(i,j-1) <= SIGNAL(i,j) &&\
				   SIGNAL(i,j+1) <= SIGNAL(i,j) && SIGNAL(i,j+2) < SIGNAL(i,j))
double calculate_Phase_Shift2(Array signal[]) {
	long num_Of_Max[2] = {0, 0};
	double phase_Shift;
	size_t i, j;
	for (i = 0; i < 2; i++) {
		for (j = 2; j < signal[i].length - 2; j++) {
			if (LOCALMAX) {
				num_Of_Max[i]++;
			}
		}
	}
	phase_Shift = FITTING_PI * (num_Of_Max[0] - num_Of_Max[1]);
	//printf("%d %d, sqt-st: %ld, %lg\n", num_Of_Max[0], num_Of_Max[1], num_Of_Max[0] - num_Of_Max[1], phase_Shift);
	return phase_Shift;
}

#undef LOCALMAX
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
/*
    for (i=0; i<2; i++) {
        printf("%d. %ld %ld\n", i, index[i][1], index[i][2]);
    }
*/	
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
	const double epsilon = 1e-14;
	double x;
	for (i = 0; i < 2; i++) {
		x = 1. - SQR(spin[i].cth);
		while (x < 0.) {
			x += epsilon;
			printf("X: %30.25lg\n", x);
		}
		spin[i].x = spin[i].chi * sqrt(1. - SQR(spin[i].cth) + epsilon) * cos(spin[i].phi);
		spin[i].y = spin[i].chi * sqrt(1. - SQR(spin[i].cth) + epsilon) * sin(spin[i].phi);
		spin[i].z = spin[i].chi * spin[i].cth;
//		puts("================================");
//		printf("%lg %lg %lg\n", spin[i].x, spin[i].y, spin[i].z);
//		printf("%30.25lg %30.25lg %30.25lg\n", spin[i].chi, spin[i].cth, spin[i].phi);
//		printf("%30.25lg %30.25lg %30.25lg\n", sqrt(x), x, spin[i].cth);
//		puts("================================");
	}
}

void make_Statistic(Statistic *stat, SimInspiralTable *params, PPNParamStruc *pparams, Params *par) {
	double *new[2];
	if (strstr(par->name, "chi")) {
		new[0] = &par->spin[0].chi;
		new[1] = &par->spin[1].chi;
	} else if (strstr(par->name, "phi")) {
		par->lower = 0.0;
		par->upper = 2. * FITTING_PI;
		par->step *= par->upper;
		new[0] = &par->spin[0].phi;
		new[1] = &par->spin[1].phi;
	} else {
		par->lower = -1.0;
		par->upper = 1.0;
		par->step *= 2.;
		new[0] = &par->spin[0].cth;
		new[1] = &par->spin[1].cth;
	}
	LALStatus status;
	CoherentGW wave[2];
	Array signal[2];
	char PNString[50];
	char filename[50];
	sprintf(filename, "%s_stat%d.txt", par->name, par->index);
	double actual[2] = {par->lower, par->lower};
	size_t i, j, s1, s2;
	double a1, a2, phi, shift;
	FILE *file = fopen(filename, "w");
	for (s1 = 0; s1 < stat->size; s1++) {
		*new[0] = actual[0];
		for (s2 = 0; s2 < stat->size; s2++) {
			*new[1] = actual[1];
			angle_To_Component(par->spin);
			params->spin1x = par->spin[0].x;
			params->spin1y = par->spin[0].y;
			params->spin1z = par->spin[0].z;
			params->spin2x = par->spin[1].x;
			params->spin2y = par->spin[1].y;
			params->spin2z = par->spin[1].z;
			for (i = 0; i < 2; i++) {
				memset(&status, 0, sizeof(LALStatus));
				memset(&wave[i], 0, sizeof(CoherentGW));
				if (!i) {
					sprintf(PNString, "SpinQuadTaylortwoPNALL");
				} else {
					sprintf(PNString, "SpinTaylortwoPN");
				}
				LALSnprintf(params->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString);
				/*
				if (s2 > 30) {
					printf("%d %d, %lg %lg\n", s1, s2, actual[0], actual[1]);fflush(stdout);
					printf("%lg %lg %lg\n", params->spin1x, params->spin1y, params->spin1z);
					printf("%lg %lg %lg\n", params->spin2x, params->spin2y, params->spin2z);
				}
				*/
				LALGenerateInspiral(&status, &wave[i], params, pparams);
				if (status.statusCode) {
					fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );fflush(stderr);
					stat->stat[s1 + stat->size * s2] = -100.;
					break;
				}
				mallocArray(&signal[i], wave[i].phi->data->length);
				for (j = 0; j < signal[i].length; j++) {
					a1  = wave[i].a->data->data[2*j];
					a2  = wave[i].a->data->data[2*j+1];
					phi     = wave[i].phi->data->data[j] - wave[i].phi->data->data[0];
					shift   = wave[i].shift->data->data[j];
					signal[i].data[j] = par->fp * (a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi)) +
										par->fc * (a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
				}
				XLALSQTPNDestroyCoherentGW(&wave[i]);
			}
			stat->stat[s1 + stat->size * s2] = calculate_Phase_Shift1(signal);
			freeArray(&signal[0]);
			freeArray(&signal[1]);
			fprintf(file, "%lg %lg %lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2]);
			if (s1 % 10 == 0 || s2 % 10 == 0) {
				printf("%lg %lg, stat= %lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2]);fflush(stdout);
			}
			actual[1] += par->step;
		}
		actual[0] += par->step;
		actual[1] = par->lower;
		fprintf(file, "\n");
	}
	fclose(file);
}

void phi_Statistic(Statistic *stat, SimInspiralTable *params, PPNParamStruc *pparams, Params *par) {
	LALStatus status;
	CoherentGW wave[2];
	Array signal[2];
	char PNString[50];
	char filename[50];
	sprintf(filename, "phi_stat%d.txt", par->index);
	double actual[2] = {0, 0};
	memset(&status, 0, sizeof(LALStatus));
	memset(&wave[0], 0, sizeof(CoherentGW));
	memset(&wave[1], 0, sizeof(CoherentGW));
	size_t i, j, s1, s2;
	double a1, a2, phi, shift;
	FILE *file = fopen(filename, "w");
	for (s1 = 0; s1 < stat->size; s1++) {
		par->spin[0].phi = actual[0];
		for (s2 = 0; s2 < stat->size; s2++) {
			par->spin[1].phi = actual[1];
			angle_To_Component(par->spin);
			params->spin1x = par->spin[0].x;
			params->spin1y = par->spin[0].y;
			params->spin1z = par->spin[0].z;
			params->spin2x = par->spin[1].x;
			params->spin2y = par->spin[1].y;
			params->spin2z = par->spin[1].z;
			for (i = 0; i < 2; i++) {
				memset(&status, 0, sizeof(LALStatus));
				memset(&wave[i], 0, sizeof(CoherentGW));
				if (!i) {
					sprintf(PNString, "SpinQuadTaylortwoPNALL");
				} else {
					sprintf(PNString, "SpinTaylortwoPN");
				}
				LALSnprintf(params->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString);
				LALGenerateInspiral(&status, &wave[i], params, pparams);
				if (status.statusCode) {
					fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );fflush(stderr);
					stat->stat[s1 + stat->size * s2] = -1.;
					break;
				}
				mallocArray(&signal[i], wave[i].phi->data->length);
				for (j = 0; j < signal[i].length; j++) {
					a1  = wave[i].a->data->data[2*j];
					a2  = wave[i].a->data->data[2*j+1];
					phi     = wave[i].phi->data->data[j] - wave[i].phi->data->data[0];
					shift   = wave[i].shift->data->data[j];
					signal[i].data[j] = par->fp * (a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi)) +
										par->fc * (a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
				}
				XLALSQTPNDestroyCoherentGW(&wave[i]);
			}
			stat->stat[s1 + stat->size * s2] = calculate_Phase_Shift2(signal);
			freeArray(&signal[0]);
			freeArray(&signal[1]);
			fprintf(file, "%lg %lg %lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2]);fflush(file);
			if (s1 % 50 == 0 && s2 % 50 == 0)
			printf("%lg %lg, stat= %lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2]);fflush(stdout);
			actual[1] += 2 * FITTING_PI * par->step;
		}
		actual[0] += 2 * FITTING_PI * par->step;
		actual[1] = par->lower;
		fprintf(file, "\n");
	}
	fclose(file);
}

void cth_Statistic(Statistic *stat, SimInspiralTable *params, PPNParamStruc *pparams, Params *par) {
	LALStatus status;
	CoherentGW wave[2];
	Array signal[2];
	char PNString[50];
	char filename[50];
	sprintf(filename, "cth_stat%d.txt", par->index);
	double actual[2] = {-1, -1};
	memset(&status, 0, sizeof(LALStatus));
	memset(&wave[0], 0, sizeof(CoherentGW));
	memset(&wave[1], 0, sizeof(CoherentGW));
	size_t i, j, s1, s2;
	double a1, a2, phi, shift;
	FILE *file = fopen(filename, "w");
	printf("%d\n", stat->size);
	for (s1 = 0; s1 < stat->size; s1++) {
		par->spin[0].cth = actual[0];
		for (s2 = 0; s2 < stat->size; s2++) {
			par->spin[1].cth = actual[1];
			angle_To_Component(par->spin);
			params->spin1x = par->spin[0].x;
			params->spin1y = par->spin[0].y;
			params->spin1z = par->spin[0].z;
			params->spin2x = par->spin[1].x;
			params->spin2y = par->spin[1].y;
			params->spin2z = par->spin[1].z;
			for (i = 0; i < 2; i++) {
				memset(&status, 0, sizeof(LALStatus));
				memset(&wave[i], 0, sizeof(CoherentGW));
				if (!i) {
					sprintf(PNString, "SpinQuadTaylortwoPNALL");
				} else {
					sprintf(PNString, "SpinTaylortwoPN");
				}
				LALSnprintf(params->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString);
				LALGenerateInspiral(&status, &wave[i], params, pparams);
				if (status.statusCode) {
					fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );
					stat->stat[s1 + stat->size * s2] = -1.;
					break;
				}
				mallocArray(&signal[i], wave[i].phi->data->length);
				for (j = 0; j < signal[i].length; j++) {
					a1  = wave[i].a->data->data[2*j];
					a2  = wave[i].a->data->data[2*j+1];
					phi     = wave[i].phi->data->data[j] - wave[i].phi->data->data[0];
					shift   = wave[i].shift->data->data[j];
					signal[i].data[j] = par->fp * (a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi)) +
										par->fc * (a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
				}
				XLALSQTPNDestroyCoherentGW(&wave[i]);
			}
			stat->stat[s1 + stat->size * s2] = calculate_Phase_Shift1(signal);
			freeArray(&signal[0]);
			freeArray(&signal[1]);
			//printf("%30.25lg %30.25lg, stat= %30.25lg, %30.25lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2], par->step);fflush(stdout);
			fprintf(file, "%lg %lg %lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2]);fflush(file);
			if (s1 % 10 == 0 && s2 % 10 == 0) {
				printf("%30.25lg %30.25lg, stat= %30.25lg, %30.25lg\n", actual[0], actual[1], stat->stat[s1 + stat->size * s2], par->step);fflush(stdout);
			}
			actual[1] += 1. * par->step;
		}
		actual[0] += 1. * par->step;
		actual[1] = par->lower;
		fprintf(file, "\n");
	}
	fclose(file);
}
