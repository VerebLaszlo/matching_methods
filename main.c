/**
 * @file main.c
 * @date 08/16/2010 05:38:16 PM
 * @author László Veréb
 * @brief 
 */

#include "fitting.h"
#include "util_math.h"

int main(int argc, char* argv[]) {
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	if (argc != 18) {
		short xx;
		for (xx = 0; xx <= argc; xx++) {
			printf("%s\n", argv[xx]);
		}
		printf(
				"                         1  2  3   4   5   6   7   8   9    10	     11      12       13 14      15   16       17\n");
		printf(
				"Correct parameter order: m1 m2 S1x S1y S1z S2x S2y S2z incl f_lower f_final distance dt PNorder Spin filename runs\n");
		return (1);
	}
    //const char *filename = argv[16];
	char PNString [50];
	sprintf(PNString, "SpinQuadTaylor%s%s", argv[14], argv[15]);
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	//	setting the parameters
	injParams.mass1 = atof(argv[1]);
	injParams.mass2 = atof(argv[2]);
	injParams.spin1x = atof(argv[3]);
	injParams.spin1y = atof(argv[4]);
	injParams.spin1z = atof(argv[5]);
	injParams.spin2x = atof(argv[6]);
	injParams.spin2y = atof(argv[7]);
	injParams.spin2z = atof(argv[8]);
	injParams.qmParameter1 = 1.;//atof(argv[9]);
	injParams.qmParameter2 = 1.;//atof(argv[10]);
	injParams.inclination = atof(argv[9]);
	injParams.f_lower = atof(argv[10]);
	injParams.f_final = atof(argv[11]);
	injParams.distance = atof(argv[12]);
	dt = ppnParams.deltaT = atof(argv[13]);
	injParams.polarization = 0;
	LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString);
	Params par;
	par.lower = 0.0;
	//par.lower = -1.0;
	par.upper = 1.;
	par.step = 0.005;
	par.theta = par.phi = par.pol = 0.;
	par.fp = 0.5 * (1. + SQR(par.theta)) * cos(par.phi) * cos(par.pol) - par.theta * sin(par.phi) * sin(par.pol);
	par.fc = 0.5 * (1. + SQR(par.theta)) * cos(par.phi) * sin(par.pol) - par.theta * sin(par.phi) * cos(par.pol);
	Statistic stat;
 	stat.size = (size_t)ceil((par.upper - par.lower) / par.step) + 1;
	stat.stat = malloc(stat.size * stat.size * sizeof(double));
	srand(10);
/*
	injParams.spin1x = -0.0318158;
	injParams.spin1y = 0.746334;
	injParams.spin1z = -0.286309;
	injParams.spin2x = 0.500275;
	injParams.spin2y = -0.0486518;
	injParams.spin2z = 0.622381;
*/
	par.spin[0].chi = par.spin[1].chi = 1.;
	par.spin[0].phi = 3.;
	par.spin[1].phi = 3.1;
	for (par.index = 0; par.index < atoi(argv[17]); par.index++) {
		//srand(time(NULL));
		par.spin[0].cth = RANDNK(-1, 1);
		par.spin[1].cth = RANDNK(-1, 1);
		///par.spin[0].phi = RANDNK(0, 2. * FITTING_PI);
		///par.spin[1].phi = RANDNK(0, 2. * FITTING_PI);
		//printf("%d: %lg %lg %lg %lg\n", par.index, par.spin[0].cth, par.spin[0].phi, par.spin[1].cth, par.spin[1].phi);
		//chi_Statistic(&stat, &injParams, &ppnParams, &par);
		phi_Statistic(&stat, &injParams, &ppnParams, &par);
		//cth_Statistic(&stat, &injParams, &ppnParams, &par);
	}
/*	FILE *file = fopen(filename, "w");
	size_t i, j;
	for (i = 0; i < chi_stat.size; i++) {
		for (j = 0; j < chi_stat.size; j++) {
			printf("%lg %lg %lg\n", par.lower + (double)i * par.step, par.lower + (double)j * par.step, chi_stat.stat[i + j * chi_stat.size]);
			fprintf(file, "%lg %lg %lg\n", par.lower + (double)i * par.step, par.lower + (double)j * par.step, chi_stat.stat[i + j * chi_stat.size]);
		}
		fprintf(file, "\n");
	}
	fclose(file);*/
	free(stat.stat);
	puts("Done!");
	return 0;
}
