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
		printf(
				"                         1  2  3   4   5   6   7   8   9    10	     11      12       13 14       15       16    17\n");
		printf(
				"Correct parameter order: m1 m2 S1x S1y S1z S2x S2y S2z incl f_lower f_final distance dt PNorder1 PNorder2 Spin1 Spin2\n");
		return (1);
	}

    const char *filename = argv[15];
	char PNString [50];
	sprintf(PNString, "SpinQuadTaylor%s%s", argv[14], argv[16]);
	filename = argv[15];
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	//	setting the parameters
	injParams.mass1 = atof(argv[1]);
	injParams.mass2 = atof(argv[2]);
	/*injParams.spin1x = atof(argv[3]);
	injParams.spin1y = atof(argv[4]);
	injParams.spin1z = atof(argv[5]);
	injParams.spin2x = atof(argv[6]);
	injParams.spin2y = atof(argv[7]);
	injParams.spin2z = atof(argv[8]);*/
	injParams.qmParameter1 = 1.;//atof(argv[9]);
	injParams.qmParameter1 = 1.;//atof(argv[10]);
	injParams.inclination = atof(argv[9]);
	injParams.f_lower = atof(argv[10]);
	injParams.f_final = atof(argv[11]);
	injParams.distance = atof(argv[12]);
	ppnParams.deltaT = atof(argv[13]);
	injParams.polarization = 0;

	Params par;
	srand(time(NULL));
	par.spin[0].cth = randnk(-1, 1);
	par.spin[1].cth = randnk(-1, 1);
	par.spin[0].phi = randnk(0, 2. * FITTING_PI);
	par.spin[1].phi = randnk(0, 2. * FITTING_PI);

	par.spin[0].m = atof(argv[1]);
	par.spin[1].m = atof(argv[2]);
	par.lower = 0.;
	par.upper = 1.;
	par.step = 0.01;
	par.theta = par.phi = par.pol = 0.;
	par.fp = 0.5 * (1. + par.theta * par.theta) * cos(par.phi) * cos(par.pol) -
		par.theta * sin(par.phi) * sin(par.pol);
	par.fc = 0.5 * (1. + par.theta * par.theta) * cos(par.phi) * sin(par.pol) -
		par.theta * sin(par.phi) * cos(par.pol);
	Statistic chi_stat = chi_Statistic(&injParams, &ppnParams, par);
	FILE *file = fopen(filename, "w");
	size_t i, j;
	for (i = 0; i < chi_stat.size; i++) {
		for (j = 0; j < chi_stat.size; j++) {
			fprintf(file, "%15.10lg %15.10lg %15.10lg\n", par.lower + (double)i * par.step, par.lower + (double)j * par.step, chi_stat.stat[i + j * chi_stat.size]);
		}
	}
	fclose(file);

	return 0;
}
