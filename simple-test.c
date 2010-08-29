#include <stdio.h>
#include "fitting.h"
#include "util_math.h"

int main () {
    int N[] = {200,200};
    int i;
    Array signal[2];

    for (i=0; i<2; i++) {
        signal[i].data = calloc(sizeof(double),N[i]);
        signal[i].length = N[i];
    }

    for (i=0; i<N[0]; i++) {
        signal[0].data[i] = sin(10*FITTING_PI/N[0]*i);
    }
	FILE *file = fopen("x.txt", "w");
    for (i=0; i<N[1]; i++) {
        //signal[1].data[i] = cos(10*FITTING_PI/N[1]*i);
        //signal[1].data[i] = sin(10*FITTING_PI/N[0]*i);
        //signal[1].data[i] = sin(9*FITTING_PI/N[0]*i);
        signal[1].data[i] = sin(7*FITTING_PI/N[0]*i);
		fprintf(file, "%d %lg %lg %lg\n", i, 10*FITTING_PI/N[0]*i, signal[0].data[i], signal[1].data[i]);
    }
	fclose(file);

    printf("Phase shift: %.4f\n",calculate_Phase_Shift(signal));
    printf("Phase shift: %.4f\n",calculate_Phase_Shift1(signal));

    return 0;
}
