#include <stdio.h>
#include "fitting.h"
#include "util_math.h"

int main () {
#define N 200
    int i;
    Array signal[2];

    signal[0].data = calloc(sizeof(double),N);
    signal[0].length = N;

    signal[1].data = calloc(sizeof(double),N);
    signal[1].length = N;

    for (i=0; i<N; i++) {
        signal[0].data[i] = sin(i);
        signal[1].data[i] = cos(i);
    }
    calculate_Phase_Shift(signal);

    return 0;
}
