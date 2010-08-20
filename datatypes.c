/**
 * @file datatypes.c
 * @date 2010-07-28 16.10.59
 * @author László Veréb
 * @brief Varios datatype managing function definitions.
 */

#include "datatypes.h"
#include <stdio.h>

Array *mallocArray(Array *array, size_t length) {
	array->length = length;
	array->data = malloc(array->length * sizeof(double));
	return array;
}

Array *callocArray(Array *array, size_t length) {
	array->length = length;
	array->data = calloc(array->length, sizeof(double));
	return array;
}

Array *reallocArray(Array *array, size_t length) {
	array->length = length;
	array->data = realloc(array, array->length * sizeof(double));
	return array;
}

void freeArray(Array *array) {
	if(array->length) {
		array->length = 0;
		free(array->data);
		array->data = NULL;
	}
}

