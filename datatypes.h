/**
 * @file datatypes.h
 * @date 2010-07-28 16.02.38
 * @author László Veréb
 * @brief Various datatype definitions and managing function declarations.
 */

#ifndef DATATYPES_H
#define DATATYPES_H

#include <stddef.h>
#include <stdlib.h>

/// Array structure
typedef struct {
	size_t length;	///< the length of the array
	double *data;	///< the data of the array
} Array;

/**		Allocates memory for the array structure.
 * @param[out]	array	: pointer to the array
 * @param[in]	length	: the length of the array
 * @return pointer to the array
 */
Array *mallocArray(Array *array, size_t length);

/**		Allocates memory for the array structure and sets it to zero.
 * @param[out]	array	: pointer to the array
 * @param[in]	length	: the length of the array
 * @return	pointer to the array
 */
Array *callocArray(Array *array, size_t length);

/**		Reallocates memory for the array structure.
 * @param[out]	array	: pointer to the array
 * @param[in]	length	: the new length of the array
 * @return pointer to the array
 */
Array *reallocArray(Array *array, size_t length);

/**		Deallocates memory of the array.
 * @param[in]	array	: pointer to the array
 */
void freeArray(Array *array);

#ifdef DEBUG
#define MallocArray(A, L)															\
	fprintf(stderr, "Allocating memory for the "#A" array with %d length.\n");		\
	mallocArray(A, L);																\
	//fprintf(stderr, "The "#A" array with %d length could not be allocated!\n", L);
	
#define CallocArray(A, L)																			\
	fprintf(stderr, "Allocating memory for the "#A" array with %d length, and setting to zero.\n");	\
	callocArray(A, L);

#define ReallocArray(A, L);															\
	fprintf(stderr, "Reallocating memory for the "#A" array with %d length.\n");	\
	reallocArray(A, L);

#define FreeArray(A);												\
	fprintf(stderr, "Deallocating memory of the "#A" array.\n");	\
	freeArray(A);

#else 
#define MallocArray(A, L) mallocArray(A, L);
#define CallocArray(A, L) callocArray(A, L);
#define ReallocArray(A, L) reallocArray(A, L);
#define FreeArray(A) freeArray(A);
#endif

#endif // DATATYPES_H
