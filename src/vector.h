#ifndef VECTOR_H
#define VECTOR_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

typedef struct VEC
{
    double * __restrict storage;
	unsigned int dim;
} VEC;

typedef struct VEC2
{
    unsigned int dim[2];
    double * __restrict storage;
} VEC2;

typedef struct VEC3
{
    unsigned int dim[3];
    double * __restrict storage;
} VEC3;

//Initialialization
VEC v_init(unsigned int dim);
VEC2 v2_init(unsigned int dim1, unsigned int dim2);
VEC3 v3_init(unsigned int dim1, unsigned int dim2, unsigned int dim3);

// Getter
double v_at(VEC v, unsigned int n);
double v2_at(VEC2 v, unsigned int m, unsigned int n);
double v3_at(VEC3 v, unsigned int m, unsigned int n, unsigned int o);

// Setter
void v_set(VEC v, unsigned int n, double value);
void v2_set(VEC2 v, unsigned int m, unsigned int n, double value);
void v3_set(VEC3 v, unsigned int m, unsigned int n, unsigned int o, double value);

// Pointer
double* v_ptr(VEC v, unsigned int n);

//Resize
void v_resize(VEC *v, unsigned int size);

//Release
void v_free(VEC *v);
void v2_free(VEC2 *v);
void v3_free(VEC3 *v);

//Slice
VEC v2_slice(VEC2 v, unsigned int n);
VEC2 v3_slice(VEC3 v, unsigned int n);
VEC v3_slice2(VEC3 v, unsigned int m, unsigned int n);

//Zeros
void v_zero(VEC v);
void v2_zero(VEC2 v);

//Copy
void v_copy(VEC src, VEC dest);
void v_copy_n(VEC src, VEC dest, unsigned int n);

//Assign
void v_assign(VEC v, const double *src, unsigned int n);

void v2_copy(VEC2 src, VEC2 dest);

#endif //VECTOR_H