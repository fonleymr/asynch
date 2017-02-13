#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <memory.h>

#include "vector.h"


//Allocates space for a one dimentional vector.
VEC v_init(unsigned int dim)
{
    VEC v;
    if (dim == 0)
        v.storage = NULL;
    else
        v.storage = (double*)calloc(dim, sizeof(double));

    v.dim = dim;
    return v;
}

//Allocates space for a vector with size entries.
double v_at(VEC v, unsigned int n)
{
    assert(v.storage != NULL);
    assert(n < v.dim);

    return v.storage[n];
}

// Set the value at the given location
void v_set(VEC v, unsigned int n, double value)
{
    assert(v.storage != NULL);
    assert(n < v.dim);

    v.storage[n] = value;
}

//Get a pointer on the n-ith element
double* v_ptr(VEC v, unsigned int n)
{
    assert(v.storage != NULL);
    assert(n < v.dim);

    return v.storage + n;
}

//Deallocates the vector v.
void v_free(VEC *v)
{
    assert(v != NULL);
    if (v->storage)
        free(v->storage);
    memset(v, 0, sizeof(VEC));
}

//Resize space for a vector with size entries.
void v_resize(VEC *v, unsigned int size)
{
    if (v->dim == size)
        return;

    assert(v != NULL);
    if (size == 0)
        v_free(v);
    else
        v->storage = (double*)realloc(v->storage, size * sizeof(double));

    v->dim = size;
}

//Assign a buffer to the vector
void v_assign(VEC v, const double *src, unsigned int n)
{
    assert(n <= v.dim);
    memcpy(v.storage, src, n * sizeof(double));
}

//Allocates space for a two dimentional vector.
VEC2 v2_init(unsigned int dim1, unsigned int dim2)
{
    VEC2 v;
    if (dim1 == 0 || dim2 == 0)
        v.storage = NULL;
    else
        v.storage = (double*)calloc(dim1 * dim2, sizeof(double));

    v.dim[0] = dim1;
    v.dim[1] = dim2;
    return v;
}

double v2_at(VEC2 v, unsigned int m, unsigned int n)
{
    assert(v.storage != NULL);
    assert(m < v.dim[0]);
    assert(n < v.dim[1]);

    return v.storage[m * v.dim[1] + n];
}

void v2_set(VEC2 v, unsigned int m, unsigned int n, double value)
{
    assert(v.storage != NULL);
    assert(m < v.dim[0]);
    assert(n < v.dim[1]);

    v.storage[m * v.dim[1] + n] = value;
}

VEC v2_slice(VEC2 v, unsigned int n)
{
    assert(n < v.dim[0]);

    VEC res;
    res.storage = v.storage + n * v.dim[1];
    res.dim = v.dim[1];
    return res;
}

//Deallocates the vector v.
void v2_free(VEC2 *v)
{
    assert(v != NULL);
    if (v->storage)
        free(v->storage);
    memset(v, 0, sizeof(VEC2));
}

//Allocates space for a three dimentional vector.
VEC3 v3_init(unsigned int dim1, unsigned int dim2, unsigned int dim3)
{
    VEC3 v;
    if (dim1 == 0 || dim2 == 0 || dim3 == 0)
        v.storage = NULL;
    else
        v.storage = (double*)calloc(dim1 * dim2 * dim3, sizeof(double));

    v.dim[0] = dim1;
    v.dim[1] = dim2;
    v.dim[2] = dim3;
    return v;
}

double v3_at(VEC3 v, unsigned int m, unsigned int n, unsigned int o)
{
    assert(v.storage != NULL);
    assert(m < v.dim[0]);
    assert(n < v.dim[1]);
    assert(o < v.dim[2]);

    return v.storage[m * v.dim[1] * v.dim[2] + n * v.dim[2] + o];
}

void v3_set(VEC3 v, unsigned int m, unsigned int n, unsigned int o, double value)
{
    assert(v.storage != NULL);
    assert(m < v.dim[0]);
    assert(n < v.dim[1]);
    assert(o < v.dim[2]);

    v.storage[m * v.dim[1] * v.dim[2] + n * v.dim[2] + o] = value;
}

VEC2 v3_slice(VEC3 v, unsigned int n)
{
    assert(n < v.dim[0]);

    VEC2 res;
    res.storage = v.storage + n * v.dim[1] * v.dim[2];
    res.dim[0] = v.dim[1];
    res.dim[1] = v.dim[2];
    return res;
}

VEC v3_slice2(VEC3 v, unsigned int m, unsigned int n)
{
    assert(m < v.dim[0]);
    assert(n < v.dim[1]);

    VEC res;
    res.storage = v.storage + m * v.dim[1] * v.dim[2] + n * v.dim[2];
    res.dim = v.dim[2];
    return res;
}

//Deallocates the vector v.
void v3_free(VEC3 *v)
{
    assert(v != NULL);
    if (v->storage)
        free(v->storage);
    memset(v, 0, sizeof(VEC3));
}

//Fill the vector with the given value
void v_zero(VEC v)
{
    memset(v.storage, 0, v.dim * sizeof(double));
}

void v2_zero(VEC2 v)
{
    memset(v.storage, 0, v.dim[0] * v.dim[1] * sizeof(double));
}

//Copies the contents of src into dest.
void v_copy(VEC src, VEC dest)
{
    assert(dest.dim >= src.dim);
    memcpy(dest.storage, src.storage, src.dim * sizeof(double));
}

//Copies the contents of src into dest.
void v_copy_n(VEC src, VEC dest, unsigned int n)
{
    assert(n <= dest.dim);
    assert(n <= src.dim);
    memcpy(dest.storage, src.storage, n * sizeof(double));
}

//Copies the contents of src into dest.
void v2_copy(VEC2 src, VEC2 dest)
{
    assert(dest.dim[0] >= src.dim[0]);
    assert(dest.dim[1] >= src.dim[1]);
    memcpy(dest.storage, src.storage, src.dim[0] * src.dim[1] * sizeof(double));
}
