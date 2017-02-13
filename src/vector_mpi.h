#ifndef VECTOR_MPI_H
#define VECTOR_MPI_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <mpi.h>

#include "vector.h"

//MPI
static inline
int MPI_Send_Vec(VEC v, int dest, int tag, MPI_Comm comm)
{
    return MPI_Send(v.storage, v.dim, MPI_DOUBLE, dest, tag, comm);
}

static inline
int MPI_Recv_Vec(VEC v, int source, int tag, MPI_Comm comm, MPI_Status* status)
{
    return MPI_Recv(v.storage, v.dim, MPI_DOUBLE, source, tag, comm, status);
}

static inline
int MPI_Pack_Vec(VEC v, void* outbuf, int outsize, int* position, MPI_Comm comm)
{
    return MPI_Pack(v.storage, v.dim, MPI_DOUBLE, outbuf, outsize, position, comm);
}

static inline
int MPI_Unpack_Vec(const void* inbuf, int insize, int* position, VEC v, MPI_Comm comm)
{
    return MPI_Unpack(inbuf, insize, position, v.storage, v.dim, MPI_DOUBLE, MPI_COMM_WORLD);
}


#endif //VECTOR_MPI_H