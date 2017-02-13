#ifndef MATHMETHODS_H
#define MATHMETHODS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vector.h"

//BLAS
void v_add(VEC u,VEC v,VEC w,unsigned int start);
void m_add(VEC2 A,VEC2 B,VEC2 C);
void mv_mlt(VEC2 A,VEC u,VEC v);
void mTv_mlt(VEC2 A,VEC u,VEC v);
void mm_mlt(VEC2 A,VEC2 B,VEC2 C);
void mTm_mlt(VEC2 A,VEC2 B,VEC2 C);
void mmT_mlt(VEC2 A,VEC2 B,VEC2 C);
void v_sub(VEC u,VEC v,VEC w,unsigned int start);
void sm_mlt(double alpha,VEC2 A,VEC2 B,unsigned int startrow,unsigned int startcol);
void dipaa(double alpha,VEC2 A,VEC2 B,unsigned int startrow,unsigned int startcol);
void diag_mm_mlt(VEC2 A,VEC2 B,VEC2 C);
void diag_mm_mltT(VEC2 A,VEC2 B,VEC2 C);
void diag_m_add(VEC2 A,VEC2 B,VEC2 C);
void diag_invert(VEC2 A,VEC2 B);
void diag_mv_mlt(VEC2 A,VEC u,VEC v);
void mTdiag_mlt(VEC2 A,VEC B,VEC2 C);
void mdiag_mlt(VEC2 A,VEC B,VEC2 C);
void diagm_mlt(VEC A,VEC2 B,VEC2 C);

void VDVT_mlt(VEC2 V,VEC D,VEC2 A);
void VDinvVT_mlt(VEC2 V,VEC D,VEC2 A,VEC temp);
void VsqrtDinvVT_mlt(VEC2 V,VEC D,VEC2 A,VEC temp);
void VsqrtDVT_mlt(VEC2 V,VEC D,VEC2 A,VEC temp);
double vTAv(VEC2 A,VEC v);

VEC lagrange_weights(unsigned short num_stages, VEC c);
//double lagrange_bary(double theta,short unsigned int s,VEC c,VEC* Q,VEC w);
void lagrange_bary(double theta,VEC c,VEC* Z,VEC w,VEC sum);
double lagrange(double theta,short unsigned int s,VEC c,VEC* Q);
double norm_inf(VEC v,VEC w,unsigned int start);
double norm_inf_u(VEC v,VEC w,unsigned int start,unsigned int t);
double vector_norminf(VEC v,unsigned int start);
void Print_Vector(VEC v);
void Print_Matrix(VEC2 A);
void Print_VectorC(VEC v);
void Print_MatrixC(VEC2 A);


void daxpy(double alpha, VEC x, VEC y, unsigned int begin);
void daxpy_u(double alpha, VEC x, VEC y, unsigned int begin, unsigned int end);
void sv_mlt(double val, VEC v, unsigned int begin);
void sv_mlt_u(double val, VEC v, unsigned int begin, unsigned int end);


#endif //MATHMETHODS_H