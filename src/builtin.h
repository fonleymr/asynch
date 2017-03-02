#ifndef DEFINETYPE_H
#define DEFINETYPE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"


void SetParamSizes(GlobalVars* globals, void* external);
void ConvertParams(double *params, unsigned int type, void* external);
void InitRoutines(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int is_dam, void* external);
void Precalculations(Link* link_i, double *global_params, double *params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void* external);
int ReadInitData(VEC global_params, VEC params, QVSData* qvs, unsigned short int is_dam, VEC y_0, unsigned int type, unsigned int diff_start, unsigned int no_init_start, void* user, void* external);
//void AssimError(unsigned int N,UnivVars* GlobalVars,ErrorData* GlobalErrors);

#endif
