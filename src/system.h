#ifndef SYSTEM_H
#define SYSTEM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vector.h"
#include "structs.h"

extern int my_rank;
extern int np;

//Initialization
void Init_List(RKSolutionList* list, VEC y0, double t0, int dim, unsigned int num_dense, unsigned short int num_stages, unsigned int list_length);

//Destructors
void Destroy_Link(Link* link_i, unsigned int list_length, int rkd_flag, Forcing* forcings, GlobalVars* GlobalVars);
void ForcingData_Free(ForcingData** forcing_buff);
void Destroy_RKMethod(RKMethod* method);
void Destroy_ErrorData(ErrorData* error);
void Destroy_List(RKSolutionList* list, unsigned int list_length);
void Destroy_UnivVars(GlobalVars* GlobalVars);

//Discontinuity list
unsigned int Insert_Discontinuity(double time, unsigned int start, unsigned int end, unsigned int* count, unsigned int size, double* array, unsigned int id);
void Insert_SendDiscontinuity(double time, unsigned int order, unsigned int* count, unsigned int size, double* array, unsigned int*order_array, unsigned int id);

//RKSolution list
RKSolutionNode* New_Step(RKSolutionList* list);
void Undo_Step(RKSolutionList* list);
void Remove_Head_Node(RKSolutionList* list);

//Workspace methods
void Create_Workspace(Workspace *workspace, unsigned int dim, unsigned short num_stages, unsigned short int max_parents);
void Destroy_Workspace(Workspace* workspace, unsigned short int s, unsigned short int max_parents);

#endif
