#ifndef RKMETHODS_H
#define RKMETHODS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdbool.h>

#include "vector.h"
#include "structs.h"

extern int my_rank;
extern int np;

double InitialStepSize(double t,Link* link_i,GlobalVars* GlobalVars,Workspace* workspace);
//void BackupParent(Link* currentp,double time,VEC* backup,UnivVars* GlobalVars);

//solver methods
//int ExplicitRKSolverLinear(Link* link_i,UnivVars* GlobalVars,int* assignments,FILE* outputfile,Workspace* workspace);
int ExplicitRKSolver(Link* link_i,GlobalVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace);
int ExplicitRKIndex1SolverDam(Link* link_i,GlobalVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace);
int ExplicitRKIndex1Solver(Link* link_i,GlobalVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace);
int ExplicitRKSolverDiscont(Link* link_i,GlobalVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace);

//Forced solution methods
int ForcedSolutionSolver(Link* link_i,GlobalVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace);

//RK methods data
void Init_RKDense3_2(RKMethod* method);
void RKDense3_2_b(double theta,VEC b);
void RKDense3_2_bderiv(double theta,VEC b);

void Init_TheRKDense4_3(RKMethod* method);
void TheRKDense4_3_b(double theta,VEC b);

void Init_DOPRI5_dense(RKMethod* method);
void DOPRI5_b(double theta,VEC b);
void DOPRI5_bderiv(double theta,VEC b);

void Init_RadauIIA3_dense(RKMethod* method);
void RadauIIA3_b(double theta,VEC b);

#endif

