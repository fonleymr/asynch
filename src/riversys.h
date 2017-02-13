#ifndef RIVERSYS_H
#define RIVERSYS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"


extern int np;
extern int my_rank;

//Routines to construct the system of equations to solve
//Link** Create_River_System_parallel(char rk_filename[],unsigned int* N,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData** my_data,int** assignments,short int** getting,RKMethod*** rk_methods,unsigned int* num_methods,UnivVars* GlobalVars,ErrorData* errors,unsigned int** save_list,unsigned int* my_save_size,unsigned int* save_size,unsigned int* peaksave_size,unsigned int*** id_to_loc,Forcing** forcings,ConnData** db_connections,Workspace** workspace,model* custom_model,void* external);


void Create_River_Network(GlobalVars* GlobalVars, Link** system, unsigned int* N, unsigned int*** id_to_loc, ConnData* db_connections);

int Load_Local_Parameters(Link* system, unsigned int N, unsigned int* my_sys, unsigned int my_N, int* assignments, short int* getting, unsigned int** id_to_loc, GlobalVars* GlobalVars, ConnData* db_connections, Model* custom_model, void* external);
int Partition_Network(Link* system, unsigned int N, GlobalVars* GlobalVars, unsigned int** my_sys, unsigned int* my_N, int** assignments, TransData** my_data, short int** getting, Model* custom_model);
int Build_RKData(Link* system, char rk_filename[], unsigned int N, unsigned int* my_sys, unsigned int my_N, int* assignments, short int* getting, GlobalVars* globals, ErrorData* error_data, RKMethod** methods, unsigned int* num_methods);
int Initialize_Model(Link* system, unsigned int N, unsigned int* my_sys, unsigned int my_N, int* assignments, short int* getting, GlobalVars* GlobalVars, Model* custom_model, void* external);
int Load_Initial_Conditions(Link* system, unsigned int N, int* assignments, short int* getting, unsigned int** id_to_loc, GlobalVars* GlobalVars, ConnData* db_connections, Model* custom_model, void* external);
int Load_Forcings(Link* system, unsigned int N, unsigned int* my_sys, unsigned int my_N, int* assignments, short int* getting, unsigned int* res_list, unsigned int res_size, unsigned int** id_to_loc, GlobalVars* GlobalVars, Forcing* forcings, ConnData* db_connections);
int Load_Dams(Link* system, unsigned int N, unsigned int* my_sys, unsigned int my_N, int* assignments, short int* getting, unsigned int** id_to_loc, GlobalVars* GlobalVars, ErrorData* GlobalError, ConnData* db_connections, unsigned int** res_list, unsigned int* res_size, unsigned int* my_res_size);
int CalculateInitialStepSizes(Link* system, unsigned int* my_sys, unsigned int my_N, GlobalVars* GlobalVars, Workspace* workspace, short int watch_forcings);
int BuildSaveLists(Link* system, unsigned int N, unsigned int* my_sys, unsigned int my_N, unsigned int* assignments, unsigned int** id_to_loc, GlobalVars* GlobalVars, unsigned int** save_list, unsigned int* save_size, unsigned int* my_save_size, unsigned int** peaksave_list, unsigned int* peaksave_size, unsigned int* my_peaksave_size, ConnData* db_connections);
int FinalizeSystem(Link* system, unsigned int N, unsigned int* my_sys, unsigned int my_N, unsigned int* assignments, short int* getting, unsigned int** id_to_loc, TransData* my_data, GlobalVars* GlobalVars, ConnData* db_connections, Workspace* workspace);

GlobalVars* Read_Global_Data(char globalfilename[], ErrorData** errors, Forcing* forcings, ConnData* db_connections, char* rkdfilename, Model* custom_model, void* external);

void ReadDBC(char* filename, ConnData* const conninfo);

int Create_SAV_Data(char filename[], Link* sys, unsigned int N, unsigned int** save_list, unsigned int* size, ConnData *conninfo, unsigned short int flag);

void ReadLineFromTextFile(FILE* globalfile, char linebuffer[], unsigned int size);

int ReadLineError(int valsread, int valswant, char message[]);

int RemoveSuffix(char* filename, const char* suffix);

int AttachParameters(char* filename, unsigned int max_size, VEC v, unsigned int string_size);

int CheckFilenameExtension(char* filename, char* extension);
int CheckWinFormat(FILE* file);
int FindPath(char* filename, char* path);
int FindFilename(char* fullpath, char* filename);

#endif

