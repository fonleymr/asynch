#if !defined(ASYNCH_MODEL_H)
#define ASYNCH_MODEL_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdbool.h>

#include <structs_fwd.h>


/// Get a model given its uid
/// 
/// \param_uid Model uid
/// \return A pointer to the model
AsynchModel const * GetModel(unsigned short model_uid);


// Right-hand side function for ODE
typedef void (DifferentialFunc) (
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned int num_parents,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    void *user,
    double *ans);                                                     

// Right-hand side function for algebraic variables
typedef void (AlgebraicFunc)(
    const double * const y_i, unsigned int num_dof,
    const double * const global_params,
    const double * const params,
    void* user,
    double * ans);

/// Jacobian of right-hand side function
typedef void (JacobianFunc)(
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned int num_parents,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    double *ans);                                    

/// RK solver
typedef int (RKSolverFunc)(
    Link* link,
    GlobalVars* globals,
    int* assignments,
    bool print_flag,
    FILE* outputfile,
    ConnData* conninfo,
    Forcing* forcings,
    Workspace* workspace);

// Function to check what "state" the state variables are in (for discontinuities)
typedef int (CheckStateFunc)(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);                                                    

// Function to check state consistency
typedef void (CheckConsistencyFunc)(
    double *y, unsigned int num_dof,
    double *params, unsigned int num_params,
    double *global_params, unsigned int num_global_params,
    void *user);                                                    

// Models function signatures
//typedef void (SetParamSizesFunc)(GlobalVars* globals, void* user);
typedef void (ConvertFunc)(double *params, unsigned int type, void* user);
typedef void (RoutinesFunc)(Link*, unsigned int, unsigned int, unsigned short int, void *user);
typedef void (PrecalculationsFunc)(Link* link_i, double *global_params, double *params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void *user);
typedef int (InitializeEqsFunc)(double *global_params, double *params, double *y_0, void *user);
typedef int* (PartitionFunc)(Link *sys, unsigned int N, Link **leaves, unsigned int num_leaves, Link ***my_sys, unsigned int *my_N, TransData *my_data, short int *getting);

typedef struct AsynchModel
{
    unsigned short uid;                 //!< Unique identifier of the model

    unsigned short dim;                 //!< Dimension of the problem at this link
    unsigned int diff_start;            //!< Starting index of differential variables in solution vectors
    unsigned int no_ini_start;          //!< Starting index of differential variables not read from disk

    unsigned int num_dense;             //!< Number of states where dense output is calculated (usually only discharge is used)
    unsigned int *dense_indices;        //!< List of indices in solution where dense output is needed

    unsigned int num_global_params;     //!< Number of global parameters
    double *global_params;              //!< List of global parameters    

    bool uses_dam;                      //!< true if this type can use dams, false else
    unsigned int num_params;            //!< The number of params at each link without a dam
    unsigned int num_dam_params_size;   //!< The number of params at each link with a dam
    unsigned int num_disk_params;       //!< Number of parameters to read from disk
    
    unsigned int area_idx;              //!< Index of upstream area (A_i) in params
    unsigned int areah_idx;             //!< Index of hillslope area (A_h) in params
    bool convertarea_flag;              //!< true if hillslope and upstream areas are converted from km^2 to m^2, false if not
    
    unsigned int min_error_tolerances;  //!< The minimum number of error tolerances needed at every link. Used for uniform error tolerances.
    
    unsigned int num_forcings;          //!< The number of forcings

    DifferentialFunc *differential;         //!< Right-hand side function for ODE
    JacobianFunc *jacobian;                 //!< jacobian of right-hand side function
    AlgebraicFunc *algebraic;               //!< Function for algebraic variables
    CheckStateFunc *check_state;            //!< Function to check what "state" the state variables are in (for discontinuities)
    CheckConsistencyFunc *check_consistency; //!< Function to check state variables
    RKSolverFunc *solver;                   //!< RK solver to use
    

    //SetParamSizesFunc *set_param_sizes;
    ConvertFunc *convert;
    RoutinesFunc *routines;
    PrecalculationsFunc *precalculations;
    InitializeEqsFunc *initialize_eqs;
    PartitionFunc *partition;
} AsynchModel;


#endif //!defined(ASYNCH_MODEL_H)
