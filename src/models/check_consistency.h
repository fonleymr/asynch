#if !defined(ASYNCH_MODEL_CHECK_CONSISTENCY_H)
#define ASYNCH_MODEL_CHECK_CONSISTENCY_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


// Shared consistency checking functions
void CheckConsistency_Nonzero_1States(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

void CheckConsistency_Nonzero_2States(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

void CheckConsistency_Nonzero_3States(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

void CheckConsistency_Nonzero_4States(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

void CheckConsistency_Model5(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

//TODO mode to specific model
void CheckConsistency_Model30(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

void CheckConsistency_Nonzero_AllStates_q(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

void CheckConsistency_Nonzero_AllStates_qs(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params);

#endif //!defined(ASYNCH_MODEL_CHECK_CONSISTENCY_H)

