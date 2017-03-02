#if !defined(ASYNCH_MODEL_CHECK_STATE_H)
#define ASYNCH_MODEL_CHECK_STATE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


//State Checks
int dam_check(
    double *y, unsigned int num_dof,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    void *user);

int dam_check2(
    double *y, unsigned int num_dof,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    void *user);

int dam_check3(
    double *y, unsigned int num_dof,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    void *user);

int dam_check_qvs(
    double *y, unsigned int num_dof,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    void *user);

#endif //!defined(ASYNCH_MODEL_CHECK_STATE_H)
