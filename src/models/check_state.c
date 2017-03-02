#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>

#include <models/check_state.h>


//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
int dam_check(
    double *y, unsigned int num_dof,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    void *user)
{
    //if (dam == 0)
    //    return 0;

    double H_spill = params[7];
    double H_max = params[8];
    double S_max = params[9];
    double alpha = params[10];
    double diam = params[11];
    double S = y[1];
    double h = H_max * pow(S / S_max, alpha);

    if (h < diam)		return 4;
    if (h <= H_spill)	return 1;
    if (h <= H_max)		return 2;
    return 3;
}