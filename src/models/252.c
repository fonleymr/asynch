#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>

#include <structs.h>
#include <rksteppers.h>
#include <models/model.h>
#include <models/check_consistency.h>


static unsigned int dense_indices[] = { 0 };

static double global_params[11];

extern const AsynchModel model_252;

// AKA ConvertParams
static void convert(double *params, unsigned int type, void* external)
{
    params[1] *= 1000;		//L_h: km -> m
    params[2] *= 1e6;		//A_h: km^2 -> m^2
}

// AKA InitRoutines
static void routines(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* user)
{
    link->dim = model_252.dim;
}

//AKA Precalculations
static void precalculations(Link* link_i, double *global_params, double *params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void* user)
{
    //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
    //The numbering is:	0   1   2    3     4   5   6   7 
    //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
    //The numbering is:        0      1        2     3   4     5        6   7  8 9  10
    double* vals = params;
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];

    double v_0 = global_params[0];
    double lambda_1 = global_params[1];
    double lambda_2 = global_params[2];
    double v_h = global_params[3];
    double k_i_factor = global_params[5];

    vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
    vals[4] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
    vals[5] = vals[4] * k_i_factor;	//[1/min] k_i
    vals[6] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
    vals[7] = A_h / 60.0;	//  c_2
}

//AKA ReadInitData
static void initialize_eqs(double *global_params, double *params, double *y_0, void *user)
{

}

// AKA Partition_System_By_Leaves
static void partition(Link* sys, unsigned int N, Link** leaves, unsigned int numleaves, unsigned int** my_sys, unsigned int* my_N, TransData* my_data, short int *getting)
{

}


//#pragma omp declare simd notinbranch uniform(num_dof, num_parents, global_params, params, forcing_values, user) 
void TopLayerHillslopeSIMD(
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned int num_parents,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    void *user,
    double *ans)
{
    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double rain = forcing_values[0];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
                                                                    //double e_pot = global_params[11];	//[m/min]
                                                                    //double e_pot = global_params[11] * (1e-3*60.0);	//[m/min]
                                                                    //double e_pot = 0.0;

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                            //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]

                            //Evaporation
    double e_p, e_t, e_s;
    double corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && corr > 1e-12)
    {
        e_p = s_p * 1e3 * e_pot / corr;
        e_t = s_t / S_L * e_pot / corr;
        e_s = s_s / (h_b - S_L) * e_pot / corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;

    //Discharge
    double discharge = -q + (q_pl + q_sl) * c_2;
#pragma novector
    for (unsigned int i = 0; i < num_parents; i++)
        discharge += y_p[i * num_dof];

    //Hillslope
    //v_set(ans, 0, discharge);
    //v_set(ans, 1, rain * c_1 - q_pl - q_pt - e_p);
    //v_set(ans, 2, q_pt - q_ts - e_t);
    //v_set(ans, 3, q_ts - q_sl - e_s);

    //ans[0] = discharge;
    ans[0] = invtau * pow(q, lambda_1) * discharge;
    ans[1] = rain * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;
}

//#pragma omp declare simd notinbranch
void CheckConsistency_Nonzero_4StatesSIMD(
    double *y, unsigned int num_dof,
    double *params, unsigned int num_params,
    double *global_params, unsigned int num_global_params,
    void *user)
{
    if (y[0] < 1e-14)	y[0] = 1e-14;
    if (y[1] < 0.0)	y[1] = 0.0;
    if (y[2] < 0.0)	y[2] = 0.0;
    if (y[3] < 0.0)	y[3] = 0.0;
}


static const AsynchModel model_252 = {

    252,            // uid

    4,              // dim
    0,              // diff_start
    4,              // no_ini_start
    
    1,              // num_dense
    dense_indices,  // dense_indices

    11,             // num_global_params
    global_params,  // global_params


    0,              // use_dam
    8,              // num_params
    0,              // num_dam_params_size
    3,              // num_disk_params

    0,              // area_idx
    2,              // areah_idx
    false,          // convertarea_flag

    4,              // min_error_tolerances

    2,              // num_forcings

    &TopLayerHillslopeSIMD,     // differential
    NULL,                   // jacobian
    NULL,                   // algebraic
    NULL,                   // check_state
    &CheckConsistency_Nonzero_4StatesSIMD,      // check_consistency
    &ExplicitRKSolver,                          // solver
    
    &convert,           // convert
    &routines,          // routines
    &precalculations,   // precalculations
    NULL,               // initialize_eqs
    NULL                // partition
};


const AsynchModel *the_model = &model_252;