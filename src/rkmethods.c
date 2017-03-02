#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "io.h"
#include "minmax.h"
#include "rkmethods.h"
#include "blas.h"

extern VEC* dump;


static double sq(double val)
{
    return val * val;
}


//Copies contents of the vectors full_k with dim entries into the vectors k with num_dense entries.
static void store_k(VEC2 full_k, VEC2 k, unsigned int s, unsigned int* dense_indices, unsigned int num_dense)
{
    unsigned int i, j;

    for (i = 0; i < s; i++)
    {
        for (j = 0; j < num_dense; j++)
            v2_set(k, i, j, v2_at(full_k, i, dense_indices[j]));
    }
}


//Calculate the weights for the Lagrange interpolation polynomial
//Assumes that 0 and each entry of c are node points
//This does not store the weight for 0
VEC lagrange_weights(unsigned short num_stages, VEC c)
{
    unsigned int i, j;
    VEC w = v_init(num_stages);

    for (i = 0; i < num_stages; i++)
    {
        w.storage[i] = 1.0 / v_at(c, i);	//For the 0 node
        for (j = 0; j < i; j++)
            w.storage[i] *= 1.0 / (v_at(c, i) - v_at(c, j));
        for (j = i + 1; j < num_stages; j++)
            w.storage[i] *= 1.0 / (v_at(c, i) - v_at(c, j));
    }

    return w;
}


//Builds a dense output RK method of order 4 at each step (order 3 for dense output)
void Init_TheRKDense4_3(RKMethod* method)
{
    method->s = 4;
    method->unique_c = 4;
    method->exp_imp = 0;
    method->A = v2_init(method->s, method->s);
    method->b = v_init(method->s);
    method->b_theta = v_init(method->s);
    method->b_theta_deriv = v_init(0);
    method->c = v_init(method->s);
    method->dense_b = &TheRKDense4_3_b;
    method->dense_bderiv = NULL;
    method->e = v_init(method->s);
    method->d = v_init(method->s);
    method->e_order = 4;
    method->e_order_ratio = 4.0 / 2.0;
    //method->e_order_ratio = 1.0/2.0;
    method->d_order = 3;
    method->d_order_ratio = 3.0 / 2.0;
    //method->d_order_ratio = 1.0/2.0;
//	method->d_max_error = 1.0/3.0;
    method->localorder = 4;

    //Build the parameters for the method
    v2_set(method->A, 1, 0, .5);
    v2_set(method->A, 2, 1, .5);
    v2_set(method->A, 3, 2, 1);

    v_set(method->b, 0, 1.0 / 6.0);
    v_set(method->b, 1, 2.0 / 6.0);
    v_set(method->b, 2, 2.0 / 6.0);
    v_set(method->b, 3, 1.0 / 6.0);

    method->dense_b(1.0, method->b_theta);

    v_set(method->c, 0, 0);
    v_set(method->c, 1, .5);
    v_set(method->c, 2, .5);
    v_set(method->c, 3, 1);

    v_set(method->e, 0, 2.0 / 3.0);
    v_set(method->e, 1, 0.0);
    v_set(method->e, 2, -4.0 / 3.0);
    v_set(method->e, 3, 2.0 / 3.0);

    v_set(method->d, 0, 2.0 / 9.0);
    v_set(method->d, 1, 0.0);
    v_set(method->d, 2, -4.0 / 9.0);
    v_set(method->d, 3, 2.0 / 9.0);

    method->w = v_init(0);
}

//The b(theta) coefficients for TheRKDense4_3()
void TheRKDense4_3_b(double theta, VEC b)
{
    double t2 = theta * theta;
    double t3 = t2 * theta;

    v_set(b, 0, theta - 3.0 / 2.0 * t2 + 2.0 / 3.0 * t3);
    v_set(b, 1, t2 - 2.0 / 3.0 * t3);
    v_set(b, 2, t2 - 2.0 / 3.0 * t3);
    v_set(b, 3, -1.0 / 2.0 * t2 + 2.0 / 3.0 * t3);
}

//Builds a dense output version of Dormand-Prince 5(4) with dense order 4
void Init_DOPRI5_dense(RKMethod* method)
{
    method->s = 7;
    method->unique_c = 6;
    method->exp_imp = 0;
    method->A = v2_init(method->s, method->s);
    method->b = v_init(method->s);
    method->b_theta = v_init(method->s);
    method->b_theta_deriv = v_init(method->s);
    method->c = v_init(method->s);
    method->dense_b = &DOPRI5_b;
    method->dense_bderiv = &DOPRI5_bderiv;
    method->e = v_init(method->s);
    method->d = v_init(method->s);
    method->e_order = 5;
    //method->e_order = 6;
    method->e_order_ratio = 6.0 / 5.0;
    //method->e_order_ratio = 1.0/5.0;
    method->d_order = 4;
    //method->d_order = 5;
    method->d_order_ratio = 5.0 / 4.0;
    //method->d_order_ratio = 1.0/4.0;
//	method->d_max_error = .6510416666666667;
    method->localorder = 5;

    //Build the parameters for the method
    v2_set(method->A, 1, 0, 1.0 / 5.0);
    v2_set(method->A, 2, 0, 3.0 / 40.0);
    v2_set(method->A, 2, 1, 9.0 / 40.0);
    v2_set(method->A, 3, 0, 44.0 / 45.0);
    v2_set(method->A, 3, 1, -56.0 / 15.0);
    v2_set(method->A, 3, 2, 32.0 / 9.0);
    v2_set(method->A, 4, 0, 19372.0 / 6561.0);
    v2_set(method->A, 4, 1, -25360.0 / 2187.0);
    v2_set(method->A, 4, 2, 64448.0 / 6561.0);
    v2_set(method->A, 4, 3, -212.0 / 729.0);
    v2_set(method->A, 5, 0, 9017.0 / 3168.0);
    v2_set(method->A, 5, 1, -355.0 / 33.0);
    v2_set(method->A, 5, 2, 46732.0 / 5247.0);
    v2_set(method->A, 5, 3, 49.0 / 176.0);
    v2_set(method->A, 5, 4, -5103.0 / 18656.0);
    v2_set(method->A, 6, 0, 35.0 / 384.0);
    v2_set(method->A, 6, 1, 0.0);
    v2_set(method->A, 6, 2, 500.0 / 1113.0);
    v2_set(method->A, 6, 3, 125.0 / 192.0);
    v2_set(method->A, 6, 4, -2187.0 / 6784.0);
    v2_set(method->A, 6, 5, 11.0 / 84.0);

    method->dense_b(1.0, method->b);
    method->dense_b(1.0, method->b_theta);
    method->dense_bderiv(1.0, method->b_theta_deriv);

    v_set(method->c, 0, 0.0);
    v_set(method->c, 1, .2);
    v_set(method->c, 2, .3);
    v_set(method->c, 3, .8);
    v_set(method->c, 4, 8.0 / 9.0);
    v_set(method->c, 5, 1.0);
    v_set(method->c, 6, 1.0);

    v_set(method->e, 0, 71.0 / 57600.0);
    v_set(method->e, 1, 0.0);
    v_set(method->e, 2, -71.0 / 16695.0);
    v_set(method->e, 3, 71.0 / 1920.0);
    v_set(method->e, 4, -17253.0 / 339200.0);
    v_set(method->e, 5, 22.0 / 525.0);
    v_set(method->e, 6, -1.0 / 40.0);

    v_set(method->d, 0, .610351562499951);
    v_set(method->d, 1, 0.0);
    v_set(method->d, 2, -2.105795148247852);
    v_set(method->d, 3, 18.310546874999346);
    v_set(method->d, 4, -25.185639003536881);
    v_set(method->d, 5, 20.749496981890658);
    v_set(method->d, 6, -12.378961267605213);

    method->w = lagrange_weights(method->unique_c, method->c);
}

//The b(theta) coefficients for DOPRI5_dense()
void DOPRI5_b(double theta, VEC b)
{
    /*
        double theta_sq = sq(theta);
        double ntheta2m3 = 3.0 - 2.0*theta;
        double thetam1_sq = sq(theta-1.0);

        v_set(b, 0, theta_sq*( ntheta2m3*(35.0/384.0) - thetam1_sq*(5.0/11282082432.0)*(2558722523.0-31403016.0*theta) ) + theta*thetam1_sq);
        v_set(b, 1, 0.0);
        v_set(b, 2, theta_sq*( ntheta2m3*(500.0/1113.0) + thetam1_sq*(100.0/32700410799.0)*(882725551.0-15701508.0*theta) ));
        v_set(b, 3, theta_sq*( ntheta2m3*(125.0/192.0) - thetam1_sq*(25.0/1880347072.0)*(443332067.0-31403016.0*theta) ));
        v_set(b, 4, theta_sq*( ntheta2m3*(-2187.0/6784.0) + thetam1_sq*(32805.0/199316789632.0)*(23143187.0-3489224.0*theta) ));
        v_set(b, 5, theta_sq*( ntheta2m3*(11.0/84.0) - thetam1_sq*(55.0/822651844.0)*(29972135.0-7076736.0*theta) ));
        v_set(b, 6, theta_sq*( (theta-1.0) + thetam1_sq*(10.0/29380423.0)*(7414447.0-829305.0*theta) ));
    */

    v_set(b, 0, sq(theta)*((3.0 - 2.0*theta)*(35.0 / 384.0) - sq(theta - 1.0)*(5.0 / 11282082432.0)*(2558722523.0 - 31403016.0*theta)) + theta*sq(theta - 1.0));
    v_set(b, 1, 0.0);
    v_set(b, 2, sq(theta)*((3.0 - 2.0*theta)*(500.0 / 1113.0) + sq(theta - 1.0)*(100.0 / 32700410799.0)*(882725551.0 - 15701508.0*theta)));
    v_set(b, 3, sq(theta)*((3.0 - 2.0*theta)*(125.0 / 192.0) - sq(theta - 1.0)*(25.0 / 1880347072.0)*(443332067.0 - 31403016.0*theta)));
    v_set(b, 4, sq(theta)*((3.0 - 2.0*theta)*(-2187.0 / 6784.0) + sq(theta - 1.0)*(32805.0 / 199316789632.0)*(23143187.0 - 3489224.0*theta)));
    v_set(b, 5, sq(theta)*((3.0 - 2.0*theta)*(11.0 / 84.0) - sq(theta - 1.0)*(55.0 / 822651844.0)*(29972135.0 - 7076736.0*theta)));
    v_set(b, 6, sq(theta)*((theta - 1.0) + sq(theta - 1.0)*(10.0 / 29380423.0)*(7414447.0 - 829305.0*theta)));

    /*
        v_set(b, 0, sq(theta)*(3.0-2.0*theta)*35.0/384.0 + theta*sq(theta-1.0) - sq(theta)*sq(theta-1.0)*5.0*(2558722523.0-31403016.0*theta)/11282082432.0);
        v_set(b, 1, 0.0);
        v_set(b, 2, sq(theta)*(3.0-2.0*theta)*500.0/1113.0 + sq(theta)*sq(theta-1.0)*100.0*(882725551.0-15701508.0*theta)/32700410799.0);
        v_set(b, 3, sq(theta)*(3.0-2.0*theta)*125.0/192.0 - sq(theta)*sq(theta-1.0)*25.0*(443332067.0-31403016.0*theta)/1880347072.0);
        v_set(b, 4, sq(theta)*(3.0-2.0*theta)*-2187.0/6784.0 + sq(theta)*sq(theta-1.0)*32805.0*(23143187.0-3489224.0*theta)/199316789632.0);
        v_set(b, 5, sq(theta)*(3.0-2.0*theta)*11.0/84.0 - sq(theta)*sq(theta-1.0)*55.0*(29972135.0-7076736.0*theta)/822651844.0);
        v_set(b, 6, sq(theta)*(theta-1.0) + sq(theta)*sq(theta-1.0)*10.0*(7414447.0-829305.0*theta)/29380423.0);
    */
}

//The b'(theta) coefficients for DOPRI5_dense()
void DOPRI5_bderiv(double theta, VEC b)
{
    const double prod1 = 6.0*theta*(1.0 - theta);
    const double prod2 = 2.0*theta*(theta - 1.0)*(2.0*theta - 1.0);
    const double prod3 = sq(theta)*sq(theta - 1.0);

    v_set(b, 0, prod1*(35.0 / 384.0) + (theta - 1.0)*(3.0*theta - 1.0) - 2.0*theta*(theta - 1.0)*(2.0*theta - 1.0)*(5.0 / 11282082432.0)*(2558722523.0 - 31403016.0*theta) + prod3*(157015080.0 / 11282082432.0));
    v_set(b, 1, 0.0);
    v_set(b, 2, prod1*(500.0 / 1113.0) + prod2*(100.0 / 32700410799.0)*(882725551.0 - 15701508.0*theta) - prod3*(1570150800.0 / 32700410799.0));
    v_set(b, 3, prod1*(125.0 / 192.0) - prod2*(25.0 / 1880347072.0)*(443332067.0 - 31403016.0*theta) + prod3*(785075400.0 / 1880347072.0));
    v_set(b, 4, prod1*(-2187.0 / 6784.0) + prod2*(32805.0 / 199316789632.0)*(23143187.0 - 3489224.0*theta) - prod3*(1144640195640.0 / 199316789632.0));
    v_set(b, 5, prod1*(11.0 / 84.0) - prod2*(55.0 / 822651844.0)*(29972135.0 - 7076736.0*theta) + prod3*(389220480.0 / 822651844.0));
    v_set(b, 6, theta*(3.0*theta - 2.0) + prod2*(10.0 / 29380423.0)*(7414447.0 - 829305.0*theta) - prod3*(8293050.0 / 29380423.0));
}

//Builds a dense output RK method of order 3 at each step (order 2 for dense output)
void Init_RKDense3_2(RKMethod* method)
{
    method->s = 3;
    method->unique_c = 3;
    method->exp_imp = 0;
    method->A = v2_init(method->s, method->s);
    method->b = v_init(method->s);
    method->b_theta = v_init(method->s);
    method->b_theta_deriv = v_init(method->s);
    method->c = v_init(method->s);
    method->dense_b = &RKDense3_2_b;
    method->dense_bderiv = &RKDense3_2_bderiv;
    method->e = v_init(method->s);
    method->d = v_init(method->s);
    method->e_order = 3;
    method->e_order_ratio = 3.0 / 2.0;
    //method->e_order_ratio = 1.0/2.0;
    method->d_order = 2;
    method->d_order_ratio = 2.0 / 2.0;
    //method->d_order_ratio = 1.0/2.0;
//	method->d_max_error = 2.0/3.0;
    method->localorder = 3;

    //Build the coefficients for the method
    v2_set(method->A, 1, 0, .5);
    v2_set(method->A, 2, 0, -1);
    v2_set(method->A, 2, 1, 2);

    method->dense_b(1.0, method->b);
    method->dense_b(1.0, method->b_theta);

    v_set(method->c, 0, 0.0);
    v_set(method->c, 1, .5);
    v_set(method->c, 2, 1.0);

    v_set(method->e, 0, 2.0 / 3.0);
    v_set(method->e, 1, -4.0 / 3.0);
    v_set(method->e, 2, 2.0 / 3.0);

    v_set(method->d, 0, 1.0 / 3.0);
    v_set(method->d, 1, -2.0 / 3.0);
    v_set(method->d, 2, 1.0 / 3.0);

    method->w = v_init(0);
}

//The b(theta) coefficients for RKDense3_2()
void RKDense3_2_b(double theta, VEC b)
{
    v_set(b, 0, theta*(1 + theta*(-3.0 / 2.0 + 2.0 / 3.0*theta)));
    v_set(b, 1, 2 * theta*theta*(1 - 2.0 / 3.0*theta));
    v_set(b, 2, theta*theta*(2.0 / 3.0*theta - .5));
}

//The b'(theta) coefficients for RKDense3_2()
void RKDense3_2_bderiv(double theta, VEC b)
{
    v_set(b, 0, 1.0 - 3.0*theta + 2.0*theta*theta);
    v_set(b, 1, 4.0*theta - 4.0*theta*theta);
    v_set(b, 2, 2.0*theta*theta - theta);
}

//Builds a dense output version of the 3 stage RadauIIA method with dense output
void Init_RadauIIA3_dense(RKMethod* method)
{
    method->s = 3;
    method->unique_c = 3;
    method->exp_imp = 1;
    method->A = v2_init(method->s, method->s);
    method->b = v_init(method->s);
    method->b_theta = v_init(method->s);
    method->b_theta_deriv = v_init(0);
    method->dense_bderiv = NULL;
    method->c = v_init(method->s);
    method->dense_b = &RadauIIA3_b;
    method->e = v_init(method->s + 1);
    method->d = v_init(method->s + 1);
    method->e_order = 4;
    method->e_order_ratio = 4.0 / 3.0;
    method->d_order = 3;
    method->d_order_ratio = 3.0 / 2.0;
    method->localorder = 5;

    //Build the parameters for the method
    v2_set(method->A, 0, 0, 0.1968154772236604);
    v2_set(method->A, 0, 1, -0.06553542585019839);
    v2_set(method->A, 0, 2, 0.02377097434822015);
    v2_set(method->A, 1, 0, 0.3944243147390873);
    v2_set(method->A, 1, 1, 0.2920734116652285);
    v2_set(method->A, 1, 2, -0.04154875212599793);
    v2_set(method->A, 2, 0, 0.3764030627004673);
    v2_set(method->A, 2, 1, 0.5124858261884216);
    v2_set(method->A, 2, 2, 0.1111111111111111);

    method->dense_b(1.0, method->b);
    method->dense_b(1.0, method->b_theta);

    v_set(method->c, 0, 0.1550510257216822);
    v_set(method->c, 1, 0.6449489742783178);
    v_set(method->c, 2, 1.0);

    v_set(method->e, 0, -2.762305454748599);
    v_set(method->e, 1, 0.379935598252729);
    v_set(method->e, 2, -0.091629609865226);
    v_set(method->e, 3, 0.2748888295956774);
    /*
        //Two trees 0
        v_set(method->d, 0, 0.8097732865210982);
        v_set(method->d, 1, -2.701318593546288);
        v_set(method->d, 2, 1.616656477429513);
        v_set(method->d, 3, 0.2748888295956774);
    */

    //Three trees 0 (Original)
    v_set(method->d, 0, -0.428298294115368);
    v_set(method->d, 1, 0.245039074384917);
    v_set(method->d, 2, -0.091629609865226);
    v_set(method->d, 3, 0.2748888295956774);

    /*
        //No d_0, two tree 0
        v_set(method->d, 0, 1.238071580636466);
        v_set(method->d, 1, -2.946357667931205);
        v_set(method->d, 2, 1.708286087294739);

    */
    method->w = lagrange_weights(method->s, method->c);
}

//The b(theta) coefficients for RadauIIA2_dense()
void RadauIIA3_b(double theta, VEC b)
{
    double t2 = theta * theta;
    double t3 = t2 * theta;

    v_set(b, 0, 1.558078204724922 * theta - 1.986947221348443 * t2 + 0.805272079323988 * t3);
    v_set(b, 1, -0.891411538058256 * theta + 3.320280554681776 * t2 - 1.916383190435099 * t3);
    v_set(b, 2, 0.3333333333333333 * theta - 1.333333333333333 * t2 + 1.111111111111111 * t3);
}

