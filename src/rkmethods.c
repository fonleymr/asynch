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


double InitialStepSize(double t, Link* link_i, GlobalVars* GlobalVars, Workspace* workspace)
{
    unsigned int start = link_i->diff_start;
    VEC y_0 = link_i->my->list.tail->y_approx;
    double t_0 = link_i->my->list.tail->t;
    VEC fy0 = workspace->temp;
    VEC fy1 = workspace->sum;
    VEC y_1 = workspace->temp2;
    VEC SC = workspace->temp3;
    unsigned int p = link_i->method->localorder;
    Link* currentp;
    RKSolutionNode* curr_node;
    double d0, d1, d2, h0, h1, largest, timediff, current_theta;
    unsigned int i, l, m, idx;
    unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    ErrorData* error = &link_i->my->error_data;

    //Build SC for this link
    for (i = 0; i < dim; i++)
        v_set(SC, i, fabs(v_at(y_0, i)) * v_at(error->reltol, i) + v_at(error->abstol, i));

    //Grab parents data
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node = currentp->my->list.head;

        VEC parent_approx = v3_slice2(workspace->temp_parent_approx, 0, i);

        if (fabs(currentp->my->list.tail->t - t) < 1e-10)	//!!!! Ugh... !!!!
            v_copy_n(currentp->my->list.tail->y_approx, parent_approx, currentp->dim);
        else
        {
            //Find the corresponding theta value and approximate solution
            while (t > curr_node->t)
                curr_node = curr_node->next;
            if (curr_node != currentp->my->list.head)
                curr_node = curr_node->prev;

            timediff = curr_node->next->t - curr_node->t;
            current_theta = (t - curr_node->t) / timediff;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            // !!!! Note: this varies with num_print. Consider doing a linear interpolation. !!!!
            for (m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                v_set(parent_approx, idx, v_at(curr_node->y_approx, idx));

                for (l = 0; l < currentp->method->s; l++)
                    v_set(parent_approx, idx, v_at(parent_approx, idx)
                        + timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node->next->k, l, m));
            }
            link_i->check_consistency(parent_approx, currentp->params, GlobalVars->global_params);
            //if(link_i->algebraic != NULL)	link_i->algebraic(temp_parent_approx[0][i],GlobalVars->global_params,currentp->params,currentp->qvs,curr_node->next->state,temp_parent_approx[0][i]);
            if (currentp->algebraic)
                currentp->algebraic(parent_approx, GlobalVars->global_params, currentp->params, currentp->qvs, curr_node->state, currentp->user, parent_approx);
        }
    }

    VEC2 parents_approx = v3_slice(workspace->temp_parent_approx, 0);
    parents_approx.dim[0] = link_i->num_parents;    //Adjust the number of parents for this specific link (temp_parent_approx as max num parents)

    //Step a
    //d0 = vector_norminf(y0,start);
    d0 = norm_inf_u(y_0, SC, start, link_i->dim);
    link_i->differential(t_0, y_0, parents_approx, GlobalVars->global_params, link_i->forcing_values, link_i->qvs, link_i->params, link_i->state, link_i->user, fy0);
    d1 = norm_inf_u(fy0, SC, start, link_i->dim);

    //Step b
    if (d0 < 1e-5 || d1 < 1e-5)	h0 = 1e-6;
    else				h0 = 0.01*(d0 / d1);

    //Step c
    //Note: This assumes the parents discharge is the same. It also assume no change in rain or state
    v_copy_n(y_0, y_1, link_i->dim);
    daxpy_u(h0, fy0, y_1, start, link_i->dim);
    link_i->check_consistency(y_1, link_i->params, GlobalVars->global_params);
    link_i->differential(t_0 + h0, y_1, parents_approx, GlobalVars->global_params, link_i->forcing_values, link_i->qvs, link_i->params, link_i->state, link_i->user, fy1);

    //Step d
    v_sub(fy1, fy0, fy1, start);
    d2 = norm_inf_u(fy1, SC, start, link_i->dim);

    //Step e
    largest = max(d1, d2);
    h1 = (largest < 1e-1) ? max(1e-6, h0*1e-3) : pow(1e-2 / largest, 1.0 / (p + 1.0));

    //Step f
    h1 = min(100.0*h0, h1);

    //Make sure returned value is not too small. This can create roundoff problems when time gets larger
    //return h1;
    //frexp(h1,&exponent);
    //if(exponent < -29)	return 1e-10;
    if (v_at(y_0, 0) < 1e-6)
        return 1e-2;
    else		
        return h1;
}


//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//Link* link_i: the link to apply a numerical method to.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolver(Link* link_i, GlobalVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int idx;
    //VEC** k;
    VEC new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    Link* currentp;
    double t_needed, timediff, current_theta;

    //Some variables to make things easier to read
    VEC y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2 A = link_i->method->A;
    VEC b = link_i->method->b;
    VEC c = link_i->method->c;
    unsigned int num_stages = link_i->method->s;
    VEC params = link_i->params;
    VEC e = link_i->method->e;
    VEC d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = GlobalVars->num_outputs;
    VEC temp = workspace->temp;
    VEC sum = workspace->sum;
    VEC* temp_k = workspace->temp_k_slices;

    //Get the approximate solutions from each parent
    for (unsigned int i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;
        for (unsigned int j = 0; j < num_stages; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j) * h, currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            //while(t_needed > curr_node[i]->t && ( fabs(curr_node[i]->t) < 1e-12 || fabs(t_needed - curr_node[i]->t)/curr_node[i]->t > 1e-12) )
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)
                curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / timediff;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            VEC parent_approx = v3_slice2(workspace->temp_parent_approx, j, i);

            for (unsigned int m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                double approx = v_at(curr_node[i]->y_approx, idx);
                for (unsigned int l = 0; l < currentp->method->s; l++)
                    approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);
                
                v_set(parent_approx, idx, approx);
            }
            link_i->check_consistency(parent_approx, params, GlobalVars->global_params);
        }
    }
    
    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    new_y = new_node->y_approx;
    /*
    int stopper = 0;
    if(link_i->ID == 0 && link_i->dim > 7)
    {
    printf("************\n");
    printf("t = %e\n",t);
    Print_Vector(y_0);
    stopper = 1;
    //Print_Vector(new_y);
    printf("************\n");
    }
    */
    //Compute the k's
    for (unsigned int i = 0; i < num_stages; i++)
    {
        v_copy_n(y_0, sum, link_i->dim);
        for (unsigned int j = 0; j < i; j++)
            //daxpy_u(h * v2_at(A, i, j), v2_slice(temp_k, j), sum, 0, link_i->dim);
            daxpy_u(h * v2_at(A, i, j), temp_k[j], sum, 0, link_i->dim);
        
        link_i->check_consistency(sum, params, GlobalVars->global_params);
        
        VEC2 y_p = v3_slice(workspace->temp_parent_approx, i);
        y_p.dim[0] = link_i->num_parents;
        
        double dt = v_at(c, i) * h;

        link_i->differential(t + dt, sum, y_p, GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp_k[i]);
    }

    //Build the solution
    v_copy_n(y_0, new_y, link_i->dim);
    for (unsigned int i = 0; i < num_stages; i++)
        //daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 0, link_i->dim);
        daxpy_u(h * v_at(b, i), temp_k[i], new_y, 0, link_i->dim);
    link_i->check_consistency(new_y, params, GlobalVars->global_params);
    new_node->state = link_i->state;

    /*
    if(stopper)
    {
    printf("************\n");
    printf("t = %e\n",t);
    //Print_Vector(y_0);
    Print_Vector(new_y);
    printf("************\n");
    getchar();
    }
    */
    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(temp_k[0], sum, link_i->dim);
    sv_mlt_u(h * v_at(e, 0), sum, 0, link_i->dim);
    for (unsigned int i = 1; i < num_stages; i++)
        //daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
        daxpy_u(h*v_at(e, i), temp_k[i], sum, 0, link_i->dim);        

    //Build SC_i
    for (unsigned int i = 0; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i))) * v_at(error->reltol, i) + v_at(error->abstol, i));

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,0);
    err_1 = norm_inf_u(sum, temp, 0, link_i->dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);

    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(temp_k[0], sum, link_i->dim);
    sv_mlt_u(h * v_at(d, 0), sum, 0, link_i->dim);
    
    for (unsigned int i = 1; i < num_stages; i++)
        //daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
        daxpy_u(h*v_at(d, i), temp_k[i], sum, 0, link_i->dim);        
    
    for (unsigned int i = 0; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i))) * v_at(error->reltol_dense, i) + v_at(error->abstol_dense, i));

    //err_d = norm_inf(sum,temp,meth->d_order_ratio,0);
    err_d = norm_inf_u(sum, temp, 0, link_i->dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h*min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);


    /*
    if(err_1 < 1.0 && err_d < 1.0)
    {
    printf("Accepted time = %e new h = %e  %e %e\n",t+h,link_i->h,err_1,err_d);
    Print_Vector(new_y);
    }
    else
    {
    printf("Rejected time = %e new h = %e  %e %e\n",t+h,link_i->h,err_1,err_d);
    Print_Vector(new_y);
    }
    */

    /*
        //Try new error control
        //This uses new idea to modify step size of the parents, not link_i
        Link* child = link_i->child;
        //Link* p;
        double err_new = 0.0;
        double step_new = step_1;
        double sum_of_errors;
        if(child != NULL)
        {
            //Works ok
            //!!!! Assumes 1d problem !!!!
            //for(i=0;i<dim;i++)	tempv_at(2, i) = 0.0;
            sum_of_errors = 0.0;
            for(i=0;i<child->num_parents;i++)	//Assumes one state is passed link to link
                sum_of_errors += child->parents[i]->error_data->abstol_densv_at(e, 0);

            Jx_simple_river(child->list->tail->y_approx,GlobalVars->global_params,child->params,temp);
            sv_mlt((t+h - child->last_t)*sum_of_errors,temp,GlobalVars->diff_start);
            //tempv_at(2, 0) *= (t+h - child->last_t) * Jx;

            err_new = norm_inf(temp,error->abstol_dense,1.0,0);
            double value_new = pow(1.0/(link_i->h * err_new),1.0/4.0);	//!!!! For Dormand & Prince !!!!
            unsigned int maxorder = 4;	// !!!! Need loop !!!!
            double largest = child->parents[0]->h;
            for(i=1;i<child->num_parents;i++)
                largest = max(largest,child->parents[i]->h);
            step_new = largest * value_new;
            //for(i=0;i<child->num_parents;i++)
            //	child->parents[i]->h = min(child->parents[i]->h,step_new);
            link_i->h = min(link_i->h,step_new);
        }
    */


    //	if(err_1 < 1.0 && err_d < 1.0 && err_new < 1.0)
    if (err_1 < 1.0 && err_d < 1.0)
        //	if(err_1 < 1.0)
    {
        /*
                //Try new error control
                //This uses new idea to modify step size of the parents, not link_i
                if(link_i->num_parents > 0)
                {
                    double err_new = norm_inf(sum,error->abstol_dense,1.0,0);
                    double value_new = pow(1.0/(link_i->h * err_new),1.0/4.0);	//!!!! For Dormand & Prince !!!!
                    unsigned int maxorder = 4;	// !!!! Need loop !!!!
                    double smallest = link_i->parents[0]->h;
                    for(i=1;i<link_i->num_parents;i++)
                        smallest = min(smallest,link_i->parents[i]->h);
                    double newh = smallest * value_new * .9;
                    for(i=0;i<link_i->num_parents;i++)
                        link_i->parents[i]->h = min(link_i->parents[i]->h,newh);
                }
        */

        //Check if a discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(workspace->temp_k, new_node->k, num_stages, dense_indices, num_dense);

        //Check if new data should be written to disk
        if (print_flag)
        {
            //while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
            {
                /*
                                //Don't write anything if using data assimilation and at a time when data is available
                                if(GlobalVars->assim_flag)
                                {
                                    if( fabs(GlobalVars->maxtime - link_i->next_save) < 1e-13 )	break;
                                    //double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                                    //if(rounded < 1e-13 && -rounded < 1e-13)		break;
                                }
                */

                if (link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);
                for (unsigned int m = 0; m < num_dense; m++)
                {
                    idx = dense_indices[m];
                    double approx = v_at(node->y_approx, idx);
                    for (unsigned int l = 0; l < link_i->method->s; l++)
                        approx += h * v_at(link_i->method->b_theta, l) * v2_at(node->next->k, l, m);

                    v_set(sum, idx, approx);
                }
                link_i->check_consistency(sum, params, GlobalVars->global_params);


                //printf("writing ID = %u t = %e, q = %e s_p = %e\n",link_i->ID,link_i->next_save,suv_at(m, 0),sum.ve[1]);
                                //Write to a file
                                //fsetpos(outputfile,&(link_i->pos));
                                //WriteStep(link_i->next_save,sum,GlobalVars,params,link_i->state,outputfile,link_i->output_user,&(link_i->pos));
                WriteStep(outputfile, link_i->ID, link_i->next_save, sum, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
                //fgetpos(outputfile,&(link_i->pos));

/*
                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                    fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));
*/
/*
                #ifdef PRINT2DATABASE
                    sprintf(conninfo->query,"%i,%u,%.12e,%.12e\n",link_i->ID,(unsigned int)(link_i->next_save * 60.0) + conninfo->time_offset,suv_at(m, 0)/link_i->Q_TM,suv_at(m, 0));
                    unsigned int length = strlen(conninfo->query);
                    (conninfo->submission_content) += length;
                    while(conninfo->submission_content > conninfo->submission_size)
                    {
                        //Allocate more space
                        (conninfo->submission_size) *= 2;
                        conninfo->submission = realloc(conninfo->submission,conninfo->submission_size);
                    }
                    strcat(conninfo->submission,conninfo->query);
                #else //Write to a file
                    fsetpos(outputfile,&(link_i->pos));
                    fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                    for(j=0;j<num_print;j++)
                        fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                    fgetpos(outputfile,&(link_i->pos));
                #endif
*/

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (link_i->peak_flag && (v_at(new_y, 0) > v_at(link_i->peak_value, 0)))
        {
            v_copy_n(new_y, link_i->peak_value, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (unsigned int j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && (fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (unsigned int i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                unsigned int l;
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_buff[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                unsigned int i;
                for (i = l + 1; i < link_i->forcing_buff[j]->nrows; i++)
                {
                    //if(link_i->forcing_buff[j]->rainfall[i][1] != forcing_buffer)
                    if (fabs(link_i->forcing_buff[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i - 1][0];
            }
        }

        //Select new step size, if forcings changed
        if (propagated)	link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);

        //Free up parents' old data
        for (unsigned int i = 0; i < link_i->num_parents; i++)
        {
            currentp = link_i->parents[i];
            while (currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(&currentp->my->list);
                currentp->current_iterations--;
                currentp->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}


//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//Link* link_i: the link to apply a numerical method to.
//VEC* sum: some space for temporary calculations. Should have same dimension as number of equations.
//VEC* temp: same as sum.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKIndex1Solver(Link* link_i, GlobalVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l, m, idx;
    //VEC** k;
    VEC new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    Link* currentp;
    double t_needed, timediff, theta;

    //Some variables to make things easier to read
    VEC y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2 A = link_i->method->A;
    VEC b = link_i->method->b;
    VEC c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC params = link_i->params;
    VEC e = link_i->method->e;
    VEC d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    const unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = GlobalVars->num_outputs;
    unsigned int* print_indices = GlobalVars->print_indices;
    VEC temp = workspace->temp;
    VEC sum = workspace->sum;
    VEC2 temp_k = workspace->temp_k;

    //Get the approximate solutions from each parent
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;
        for (j = 0; j < s; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j) * h, currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            theta = (t_needed - curr_node[i]->t) / timediff;
            currentp->method->dense_b(theta, currentp->method->b_theta);

            VEC parent_approx = v3_slice2(workspace->temp_parent_approx, j, i);

            for (m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                double approx = v_at(curr_node[i]->y_approx, idx);
                for (l = 0; l < currentp->method->s; l++)
                    approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                v_set(parent_approx, idx, approx);
            }

            //Build the algebraic variable
            currentp->check_consistency(parent_approx, currentp->params, GlobalVars->global_params);
            //link_i->algebraic(temp_parent_approx[j][i],GlobalVars->global_params,currentp->params,currentp->qvs,curr_node[i]->next->state,temp_parent_approx[j][i]);
            if (currentp->algebraic)
                currentp->algebraic(parent_approx, GlobalVars->global_params, currentp->params, currentp->qvs, curr_node[i]->state, currentp->user, parent_approx);
        }
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    new_y = new_node->y_approx;

    //Compute the k's
    for (i = 0; i < s; i++)
    {
        VEC2 parent_approx = v3_slice(workspace->temp_parent_approx, i);

        v_copy_n(y_0, sum, link_i->dim);
        for (j = 0; j < i; j++)
            daxpy_u(h * v2_at(A, i, j), v2_slice(temp_k, j), sum, 1, link_i->dim);
        link_i->check_consistency(sum, params, GlobalVars->global_params);
        link_i->differential(t + v_at(c, i) * h, sum, parent_approx, GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, v2_slice(temp_k, i));
    }

    //Build the solution
    v_copy_n(y_0, new_y, link_i->dim);
    for (i = 0; i < s; i++)
        daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 1, link_i->dim);
    //new_node->state = 0;
    new_node->state = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);
    link_i->check_consistency(new_y, params, GlobalVars->global_params);
    link_i->algebraic(new_y, GlobalVars->global_params, params, link_i->qvs, new_node->state, link_i->user, new_y);

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    sv_mlt_u(h * v_at(e, 0), sum, 1, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 1, link_i->dim);

    //Build SC_i
    for (i = 1; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * v_at(error->reltol, i) + v_at(error->abstol, i)));

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,1);
    err_1 = norm_inf_u(sum, temp, 1, link_i->dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);

    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    sv_mlt_u(h * v_at(d, 0), sum, 1, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 1, link_i->dim);

    for (i = 1; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * v_at(error->reltol_dense, i) + v_at(error->abstol_dense, i)));

    err_d = norm_inf_u(sum, temp, 1, link_i->dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h*min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);
    //link_i->h = step_1;

    if (err_1 < 1.0 && err_d < 1.0)
        //if(err_1 < 1.0)
    {
        //Check if a discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        //new_node->state = 0;
        store_k(temp_k, new_node->k, s, dense_indices, num_dense);

        //Check if new data should be written to disk
        if (print_flag)
        {
            //while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
            {
                /*
                                //Don't write anything if using data assimilation and at a time when data is available
                                if(GlobalVars->assim_flag)
                                {
                                    double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                                    if(rounded < 1e-13 && -rounded < 1e-13)		break;
                                }
                */

                if (link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(theta, link_i->method->b_theta);

                //v_copy(node->y_approx,sum);
                //for(l=0;l<link_i->method->s;l++)
                //	daxpy(h * link_i->method->b_theta.ve[l],node->next->k[l],sum,1);
                for (m = 0; m < num_dense; m++)
                {
                    idx = dense_indices[m];
                    double approx = v_at(node->y_approx, idx);
                    for (l = 0; l < link_i->method->s; l++)
                        approx += h * v_at(link_i->method->b_theta, l) * v2_at(node->next->k, l, m);

                    v_set(sum, idx, approx);
                }
                link_i->algebraic(sum, GlobalVars->global_params, link_i->params, link_i->qvs, link_i->my->list.tail->state, link_i->user, sum);
                link_i->check_consistency(sum, params, GlobalVars->global_params);

                //Write to a file
                //fsetpos(outputfile,&(link_i->pos));
                //WriteStep(link_i->next_save,sum,GlobalVars,params,link_i->state,outputfile,link_i->output_user,&(link_i->pos));
                WriteStep(outputfile, link_i->ID, link_i->next_save, sum, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
                //fgetpos(outputfile,&(link_i->pos));
/*
                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                    fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));
*/

/*
                #ifdef PRINT2DATABASE
                    sprintf(conninfo->query,"%i,%u,%.12e,%.12e\n",link_i->ID,(unsigned int)(link_i->next_save * 60.0) + conninfo->time_offset,suv_at(m, 0)/link_i->Q_TM,suv_at(m, 0));
                    unsigned int length = strlen(conninfo->query);
                    (conninfo->submission_content) += length;
                    while(conninfo->submission_content > conninfo->submission_size)
                    {
                        //Allocate more space
                        (conninfo->submission_size) *= 2;
                        conninfo->submission = realloc(conninfo->submission,conninfo->submission_size);
                    }
                    strcat(conninfo->submission,conninfo->query);
                #else //Write to a file
                    fsetpos(outputfile,&(link_i->pos));
                    fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                    for(j=0;j<num_print;j++)
                        fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                    fgetpos(outputfile,&(link_i->pos));
                #endif
*/

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (link_i->peak_flag && (v_at(new_y, 0) > v_at(link_i->peak_value, 0)))
        {
            v_copy_n(new_y, link_i->peak_value, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8)
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_buff[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_buff[j]->nrows; i++)
                {
                    //if(link_i->forcing_buff[j]->rainfall[i][1] != forcing_buffer)
                    if (fabs(link_i->forcing_buff[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i - 1][0];
            }
        }

        //Select new step size, if forcings changed
        if (propagated)	link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);

        //Free up parents' old data
        for (i = 0; i < link_i->num_parents; i++)
        {
            currentp = link_i->parents[i];
            while (currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(&currentp->my->list);
                currentp->current_iterations--;
                currentp->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}


//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//This method should be used for links that have a dam.
//Link* link_i: the link to apply a numerical method to.
//VEC* sum: some space for temporary calculations. Should have same dimension as number of equations.
//VEC* temp: same as sum.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKIndex1SolverDam(Link* link_i, GlobalVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l, m, idx;
    VEC new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    Link* currentp;
    double t_needed, timediff, current_theta;

    //Some variables to make things easier to read
    VEC y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2 A = link_i->method->A;
    VEC b = link_i->method->b;
    VEC c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC params = link_i->params;
    VEC e = link_i->method->e;
    VEC d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    const unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = GlobalVars->num_outputs;
    unsigned int* print_indices = GlobalVars->print_indices;
    VEC temp = workspace->temp;
    VEC sum = workspace->sum;
    VEC2 temp_k = workspace->temp_k;

    //Get the approximate solutions from each parent
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;
        for (j = 0; j < s; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j) * h, currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / timediff;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            VEC parent_approx = v3_slice2(workspace->temp_parent_approx, j, i);

            //v_copy(curr_node[i]->y_approx,temp_parent_approx[j][i]);
            //for(l=0;l<currentp->method->s;l++)
            //	daxpy(timediff*currentp->method->b_theta.ve[l],curr_node[i]->next->k[l],temp_parent_approx[j][i],1);
            for (m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                double approx = v_at(curr_node[i]->y_approx, idx);
                for (l = 0; l < currentp->method->s; l++)
                    approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                v_set(parent_approx, idx, approx);
            }

            //Build the algebraic variable
            currentp->check_consistency(parent_approx, currentp->params, GlobalVars->global_params);
            //link_i->algebraic(temp_parent_approx[j][i],GlobalVars->global_params,currentp->params,currentp->qvs,curr_node[i]->next->state,temp_parent_approx[j][i]);
            if (currentp->algebraic)
                currentp->algebraic(parent_approx, GlobalVars->global_params, currentp->params, currentp->qvs, curr_node[i]->state, currentp->user, parent_approx);
        }
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    new_y = new_node->y_approx;

    //Compute the k's
    for (i = 0; i < s; i++)
    {
        v_copy_n(y_0, sum, link_i->dim);
        for (j = 0; j < i; j++)
            daxpy_u(h * v2_at(A, i, j), v2_slice(temp_k, j), sum, 1, link_i->dim);
        link_i->check_consistency(sum, params, GlobalVars->global_params);
        link_i->differential(t + v_at(c, i) * h, sum, v3_slice(workspace->temp_parent_approx, i), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, v2_slice(temp_k, i));
    }

    //Build the solution
    v_copy_n(y_0, new_y, link_i->dim);
    for (i = 0; i < s; i++)	daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 1, link_i->dim);
    link_i->check_consistency(new_y, params, GlobalVars->global_params);
    new_node->state = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);
    link_i->algebraic(new_y, GlobalVars->global_params, params, link_i->qvs, new_node->state, link_i->user, new_y);

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    sv_mlt_u(h * v_at(e, 0), sum, 1, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 1, link_i->dim);

    //Build SC_i
    for (i = 1; i < dim; i++)        
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * v_at(error->reltol, i) + v_at(error->abstol, i)));

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,1);
    err_1 = norm_inf_u(sum, temp, 1, link_i->dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);


    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    sv_mlt_u(h * v_at(d, 0), sum, 1, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 1, link_i->dim);

    for (i = 1; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * v_at(error->reltol_dense, i) + v_at(error->abstol_dense, i)));

    //err_d = norm_inf(sum,temp,meth->d_order_ratio,1);
    err_d = norm_inf_u(sum, temp, 1, link_i->dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h*min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);
    //link_i->h = step_1;

    if (err_1 < 1.0 && err_d < 1.0)
        //if(err_1 < 1.0)
    {
        //Check for issues with state discontinuities
        if (link_i->rejected == 2)
        {
            int old_state = link_i->state;
            link_i->state = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);

            if (old_state != link_i->state)
            {
                //Pass the time to the next link
                double xh = t + h;
                Link* next = link_i->child;
                Link* prev = link_i;

                for (i = 1; i < GlobalVars->max_localorder && next != NULL; i++)
                {
                    if (assignments[next->location] == my_rank && i < next->method->localorder)
                    {
                        //Insert the time into the discontinuity list
                        next->discont_end = Insert_Discontinuity(xh, next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                    }
                    else if (assignments[next->location] != my_rank)
                    {
                        //Store the time to send to another process
                        Insert_SendDiscontinuity(xh, i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                        break;
                    }

                    prev = next;
                    next = next->child;
                }

                link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
            }
        }
        else
        {
            //Check for discontinuities
            unsigned int sign_l, sign_h, sign_m;
            //sign_l = link_i->state_check(y_0,GlobalVars->global_params,params,link_i->qvs,link_i->dam);
            sign_l = link_i->state;
            sign_h = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);
            new_node->state = sign_h;
            if (sign_l != sign_h)
            {
                double xl = t;
                double xh = t + h;
                double xm = (xl + xh) / 2.0;

                //Form the derivative of the interpolant
                meth->dense_bderiv(1.0, meth->b_theta_deriv);
                v_zero(sum);
                for (i = 0; i < s; i++)
                    daxpy_u(h * v_at(meth->b_theta_deriv, i), v2_slice(temp_k, i), sum, 0, link_i->dim);

                //Get the approximate solutions from each parent
                for (i = 0; i < link_i->num_parents; i++)
                {
                    currentp = link_i->parents[i];
                    curr_node[i] = currentp->my->list.head;

                    //Find the needed value of t and corresponding node for y_p
                    //Assuming everything needed is already calculated
                    t_needed = min(t + h, currentp->last_t);
                    //t_needed = t + h;

                    //Find the corresponding theta value and approximate solution
                    while (t_needed > curr_node[i]->t)
                        curr_node[i] = curr_node[i]->next;
                    if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

                    timediff = curr_node[i]->next->t - curr_node[i]->t;
                    current_theta = (t_needed - curr_node[i]->t) / timediff;
                    currentp->method->dense_b(current_theta, currentp->method->b_theta);

                    VEC parent_approx = v3_slice2(workspace->temp_parent_approx, 0, i);

                    for (m = 0; m < currentp->num_dense; m++)
                    {
                        idx = currentp->dense_indices[m];
                        double approx = v_at(curr_node[i]->y_approx, idx);
                        for (l = 0; l < currentp->method->s; l++)
                            approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                        v_set(parent_approx, idx, approx);
                    }

                    if (currentp->algebraic)
                        currentp->algebraic(parent_approx, GlobalVars->global_params, currentp->params, currentp->qvs, curr_node[i]->next->state, currentp->user, parent_approx);
                }

                //Exact derivative at time t + h
                link_i->differential(t + h, new_y, v3_slice(workspace->temp_parent_approx, 0), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp);

                //Form the defect for a stopping criteria
                v_sub(sum, temp, temp, 1);
                double norm_defect = vector_norminf(temp, 1);
                double prev_error;
                double curr_error = (xh - xl)*norm_defect;

                //Locate the discontinuity
                do
                {
                    prev_error = curr_error;

                    //Build the solution at time xm
                    v_copy_n(y_0, sum, link_i->dim);
                    current_theta = (xm - t) / h;
                    meth->dense_b(current_theta, meth->b_theta);
                    for (i = 0; i < s; i++)
                        daxpy_u(h * v_at(meth->b_theta, i), v2_slice(temp_k, i), sum, 1, link_i->dim);
                    
                    link_i->algebraic(sum, GlobalVars->global_params, link_i->params, link_i->qvs, new_node->state, link_i->user, sum);
                    sign_m = link_i->state_check(sum, GlobalVars->global_params, params, link_i->qvs, link_i->dam);

                    if (sign_l != sign_m)
                    {
                        sign_h = sign_m;
                        xh = xm;
                    }
                    else if (sign_m != sign_h)
                    {
                        sign_l = sign_m;
                        xl = xm;
                    }
                    else
                        printf("Uh oh....\n");

                    xm = (xl + xh) / 2.0;
                    curr_error = (xh - xl)*norm_defect;

                } while (curr_error > GlobalVars->discont_tol  &&  fabs(curr_error - prev_error) > 1e-13);

                //Prepare to redo the step
                if (curr_error <= GlobalVars->discont_tol)
                {
                    //Prepare to step on the time of the discontinuity next iteration
                    link_i->h = xh - t;
                    Undo_Step(&link_i->my->list);

                    return 2;
                }
                else
                {
                    //printf("[%i]: Notice: Discontinuity error has converged. Rejecting current step.\n",my_rank);
                    //printf("[%i]: ID = %u t = %f old h = %e current error = %.16e  previous error = %.16e\n",my_rank,link_i->ID,t,h,curr_error,prev_error);

                    //Cut the step size and redo the iteration
                    link_i->h = h * 0.5;
                    Undo_Step(&link_i->my->list);

                    return 0;
                }
            }
        }

        //Check if a propagated discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(temp_k, new_node->k, s, dense_indices, num_dense);

        //Check if new data should be written to disk
        if (print_flag)
        {
            //while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-8))
            {
                /*
                                //Don't write anything if using data assimilation and at a time when data is available
                                if(GlobalVars->assim_flag)
                                {
                                    double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                                    if(rounded < 1e-13 && -rounded < 1e-13)		break;
                                }
                */
                if (link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);

                //v_copy(node->y_approx,sum);
                //for(l=0;l<link_i->method->s;l++)
                //	daxpy(h * link_i->method->b_theta.ve[l],node->next->k[l],sum,1);
                for (m = 0; m < num_dense; m++)
                {
                    idx = dense_indices[m];
                    double approx = v_at(node->y_approx, idx);
                    for (l = 0; l < link_i->method->s; l++)
                        approx += h * v_at(link_i->method->b_theta, l) * v2_at(node->next->k, l, m);

                    v_set(sum, idx, approx);
                }
                link_i->algebraic(sum, GlobalVars->global_params, link_i->params, link_i->qvs, link_i->my->list.tail->state, link_i->user, sum);
                link_i->check_consistency(sum, params, GlobalVars->global_params);

                //Write to a file
                //fsetpos(outputfile,&(link_i->pos));
                //WriteStep(link_i->next_save,sum,GlobalVars,params,link_i->state,outputfile,link_i->output_user,&(link_i->pos));
                WriteStep(outputfile, link_i->ID, link_i->next_save, sum, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
                //fgetpos(outputfile,&(link_i->pos));
/*
                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                    fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));
*/
/*
                #ifdef PRINT2DATABASE
                    sprintf(conninfo->query,"%i,%u,%.12e,%.12e\n",link_i->ID,(unsigned int)(link_i->next_save * 60.0) + conninfo->time_offset,suv_at(m, 0)/link_i->Q_TM,suv_at(m, 0));
                    unsigned int length = strlen(conninfo->query);
                    (conninfo->submission_content) += length;
                    while(conninfo->submission_content > conninfo->submission_size)
                    {
                        //Allocate more space
                        (conninfo->submission_size) *= 2;
                        conninfo->submission = realloc(conninfo->submission,conninfo->submission_size);
                    }
                    strcat(conninfo->submission,conninfo->query);
                #else //Write to a file
                    fsetpos(outputfile,&(link_i->pos));
                    fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                    for(j=0;j<num_print;j++)
                        fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                    fgetpos(outputfile,&(link_i->pos));
                #endif
*/

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (link_i->peak_flag && (v_at(new_y, 0) > v_at(link_i->peak_value, 0)))
        {
            v_copy_n(new_y, link_i->peak_value, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8)
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_buff[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_buff[j]->nrows; i++)
                {
                    //if(link_i->forcing_buff[j]->rainfall[i][1] != forcing_buffer)
                    if (fabs(link_i->forcing_buff[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i - 1][0];
            }
        }

        //Select new step size, if forcings changed
        if (propagated)	link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);

        //Free up parents' old data
        for (i = 0; i < link_i->num_parents; i++)
        {
            currentp = link_i->parents[i];
            while (currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(&currentp->my->list);
                currentp->current_iterations--;
                currentp->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}



//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//This method should be used for links that have a dam.
//Link* link_i: the link to apply a numerical method to.
//VEC* sum: some space for temporary calculations. Should have same dimension as number of equations.
//VEC* temp: same as sum.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolverDiscont(Link* link_i, GlobalVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l, m, idx;
    VEC new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    Link* currentp;
    double t_needed, timediff, current_theta;

    //Some variables to make things easier to read
    VEC y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2 A = link_i->method->A;
    VEC b = link_i->method->b;
    VEC c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC params = link_i->params;
    VEC e = link_i->method->e;
    VEC d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    const unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = GlobalVars->num_outputs;
    unsigned int* print_indices = GlobalVars->print_indices;
    VEC temp = workspace->temp;
    VEC sum = workspace->sum;
    VEC2 temp_k = workspace->temp_k;

    //Get the approximate solutions from each parent
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;
        for (j = 0; j < s; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j) * h, currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / timediff;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            VEC parent_approx = v3_slice2(workspace->temp_parent_approx, j, i);

            //v_copy(curr_node[i]->y_approx,temp_parent_approx[j][i]);
            //for(l=0;l<currentp->method->s;l++)
            //	daxpy(timediff*currentp->method->b_theta.ve[l],curr_node[i]->next->k[l],temp_parent_approx[j][i],1);
            for (m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                double approx = v_at(curr_node[i]->y_approx, idx);
                for (l = 0; l < currentp->method->s; l++)
                    approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                v_set(parent_approx, idx, approx);
            }

            //Build the algebraic variable
            currentp->check_consistency(parent_approx, currentp->params, GlobalVars->global_params);
        }
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    new_y = new_node->y_approx;

    //Compute the k's
    for (i = 0; i < s; i++)
    {
        v_copy_n(y_0, sum, link_i->dim);
        for (j = 0; j < i; j++)
            daxpy_u(h * v2_at(A, i, j), v2_slice(temp_k, j), sum, 0, link_i->dim);
        link_i->check_consistency(sum, params, GlobalVars->global_params);
        link_i->differential(t + v_at(c, i) * h, sum, v3_slice(workspace->temp_parent_approx, i), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, v2_slice(temp_k, i));
    }

    //Build the solution
    v_copy_n(y_0, new_y, link_i->dim);
    for (i = 0; i < s; i++)	daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 0, link_i->dim);
    link_i->check_consistency(new_y, params, GlobalVars->global_params);
    new_node->state = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    sv_mlt_u(h * v_at(e, 0), sum, 0, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 0, link_i->dim);

    //Build SC_i
    for (i = 0; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * v_at(error->reltol, i) + v_at(error->abstol, i)));

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,1);
    err_1 = norm_inf_u(sum, temp, 0, link_i->dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);


    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    sv_mlt_u(h * v_at(d, 0), sum, 0, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
    for (i = 0; i < dim; i++)        
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * v_at(error->reltol_dense, i) + v_at(error->abstol_dense, i)));

    //err_d = norm_inf(sum,temp,meth->d_order_ratio,1);
    err_d = norm_inf_u(sum, temp, 0, link_i->dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h*min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);
    //link_i->h = step_1;

    if (err_1 < 1.0 && err_d < 1.0)
        //if(err_1 < 1.0)
    {
        //Check for issues with state discontinuities
        if (link_i->rejected == 2)
        {
            int old_state = link_i->state;
            link_i->state = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);

            if (old_state != link_i->state)
            {
                //Pass the time to the next link
                double xh = t + h;
                Link* next = link_i->child;
                Link* prev = link_i;

                for (i = 1; i < GlobalVars->max_localorder && next != NULL; i++)
                {
                    if (assignments[next->location] == my_rank && i < next->method->localorder)
                    {
                        //Insert the time into the discontinuity list
                        next->discont_end = Insert_Discontinuity(xh, next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                    }
                    else if (assignments[next->location] != my_rank)
                    {
                        //Store the time to send to another process
                        Insert_SendDiscontinuity(xh, i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                        break;
                    }

                    prev = next;
                    next = next->child;
                }

                link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
            }
        }
        else
        {
            //Check for discontinuities
            unsigned int sign_l, sign_h, sign_m;
            //sign_l = link_i->state_check(y_0,GlobalVars->global_params,params,link_i->qvs,link_i->dam);
            sign_l = link_i->state;
            sign_h = link_i->state_check(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->dam);
            new_node->state = sign_h;
            if (sign_l != sign_h)
            {
                double xl = t;
                double xh = t + h;
                double xm = (xl + xh) / 2.0;

                //Form the derivative of the interpolant
                meth->dense_bderiv(1.0, meth->b_theta_deriv);
                v_zero(sum);
                for (i = 0; i < s; i++)
                    daxpy_u(h * v_at(meth->b_theta_deriv, i), v2_slice(temp_k, i), sum, 0, link_i->dim);

                //Get the approximate solutions from each parent
                for (i = 0; i < link_i->num_parents; i++)
                {
                    currentp = link_i->parents[i];
                    curr_node[i] = currentp->my->list.head;

                    //Find the needed value of t and corresponding node for y_p
                    //Assuming everything needed is already calculated
                    t_needed = min(t + h, currentp->last_t);
                    //t_needed = t + h;

                    //Find the corresponding theta value and approximate solution
                    while (t_needed > curr_node[i]->t)
                        curr_node[i] = curr_node[i]->next;
                    if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

                    timediff = curr_node[i]->next->t - curr_node[i]->t;
                    current_theta = (t_needed - curr_node[i]->t) / timediff;
                    currentp->method->dense_b(current_theta, currentp->method->b_theta);

                    VEC parent_approx = v3_slice2(workspace->temp_parent_approx, 0, i);

                    for (m = 0; m < currentp->num_dense; m++)
                    {
                        idx = currentp->dense_indices[m];
                        double approx  = v_at(curr_node[i]->y_approx, idx);
                        for (l = 0; l < currentp->method->s; l++)
                            approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                        v_set(parent_approx, idx, approx);
                    }
                }

                //Exact derivative at time t + h
                link_i->differential(t + h, new_y, v3_slice(workspace->temp_parent_approx, 0), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp);

                //Form the defect for a stopping criteria
                v_sub(sum, temp, temp, 0);
                double norm_defect = vector_norminf(temp, 0);
                double prev_error;
                double curr_error = (xh - xl)*norm_defect;

                //Locate the discontinuity
                do
                {
                    prev_error = curr_error;

                    //Build the solution at time xm
                    v_copy_n(y_0, sum, link_i->dim);
                    current_theta = (xm - t) / h;
                    meth->dense_b(current_theta, meth->b_theta);
                    for (i = 0; i < s; i++)
                        daxpy_u(h * v_at(meth->b_theta, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
                    sign_m = link_i->state_check(sum, GlobalVars->global_params, params, link_i->qvs, link_i->dam);

                    if (sign_l != sign_m)
                    {
                        sign_h = sign_m;
                        xh = xm;
                    }
                    else if (sign_m != sign_h)
                    {
                        sign_l = sign_m;
                        xl = xm;
                    }
                    else
                        printf("Uh oh....\n");

                    xm = (xl + xh) / 2.0;
                    curr_error = (xh - xl)*norm_defect;

                } while (curr_error > GlobalVars->discont_tol  &&  fabs(curr_error - prev_error) > 1e-13);

                //Prepare to redo the step
                if (curr_error <= GlobalVars->discont_tol)
                {
                    //Prepare to step on the time of the discontinuity next iteration
                    link_i->h = xh - t;
                    Undo_Step(&link_i->my->list);

                    return 2;
                }
                else
                {
                    //printf("[%i]: Notice: Discontinuity error has converged. Rejecting current step.\n",my_rank);
                    //printf("[%i]: ID = %u t = %f old h = %e current error = %.16e  previous error = %.16e\n",my_rank,link_i->ID,t,h,curr_error,prev_error);

                    //Cut the step size and redo the iteration
                    link_i->h = h * 0.5;
                    Undo_Step(&link_i->my->list);

                    return 0;
                }
            }
        }

        //Check if a propagated discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(temp_k, new_node->k, s, dense_indices, num_dense);

        //Check if new data should be written to disk
        if (print_flag)
        {
            //while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
            {
                /*
                                //Don't write anything if using data assimilation and at a time when data is available
                                if(GlobalVars->assim_flag)
                                {
                                    double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                                    if(rounded < 1e-13 && -rounded < 1e-13)		break;
                                }
                */
                if (link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);

                //v_copy(node->y_approx,sum);
                //for(l=0;l<link_i->method->s;l++)
                //	daxpy(h * link_i->method->b_theta.ve[l],node->next->k[l],sum,0);
                for (m = 0; m < num_dense; m++)
                {
                    idx = dense_indices[m];
                    double approx = v_at(node->y_approx, idx);
                    for (l = 0; l < link_i->method->s; l++)
                        approx += h * v_at(link_i->method->b_theta, l) * v2_at(node->next->k, l, m);

                    v_set(sum, idx, approx);
                }
                link_i->check_consistency(sum, params, GlobalVars->global_params);

                //Write to a file
                //fsetpos(outputfile,&(link_i->pos));
                //WriteStep(link_i->next_save,sum,GlobalVars,params,link_i->state,outputfile,link_i->output_user,&(link_i->pos));
                WriteStep(outputfile, link_i->ID, link_i->next_save, sum, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
                //fgetpos(outputfile,&(link_i->pos));
/*
                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                    fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));
*/
/*
                #ifdef PRINT2DATABASE
                    sprintf(conninfo->query,"%i,%u,%.12e,%.12e\n",link_i->ID,(unsigned int)(link_i->next_save * 60.0) + conninfo->time_offset,suv_at(m, 0)/link_i->Q_TM,suv_at(m, 0));
                    unsigned int length = strlen(conninfo->query);
                    (conninfo->submission_content) += length;
                    while(conninfo->submission_content > conninfo->submission_size)
                    {
                        //Allocate more space
                        (conninfo->submission_size) *= 2;
                        conninfo->submission = realloc(conninfo->submission,conninfo->submission_size);
                    }
                    strcat(conninfo->submission,conninfo->query);
                #else //Write to a file
                    fsetpos(outputfile,&(link_i->pos));
                    fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                    for(j=0;j<num_print;j++)
                        fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                    fgetpos(outputfile,&(link_i->pos));
                #endif
*/
                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (link_i->peak_flag && (v_at(new_y, 0) > v_at(link_i->peak_value, 0)))
        {
            v_copy_n(new_y, link_i->peak_value, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8)
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_buff[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_buff[j]->nrows; i++)
                {
                    //if(link_i->forcing_buff[j]->rainfall[i][1] != forcing_buffer)
                    if (fabs(link_i->forcing_buff[j]->data[i][1] - forcing_buffer) < 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i - 1][0];
            }
        }

        //Select new step size, if forcings changed
        if (propagated)	link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);

        //Free up parents' old data
        for (i = 0; i < link_i->num_parents; i++)
        {
            currentp = link_i->parents[i];
            while (currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(&currentp->my->list);
                currentp->current_iterations--;
                currentp->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}


#if defined(ASYNCH_HAVE_IMPLICIT_SOLVER)

//Computes one step of an implicit RK method to solve the ODE at a link. Assumes parents have enough computed solutions.
//Link* link_i: the link to apply a numerical method to.
//VEC* sum: some space for temporary calculations. Should have same dimension as number of equations.
//VEC* temp: same as sum.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int RadauRKSolver(Link* link_i,UnivVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace)
{
    unsigned int i,j,l,m;
    VEC** k;
    VEC* new_y;
    RKSolutionNode *curr_node[link_i->num_parents],*node,*new_node;
    Link* currentp;
    double t_needed,timediff,current_theta;

    //Some variables to make things easier to read
    VEC* y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2* A = link_i->method->A;
    VEC* b = link_i->method->b;
    VEC* c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC* params = link_i->params;
    VEC* e = link_i->method->e;
    VEC* d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = link_i->error_data;
    const unsigned int dim = link_i->dim;
    VEC* temp = workspace->temp;
    VEC* sum = workspace->sum;
    VEC*** temp_parent_approx = workspace->temp_parent_approx;

    //Get the approximate solutions from each parent
    for(i=0;i<link_i->num_parents;i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;

        //Find the needed value of t_0 and corresponding node for y_p for the extra stage
        t_needed = t;

        //Find the corresponding theta value and approximate solution
        while(t_needed > curr_node[i]->t)
            curr_node[i] = curr_node[i]->next;
        if(curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

        timediff = curr_node[i]->next->t - curr_node[i]->t;
        current_theta = (t_needed-curr_node[i]->t)/timediff;
        currentp->method->dense_b(current_theta,currentp->method->b_theta);

        v_copy(curr_node[i]->y_approx,temp_parent_approx[s][i]);
        for(l=0;l<currentp->method->s;l++)
            daxpy(timediff*currentp->method->b_theta.ve[l],curr_node[i]->next->k[l],temp_parent_approx[s][i],0);


        //currentp = link_i->parents[i];
        //curr_node[i] = currentp->my->list.head;
        for(j=0;j<s;j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j)*h,currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while(t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if(curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed-curr_node[i]->t)/timediff;
            currentp->method->dense_b(current_theta,currentp->method->b_theta);

            v_copy(curr_node[i]->y_approx,temp_parent_approx[j][i]);
            for(l=0;l<currentp->method->s;l++)
                daxpy(timediff*currentp->method->b_theta.ve[l],curr_node[i]->next->k[l],temp_parent_approx[j][i],0);
        }
    }

    //Do the RK method to get the next approximation

    //Allocate space. !!!! Do this elsewhere !!!!
    int info;
    int* ipiv = workspace->ipiv;
    VEC* RHS = workspace->RHS;
//	VEC2* CoefMat = workspace->CoefMat;
    VEC2* CoefMat = link_i->CoefMat;
//	VEC2* JMatrix = workspace->JMatrix;
    VEC2* tempmat = workspace->JMatrix;
    VEC** Z_i = workspace->Z_i;
    VEC** old_Z = link_i->Z_i;
    VEC* err = workspace->err;
    unsigned int iters = 10;
    double normDelta = -1.0;
    double normDeltam1 = -1.0;
    double Theta,eta,error_est;
    unsigned short int abort = 0;
    unsigned short int newt = 0;
    unsigned short int newtmax = 20;

    //Construct a guess for the nonlinear solver
    if(link_i->h_old < 0.0)
        for(i=0;i<s;i++)	for(j=0;j<dim;j++)	Z_i[i].ve[j] = 0.0;
    else
    {
        for(i=0;i<s;i++)
        {
            lagrange_bary(1.0 + h/link_i->h_old*v_at(c, i),c,old_Z,meth->w,Z_i[i]);
            v_add(Z_i[i],link_i->sol_diff,Z_i[i],0);
            if(Z_i[i].ve[0] + y_v_at(0, 0) < 0.0)	Z_i[i].ve[0] = 0.0 - y_v_at(0, 0);
            if(dim > 1 && Z_i[i].ve[1] + y_0.ve[1] < 0.0)	Z_i[i].ve[1] = -y_0.ve[1];
            if(dim > 2 && Z_i[i].ve[2] + y_0.ve[2] < 0.0)	Z_i[i].ve[2] = -y_0.ve[2];
            if(dim > 2 && Z_i[i].ve[2] + y_0.ve[2] > 1.0)	Z_i[i].ve[2] = 1.0 - y_0.ve[2] - 0.0;
            if(dim > 3 && Z_i[i].ve[3] + y_0.ve[3] < 0.0)	Z_i[i].ve[3] = -y_0.ve[3];
            if(dim > 3 && Z_i[i].ve[3] + y_0.ve[3] > 1.0 - Z_i[i].ve[2] - y_0.ve[2])	Z_i[i].ve[3] = 1.0 - Z_i[i].ve[2] - y_0.ve[2] - y_0.ve[3];
        }
    }

    //Build the coefficient matrix, if necessary
    if(link_i->compute_J == 1)
        link_i->jacobian(t,y_0,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,params,link_i->JMatrix);

    if(link_i->compute_LU == 1)
    {
        for(i=0;i<s;i++) for(j=0;j<s;j++)
        {
            for(l=0;l<dim;l++) for(m=0;m<dim;m++)
                CoefMat.me[i*dim + l][j*dim + m] = -h*A.me[i][j]*link_i->JMatrix.me[l][m];
        }
        for(i=0;i<s*dim;i++)	CoefMat.me[i][i] += 1.0;
        info = clapack_dgetrf(CblasRowMajor,s*dim,s*dim,CoefMat->array,CoefMat->m,ipiv);
    }

    //Perform the simplified Newton iterations
    do
    {
        //Initialize RHS
        for(j=0;j<s;j++) for(l=0;l<dim;l++)
            RHS.ve[l+j*dim] = -Z_i[j].ve[l];

        //Evaluate the guess in the system created by the numerical method
        for(j=0;j<s;j++)
        {
            v_add(y_0,Z_i[j],sum,0);
            if(suv_at(m, 0) < 0.0)	suv_at(m, 0) = 0.0;	//A negative discharge is clearly a poor approximation
            if(dim > 1 && sum.ve[1] < 0.0)	sum.ve[1] = 0.0;
            if(dim > 2 && sum.ve[2] < 0.0)	sum.ve[2] = 0.0;
            if(dim > 2 && sum.ve[2] > 1.0)	sum.ve[2] = 1.0;
            if(dim > 3 && sum.ve[3] < 0.0)	sum.ve[3] = 0.0;
            if(dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];
            link_i->f(t+v_at(c, j)*h,sum,temp_parent_approx[j],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->state,link_i->user,temp);
            for(l=0;l<s;l++)
                daxpy(h*A.me[l][j],temp,RHS,l*dim);
        }

        //Make sure RHS has valid values
        for(j=0;j<s;j++)
        {
            if(RHS.ve[j*dim] + (y_v_at(0, 0) + Z_i[j].ve[0]) < 0.0)	RHS.ve[j*dim] = -(y_v_at(0, 0) + Z_i[j].ve[0]);
            if(dim > 1 && RHS.ve[j*dim + 1] + (y_0.ve[1] + Z_i[j].ve[1]) < 0.0)	RHS.ve[j*dim + 1] = -(y_0.ve[1] + Z_i[j].ve[1]);
            if(dim > 2 && RHS.ve[j*dim + 2] + (y_0.ve[2] + Z_i[j].ve[2]) < 0.0)	RHS.ve[j*dim + 2] = -(y_0.ve[2] + Z_i[j].ve[2]);
            if(dim > 2 && RHS.ve[j*dim + 2] + (y_0.ve[2] + Z_i[j].ve[2]) > 1.0)	RHS.ve[j*dim + 2] = 1.0 -(y_0.ve[2] + Z_i[j].ve[2]);
            if(dim > 3 && RHS.ve[j*dim + 3] + (y_0.ve[3] + Z_i[j].ve[3]) < 0.0)	RHS.ve[j*dim + 3] = -(y_0.ve[3] + Z_i[j].ve[3]);
            if(dim > 3 && RHS.ve[j*dim + 3] + (y_0.ve[3] + Z_i[j].ve[3]) > 1.0 - RHS.ve[j*dim + 2] - (y_0.ve[2] + Z_i[j].ve[2]))	RHS.ve[j*dim + 3] = 1.0 - RHS.ve[j*dim + 2] - (y_0.ve[2] + Z_i[j].ve[2]) - (y_0.ve[3] + Z_i[j].ve[3]);
        }

        //Solve the system
        info = clapack_dgetrs(CblasRowMajor,111,s*dim,1,CoefMat->array,CoefMat->m,ipiv,RHS.ve,s*dim);

        //Form the solution to the iteration
        for(j=0;j<s;j++)	v_add(RHS,Z_i[j],Z_i[j],j*dim);

        //Make sure Z_i has valid values
        for(j=0;j<s;j++)
        {
            if(Z_i[j].ve[0] + y_v_at(0, 0) < 0.0)	Z_i[j].ve[0] = 0.0 - y_v_at(0, 0);
            if(dim > 1 && Z_i[j].ve[1] + y_0.ve[1] < 0.0)	Z_i[j].ve[1] = -y_0.ve[1];
            if(dim > 2 && Z_i[j].ve[2] + y_0.ve[2] < 0.0)	Z_i[j].ve[2] = -y_0.ve[2];
            if(dim > 2 && Z_i[j].ve[2] + y_0.ve[2] > 1.0)	Z_i[j].ve[2] = 1.0 - y_0.ve[2] - 0.0;
            if(dim > 3 && Z_i[j].ve[3] + y_0.ve[3] < 0.0)	Z_i[j].ve[3] = -y_0.ve[3];
            if(dim > 3 && Z_i[j].ve[3] + y_0.ve[3] > 1.0 - Z_i[j].ve[2] - y_0.ve[2])	Z_i[j].ve[3] = 1.0 - Z_i[j].ve[2] - y_0.ve[2] - y_0.ve[3];
        }

        //Check if solution is accurate enough
        normDeltam1 = normDelta;
        normDelta = vector_norminf(RHS);
        Theta = normDelta / normDeltam1;
        if(normDeltam1 < 0.0)	eta = link_i->last_eta;
        else			eta = Theta / (1.0 - Theta);

        //Check that the error estimation is roughly correct
        //error_est = pow(Theta,imax-i) / (1.0 - Theta) * normDelta;

        //Check abort conditions
        newt++;
        if(Theta >= 1.0 || newt == newtmax)	abort = 1;
        else			abort = 0;
    } while( eta*normDelta > 1e-10 && abort == 0 );

    //Check that the nonlinear solver converged
    if(abort)
    {
        link_i->h = link_i->h / 2;
        if(link_i->ID == 10)	printf("Solver diverged. ID: %u Time: %f New Step: %.16f\n",link_i->ID,link_i->last_t,link_i->h);
        link_i->compute_J = 1;
        link_i->compute_LU = 1;
        return 0;
    }
    else
        link_i->last_eta = pow(eta,.8);

    //Check if jacobian should be recomputed
    if(normDeltam1 < 0.0 || Theta < 1e-3)
        link_i->compute_J = 0;
    else
        link_i->compute_J = 1;

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    k = new_node->k;
    new_y = new_node->y_approx;

    //Build the solution
    v_add(y_0,Z_i[s-1],new_y,0);
    if(new_v_at(y, 0) < 0.0)	new_v_at(y, 0) = 1e-15;	//A negative discharge is clearly a poor approximation
    if(dim > 1 && new_y.ve[1] < 0.0)	new_y.ve[1] = 0.0;
    if(dim > 2 && new_y.ve[2] < 0.0)	new_y.ve[2] = 0.0;
    if(dim > 2 && new_y.ve[2] > 1.0)	new_y.ve[2] = 1.0 - 1e-15;
    if(dim > 3 && new_y.ve[3] < 0.0)	new_y.ve[3] = 1e-15;
    if(dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];

    //Error estimation and step size selection

    //Compute LU factorization of (I - h e[s] JMatrix)
    dipaa(-h*e.ve[s],link_i->JMatrix,tempmat,0,0);
    info = clapack_dgetrf(CblasRowMajor,dim,dim,tempmat->array,dim,ipiv);

    //Build err for y_1
    v_copy(Z_i[0],sum);
    sv_mlt(v_at(e, 0),sum,0);
    for(i=1;i<s;i++)	daxpy(v_at(e, i),Z_i[i],sum,0);
    v_copy(sum,err);
    link_i->f(t,y_0,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->state,link_i->user,temp);	//Compute the extra k
    daxpy(h*e.ve[s],temp,err,0);	//Extra k
    info = clapack_dgetrs(CblasRowMajor,111,dim,1,tempmat->array,dim,ipiv,err.ve,dim);

    if(link_i->rejected == 0)	//If the previous step was rejected
    {
        v_add(y_0,err,temp,0);
        v_copy(sum,err);
        link_i->f(t,temp,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->state,link_i->user,temp);
        daxpy(h*e.ve[s],temp,err,0);
        info = clapack_dgetrs(CblasRowMajor,111,dim,1,tempmat->array,dim,ipiv,err.ve,dim);
    }

    //Build SC_i for y_1
    for(i=0;i<dim;i++)
        temv_at(p, i) = max(fabs(new_v_at(y, i)),fabs(y_v_at(0, i))) * error->reltov_at(l, i) + error->abstov_at(l, i);

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,0);
    double err_1 = norm_inf(err,temp,1.0,0);
    double value_1 = pow(1.0/err_1,1.0/meth->e_order);

    //Build err for dense output
//	v_copy(Z_i[0],sum);
//	sv_mlt(v_at(d, 0),sum,0);
//	for(i=1;i<s;i++)	daxpy(v_at(d, i),Z_i[i],sum,0);
//	v_copy(sum,err);
//	link_i->f(t,y_0,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->user,temp);	//Compute the extra k
//	daxpy(h*d.ve[s],temp,err,0);	//Extra k
//	info = clapack_dgetrs(CblasRowMajor,111,dim,1,JMatrix->array,dim,ipiv,err.ve,dim);	//This assumes d.ve[s] == e.ve[s]

//	if(link_i->rejected == 0)	//If the previous step was rejected
//	{
//		v_add(y_0,err,temp,0);
//		v_copy(sum,err);
//		link_i->f(t,temp,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->user,temp);
//		daxpy(h*d.ve[s],temp,err,0);
//		info = clapack_dgetrs(CblasRowMajor,111,dim,1,JMatrix->array,dim,ipiv,err.ve,dim);
//	}

    //Build SC_i for dense
//	for(i=0;i<dim;i++)
//		temv_at(p, i) = max(fabs(new_v_at(y, i)),fabs(y_v_at(0, i))) * error->reltol_densv_at(e, i) + error->abstol_densv_at(e, i);

//	double err_d = norm_inf(err,temp,1.0,0);
//	double value_d = pow(1.0/err_d,1.0/meth->d_order);

    //Determine a new step size for the next step
    //double step_1 = h*min(error->facmax,max(error->facmin,error->fac * value_1));
    double step_1 = error->fac * (2*newtmax + 1)/(2*newtmax + newt) * h * value_1;
    double step_log;
    if(link_i->h_old > 0.0)
        step_log = error->fac * h * value_1 * h/link_i->h_old * value_1 / link_i->value_old;
    else
        step_log = 1e5;
//	double step_d = h*min(error->facmax,max(error->facmin,error->fac * value_d));
//	link_i->h = min(step_1,step_d);
//	link_i->h = step_1;

//	if(err_1 < 1.0 && err_d < 1.0)
    if(err_1 < 1.0)
    {
        //Check if an LU factorization should take place in the next step
        if( link_i->compute_J == 1 || h > step_1 || step_1 > 1.2*h )
        {
            link_i->h = min(step_1,step_log);
            link_i->compute_LU = 1;
        }
        else
            link_i->compute_LU = 0;

        //Save the new data
        for(j=0;j<s;j++)
        {
            v_add(y_0,Z_i[j],sum,0);
            if(suv_at(m, 0) < 0.0)	suv_at(m, 0) = 1e-10;	//A negative discharge is clearly a poor approximation
            if(dim > 1 && sum.ve[1] < 0.0)	sum.ve[1] = 0.0;
            if(dim > 2 && sum.ve[2] < 0.0)	sum.ve[2] = 0.0;
            if(dim > 2 && sum.ve[2] > 1.0)	sum.ve[2] = 1.0 - 1e-10;
            if(dim > 3 && sum.ve[3] < 0.0)	sum.ve[3] = 1e-10;
            if(dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];
            link_i->f(t+v_at(c, j)*h,sum,temp_parent_approx[j],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->state,link_i->user,k[j]);
        }
        link_i->last_t = t + h;
        link_i->current_iterations++;
        link_i->h_old = h;
        link_i->value_old = value_1;
        v_sub(y_0,new_y,link_i->sol_diff,0);
        for(i=0;i<s;i++)	v_copy(Z_i[i],old_Z[i]);

        //Check if new data should be written to disk
        if(print_flag)
        {
            while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            {
                //Don't write anything if using data assimilation and at a time when data is available
                if(GlobalVars->assim_flag)
                {
                    double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                    if(rounded < 1e-13 && -rounded < 1e-13)		break;
                }

                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t)/h;
                link_i->method->dense_b(current_theta,link_i->method->b_theta);

                v_copy(node->y_approx,sum);
                for(l=0;l<link_i->method->s;l++)
                    daxpy(h * link_i->method->b_theta.ve[l],node->next->k[l],sum,0);
                if(suv_at(m, 0) < 0.0)	suv_at(m, 0) = 0.0;
                if(dim > 1 && sum.ve[1] < 0.0)	sum.ve[1] = 0.0;
                if(dim > 2 && new_y.ve[2] < 0.0)	new_y.ve[2] = 0.0;
                if(dim > 2 && new_y.ve[2] > 1.0)	new_y.ve[2] = 1.0;
                if(dim > 3 && new_y.ve[3] < 0.0)	new_y.ve[3] = 0.0;
                if(dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];

                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                    fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if(new_v_at(y, 0) > link_i->peak_valuv_at(e, 0))
        {
            v_copy(new_y,link_i->peak_value);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for(j=0;j<GlobalVars->num_forcings;j++)
        {
            if( link_i->forcing_buff[j] && fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-14 )
            {
                //Propagate the discontinuity to downstream links
                if(!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for(i=0;i<GlobalVars->max_localorder && next != NULL;i++)
                    {
                        if(assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j],next->discont_start,next->discont_end,&(next->discont_count),GlobalVars->discont_size,next->discont,next->ID);
                        }
                        else if(next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j],i,&(prev->discont_send_count),GlobalVars->discont_size,prev->discont_send,prev->discont_order_send,prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for(l=link_i->forcing_indices[j]+1;l<link_i->forcing_buff[j]->n_times;l++)
                    if( fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->rainfall[l][0]) < 1e-8 )	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->rainfall[l][1];
                link_i->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                for(i=l+1;i<link_i->forcing_buff[j]->n_times;i++)
                {
                    if(link_i->forcing_buff[j]->rainfall[i][1] != forcing_buffer)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->rainfall[i][0];
                        break;
                    }
                }
                if(i == link_i->forcing_buff[j]->n_times)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->rainfall[i-1][0];
            }
        }

        //Select new step size, if forcings changed
        if(propagated)	link_i->h = InitialStepSize(link_i->last_t,link_i,GlobalVars,workspace);

        //Free up parents' old data
        for(i=0;i<link_i->num_parents;i++)
        {
            currentp = link_i->parents[i];
            while(currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(currentp->list);
                currentp->current_iterations--;
                currentp->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        link_i->h = min(step_1,step_log);
        link_i->compute_J = 1;
        link_i->compute_LU = 1;

        return 0;
    }

}

#endif // defined(ASYNCH_HAVE_IMPLICIT_SOLVER)


//Computes a solution for a location with forced system states.
//Computes a solution at either the last time of the upstream link, or the time when a change in the system state occurs.
//Should return 1.
int ForcedSolutionSolver(Link* link_i, GlobalVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l;
    VEC new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *new_node;
    Link* currentp;
    double t_needed;
    short int change_value = 0;

    //Some variables to make things easier to read
    VEC y_0 = link_i->my->list.tail->y_approx;
    double t = link_i->my->list.tail->t;
    double h = link_i->h;
    unsigned int s = link_i->method->s;
    VEC params = link_i->params;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    const unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = GlobalVars->num_outputs;
    VEC2 temp_k = workspace->temp_k;

    //Find the next time to step on
    t_needed = GlobalVars->maxtime;
    for (i = 0; i < GlobalVars->num_forcings; i++)
    {
        if (forcings[i].active && link_i->forcing_buff[i])
        {
            if (t_needed > link_i->forcing_change_times[i])
            {
                t_needed = link_i->forcing_change_times[i];
                change_value = 1;
            }
        }
    }

    for (i = 0; i < link_i->num_parents; i++)
    {
        if (t_needed > link_i->parents[i]->last_t)
        {
            t_needed = link_i->parents[i]->last_t;
            change_value = 0;
        }
    }

    h = t_needed - t;

    //Setup the current nodes at each parent (for deleting data)
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;

        //Find the corresponding theta value and approximate solution
        while (t_needed > curr_node[i]->t)
            curr_node[i] = curr_node[i]->next;

        if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;
    }

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    new_y = new_node->y_approx;

    //Compute the k's
    v2_zero(temp_k);

    //Build the solution
    if (change_value)
    {
        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && (fabs(new_node->t - link_i->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_buff[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_buff[j]->nrows; i++)
                {
                    if (fabs(link_i->forcing_buff[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i - 1][0];
            }
        }
    }
    link_i->differential(t + h, y_0, v2_init(0, 0), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, new_y);
    if (link_i->state_check)
        new_node->state = link_i->state_check(new_y, GlobalVars->global_params, link_i->params, link_i->qvs, link_i->dam);

    //Set stepsize
    link_i->h = h;

    //Ignore propagated discontinuities
    link_i->discont_count = 0;
    link_i->discont_start = 0;
    link_i->discont_end = GlobalVars->discont_size - 1;

    //Save the new data
    link_i->last_t = t + h;
    link_i->current_iterations++;
    store_k(temp_k, new_node->k, s, dense_indices, num_dense);

    //Check if new data should be written to disk
    if (print_flag)
    {
        while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
        {
            if (link_i->disk_iterations == link_i->expected_file_vals)
            {
                printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                break;
            }
            (link_i->disk_iterations)++;

            //Write to a file
            if (change_value && fabs((link_i->next_save - link_i->last_t) / link_i->next_save) < 1e-12)
                WriteStep(outputfile, link_i->ID, link_i->next_save, new_y, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
            else
                WriteStep(outputfile, link_i->ID, link_i->next_save, y_0, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
            link_i->next_save += link_i->print_time;
        }
    }

    //Check if this is a max discharge
    if (link_i->peak_flag && (v_at(new_y, 0) > v_at(link_i->peak_value, 0)))
    {
        v_copy_n(new_y, link_i->peak_value, link_i->dim);
        link_i->peak_time = link_i->last_t;
    }


    //Check if the newest step is on a change in rainfall
    if (!change_value)
    {
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && (fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_buff[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_buff[j]->nrows; i++)
                {
                    if (fabs(link_i->forcing_buff[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->data[i - 1][0];
            }
        }
    }

    //Free up parents' old data
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        while (currentp->my->list.head != curr_node[i])
        {
            Remove_Head_Node(&currentp->my->list);
            currentp->current_iterations--;
            currentp->iters_removed++;
        }
    }

    //if(link_i->ID == 2)
    //printf("%f %f %u\n",link_i->last_t,link_i->forcing_values[2],link_i->forcing_indices[2]);

    return 1;
}


#if defined(ASYNCH_HAVE_ASSIM_SOLVER)

//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//This is for use with data assimilation. It assumes the vectors are a lot smaller than they really are.
//Link* link_i: the link to apply a numerical method to.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolver_DataAssim(Link* link_i,UnivVars* GlobalVars,int* assignments,bool print_flag,FILE* outputfile,ConnData* conninfo,Forcing* forcings,Workspace* workspace)
{
    unsigned int i,j,l,m,idx;
    VEC* new_y;
    RKSolutionNode *curr_node[link_i->num_parents],*node,*new_node;
    Link* currentp;
    double t_needed,timediff,current_theta;

    //Some variables to make things easier to read
    VEC* y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2* A = link_i->method->A;
    VEC* b = link_i->method->b;
    VEC* c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC* params = link_i->params;
    VEC* e = link_i->method->e;
    VEC* d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = link_i->error_data;
    //const unsigned int dim = GlobalVars->dim;
    //unsigned int num_dense = GlobalVars->num_dense;
    unsigned int* dense_indices = GlobalVars->dense_indices;
    unsigned int num_print = GlobalVars->num_print;
    VEC* temp = workspace->temp;
    VEC* sum = workspace->sum;
    VEC*** temp_parent_approx = workspace->temp_parent_approx;
    VEC** temp_k = workspace->temp_k;
    unsigned int problem_dim = GlobalVars->problem_dim;

    //Calculate total number of upstream links
    unsigned int practical_dim,num_dense_parents[link_i->num_parents],num_upstream_links = 0;
    for(i=0;i<link_i->num_parents;i++)
    {
        num_upstream_links += link_i->numupstream[i];
        num_dense_parents[i] = (link_i->numupstream[i]*problem_dim + 1) * problem_dim;	//!!!! This is actually an upper bound. The extra problem_dim assumes all states are needed. !!!!
    }
    practical_dim = 2*problem_dim + (problem_dim - 1)*(problem_dim - 1) + problem_dim * num_upstream_links;

    //Get the approximate solutions from each parent
    for(i=0;i<link_i->num_parents;i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;

        for(j=0;j<s;j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j)*h,currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while(t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if(curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed-curr_node[i]->t)/timediff;
            currentp->method->dense_b(current_theta,currentp->method->b_theta);

            //for(m=0;m<num_dense;m++)
            for(m=0;m<num_dense_parents[i];m++)
            {
                idx = dense_indices[m];
                temp_parent_approx[j][i].ve[idx] = curr_node[i]->y_approx.ve[idx];
                for(l=0;l<currentp->method->s;l++)
                    temp_parent_approx[j][i].ve[idx] += timediff*currentp->method->b_theta.ve[l] * curr_node[i]->next->k[l].ve[m];
            }
            link_i->check_consistency(temp_parent_approx[j][i],params,GlobalVars->global_params);
//printf("ID = %u parent = %u num_dense_parents = %u\n",link_i->ID,currentp->ID,num_dense_parents[i]);
//Print_Vector(temp_parent_approx[j][i]);
        }
//getchar();
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    new_y = new_node->y_approx;

    //Compute the k's
    for(i=0;i<s;i++)
    {
        v_copy_n(y_0,sum,practical_dim);
        for(j=0;j<i;j++)
            daxpy_u(h*A.me[i][j],v2_slice(temp_k, j),sum,0,practical_dim);
        link_i->check_consistency(sum,params,GlobalVars->global_params);
        link_i->f(t + v_at(c, i) * h,sum,temp_parent_approx[i],link_i->num_parents,GlobalVars->global_params,link_i->forcing_values,link_i->qvs,params,link_i->state,link_i->user,v2_slice(temp_k, i));
    }

    //Build the solution
    v_copy_n(y_0,new_y,practical_dim);
    for(i=0;i<s;i++)	daxpy_u(h*v_at(b, i),v2_slice(temp_k, i),new_y,0,practical_dim);
    link_i->check_consistency(new_y,params,GlobalVars->global_params);

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(v2_slice(temp_k, 0),sum,practical_dim);
    sv_mlt_u(h * v_at(e, 0),sum,0,practical_dim);
    for(i=1;i<s;i++)	daxpy_u(h*v_at(e, i),v2_slice(temp_k, i),sum,0,practical_dim);

    //Build SC_i
    //for(i=0;i<dim;i++)
    for(i=0;i<practical_dim;i++)
        temv_at(p, i) = max(fabs(new_v_at(y, i)),fabs(y_v_at(0, i))) * error->reltov_at(l, i) + error->abstov_at(l, i);

    err_1 = norm_inf_u(sum,temp,0,practical_dim);
    double value_1 = pow(1.0/err_1,1.0/meth->e_order);

    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(v2_slice(temp_k, 0),sum,practical_dim);
    sv_mlt_u(h * v_at(d, 0),sum,0,practical_dim);
    for(i=1;i<s;i++)	daxpy_u(h*v_at(d, i),v2_slice(temp_k, i),sum,0,practical_dim);

    //for(i=0;i<dim;i++)
    for(i=0;i<practical_dim;i++)
        temv_at(p, i) = max(fabs(new_v_at(y, i)),fabs(y_v_at(0, i))) * error->reltol_densv_at(e, i) + error->abstol_densv_at(e, i);

    err_d = norm_inf_u(sum,temp,0,practical_dim);
    double value_d = pow(1.0/err_d,1.0/meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax,max(error->facmin,error->fac * value_1));
    double step_d = h*min(error->facmax,max(error->facmin,error->fac * value_d));
    link_i->h = min(step_1,step_d);

//printf("ID = %u t_0 = %e 1 = %e d = %e\n",link_i->ID,t,step_1,step_d);
//Print_Vector(new_y);
//getchar();

    if(err_1 < 1.0 && err_d < 1.0)
    {
        //Check if a discontinuity has been stepped on
        if(link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t,link_i,GlobalVars,workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(temp_k,new_node->k,s,dense_indices,num_dense);

        //Check if new data should be written to disk
        if(print_flag)
        {
            while( t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t)/link_i->next_save < 1e-12) )
            {
                if(link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n",my_rank,link_i->ID,link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t)/h;
                link_i->method->dense_b(current_theta,link_i->method->b_theta);
                //for(m=0;m<num_dense;m++)
                for(m=0;m<problem_dim;m++)	//!!!! This is wasteful if not all states are printed. Also forbids printing variational eq states. !!!!
                {
                    idx = dense_indices[m];
                    sum.ve[idx] = node->y_approx.ve[idx];
                    for(l=0;l<link_i->method->s;l++)
                        sum.ve[idx] += h*link_i->method->b_theta.ve[l] * node->next->k[l].ve[m];
                }
                link_i->check_consistency(sum,params,GlobalVars->global_params);

                //Write to a file
                WriteStep(link_i->next_save,sum,GlobalVars,params,link_i->state,outputfile,link_i->output_user,&(link_i->pos_offset));

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if( link_i->peak_flag && (new_v_at(y, 0) > link_i->peak_valuv_at(e, 0)) )
        {
            v_copy_n(new_y,link_i->peak_value,practical_dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for(j=0;j<GlobalVars->num_forcings;j++)
        {
            if( forcings[j].active && link_i->forcing_buff[j] && (fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8) )
            {
                //Propagate the discontinuity to downstream links
                if(!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for(i=0;i<GlobalVars->max_localorder && next != NULL;i++)
                    {
                        if(assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j],next->discont_start,next->discont_end,&(next->discont_count),GlobalVars->discont_size,next->discont,next->ID);
                        }
                        else if(next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j],i,&(prev->discont_send_count),GlobalVars->discont_size,prev->discont_send,prev->discont_order_send,prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                for(l=link_i->forcing_indices[j]+1;l<link_i->forcing_buff[j]->n_times;l++)
                    if( fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->rainfall[l][0]) < 1e-8 )	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->rainfall[l][1];
                link_i->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                for(i=l+1;i<link_i->forcing_buff[j]->n_times;i++)
                {
                    if( fabs(link_i->forcing_buff[j]->rainfall[i][1] - forcing_buffer) > 1e-8 )
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->rainfall[i][0];
                        break;
                    }
                }
                if(i == link_i->forcing_buff[j]->n_times)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->rainfall[i-1][0];
            }
        }

        //Select new step size, if forcings changed
        if(propagated)	link_i->h = InitialStepSize(link_i->last_t,link_i,GlobalVars,workspace);

        //Free up parents' old data
        for(i=0;i<link_i->num_parents;i++)
        {
            currentp = link_i->parents[i];
            while(currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(currentp->list);
                currentp->current_iterations--;
                currentp->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}

#endif // defined(ASYNCH_HAVE_ASSIM_SOLVER)
