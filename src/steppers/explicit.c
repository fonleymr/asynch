#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <memory.h>

#include <minmax.h>
#include <system.h>
#include <blas.h>
#include <io.h>
#include <rksteppers.h>


extern AsynchModel *the_model;



void TopLayerHillslopeSIMD(
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned int num_parents,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    void *user,
    double *ans);


//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//Link* link_i: the link to apply a numerical method to.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolver(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int idx;
    //VEC** k;
    
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    Link* currentp;
    double t_needed, current_theta;

    //Some variables to make things easier to read
    double *y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    const double * const A = link_i->method->A;
    double *b = link_i->method->b;
    const double * const c = link_i->method->c;
    const double * const e = link_i->method->e;
    const double * const d = link_i->method->d;
    unsigned int num_stages = link_i->method->num_stages;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = globals->num_outputs;
    double *temp = workspace->temp;
    double *sum = workspace->sum;
    double **temp_k = workspace->temp_k_slices;

    //Get the approximate solutions from each parent
#pragma loop_count min(0), max(8)
    for (unsigned int i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;
#pragma loop_count min(5), max(7)
        for (unsigned int j = 0; j < num_stages; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + c[j] * h, currentp->last_t);

            //Find the corresponding theta value and approximate solution
            //while(t_needed > curr_node[i]->t && ( fabs(curr_node[i]->t) < 1e-12 || fabs(t_needed - curr_node[i]->t)/curr_node[i]->t > 1e-12) )
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)
                curr_node[i] = curr_node[i]->prev;

            double dt = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / dt;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            //[num_stages][max_parents][dim]
            double *parent_approx = workspace->stages_parents_approx + j * globals->max_parents * dim + i * dim;

            //for (unsigned int m = 0; m < currentp->num_dense; m++)
            //{
            //    idx = currentp->dense_indices[m];
            //    double approx = v_at(curr_node[i]->y_approx, idx);
            //    for (unsigned int l = 0; l < currentp->method->s; l++)
            //        approx += dt * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

            //    v_set(parent_approx, idx, approx);
            //}
#pragma loop_count max(1)
            for (unsigned int m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                double approx = curr_node[i]->y_approx[idx];
                //for (unsigned int l = 0; l < currentp->method->num_stages; l++)
#pragma loop_count min(5), max(7)
                for (int l = 0; l < currentp->method->num_stages; l++)
                    approx += dt * currentp->method->b_theta[l] *
                        curr_node[i]->next->k[l * currentp->num_dense + m];

                parent_approx[idx] = approx;
            }

            the_model->check_consistency(parent_approx, currentp->dim, currentp->params, the_model->num_params, globals->global_params, the_model->num_global_params, currentp->user);
        }
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    double *new_y = new_node->y_approx;
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
#pragma loop_count min(5), max(7)
    for (unsigned int i = 0; i < num_stages; i++)
    {
        //v_copy_n(y_0, sum, link_i->dim);
        memcpy(sum, y_0, link_i->dim * sizeof(double));
        for (unsigned int j = 0; j < i; j++)
        {
            //daxpy_u(h * v2_at(A, i, j), v2_slice(temp_k, j), sum, 0, link_i->dim);
            //daxpy_u(h * v2_at(A, i, j), temp_k[j], sum, 0, link_i->dim);

            double alpha = h * A[i * num_stages + j];
            daxpy(alpha, temp_k[j], sum, 0, link_i->dim);
        }

        //link_i->check_consistency(sum, params, globals->global_params);
        the_model->check_consistency(sum, link_i->dim, link_i->params, the_model->num_params, globals->global_params, the_model->num_global_params, link_i->user);

        //[num_stages][max_parents][dim]
        double *y_p = workspace->stages_parents_approx + i * globals->max_parents * link_i->dim;

        double dt = c[i] * h;

        the_model->differential(
            t + dt,
            sum, link_i->dim,
            y_p, link_i->num_parents,
            globals->global_params,
            link_i->params,
            link_i->my->forcing_values,            
            link_i->user,
            temp_k[i]);
    }

    //Build the solution
    dcopy(y_0, new_y, 0, link_i->dim);
    for (unsigned int i = 0; i < num_stages; i++)
        //daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 0, link_i->dim);
        daxpy(h * b[i], temp_k[i], new_y, 0, link_i->dim);

    the_model->check_consistency(new_y, link_i->dim, link_i->params, the_model->num_params, globals->global_params, the_model->num_global_params, link_i->user);
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
    dcopy(temp_k[0], sum, 0, link_i->dim);
    dscal(h * e[0], sum, 0, link_i->dim);
    for (unsigned int i = 1; i < num_stages; i++)
        //daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
        daxpy(h * e[i], temp_k[i], sum, 0, link_i->dim);

    //Build SC_i
    for (unsigned int i = 0; i < dim; i++)
        temp[i] = max(fabs(new_y[i]), fabs(y_0[i])) * error->reltol[i] + error->abstol[i];

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,0);
    err_1 = nrminf(sum, temp, 0, link_i->dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);

    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    dcopy(temp_k[0], sum, 0, link_i->dim);
    dscal(h * d[0], sum, 0, link_i->dim);

    for (unsigned int i = 1; i < num_stages; i++)
        //daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
        daxpy(h * d[i], temp_k[i], sum, 0, link_i->dim);

    for (unsigned int i = 0; i < dim; i++)
        temp[i] = max(fabs(new_y[i]), fabs(y_0[i])) * error->reltol_dense[i] + error->abstol_dense[i];

    //err_d = norm_inf(sum,temp,meth->d_order_ratio,0);
    err_d = nrminf(sum, temp, 0, link_i->dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h * min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h * min(error->facmax, max(error->facmin, error->fac * value_d));
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

    Jx_simple_river(child->list->tail->y_approx,globals->global_params,child->params,temp);
    sv_mlt((t+h - child->last_t)*sum_of_errors,temp,globals->diff_start);
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
            link_i->discont_start = (link_i->discont_start + 1) % globals->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, globals, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(workspace->temp_k, link_i->dim, new_node->k, num_stages, dense_indices, num_dense);

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
                    double approx = node->y_approx[idx];
                    for (unsigned int l = 0; l < link_i->method->num_stages; l++)
                        approx += h * link_i->method->b_theta[l] * node->next->k[l * num_dense + m];

                    sum[idx] = approx;
                }

                the_model->check_consistency(sum, link_i->dim, link_i->params, the_model->num_params, globals->global_params, the_model->num_global_params, link_i->user);

                WriteStep(globals->outputs, globals->num_outputs, outputfile, link_i->ID, link_i->next_save, sum, link_i->dim, &(link_i->pos_offset));

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a peak value
        if (link_i->peak_flag && (new_y[0] > link_i->peak_value[0]))
        {
            dcopy(new_y, link_i->peak_value, 0, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (unsigned int j = 0; j < globals->num_forcings; j++)
        {
            if (forcings[j].active && (link_i->my->forcing_data[j].num_points > 0) && (fabs(link_i->last_t - link_i->my->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (unsigned int i = 0; i < globals->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->my->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), globals->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->my->forcing_change_times[j], i, &(prev->discont_send_count), globals->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->my->forcing_data[j].n_times;l++)
                unsigned int l;
                for (l = link_i->my->forcing_indices[j] + 1; l < link_i->my->forcing_data[j].num_points; l++)
                    if (fabs(link_i->my->forcing_change_times[j] - link_i->my->forcing_data[j].data[l].time) < 1e-8)
                        break;
                link_i->my->forcing_indices[j] = l;

                double forcing_buffer = link_i->my->forcing_data[j].data[l].value;
                link_i->my->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                unsigned int i;
                for (i = l + 1; i < link_i->my->forcing_data[j].num_points; i++)
                {
                    //if(link_i->my->forcing_data[j].rainfall[i].value != forcing_buffer)
                    if (fabs(link_i->my->forcing_data[j].data[i].value - forcing_buffer) > 1e-8)
                    {
                        link_i->my->forcing_change_times[j] = link_i->my->forcing_data[j].data[i].time;
                        break;
                    }
                }
                if (i == link_i->my->forcing_data[j].num_points)
                    link_i->my->forcing_change_times[j] = link_i->my->forcing_data[j].data[i - 1].time;
            }
        }

        //Select new step size, if forcings changed
        if (propagated)
            link_i->h = InitialStepSize(link_i->last_t, link_i, globals, workspace);

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
