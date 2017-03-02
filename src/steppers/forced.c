#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>

#include <minmax.h>
#include <system.h>
#include <blas.h>
#include <io.h>
#include <rksteppers.h>

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
        if (forcings[i].active && link_i->forcing_data[i])
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
            if (forcings[j].active && link_i->forcing_data[j] && (fabs(new_node->t - link_i->forcing_change_times[j]) < 1e-8))
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
                //for(l=1;l<link_i->forcing_data[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_data[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_data[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_data[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_data[j]->nrows; i++)
                {
                    if (fabs(link_i->forcing_data[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_data[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_data[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_data[j]->data[i - 1][0];
            }
        }
    }
    link_i->differential(t + h, y_0, v2_init(0, 0), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, new_y);
    if (link_i->check_state)
        new_node->state = link_i->check_state(new_y, GlobalVars->global_params, link_i->params, link_i->qvs, link_i->is_dam);

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
                WriteStep(global->outputs, outputfile, link_i->ID, link_i->next_save, new_y, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
            else
                WriteStep(global->outputs, outputfile, link_i->ID, link_i->next_save, y_0, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
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
            if (forcings[j].active && link_i->forcing_data[j] && (fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8))
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
                //for(l=1;l<link_i->forcing_data[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_data[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_data[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_data[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_data[j]->nrows; i++)
                {
                    if (fabs(link_i->forcing_data[j]->data[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_data[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_data[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_data[j]->data[i - 1][0];
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