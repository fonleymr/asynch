#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "blas.h"
#include "system.h"

//Frees link.
//Link* link: link to be freed.
//unsigned int list_length: the length of the list stored with link.
//int rkd_flag: should be 1 if an .rkd file was used, 0 if not.
void Destroy_Link(Link* link, unsigned int list_length, int rkd_flag, Forcing* forcings, GlobalVars* GlobalVars)
{
    unsigned int i;
    assert(link != NULL);

    v_free(&link->params);

    if (link->method != NULL)
    {
        v_free(&link->forcing_values);
        free(link->forcing_indices);
        free(link->forcing_change_times);
        if (rkd_flag)
            Destroy_ErrorData(&link->my->error_data);
        Destroy_List(&link->my->list, list_length);
        //v_free(link->params);
        v_free(&link->peak_value);
        if (link->discont != NULL)
            free(link->discont);
        if (link->discont_send != NULL)
        {
            free(link->discont_send);
            free(link->discont_order_send);
        }
        for (i = 0; i < GlobalVars->num_forcings; i++)
        {
            if (link->forcing_buff && forcings[i].flag != 4 && forcings[i].flag != 7 && link->forcing_buff[i] != NULL)
                ForcingData_Free(&(link->forcing_buff[i]));
        }
        free(link->forcing_buff);
        if (link->qvs != NULL)
        {
            /*
                        for(i=0;i<link->qvs->n_values;i++)
                            free(link->qvs->points[i]);
            */
            free(link->qvs->points);
            free(link->qvs->points_array);
            free(link->qvs);
        }
#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
        m_free(&link->JMatrix);
        m_free(&link->CoefMat);
        v_free(&link->sol_diff);
        if (link->Z_i != NULL)
        {
            for (i = 0; i < GlobalVars->max_s; i++)
                v_free(&link->Z_i[i]);
            free(link->Z_i);
        }
#endif
        /*
                if(GlobalVars->template_flag && link->equations != NULL)
                {
                    //mupRelease(link->equations->parser);
                    printf("Warning: Not freeing parser in system.c, Destroy_Link\n");
                    v_free(&link->equations->variable_values);
                    free(link->equations);
                }
        */
    }

    if (link->dense_indices)
        free(link->dense_indices);

    free(link->parents);
}

//Frees rain
//ForcingData** forcing_buff: forcing data to be freed
void ForcingData_Free(ForcingData** forcing_buff)
{
    unsigned int i;
    if (forcing_buff && *forcing_buff)
    {
        for (i = 0; i < (*forcing_buff)->nrows; i++)
            free((*forcing_buff)->data[i]);
        free((*forcing_buff)->data);
        free(*forcing_buff);
    }
}


//Creates a list to hold the data for an ODE.
//VEC* y0: the initial data.
//double t0: the initial time.
//int dim: the dimension of the ODE.
//int num_stages: the number of stages in the RKMethod.
//unsigned int list_length: the maximum number of steps to store in the list.
void Init_List(RKSolutionList* list, VEC y0, double t0, int dim, unsigned int num_dense, unsigned short int num_stages, unsigned int list_length)
{
    assert(list_length > 0);
    if (list_length < 2)
        printf("Warning in Create_List: list_length is %u.\n", list_length);

    memset(list, 0, sizeof(RKSolutionList));
    list->nodes = (RKSolutionNode*)calloc(list_length, sizeof(RKSolutionNode));
    //for(i=0;i<list_length;i++)
    //  list->nodes[i] = (RKSolutionNode*) malloc(sizeof(RKSolutionNode));

    //Set the next and prev ptrs for each node
    list->nodes[0].next = &list->nodes[1];
    list->nodes[0].prev = &list->nodes[list_length - 1];

    for (unsigned int i = 1; i < list_length - 1; i++)
    {
        list->nodes[i].next = &list->nodes[i + 1];
        list->nodes[i].prev = &list->nodes[i - 1];
    }

    list->nodes[list_length - 1].next = &list->nodes[0];
    list->nodes[list_length - 1].prev = &list->nodes[list_length - 2];

    //Allocate space for all the vectors
    list->y_storage = v2_init(list_length, dim);
    list->k_storage = v3_init(list_length, num_stages, num_dense);

    for (unsigned int i = 0; i < list_length; i++)
    {
        list->nodes[i].y_approx = v2_slice(list->y_storage, i);
        list->nodes[i].k = v3_slice(list->k_storage, i);



        //list->nodes[i].k = list->k_vectors + i * num_stages;
        //for (j = 0; j < num_stages; j++)
        //{
        //    list->nodes[i].k[j].storage = list->k_storage + i * j * num_dense;
        //    list->nodes[i].k[j].dim = num_dense;
        //}
    }


    //for(i=0;i<list_length;i++)
    //{
    //	list->nodes[i].y_approx = v_init(dim);
    //	list->nodes[i].k = (VEC*) malloc(s * sizeof(VEC));
    //	for(j=0;j<s;j++)
 //           list->nodes[i].k[j] = v_init(num_dense);
    //}

    //Set remaining fields
    list->head = &list->nodes[0];
    list->tail = &list->nodes[0];
    list->num_stages = num_stages;

    //Store the initial step
    list->head->t = t0;
    v_copy(y0, list->head->y_approx);
}

//Frees the data list.
void Destroy_List(RKSolutionList* list, unsigned int list_length)
{
    free(list->nodes);
    v2_free(&list->y_storage);
    v3_free(&list->k_storage);
}

//Removes the first node in list.
void Remove_Head_Node(RKSolutionList* list)
{
    if (list->head == list->tail)	//Keep the tail at head if the list has no data
        list->tail = list->tail->next;

    list->head = list->head->next;
}

//Adds a new step to list, after tail. The new node is returned.
RKSolutionNode* New_Step(RKSolutionList* list)
{
    list->tail = list->tail->next;
    return list->tail;
}

//This method undoes the last step taken.
//Used for rejecting steps.
void Undo_Step(RKSolutionList* list)
{
    list->tail = list->tail->prev;
}

//Frees an RKMethod
void Destroy_RKMethod(RKMethod* method)
{
    assert(method != NULL);
    v2_free(&method->A);
    v_free(&method->b);
    v_free(&method->b_theta);
    v_free(&method->b_theta_deriv);
    v_free(&method->c);
    v_free(&method->e);
    v_free(&method->d);
    v_free(&method->w);
}

//Frees an ErrorData
void Destroy_ErrorData(ErrorData* error)
{
    assert(error != NULL);
    v_free(&error->abstol);
    v_free(&error->reltol);
    v_free(&error->abstol_dense);
    v_free(&error->reltol_dense);
}

//Allocates workspace for RK solvers
void Create_Workspace(Workspace *workspace, unsigned int dim, unsigned short num_stages, unsigned short int max_parents)
{
    memset(workspace, 0, sizeof(Workspace));

    workspace->temp = v_init(dim);
    workspace->sum = v_init(dim);
    workspace->temp2 = v_init(dim);
    workspace->temp3 = v_init(dim);

    //workspace->temp_parent_approx = (VEC**)malloc(num_stages * sizeof(VEC*));
    //for (unsigned int i = 0; i < num_stages; i++)
    //{
    //    workspace->temp_parent_approx[i] = malloc(max_parents * sizeof(VEC));
    //    for (unsigned int j = 0; j < max_parents; j++)
    //        workspace->temp_parent_approx[i][j] = v_init(dim);
    //}

    workspace->temp_parent_approx = v3_init(num_stages, max_parents, dim);

    //workspace->temp_k = (VEC*)malloc(num_stages * sizeof(VEC));
    //for (unsigned int i = 0; i < num_stages; i++)
    //    workspace->temp_k[i] = v_init(dim);

    workspace->temp_k = v2_init(num_stages, dim);

    for (unsigned int i = 0; i < num_stages; i++)
        workspace->temp_k_slices[i] = v2_slice(workspace->temp_k, i);

#if defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
    workspace->ipiv = (int*)malloc(s*dim * sizeof(int));
    workspace->rhs = v_init(s*dim);
    //workspace->CoefMat = m_get(s*dim,s*dim);
    workspace->JMatrix = m_get(dim, dim);
    workspace->Z_i = (VEC*)malloc(s * sizeof(VEC));
    for (i = 0; i < s; i++)	workspace->Z_i[i] = v_init(dim);
    workspace->err = v_init(dim);
#endif // defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
}

//Deallocates workspace for RK solvers
void Destroy_Workspace(Workspace* workspace, unsigned short int num_stages, unsigned short int max_parents)
{
    v_free(&workspace->temp);
    v_free(&workspace->sum);
    v_free(&workspace->temp2);
    v_free(&workspace->temp3);

    //for (unsigned int i = 0; i < num_stages; i++)
    //{
    //    for (unsigned int j = 0; j < max_parents; j++)
    //        v_free(&workspace->temp_parent_approx[i][j]);
    //    free(workspace->temp_parent_approx[i]);
    //}
    v3_free(&workspace->temp_parent_approx);

    //for (unsigned int i = 0; i < num_stages; i++)
    //    v_free(&workspace->temp_k[i]);
    v2_free(&workspace->temp_k);
    
#if defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
    free(workspace->ipiv);
    v_free(&workspace->rhs);
    //m_free(&workspace->CoefMat);
    m_free(&workspace->JMatrix);
    for (i = 0; i < s; i++)
        v_free(&workspace->Z_i[i]);
    free(workspace->Z_i);
    v_free(&workspace->err);
#endif // defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
}

//Deallocates UnivVars
void Destroy_UnivVars(GlobalVars* GlobalVars)
{
    unsigned int i;
    free(GlobalVars->peakflow_function_name);
    free(GlobalVars->output_types);
    free(GlobalVars->output_sizes);
    for (i = 0; i < GlobalVars->num_outputs; i++)
        free(GlobalVars->output_names[i]);

    free(GlobalVars->output_names);
    free(GlobalVars->outputs);
    if (GlobalVars->rsv_filename)
        free(GlobalVars->rsv_filename);
    if (GlobalVars->rvr_filename)
        free(GlobalVars->rvr_filename);
    if (GlobalVars->prm_filename)	free(GlobalVars->prm_filename);
    if (GlobalVars->init_filename)	free(GlobalVars->init_filename);
    if (GlobalVars->rain_filename)	free(GlobalVars->rain_filename);
    if (GlobalVars->dam_filename)	free(GlobalVars->dam_filename);
    if (GlobalVars->hydrosave_filename)	free(GlobalVars->hydrosave_filename);
    if (GlobalVars->peaksave_filename)	free(GlobalVars->peaksave_filename);
    if (GlobalVars->temp_filename)	free(GlobalVars->temp_filename);
    if (GlobalVars->peakfilename)	free(GlobalVars->peakfilename);
    if (GlobalVars->hydros_loc_filename)	free(GlobalVars->hydros_loc_filename);
    if (GlobalVars->peaks_loc_filename)	free(GlobalVars->peaks_loc_filename);
    if (GlobalVars->hydro_table)	free(GlobalVars->hydro_table);
    if (GlobalVars->peak_table)	free(GlobalVars->peak_table);
    if (GlobalVars->dump_table)	free(GlobalVars->dump_table);
    if (GlobalVars->dump_loc_filename)	free(GlobalVars->dump_loc_filename);
    v_free(&GlobalVars->global_params);
    free(GlobalVars->print_indices);
    free(GlobalVars);
}

//Inserts a time into the list of discontinuities
unsigned int Insert_Discontinuity(double time, unsigned int start, unsigned int end, unsigned int* count, unsigned int size, double* array, unsigned int id)
{
    if (*count == 0)
    {
        array[start] = time;
        (*count)++;
        return (end + 1) % size;
    }
    else if ((end + 1) % size == start)
    {
        //printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer. id = %u time = %f\n",id,time);
        return end;
    }

    unsigned int curr = end;
    unsigned int toofar = (start == 0) ? size - 1 : start - 1;

    while (curr != toofar && time < array[curr])
        curr = (curr == 0) ? size - 1 : curr - 1;
    if (curr != toofar && time == array[curr])	return end;

    curr = end;
    while (curr != toofar && time < array[curr])
    {
        array[(curr + 1) % size] = array[curr];
        curr = (curr == 0) ? size - 1 : curr - 1;
    }

    curr = (curr + 1) % size;
    end = (end + 1) % size;
    array[curr] = time;
    (*count)++;
    /*
        if(end == start)
        {
            printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer.\n");
        }
    */
    return end;
}


//Inserts a time into the list of discontinuities to be send to another processes
void Insert_SendDiscontinuity(double time, unsigned int order, unsigned int* count, unsigned int size, double* array, unsigned int* order_array, unsigned int id)
{
    if (*count == 0)
    {
        array[0] = time;
        order_array[0] = order;
        *count = 1;
        return;
    }
    else if (*count + 1 == size)
    {
        //printf("Warning: A discontinuity is being ignored in send. Increase the size of the discontinuity buffer. id = %i time = %f\n",id,time);
        return;
    }

    unsigned int curr = *count - 1;

    while (curr < size && time < array[curr])	curr--;
    if (curr < size && time == array[curr])
    {
        if (order_array[curr] < order)	order_array[curr] = order;
        return;
    }

    curr = *count - 1;
    while (curr < size && time < array[curr])
    {
        array[curr + 1] = array[curr];
        order_array[curr + 1] = order_array[curr];
        curr--;
    }

    curr++;
    array[curr] = time;
    order_array[curr] = order;
    (*count)++;

    /*
        if(*count == size)
        {
            printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer.\n");
        }
    */

}

