#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif


#include "consistency.h"


void CheckConsistency_Nonzero_1States(VEC y,VEC params,VEC global_params)
{
	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
}

void CheckConsistency_Nonzero_2States(VEC y,VEC params,VEC global_params)
{
	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	if(y.storage[1] < 0.0)	y.storage[1] = 0.0;
}

void CheckConsistency_Nonzero_3States(VEC y,VEC params,VEC global_params)
{
	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	if(y.storage[1] < 0.0)	y.storage[1] = 0.0;
	if(y.storage[2] < 0.0)	y.storage[2] = 0.0;
}

void CheckConsistency_Nonzero_4States(VEC y,VEC params,VEC global_params)
{
	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	if(y.storage[1] < 0.0)	y.storage[1] = 0.0;
	if(y.storage[2] < 0.0)	y.storage[2] = 0.0;
	if(y.storage[3] < 0.0)	y.storage[3] = 0.0;
}

void CheckConsistency_Model5(VEC y,VEC params,VEC global_params)
{
	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	if(y.storage[1] < 0.0)	y.storage[1] = 0.0;
	if(y.storage[2] < 0.0)	y.storage[2] = 0.0;
	if(y.storage[3] < 0.0)	y.storage[3] = 0.0;
	if(y.storage[3] > 1.0 - y.storage[2])	y.storage[3] = 1.0 - y.storage[2];
}

void CheckConsistency_Model30(VEC y,VEC params,VEC global_params)
{
	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	if(y.storage[1] < 0.0)	y.storage[1] = 0.0;
	if(y.storage[2] < 0.0)	y.storage[2] = 0.0;
	if(y.storage[2] > params.storage[4])	y.storage[2] = v_at(params, 4);
	if(y.storage[3] < 0.0)	y.storage[3] = 0.0;
	else if(y.storage[3] > 1.0)	y.storage[3] = 1.0;
}

void CheckConsistency_Nonzero_AllStates_q(VEC y,VEC params,VEC global_params)
{
	unsigned int i;

	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	for(i=1;i<y.dim;i++)
		if(y.storage[i] < 0.0)	y.storage[i] = 0.0;
}

void CheckConsistency_Nonzero_AllStates_qs(VEC y,VEC params,VEC global_params)
{
	unsigned int i;

	if(y.storage[0] < 1e-14)	y.storage[0] = 1e-14;
	if(y.storage[1] < 1e-14)	y.storage[1] = 1e-14;
	for(i=2;i<y.dim;i++)
		if(y.storage[i] < 0.0)	y.storage[i] = 0.0;
}


