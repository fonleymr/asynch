#if !defined(CONSISTENCY_H)
#define CONSISTENCY_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vector.h"

//Consistency
void CheckConsistency_Nonzero_1States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_2States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_3States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_4States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Model5(VEC y,VEC params,VEC global_params);
void CheckConsistency_Model30(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_AllStates_q(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_AllStates_qs(VEC y,VEC params,VEC global_params);

#endif //CONSISTENCY_H

