#if !defined(ASYNCH_MINMAX_H)
#define ASYNCH_MINMAX_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#define min(val1,val2)	(((val1) > (val2)) ? (val2) : (val1))
#define max(val1,val2)	(((val1) < (val2)) ? (val2) : (val1))
#endif

#endif //ASYNCH_MINMAX_H