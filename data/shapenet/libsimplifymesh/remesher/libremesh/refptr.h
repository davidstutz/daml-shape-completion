#include "defines.h"

#if REMESHER_PARALLELIZATION
# include "refptr_mt.h"
#else
# include "refptr_st.h"
#endif
