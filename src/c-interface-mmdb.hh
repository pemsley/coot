

#include "mmdb_manager.h"

#ifdef __cplusplus
#ifdef USE_GUILE
#include <libguile.h>	
#endif // USE_GUILE
#endif 

// return 0 on failure
CMMDBManager *
mmdb_manager_from_scheme_expression(SCM molecule_expression);
SCM display_scm(SCM o);

