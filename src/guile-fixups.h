#ifdef USE_GUILE
#include <guile/gh.h>
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    
#define scm_to_int gh_scm2int
#define scm_to_locale_string SCM_STRING_CHARS
#define scm_from_locale_string scm_makfrom0str
#define scm_to_double  gh_scm2double
#define scm_is_true gh_scm2bool
#define scm_is_string gh_string_p
#define scm_car SCM_CAR
#endif // SCM version
#endif // USE_GUILE

