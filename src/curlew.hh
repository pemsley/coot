/*! \file
  \brief Coot Scripting Interface - Curlew
*/
#include <string>

/*  ----------------------------------------------------------------------- */
/*                  curlew                                                  */
/*  ----------------------------------------------------------------------- */
/*! \name curlew */
/* \{ */
/*! \brief activate the curlew dialog */
// void curlew();

/*! \brief register an extension */
void register_extension(const std::string &name, const std::string &version);

std::string version_for_extension(const std::string &name);

void remove_file_curlew_menu_item_maybe();

/* \} */
