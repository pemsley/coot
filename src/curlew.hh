
/*! \file
  \brief Coot Scripting Interface - Curlew
*/

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

/* \} */
