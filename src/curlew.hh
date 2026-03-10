/*! \file
  \brief Coot Scripting Interface - Curlew
*/
#include <string>
#include <vector>
#include <utility>

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

/*! \brief get the list of available curlew extensions
    @return a vector of pairs of [name, file-name] */
std::vector<std::pair<std::string, std::string> > curlew_get_extension_list();

/*! \brief download and install a curlew extension
    @param extension_file_name the file-name of the extension (e.g. "coot_goodsell_menu.py")
    @return 1 on success, 0 on failure */
int curlew_download_and_install_extension(const std::string &extension_file_name);

/* \} */
