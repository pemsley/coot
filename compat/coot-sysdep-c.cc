// vim: sw=4 ts=4 sts=4 expandtab

#include "coot-sysdep-c.h"
#include "coot-sysdep.h"

#include <cstring>
#include <cstdlib>

void coot_c_free_found(Coot_C_Found *found) {
    if (!found)
        return;

    for (size_t idx = 0; idx < found->n_found; idx++)
        std::free(found->found[idx]);

    std::free(found);
}

Coot_C_Found coot_c_gather_files_by_patterns(const char *dir_path, const char **patterns, size_t n_patterns) {
    std::vector<std::string> cxx_patterns;
    for (size_t idx = 0; idx < n_patterns; idx++)
        cxx_patterns.push_back(patterns[idx]);

    const auto cxx_found = coot::gather_files_by_patterns(dir_path, cxx_patterns);
    char **found = static_cast<char **>(std::calloc(cxx_found.size(), sizeof(char *)));
    for (size_t idx = 0; idx < cxx_found.size(); idx++) {
        auto str = cxx_found[idx].c_str();
        auto len = std::strlen(str);
        found[idx] = static_cast<char *>(std::malloc(len + 1));
        std::strcpy(found[idx], str);
    }

    Coot_C_Found ret { found, cxx_found.size() };
    return ret;
}
