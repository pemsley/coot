// vim: sw=4 ts=4 sts=4 expandtab :

#ifndef _COMPAT_COOT_WIN32_H
#define _COMPAT_COOT_WIN32_H

#define COOT_ENABLE_WINAPI_SUSPENSION

#include <string>

namespace coot {
    std::wstring local_to_wide_string(const std::string &str, bool *ok = nullptr);
    std::string wide_string_to_local(const std::wstring &w_str, bool *ok = nullptr);
    void sleep(unsigned int secs);
    void usleep(unsigned int usecs);
} // namespace coot

#endif // _COMPAT_COOT_WIN32_H
