// vim: sw=4 ts=4 sts=4 expandtab :

#include "coot-sysdep.h"

#include <fcntl.h>
#include <glob.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

#include <array>
#include <cstring>

#define COOT_PATH_MAX 4096UL

namespace coot {

std::string current_working_dir() {
    std::array<char, COOT_PATH_MAX> bytes{};
    if (getcwd(bytes.data(), bytes.size()))
        return std::string{bytes.data()};
    return {};
}

std::vector<std::string> gather_files_by_patterns(const std::string &dir_path, const std::vector<std::string> &patterns) {
    if (patterns.empty())
        return {};

    glob_t gb;
    int gb_flags = 0;
    for (const auto &p : patterns) {
        std::string path = dir_path + "/" + p;
        glob(path.c_str(), gb_flags, nullptr, &gb);
        gb_flags = GLOB_APPEND;
    }

    std::vector<std::string> found;
    size_t count = gb.gl_pathc;
    for (char **p = gb.gl_pathv; count ; p++, count--)
        found.push_back(*p);

    globfree(&gb);

    return found;
}

std::string get_fixed_font() {
    return "Sans 9";
}

bool is_dir(const std::string &file_path) {
    struct stat buf;
    if (stat(file_path.c_str(), &buf) == -1)
        return false;

    return S_ISDIR(buf.st_mode);
}

bool is_link(const std::string &file_path) {
    struct stat buf;
    if (stat(file_path.c_str(), &buf) == -1)
        return false;

    return S_ISLNK(buf.st_mode);
}

bool is_regular_file(const std::string &file_path) {
    struct stat buf;
    stat(file_path.c_str(), &buf);

    return S_ISREG(buf.st_mode);
}

bool rename(const char *old_file_path, const char *new_file_path, std::string &error_message) {
    int ret = ::rename(old_file_path, new_file_path);
    if (ret == 0)
        return true;

    error_message = std::string{std::strerror(ret)};

    return false;
}

void set_os_error_mode() {
    // NOOP
}

void sleep(unsigned int secs) {
    ::sleep(secs);
}

void usleep(useconds_t usecs) {
    ::usleep(usecs);
}

} // namespace coot
