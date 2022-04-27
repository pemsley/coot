// vim: sw=4 ts=4 sts=4 expandtab :

#include "helpers.h"

#include <algorithm>
#include <cctype>

static
std::pair<std::string, std::string> split_uri(const std::string &uri) {
    static const std::string DELIM{"://"};
    static const auto DELIM_LEN = DELIM.length();

    auto pos = uri.find_first_of("://");
    if (pos == std::string::npos)
        return { "", "" };

    auto scheme = uri.substr(0, pos);
    auto tail = uri.substr(pos + DELIM_LEN);

    std::transform(
        scheme.begin(),
        scheme.end(),
        scheme.begin(),
        [](unsigned char ch) { return std::tolower(ch); }
    );

    pos = tail.find_first_of('?');
    if (pos != std::string::npos)
        tail = tail.substr(0, pos);

    pos = tail.find_first_of('#');
    if (pos != std::string::npos)
        tail = tail.substr(0, pos);

    return { scheme, tail };
}

namespace coot {

std::string uri_to_file_name(const std::string &uri) {
    auto splitted = split_uri(uri);
    const auto &scheme = splitted.first;
    const auto &path = splitted.second;

    if (scheme.empty())
        return ""; // Argument is not a valid URI
    if (path.empty())
        return "";

    if (scheme == "file" && path[0] == '/')
        return path.substr(1);
    return path;
}

} // namespace coot
