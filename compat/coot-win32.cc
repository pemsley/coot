// vim: sw=4 ts=4 sts=4 expandtab :

#include "coot-sysdep.h"

#include <windows.h>

#include <algorithm>
#include <cwctype>
#include <array>
#include <string>

#define W_PATH_BUF_SIZE 32768UL
#define W_PATH_PREFIX L"\\\\?\\"

static
std::wstring absolute_path(const std::wstring &path) {
    std::array<wchar_t, W_PATH_BUF_SIZE> buf{};
    LPWSTR file_name = nullptr;

    if (GetFullPathNameW(path.c_str(), buf.size(), buf.data(), &file_name) == 0)
        return {};

    return std::wstring{buf.data()};
}

static
std::string get_error_string(DWORD error_code) {
    LPSTR buf = NULL;

    if (FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        error_code,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_NEUTRAL),
        (LPSTR) &buf,
        0,
        NULL
    ))
        return buf;
    else
        return "Unknown error";
}

static
std::wstring windowsize_path(const std::wstring &path) {
    std::wstring win_path{path};
    std::transform(win_path.begin(), win_path.end(), win_path.begin(),
                   [](wchar_t ch) {
                       if (ch == L'/')
                           return L'\\';
                        return ch;
                    });
    /* W_PATH_PREFIX increases the maximum path length limit to 32767 characters when used with Unicode
       variants of the respective WinAPI functions */
    if (win_path.length() > 4 && win_path.substr(0, 4) != W_PATH_PREFIX)
        return W_PATH_PREFIX + win_path;
    return win_path;
}

static
DWORD get_file_attributes(const std::string &file_path) {
    bool ok;
    std::wstring w_file_path = windowsize_path(coot::local_to_wide_string(file_path, &ok));
    if (!ok)
        return INVALID_FILE_ATTRIBUTES;

    return GetFileAttributesW(w_file_path.c_str());
}

static
bool is_regular_file_attrs(DWORD attrs) {
    const DWORD REG_FILE_MASK = ~(FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_HIDDEN | FILE_ATTRIBUTE_DEVICE | FILE_ATTRIBUTE_OFFLINE);
    if (attrs == FILE_ATTRIBUTE_NORMAL)
        return true;
    return attrs & REG_FILE_MASK;
}

static
bool paths_equal(const std::wstring &first, const std::wstring &second) {
    if (first.length() != second.length())
        return false;

    for (size_t idx = 0; idx < first.length(); idx++) {
        wchar_t a = std::towlower(first[idx]);
        wchar_t b = std::towlower(second[idx]);
        if (a != b)
            return false;
    }

    return true;
}

namespace coot {

std::string current_working_dir() {
    // ret DOES include the null-terminator
    DWORD ret = GetCurrentDirectoryW(0, NULL);
    if (ret < 1)
        return {};

    std::vector<wchar_t> w_bytes(ret);
    DWORD ret2 = GetCurrentDirectoryW(ret, w_bytes.data());
    // ret2 does NOT include the null-terminator when GetCurrentDirectoryW() is called like this
    if (ret2 < 1 || ret2 != ret - 1)
        return {};

    return wide_string_to_local(w_bytes.data());
}

std::vector<std::string> gather_files_by_patterns(const std::string &dir_path, const std::vector<std::string> &patterns) {
    bool ok;

    std::wstring w_dir_path = windowsize_path(local_to_wide_string(dir_path, &ok));
    if (!ok)
        return {};

    std::vector<std::wstring> w_patterns;
    for (const auto &p : patterns) {
        auto w_p = local_to_wide_string(p, &ok);
        if (!ok)
            return {};
        w_patterns.push_back(std::move(w_p));
    }

    std::vector<std::string> found;
    WIN32_FIND_DATAW find_data;
    for (const auto &p : w_patterns) {
        std::wstring glob = w_dir_path + L"\\" + p;

        HANDLE h = FindFirstFileW(glob.c_str(), &find_data);
        if (h == INVALID_HANDLE_VALUE)
            continue;

        do {
            if (!is_regular_file_attrs(find_data.dwFileAttributes))
                continue;

            std::string file_path = dir_path + "\\" + wide_string_to_local(find_data.cFileName, &ok);
            if (ok)
                found.push_back(std::move(file_path));
        } while (FindNextFileW(h, &find_data) != 0);

        FindClose(h);
    }

    return found;
}

std::string get_fixed_font() {
    return "monospace";
}

bool is_dir(const std::string &file_path) {
    DWORD attrs = get_file_attributes(file_path);
    if (attrs == INVALID_FILE_ATTRIBUTES)
        return false;

    return attrs & FILE_ATTRIBUTE_DIRECTORY;
}

bool is_link(const std::string &file_path) {
    std::wstring w_file_path = windowsize_path(local_to_wide_string(file_path));
    w_file_path = absolute_path(w_file_path);

    HANDLE h = CreateFileW(
        w_file_path.c_str(),
        GENERIC_READ,
        FILE_SHARE_WRITE,
        NULL,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL
    );
    if (h == INVALID_HANDLE_VALUE)
        return false;

    std::array<wchar_t, W_PATH_BUF_SIZE> buf{};
    DWORD ret = GetFinalPathNameByHandleW(h, buf.data(), buf.size(), FILE_NAME_NORMALIZED);
    CloseHandle(h);

    if (ret == 0)
        return false;

    std::wstring w_final_path{buf.data()};
    return !paths_equal(w_file_path, w_final_path);
}

bool is_regular_file(const std::string & file_path) {
    DWORD attrs = get_file_attributes(file_path);
    if (attrs == INVALID_FILE_ATTRIBUTES)
        return false;

    return is_regular_file_attrs(attrs);
}

std::wstring local_to_wide_string(const std::string &str, bool *ok) {
    int ret = MultiByteToWideChar(
        CP_ACP,
        MB_ERR_INVALID_CHARS,
        str.c_str(),
        -1,
        nullptr,
        0
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    std::vector<wchar_t> w_bytes(ret);
    ret = MultiByteToWideChar(
        CP_ACP,
        MB_ERR_INVALID_CHARS,
        str.c_str(),
        -1,
        w_bytes.data(),
        ret
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    if (ok) *ok = true;
    return std::wstring(w_bytes.data());
}

bool rename(const char *old_file_path, const char *new_file_path, std::string &error_message) {
    DWORD attrs = get_file_attributes(old_file_path);
    if (attrs == INVALID_FILE_ATTRIBUTES) {
        error_message = "File does not exist";
        return false; // File to rename does not exist
    }

    bool ok;
    std::wstring w_old_file_path = local_to_wide_string(old_file_path, &ok);
    if (!ok) {
        error_message = "Failed to convert old_file_path to Unicode";
        return false;
    }
    std::wstring w_new_file_path = local_to_wide_string(new_file_path, &ok);
    if (!ok) {
        error_message = "Failed to convert new_file_path to Unicode";
        return false;
    }

    w_old_file_path = windowsize_path(w_old_file_path);
    w_new_file_path = windowsize_path(w_new_file_path);

    BOOL ret = MoveFileExW(
        w_old_file_path.c_str(),
        w_new_file_path.c_str(),
        MOVEFILE_COPY_ALLOWED | MOVEFILE_REPLACE_EXISTING
    );
    if (ret)
        return true;

    error_message = get_error_string(GetLastError());
    return false;
}

void set_os_error_mode() {
    SetErrorMode(SetErrorMode(SEM_NOGPFAULTERRORBOX) | SEM_NOGPFAULTERRORBOX);
}

std::string wide_string_to_local(const std::wstring &w_str, bool *ok) {
    int ret = WideCharToMultiByte(
        CP_ACP,
        WC_NO_BEST_FIT_CHARS,
        w_str.c_str(),
        -1,
        nullptr,
        0,
        nullptr,
        nullptr
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    std::vector<char> bytes(ret);
    ret = WideCharToMultiByte(
        CP_ACP,
        WC_NO_BEST_FIT_CHARS,
        w_str.c_str(),
        -1,
        bytes.data(),
        ret,
        nullptr,
        nullptr
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    if (ok) *ok = true;
    return std::string(bytes.data());
}

void sleep(unsigned int secs) {
    Sleep(1000 * secs);
}

void usleep(unsigned int usecs) {
    Sleep(usecs / 1000);
}

} // namespace coot
