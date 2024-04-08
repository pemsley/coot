/*
 * src/curl-utils.hh
 *
 * Copyright 2023 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef CURL_UTILS_HH
#define CURL_UTILS_HH

#ifdef USE_LIBCURL

#include <Python.h>

#include <string>
#include <optional>
#include <curl/curl.h>
#include "gtk-utils.hh"

// return 0 on success.
int coot_get_url(const std::string &url, const std::string &file_name);

#ifdef SWIG
// don't look
#else
int coot_get_url_with_notifier(const std::string &url, const std::string &file_name, std::optional<ProgressNotifier> notifier = std::nullopt);
int coot_get_url_and_activate_curl_hook(const std::string &url, const std::string &file_name, short int do_hook_flag, std::optional<ProgressNotifier> notifier = std::nullopt);
#endif

#ifdef USE_GUILE
// this handles URLs that are strings, not binaries.
SCM coot_get_url_as_string(const char *url);

SCM curl_progress_info(const char *file_name);
// for the callback of the update binary progress bar.  How much done
// is the file that I am downloading?
SCM curl_progress_info(const char *file_name);
#endif /* USE_GUILE */

#ifdef USE_PYTHON
// this handles URLs that are strings, not binaries.
PyObject *coot_get_url_as_string_py(const char *url);
// for the callback of the update binary progress bar.  How much done
// is the file that I am downloading? Not absolutely required for python
PyObject *curl_progress_info_py(const char *file_name);
#endif /* USE_PYTHON */
// internal use
size_t write_coot_curl_data(void *buffer, size_t size, size_t nmemb, void *userp);
// internal use
size_t write_coot_curl_data_to_file(void *buffer, size_t size, size_t nmemb, void *userp);

void *wrapped_curl_easy_perform(void *data);
#endif

#endif // CURL_UTILS_HH
