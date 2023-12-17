/* lhasa/glog_replacement.hpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
#ifndef LHASA_GLOG_REPLACEMENT_HPP
#define LHASA_GLOG_REPLACEMENT_HPP

namespace lhasa::impl {
    enum class LogLevel {
        Debug,
        Info,
        Message,
        Warning,
        Critical,
        Error
        // This is all we need for now.
    };
    void glog_replacement(LogLevel level, const char* format, ...) __attribute__((__format__ (__printf__, 2, 3)));
};

#define g_error(...)    lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Error,    __VA_ARGS__)          
#define g_message(...)  lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Message,  __VA_ARGS__)
#define g_critical(...) lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Critical, __VA_ARGS__)
#define g_warning(...)  lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Warning,  __VA_ARGS__)
// Just ignore the "once"
#define g_warning_once(...)  lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Warning,  __VA_ARGS__)
#define g_info(...)     lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Info,     __VA_ARGS__)
#define g_debug(...)    lhasa::impl::glog_replacement(lhasa::impl::LogLevel::Debug,    __VA_ARGS__)

#endif // LHASA_GLOG_REPLACEMENT_HPP