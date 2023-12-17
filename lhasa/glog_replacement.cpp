/* lhasa/glog_replacement.cpp
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
#include "glog_replacement.hpp"
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <stdexcept>

void lhasa::impl::glog_replacement(lhasa::impl::LogLevel level, const char* format, ...) {
    va_list args;
    va_start(args, format);
    std::string m_format;
    switch(level) {
        case LogLevel::Debug:{
            m_format += "[Debug] - ";
            break;
        }
        case LogLevel::Info:{
            m_format += "[Info] - ";
            break;
        }
        case LogLevel::Message:{
            m_format += "[Message] - ";
            break;
        }
        case LogLevel::Warning:{
            m_format += "[Warning] - ";
            break;
        }
        case LogLevel::Critical:{
            m_format += "[Critical] - ";
            break;
        }
        default:
        case LogLevel::Error:{
            m_format += "[Error] - ";
            break;
        }
    }
    m_format.append(format);
    m_format.append("\n");
    vprintf(m_format.c_str(), args);
    va_end(args);

    if(level == LogLevel::Error){
        // todo: this might need improvement
        throw std::runtime_error("g_error() called. Fatal error.");
    }
}