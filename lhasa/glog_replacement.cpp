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