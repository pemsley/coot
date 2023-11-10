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