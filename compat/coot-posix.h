// vim: sw=4 ts=4 sts=4 expandtab :

#ifndef _COMPAT_COOT_POSIX_H
#define _COMPAT_COOT_POSIX_H

#include <unistd.h>

namespace coot {
    void sleep(unsigned int secs);
    void usleep(useconds_t usecs);
}

#endif // _COMPAT_COOT_POSIX_H
