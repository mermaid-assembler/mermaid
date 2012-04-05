#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <csignal>

#include "utils.h"

void _panic(const char *file, int line, const char *fmt,...)
{
	va_list ap;

	va_start(ap, fmt);

	fprintf(stderr, "panic at %s:%d: ", file, line);
	vfprintf(stderr, fmt, ap);
	fprintf(stderr, "\n");

        raise(SIGABRT);

        /* To shut gcc up. */
        exit(1);
}
