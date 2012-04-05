#ifndef _UTILS_H_
#define _UTILS_H_

void _panic(const char* file, int line, const char* fmt, ...) __attribute__((noreturn));

#define panic(...) _panic(__FILE__, __LINE__, __VA_ARGS__)

#endif /* _UTILS_H_ */
