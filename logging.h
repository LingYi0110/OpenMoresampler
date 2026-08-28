#ifndef MORESAMPLER_LOGGING_H
#define MORESAMPLER_LOGGING_H

#include <stdio.h>

extern FILE *log_fp;
extern wchar_t *log_filename;

int log_tee_printf(const char *format, ...);
int log_file_printf(const char *format, ...);
void log_printf_w(const wchar_t *str);
void log_close(void);

#endif // MORESAMPLER_LOGGING_H
