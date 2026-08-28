#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <wchar.h>

#include "logging.h"


FILE *log_fp = NULL;
wchar_t *log_filename = NULL;

// Print to both the log file (if open) and stdout. Returns the stdout count.
int log_tee_printf(const char *format, ...) {
    va_list args;
    va_start(args, format);

    if (log_fp) {
        vfprintf(log_fp, format, args);
        fflush(log_fp);
    }
    return vprintf(format, args);
}

// Print to the log file only; nothing goes to stdout.
int log_file_printf(const char *format, ...) {
    va_list args;
    va_start(args, format);

    if (log_fp) {
        vfprintf(log_fp, format, args);
        return fflush(log_fp);
    }
    return 0;
}

// Write a wide string verbatim to the log file only (raw bytes, no formatting).
void log_printf_w(const wchar_t *str) {
    if (log_fp) {
        fwrite(str, sizeof(wchar_t), wcslen(str), log_fp);
        fflush(log_fp);
    }
}

void log_close(void) {
    if (log_filename) {
        free(log_filename);
    }
    if (log_fp) {
        fclose(log_fp);
    }
}
