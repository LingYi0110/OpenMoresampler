#ifndef MORESAMPLER_UTIL_H
#define MORESAMPLER_UTIL_H

#include <stdio.h>
#include <wchar.h>
#include <windows.h>

wchar_t *get_argv_w(int argc, char **argv, int index);
wchar_t *path_get_filename(wchar_t *path);
wchar_t *wstr_concat(const wchar_t *a, const wchar_t *b);
wchar_t *mbs_to_wchar(const char *mbs);
char *wchar_to_mbs(const wchar_t *wstr);
WCHAR *fgets_utf8_to_wchar(FILE *stream);
wchar_t *get_cwd_w(void);
BOOL path_is_directory(const wchar_t *path);

wchar_t *wstr_dup(const wchar_t *s);
wchar_t *wstr_ndup(const wchar_t *s, size_t n);
wchar_t *utf8_to_wchar(const char *s);
char *wchar_to_utf8(const wchar_t *s);
wchar_t *path_dir_prefix(const wchar_t *path);
wchar_t *path_basename(const wchar_t *path);
int count_lines(const wchar_t *s);
wchar_t *read_file_acp(const wchar_t *path);
wchar_t **list_files(const wchar_t *dir, const wchar_t *suffix, size_t *count);

#endif /* MORESAMPLER_UTIL_H */
