#include "utils.h"

#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <shlwapi.h>

wchar_t *get_argv_w(int argc, char **argv, int index)
{
    wchar_t *dst = malloc(0x800 * sizeof(wchar_t));
    int num_args = argc;
    LPWSTR *wargv = CommandLineToArgvW(GetCommandLineW(), &num_args);
    return wcscpy(dst, wargv[index]);
}

wchar_t *path_get_filename(wchar_t *path)
{
    signed int i = (signed int)wcslen(path) - 1;
    if (i < 0)
        return path;

    wchar_t c = path[i];
    if (c == L'/' || c == L'\\')
        return &path[i + 1];

    signed int j = i;
    while (--j != -1) {
        c = path[j];
        i = j;
        if (c == L'/' || c == L'\\')
            return &path[i + 1];
    }
    return path;
}

wchar_t *wstr_concat(const wchar_t *a, const wchar_t *b)
{
    size_t na = wcslen(a);
    size_t nb = wcslen(b);
    wchar_t *out = malloc(sizeof(wchar_t) * (na + nb + 1));
    wcscpy(out, a);
    wcscat(out, b);
    return out;
}

wchar_t *mbs_to_wchar(const char *mbs)
{
    if (!mbs)
        return NULL;
    int n = MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, mbs, -1, NULL, 0);
    wchar_t *out = calloc((size_t)n, sizeof(wchar_t));
    MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, mbs, -1, out, n);
    return out;
}

char *wchar_to_mbs(const wchar_t *wstr)
{
    if (!wstr)
        return NULL;
    int n = WideCharToMultiByte(CP_ACP, WC_NO_BEST_FIT_CHARS, wstr, -1,
                                NULL, 0, NULL, NULL);
    char *out = malloc((size_t)n);
    WideCharToMultiByte(CP_ACP, WC_NO_BEST_FIT_CHARS, wstr, -1,
                        out, n, NULL, NULL);
    return out;
}

WCHAR *fgets_utf8_to_wchar(FILE *stream)
{
    WCHAR *wbuf = calloc(0x2000, sizeof(WCHAR));
    char *mbuf = calloc(0x2000, 1);
    if (fgets(mbuf, 0x2000, stream) != NULL) {
        MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, mbuf, -1, wbuf, 0x2000);
    } else {
        free(wbuf);
        wbuf = NULL;
    }
    free(mbuf);
    return wbuf;
}

wchar_t *get_cwd_w(void)
{
    wchar_t *buf = calloc(0x800, sizeof(wchar_t));
    _wgetcwd(buf, 2048);
    return buf;
}

BOOL path_is_directory(const wchar_t *path)
{
    return PathIsDirectoryW(path);
}

wchar_t *wstr_dup(const wchar_t *s)
{
    size_t n;
    wchar_t *d;

    if (s == NULL)
        return NULL;
    n = wcslen(s);
    d = (wchar_t *)malloc((n + 1) * sizeof(wchar_t));
    wcscpy(d, s);
    d[n] = 0;
    return d;
}

wchar_t *wstr_ndup(const wchar_t *s, size_t n)
{
    wchar_t *d = (wchar_t *)malloc((n + 1) * sizeof(wchar_t));
    wcsncpy(d, s, n);
    d[n] = 0;
    return d;
}

wchar_t *utf8_to_wchar(const char *s)
{
    int n;
    wchar_t *d;

    if (s == NULL)
        return NULL;
    n = MultiByteToWideChar(65001, 0, s, -1, NULL, 0);
    d = (wchar_t *)calloc((size_t)n, sizeof(wchar_t));
    MultiByteToWideChar(65001, 0, s, -1, d, n);
    return d;
}

char *wchar_to_utf8(const wchar_t *s)
{
    int n;
    char *d;

    if (s == NULL)
        return NULL;
    n = WideCharToMultiByte(65001, 0, s, -1, NULL, 0, NULL, NULL);
    d = (char *)malloc((size_t)n);
    WideCharToMultiByte(65001, 0, s, -1, d, n, NULL, NULL);
    return d;
}

wchar_t *path_dir_prefix(const wchar_t *path)
{
    size_t n = wcslen(path);
    wchar_t *d = (wchar_t *)calloc(n + 2, sizeof(wchar_t));
    wchar_t last;

    wcscpy(d, path);
    if (n == 0)
        return d;
    last = d[n - 1];
    if (last != L'/' && last != L'\\')
        d[n] = L'/';
    return d;
}

wchar_t *path_basename(const wchar_t *path)
{
    signed int i = (signed int)wcslen(path) - 1;
    const wchar_t *src;
    size_t n;
    wchar_t *d;

    if (i >= 0) {
        wchar_t c = path[i];
        if (c != L'\\' && c != L'/') {
            while (--i != -1) {
                c = path[i];
                if (c == L'\\' || c == L'/')
                    break;
            }
        }
    }
    src = &path[i + 1];
    n = wcslen(src);
    d = (wchar_t *)malloc((n + 1) * sizeof(wchar_t));
    wcscpy(d, src);
    d[n] = 0;
    return d;
}

int count_lines(const wchar_t *s)
{
    signed int n = (signed int)wcslen(s);
    int lines = 0;
    signed int i = 0;
    wchar_t c;

    if (n <= 0)
        return 0;
    c = *s;
    while (1) {
        while (c != 13 && c != 10) {
            if (n <= ++i)
                return lines + 1;
            c = s[i];
        }
        if (n <= i)
            return lines + 1;
        while (n != ++i) {
            c = s[i];
            if (c != 10 && c != 13) {
                ++lines;
                break;
            }
        }
        if (i >= n)
            return lines + 1;
    }
}

wchar_t *read_file_acp(const wchar_t *path)
{
    FILE *fp;
    long n;
    void *raw;
    wchar_t *out;
    int nw;

    fp = _wfopen(path, L"rb");
    if (fp == NULL)
        return NULL;
    fseek(fp, 0, SEEK_END);
    n = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    raw = calloc((size_t)n + 1, 1);
    fread(raw, 1, (size_t)n, fp);
    nw = MultiByteToWideChar(0, 1, (LPCCH)raw, -1, NULL, 0);
    out = (wchar_t *)calloc((size_t)nw, sizeof(wchar_t));
    MultiByteToWideChar(0, 1, (LPCCH)raw, -1, out, nw);
    free(raw);
    fclose(fp);
    return out;
}

wchar_t **list_files(const wchar_t *dir, const wchar_t *suffix, size_t *count)
{
    wchar_t *pref;
    wchar_t *pat;
    HANDLE h;
    WIN32_FIND_DATAW fd;
    wchar_t **raw;
    wchar_t **out;
    size_t nraw;
    size_t nout;
    size_t i;
    signed int slen;

    *count = 0;
    pref = path_dir_prefix(dir);
    pat = wstr_concat(pref, L"*");
    free(pref);
    h = FindFirstFileW(pat, &fd);
    if (h == INVALID_HANDLE_VALUE) {
        free(pat);
        return NULL;
    }
    nraw = 0;
    raw = NULL;
    do {
        size_t ln = wcslen(fd.cFileName);
        wchar_t *name = (wchar_t *)malloc((ln + 1) * sizeof(wchar_t));
        wcscpy(name, fd.cFileName);
        name[ln] = 0;
        raw = (wchar_t **)realloc(raw, (nraw + 1) * sizeof(wchar_t *));
        raw[nraw++] = name;
    } while (FindNextFileW(h, &fd));
    FindClose(h);
    free(pat);

    slen = (signed int)wcslen(suffix);
    out = (wchar_t **)calloc(nraw, sizeof(wchar_t *));
    nout = 0;
    for (i = 0; i < nraw; i++) {
        signed int ln = (signed int)wcslen(raw[i]);
        if (slen < ln && wcscmp(&raw[i][ln - slen], suffix) == 0)
            out[nout++] = raw[i];
        else
            free(raw[i]);
    }
    free(raw);
    *count = nout;
    return out;
}
