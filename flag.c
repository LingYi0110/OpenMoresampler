#include "flag.h"

#include <locale.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"

int g_flag_capacity = 0;
int g_flag_count = 0;
flag_table **g_flag_tables = NULL;

void flag_init(void)
{
    setlocale(LC_ALL, "ja_JP.utf8");
    g_flag_count = 0;
    g_flag_tables = calloc(10, sizeof(flag_table *));
    g_flag_capacity = 10;
}

flag_table *flag_create(void)
{
    flag_table *m = malloc(sizeof(flag_table));
    m->capacity = FLG_COUNT;
    m->values = calloc(FLG_COUNT, sizeof(float));
    m->set_flags = calloc(FLG_COUNT, sizeof(int));
    return m;
}

void flag_apply(flag_table *dst, flag_table *src, float scale)
{
    int n = src->capacity;
    if (n <= 0)
        return;

    if (scale == 100000.0f) {
        // Assign source values without scaling.
        for (int i = 0; i < n; i++) {
            if (src->set_flags[i] != 0 && dst->set_flags[i] == 0) {
                dst->values[i] = src->values[i];
                dst->set_flags[i] = 2;
            }
        }
        return;
    }

    for (int i = 0; i < n; i++) {
        if (src->set_flags[i] == 0)
            continue;
        if (dst->set_flags[i] == 1 || dst->set_flags[i] == 0) {
            dst->values[i] = src->values[i] * scale * 0.01f + dst->values[i];
            if (dst->set_flags[i] == 0)
                dst->set_flags[i] = 1;
        }
    }
}

// Parse a UTAU flag string and apply recognized flag values.
float flag_parse(float default_value, flag_table *tbl, const char *str)
{
    if (!str || !tbl)
        return default_value;

    signed int len = (signed int)strlen(str);
    if (len <= 0)
        return default_value;

    int i = 0;
    while (i < len) {
        char tag[4] = {0, 0, 0, 0};
        int meta_idx = -1;

        char c = str[i];
        if (c == 'M' || c == ':') {
            // Continue with tag parsing below.
        }

        if (c == '_') {
            if (++i >= len) {
                log_tee_printf("Error: invalid flag sequence '%s'.\n", str);
                return default_value;
            }
            tag[0] = '_';
            tag[1] = str[i];
            i++;
        } else if (c == '/') {
            if (++i >= len) {
                log_tee_printf("Error: invalid flag sequence '%s'.\n", str);
                return default_value;
            }
            // Skip slash form.
            continue;
        } else {
            tag[0] = str[i++];
            if (tag[0] == 'M' && i < len && (unsigned)(str[i] - '0') <= 9) {
                // Parse a meta flag M<digits>[.].
                int j = i;
                while (j < len && (unsigned)(str[j] - '0') <= 9)
                    j++;
                char tmp[32];
                int n = j - i;
                if (n >= (int)sizeof(tmp))
                    n = (int)sizeof(tmp) - 1;
                memcpy(tmp, str + i, (size_t)n);
                tmp[n] = 0;
                meta_idx = atoi(tmp);
                i = j;
                if (i < len && str[i] == '.')
                    i++;
                tag[0] = 'M';
                tag[1] = 0;
            } else if (i < len) {
                // Parse two-letter tags like Mt and Mb.
                char c2 = str[i];
                if (c2 != '+' && c2 != '-' &&
                    !((c2 >= '0' && c2 <= '9') || c2 == '.')) {
                    tag[1] = c2;
                    i++;
                }
            }
        }

        // Parse an optional signed number.
        int num_start = i;
        if (i < len && (str[i] == '+' || str[i] == '-' ||
                        (str[i] >= '0' && str[i] <= '9') || str[i] == '.')) {
            if (str[i] == '+' || str[i] == '-')
                i++;
            while (i < len &&
                   ((str[i] >= '0' && str[i] <= '9') || str[i] == '.'))
                i++;
        }
        char numbuf[64];
        int nlen = i - num_start;
        if (nlen <= 0) {
            numbuf[0] = 0;
        } else {
            if (nlen >= (int)sizeof(numbuf))
                nlen = (int)sizeof(numbuf) - 1;
            memcpy(numbuf, str + num_start, (size_t)nlen);
            numbuf[nlen] = 0;
        }
        double v = numbuf[0] ? atof(numbuf) : 0.0;
        if (v > 10000.0)
            v = 10000.0;
        if (v < -10000.0)
            v = -10000.0;

        if (meta_idx >= 0) {
            if (meta_idx < g_flag_capacity && g_flag_tables[meta_idx]) {
                float scale = numbuf[0] ? (float)v : 100000.0f;
                flag_apply(tbl, g_flag_tables[meta_idx], scale);
            } else {
                log_tee_printf("Error: undefined meta flag M%d. Skipped.\n",
                               meta_idx);
            }
            continue;
        }

        // Match the built-in flag table.
        static const char *names[] = {
            "g",  "t",  "Y",  "e",  "b",  "A",  "P",  "u",
            "Mt", "Mb", "Md", "Mm", "Ms", "ME", "Mo", "Mr",
            "MC", "MG", "MD", "_m", "Me", "_f", "_o", "_s", "Mp"
        };
        int nflags = (int)(sizeof(names) / sizeof(names[0]));
        int matched = -1;
        for (int k = 0; k < nflags; k++) {
            if (strcmp(tag, names[k]) == 0) {
                matched = k;
                break;
            }
        }
        if (matched >= 0 && matched < tbl->capacity) {
            if (tbl->set_flags[matched] == 0) {
                tbl->values[matched] = (float)v;
                tbl->set_flags[matched] = 2;
            }
        }
    }
    return default_value;
}
