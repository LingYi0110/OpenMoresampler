#ifndef MORESAMPLER_FLAG_H
#define MORESAMPLER_FLAG_H

// Indices into the flag table; order matches flag_parse names[].
enum {
    FLG_g = 0,
    FLG_t,
    FLG_Y,
    FLG_e,
    FLG_b,
    FLG_A,
    FLG_P,
    FLG_u,
    FLG_Mt,
    FLG_Mb,
    FLG_Md,
    FLG_Mm,
    FLG_Ms,
    FLG_ME,
    FLG_Mo,
    FLG_Mr,
    FLG_MC,
    FLG_MG,
    FLG_MD,
    FLG__m,
    FLG_Me,
    FLG__f,
    FLG__o,
    FLG__s,
    FLG_Mp,
    FLG_COUNT
};

typedef struct flag_table {
    float *values;
    int *set_flags;
    int capacity;
} flag_table;

extern int g_flag_capacity;
extern int g_flag_count;
extern flag_table **g_flag_tables;

void flag_init(void);
flag_table *flag_create(void);
void flag_apply(flag_table *dst, flag_table *src, float scale);

// Parse a UTAU flag string, including meta flags.
void flag_parse(flag_table *tbl, const char *str);

#endif // MORESAMPLER_FLAG_H
