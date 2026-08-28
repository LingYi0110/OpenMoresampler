#ifndef MORESAMPLER_XCHG_H
#define MORESAMPLER_XCHG_H

#include <stdio.h>
#include <wchar.h>

#define XCHG_FRAME_SHIFT 0.0058049886f

enum xchg_node_type {
    XCHG_TYPE_DICT = 1,
    XCHG_TYPE_LIST = 2,
    XCHG_TYPE_FLOAT = 3,
    XCHG_TYPE_FLOAT_ARRAY = 5,
    XCHG_TYPE_BYTE_ARRAY = 6,
    XCHG_TYPE_OFFSET_ARRAY = 7,
};

typedef struct xchg_node {
    int type;
    void **children;
    void **keys;
    void *float_data;
    float float_value;
    void *byte_data;
    int rows;
    int columns;
} xchg_node;

xchg_node *xchg_node_create(int node_type, int rows, int columns);
void xchg_node_destroy(xchg_node *node);
void xchg_node_destroy_impl(xchg_node *node, int shallow);
void xchg_node_destroy_shallow(xchg_node *node);
char *xchg_strdup(const char *str);

void xchg_node_serialize_counted(const char *data_name, xchg_node *node, FILE *stream);
void xchg_node_serialize_terminated(const char *data_name, xchg_node *node, FILE *stream);
xchg_node *xchg_node_deserialize_counted(FILE *stream, float *params);
int xchg_node_deserialize_terminated_into(xchg_node *parent, FILE *stream);

int xchg_export_counted(xchg_node *node, const wchar_t *filename);
int xchg_export_terminated(xchg_node *node, const wchar_t *filename);
xchg_node *xchg_import(const wchar_t *filename, float start, float end);

void *xchg_node_set_field(xchg_node *node, void *value, char *key, int index);
xchg_node *xchg_node_get_field(xchg_node *node, char *key);
int xchg_node_has_fields(xchg_node *node, const char *first, ...);

int xchg_file_lock(FILE *stream);

#endif
