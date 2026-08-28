#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <windows.h>
#include <io.h>

#include "xchg.h"

char *xchg_strdup(const char *str)
{
    size_t size;
    char *copy;

    size = strlen(str) + 1;
    copy = malloc(size);
    memcpy(copy, str, size);
    return copy;
}

xchg_node *xchg_node_create(int node_type, int rows, int columns)
{
    xchg_node *node = malloc(sizeof(xchg_node));
    node->type = node_type;
    node->children = NULL;
    node->keys = NULL;
    node->float_data = NULL;
    node->byte_data = NULL;
    node->float_value = 0.0f;
    node->rows = rows;
    node->columns = columns;

    if ((node_type == XCHG_TYPE_DICT || node_type == XCHG_TYPE_LIST ||
         node_type == XCHG_TYPE_OFFSET_ARRAY) &&
        (node->children = calloc(rows * columns, sizeof(void *)),
         node_type == XCHG_TYPE_DICT)) {
        node->keys = calloc(rows * columns, sizeof(void *));
        return node;
    } else if (node_type == XCHG_TYPE_FLOAT_ARRAY) {
        node->float_data = calloc(rows * columns, sizeof(float));
        return node;
    } else {
        if (node_type == XCHG_TYPE_BYTE_ARRAY)
            node->byte_data = malloc(rows * columns);
        return node;
    }
}

void xchg_node_destroy_impl(xchg_node *node, int shallow)
{
    if (!node)
        return;

    int t = node->type;
    if ((t == XCHG_TYPE_DICT || t == XCHG_TYPE_LIST ||
         t == XCHG_TYPE_OFFSET_ARRAY) && !shallow) {
        int count = node->rows * node->columns;
        if (count > 0) {
            for (int i = 0; i < count; i++) {
                if (node->children[i])
                    xchg_node_destroy_impl(node->children[i], 0);
            }
        }
    }
    free(node->children);

    if (node->type == XCHG_TYPE_DICT) {
        int count = node->rows * node->columns;
        if (count > 0) {
            for (int i = 0; i < count; i++) {
                if (node->keys[i])
                    free(node->keys[i]);
            }
        }
        free(node->keys);
    }

    if (node->type == XCHG_TYPE_FLOAT_ARRAY)
        free(node->float_data);

    if (node->type == XCHG_TYPE_BYTE_ARRAY)
        free(node->byte_data);

    free(node);
}

void xchg_node_destroy(xchg_node *node)
{
    xchg_node_destroy_impl(node, 0);
}

void xchg_node_destroy_shallow(xchg_node *node)
{
    int t;
    int count;
    int i;
    void **keys;

    if (node == NULL)
        return;
    t = node->type;
    if (((unsigned int)(t - 1) <= 1 || t == XCHG_TYPE_OFFSET_ARRAY)
        && (free(node->children), t = node->type, t == XCHG_TYPE_DICT)) {
        count = node->rows * node->columns;
        keys = node->keys;
        if (count > 0) {
            for (i = 0; i < count; i++) {
                if (keys[i] != NULL) {
                    free(keys[i]);
                    count = node->rows * node->columns;
                    keys = node->keys;
                }
            }
        }
        free(keys);
        t = node->type;
        if (t == XCHG_TYPE_FLOAT_ARRAY)
            goto free_fdata;
    } else if (t == XCHG_TYPE_FLOAT_ARRAY) {
free_fdata:
        free(node->float_data);
        if (node->type == XCHG_TYPE_BYTE_ARRAY)
            goto free_bdata;
        goto free_node;
    }
    if (t == XCHG_TYPE_BYTE_ARRAY) {
free_bdata:
        free(node->byte_data);
        free(node);
        return;
    }
free_node:
    free(node);
}

void *xchg_node_set_field(xchg_node *node, void *value, char *key, int index)
{
    node->children[index] = value;
    size_t key_size = strlen(key) + 1;
    void *key_copy = malloc(key_size);
    memcpy(key_copy, key, key_size);
    node->keys[index] = key_copy;
    return key_copy;
}


int xchg_node_has_fields(xchg_node *node, const char *first, ...)
{
    va_list ap;
    const char *key;
    int count;
    int i;

    if (!first)
        return 1;
    if (!node)
        return 0;
    va_start(ap, first);
    key = first;
    count = node->rows * node->columns;
    while (key) {
        int found = 0;
        if (count <= 0) {
            va_end(ap);
            return 0;
        }
        for (i = 0; i < count; i++) {
            if (node->keys && node->keys[i] &&
                strcmp(key, (char *)node->keys[i]) == 0) {
                if (!node->children || !node->children[i]) {
                    va_end(ap);
                    return 0;
                }
                found = 1;
                break;
            }
        }
        if (!found) {
            va_end(ap);
            return 0;
        }
        key = va_arg(ap, const char *);
    }
    va_end(ap);
    return 1;
}

xchg_node *xchg_node_get_field(xchg_node *node, char *key)
{
    int count = node->rows * node->columns;
    if (count <= 0)
        return NULL;
    for (int i = 0; i < count; i++) {
        if (strcmp(key, (char *)node->keys[i]) == 0)
            return node->children[i];
    }
    return NULL;
}

void xchg_node_serialize_counted(const char *data_name, xchg_node *node, FILE *stream)
{
    int name_len = strlen(data_name);
    if (name_len > 0) {
        fputc(name_len, stream);
        fwrite(data_name, 1, name_len, stream);
    }

    int buffer;
    float float_buffer;

    switch (node->type) {
    case XCHG_TYPE_DICT: {
        fputc(XCHG_TYPE_DICT, stream);
        buffer = node->rows * node->columns;
        fwrite(&buffer, 4, 1, stream);
        int count = node->rows * node->columns;
        for (int i = 0; i < count; i++) {
            if (node->keys[i])
                xchg_node_serialize_counted((char *)node->keys[i],
                                 node->children[i], stream);
        }
        break;
    }
    case XCHG_TYPE_LIST: {
        fputc(XCHG_TYPE_LIST, stream);
        buffer = node->rows;
        fwrite(&buffer, 4, 1, stream);
        buffer = node->columns;
        fwrite(&buffer, 4, 1, stream);
        int count = node->rows * node->columns;
        for (int i = 0; i < count; i++)
            xchg_node_serialize_counted("", node->children[i], stream);
        break;
    }
    case XCHG_TYPE_FLOAT: {
        fputc(XCHG_TYPE_FLOAT, stream);
        float_buffer = node->float_value;
        fwrite(&float_buffer, 4, 1, stream);
        break;
    }
    case XCHG_TYPE_FLOAT_ARRAY: {
        fputc(XCHG_TYPE_FLOAT_ARRAY, stream);
        buffer = node->rows;
        fwrite(&buffer, 4, 1, stream);
        buffer = node->columns;
        fwrite(&buffer, 4, 1, stream);
        int count = node->rows * node->columns;
        float *buf = calloc(count, sizeof(float));
        if (count > 0) {
            float *src = node->float_data;
            for (int i = 0; i < count; i++)
                buf[i] = src[i];
        }
        fwrite(buf, 4, count, stream);
        free(buf);
        break;
    }
    case XCHG_TYPE_BYTE_ARRAY: {
        fputc(XCHG_TYPE_BYTE_ARRAY, stream);
        buffer = node->rows * node->columns;
        fwrite(&buffer, 4, 1, stream);
        fwrite(node->byte_data, 1, node->rows * node->columns, stream);
        break;
    }
    case XCHG_TYPE_OFFSET_ARRAY: {
        fputc(XCHG_TYPE_OFFSET_ARRAY, stream);
        buffer = node->rows;
        fwrite(&buffer, 4, 1, stream);
        int count = node->rows;
        int table_size = count + 1;
        int *offsets = calloc(table_size, sizeof(int));
        fwrite(offsets, 4, table_size, stream);
        for (int i = 0; i < count; i++) {
            offsets[i] = ftell(stream);
            xchg_node_serialize_counted("", node->children[i], stream);
        }
        offsets[count] = ftell(stream);
        long end_pos = offsets[count];
        long first_pos = offsets[0];
        long rel = first_pos - end_pos - (4 * table_size);
        fseek(stream, rel, SEEK_CUR);
        fwrite(offsets, 4, table_size, stream);
        fseek(stream, end_pos - first_pos, SEEK_CUR);
        free(offsets);
        break;
    }
    default:
        fwrite("Object type not supported.\n", 1, 27, stderr);
        break;
    }
}

void xchg_node_serialize_terminated(const char *data_name, xchg_node *node, FILE *stream)
{
    size_t name_len = strlen(data_name);
    fwrite(data_name, 1, name_len + 1, stream);

    int buffer;
    float float_buffer;

    switch (node->type) {
    case XCHG_TYPE_DICT: {
        fputc(XCHG_TYPE_DICT, stream);
        int count = node->rows * node->columns;
        for (int i = 0; i < count; i++) {
            if (node->keys[i])
                xchg_node_serialize_terminated((char *)node->keys[i],
                                    node->children[i], stream);
        }
        break;
    }
    case XCHG_TYPE_LIST: {
        fputc(XCHG_TYPE_LIST, stream);
        buffer = node->rows;
        fwrite(&buffer, 4, 1, stream);
        buffer = node->columns;
        fwrite(&buffer, 4, 1, stream);
        int count = node->rows * node->columns;
        for (int i = 0; i < count; i++)
            xchg_node_serialize_terminated(" ", node->children[i], stream);
        break;
    }
    case XCHG_TYPE_FLOAT: {
        fputc(XCHG_TYPE_FLOAT, stream);
        float_buffer = node->float_value;
        fwrite(&float_buffer, 4, 1, stream);
        break;
    }
    case XCHG_TYPE_FLOAT_ARRAY: {
        fputc(XCHG_TYPE_FLOAT_ARRAY, stream);
        buffer = node->rows;
        fwrite(&buffer, 4, 1, stream);
        buffer = node->columns;
        fwrite(&buffer, 4, 1, stream);
        int count = node->rows * node->columns;
        float *buf = calloc(count, sizeof(float));
        if (count > 0) {
            float *src = node->float_data;
            for (int i = 0; i < count; i++)
                buf[i] = src[i];
        }
        fwrite(buf, 4, count, stream);
        free(buf);
        break;
    }
    default:
        fwrite("Object type not supported.\n", 1, 27, stderr);
        break;
    }
    fputc(0, stream);
}

xchg_node *xchg_node_deserialize_counted(FILE *stream, float *params)
{
    int type = fgetc(stream);
    int buffer;

    switch (type) {
    case XCHG_TYPE_DICT: {
        fread(&buffer, 4, 1, stream);
        xchg_node *node = xchg_node_create(XCHG_TYPE_DICT, buffer, 1);
        if (node->rows <= 0)
            return node;
        for (int i = 0; i < node->rows; i++) {
            int key_len = fgetc(stream);
            char key_buf[284];
            fread(key_buf, 1, key_len, stream);
            key_buf[key_len] = 0;
            xchg_node *child = xchg_node_deserialize_counted(stream, params);
            if (!child) {
                xchg_node_destroy_impl(node, 0);
                return NULL;
            }
            if (key_len == 2 && key_buf[0] == '_' && key_buf[1] == 'c')
                *params = child->float_value;
            node->children[i] = child;
            size_t name_size = strlen(key_buf) + 1;
            void *key_copy = malloc(name_size);
            memcpy(key_copy, key_buf, name_size);
            node->keys[i] = key_copy;
        }
        return node;
    }
    case XCHG_TYPE_LIST: {
        int dim1, dim2;
        fread(&dim1, 4, 1, stream);
        fread(&dim2, 4, 1, stream);
        xchg_node *node = xchg_node_create(XCHG_TYPE_LIST, dim1, dim2);
        int count = dim1 * dim2;
        if (count <= 0)
            return node;
        for (int i = 0; i < count; i++) {
            node->children[i] = xchg_node_deserialize_counted(stream, params);
            if (!node->children[i]) {
                xchg_node_destroy_impl(node, 0);
                return NULL;
            }
        }
        return node;
    }
    case XCHG_TYPE_FLOAT: {
        xchg_node *node = malloc(sizeof(xchg_node));
        node->type = XCHG_TYPE_FLOAT;
        node->children = NULL;
        node->keys = NULL;
        node->float_data = NULL;
        node->byte_data = NULL;
        node->float_value = 0.0f;
        node->rows = 0;
        node->columns = 0;
        fread(&node->float_value, 4, 1, stream);
        return node;
    }
    case XCHG_TYPE_FLOAT_ARRAY: {
        int dim1, dim2;
        fread(&dim1, 4, 1, stream);
        fread(&dim2, 4, 1, stream);
        xchg_node *node = xchg_node_create(XCHG_TYPE_FLOAT_ARRAY, dim1, dim2);
        int count = dim1 * dim2;
        void *buf = calloc(count, 4);
        fread(buf, 4, count, stream);
        if (count > 0) {
            float *dst = node->float_data;
            float *src = buf;
            for (int i = 0; i < count; i++)
                dst[i] = src[i];
        }
        free(buf);
        return node;
    }
    case XCHG_TYPE_BYTE_ARRAY: {
        fread(&buffer, 4, 1, stream);
        int size = buffer;
        xchg_node *node = xchg_node_create(XCHG_TYPE_BYTE_ARRAY, size, 1);
        fread(node->byte_data, 1, size, stream);
        return node;
    }
    case XCHG_TYPE_OFFSET_ARRAY: {
        fread(&buffer, 4, 1, stream);
        int count = buffer;
        float frame_shift = params[0];
        float start_time = params[1];
        float end_time = params[2];
        int start_frame = (int)rint(start_time / frame_shift);
        if (end_time <= 0.0f) {
            end_time = end_time + (float)count * frame_shift;
            params[2] = end_time;
        }
        int table_size = count + 1;
        int read_count = count - start_frame;
        int total_request = (int)rint((end_time - start_time) / frame_shift) + 5;
        if (count >= start_frame + total_request)
            read_count = total_request;
        int *offsets = malloc(4 * table_size);
        fread(offsets, 4, table_size, stream);
        xchg_node *node = xchg_node_create(XCHG_TYPE_OFFSET_ARRAY, count, 1);
        fseek(stream, offsets[start_frame] - offsets[0], SEEK_CUR);
        if (count > 0)
            memset(node->children, 0, count * sizeof(void *));
        for (int i = 0; i < read_count; i++) {
            node->children[start_frame + i] = xchg_node_deserialize_counted(stream, params);
            if (!node->children[start_frame + i]) {
                xchg_node_destroy_impl(node, 0);
                free(offsets);
                return NULL;
            }
        }
        fseek(stream, offsets[count], SEEK_SET);
        free(offsets);
        return node;
    }
    default:
        fwrite("Object type not supported.\n", 1, 27, stderr);
        return NULL;
    }
}

int xchg_node_deserialize_terminated_into(xchg_node *parent, FILE *stream)
{
    char *name = malloc(32);
    int name_cap = 32;
    int name_len = 0;
    int ch = fgetc(stream);
    if (ch != 0) {
        name[0] = ch;
        name_len = 1;
        while (1) {
            ch = fgetc(stream);
            if (ch == 0)
                break;
            name[name_len++] = ch;
            if (name_cap < name_len) {
                name_cap += 16;
                name = realloc(name, name_cap);
            }
        }
    }
    name[name_len] = 0;

    int type = fgetc(stream);
    parent->children = realloc(parent->children,
            sizeof(void *) * (parent->rows + 1));
    parent->keys = realloc(parent->keys,
            sizeof(void *) * (parent->rows + 1));
    int idx = parent->rows;
    parent->keys[idx] = name;
    parent->rows = idx + 1;

    switch (type) {
    case XCHG_TYPE_DICT: {
        xchg_node *child = malloc(sizeof(xchg_node));
        child->type = XCHG_TYPE_DICT;
        child->float_data = NULL;
        child->byte_data = NULL;
        child->float_value = 0.0f;
        child->rows = 0;
        child->columns = 1;
        child->children = calloc(0, sizeof(void *));
        child->keys = calloc(0, sizeof(void *));
        parent->children[idx] = child;
        int next = fgetc(stream);
        if (next == 0) {
            fseek(stream, -1, SEEK_CUR);
        } else {
            while (1) {
                fseek(stream, -1, SEEK_CUR);
                if (xchg_node_deserialize_terminated_into(child, stream) == 0)
                    return 0;
                next = fgetc(stream);
                if (next == 0) {
                    fseek(stream, -1, SEEK_CUR);
                    break;
                }
            }
        }
        break;
    }
    case XCHG_TYPE_LIST: {
        int dim1, dim2;
        fread(&dim1, 4, 1, stream);
        fread(&dim2, 4, 1, stream);
        xchg_node *child = xchg_node_create(XCHG_TYPE_LIST, dim1, dim2);
        parent->children[idx] = child;
        int count = dim1 * dim2;
        for (int i = 0; i < count; i++) {
            xchg_node *elem = malloc(sizeof(xchg_node));
            elem->type = XCHG_TYPE_DICT;
            elem->float_data = NULL;
            elem->byte_data = NULL;
            elem->float_value = 0.0f;
            elem->rows = 0;
            elem->columns = 1;
            elem->children = calloc(0, sizeof(void *));
            elem->keys = calloc(0, sizeof(void *));
            if (xchg_node_deserialize_terminated_into(elem, stream) == 0)
                return 0;
            void *value = NULL;
            int elem_count = elem->rows * elem->columns;
            for (int j = 0; j < elem_count; j++) {
                if (strcmp(" ", (char *)elem->keys[j]) == 0) {
                    value = elem->children[j];
                    break;
                }
            }
            ((void **)child->children)[i] = value;
            xchg_node_destroy_impl(elem, 0);
        }
        break;
    }
    case XCHG_TYPE_FLOAT: {
        xchg_node *child = malloc(sizeof(xchg_node));
        child->type = XCHG_TYPE_FLOAT;
        child->children = NULL;
        child->keys = NULL;
        child->float_data = NULL;
        child->byte_data = NULL;
        child->float_value = 0.0f;
        child->rows = 0;
        child->columns = 1;
        parent->children[idx] = child;
        fread(&child->float_value, 4, 1, stream);
        break;
    }
    case XCHG_TYPE_FLOAT_ARRAY: {
        int dim1, dim2;
        fread(&dim1, 4, 1, stream);
        fread(&dim2, 4, 1, stream);
        xchg_node *child = xchg_node_create(XCHG_TYPE_FLOAT_ARRAY, dim1, dim2);
        parent->children[idx] = child;
        int count = dim1 * dim2;
        void *buf = calloc(count, 4);
        fread(buf, 4, count, stream);
        if (count > 0) {
            float *dst = child->float_data;
            float *src = buf;
            for (int i = 0; i < count; i++)
                dst[i] = src[i];
        }
        free(buf);
        break;
    }
    default:
        fwrite("Object type not supported.\n", 1, 27, stderr);
        return 0;
    }
    fseek(stream, 1, SEEK_CUR);
    return 1;
}

int xchg_export_counted(xchg_node *node, const wchar_t *filename)
{
    FILE *stream = _wfopen(filename, L"wb");
    if (!stream)
        return 0;
    xchg_file_lock(stream);
    xchg_node_serialize_counted("data", node, stream);
    fclose(stream);
    return 1;
}

int xchg_export_terminated(xchg_node *node, const wchar_t *filename)
{
    FILE *stream = _wfopen(filename, L"wb");
    if (!stream)
        return 0;
    xchg_file_lock(stream);
    xchg_node_serialize_terminated("data", node, stream);
    fclose(stream);
    return 1;
}

xchg_node *xchg_import(const wchar_t *filename, float start, float end)
{
    FILE *stream = _wfopen(filename, L"rb");
    if (!stream)
        return NULL;
    xchg_file_lock(stream);

    unsigned char header[4];
    fread(header, 1, 4, stream);
    fseek(stream, 0, SEEK_SET);

    xchg_node *result = NULL;

    if (header[0] == 'd' && header[1] == 'a' &&
        header[2] == 't' && header[3] == 'a') {
        xchg_node *wrapper = malloc(sizeof(xchg_node));
        wrapper->type = XCHG_TYPE_DICT;
        wrapper->float_data = NULL;
        wrapper->byte_data = NULL;
        wrapper->float_value = 0.0f;
        wrapper->rows = 0;
        wrapper->columns = 1;
        wrapper->children = calloc(0, sizeof(void *));
        wrapper->keys = calloc(0, sizeof(void *));
        if (xchg_node_deserialize_terminated_into(wrapper, stream) == 0) {
            fclose(stream);
            xchg_node_destroy_impl(wrapper, 0);
            return NULL;
        }
        result = wrapper;
    } else {
        int name_len = fgetc(stream);
        char name_buf[284];
        fread(name_buf, 1, name_len, stream);
        name_buf[name_len] = 0;
        float params[3];
        params[0] = XCHG_FRAME_SHIFT;
        params[1] = start;
        params[2] = end;
        xchg_node *child = xchg_node_deserialize_counted(stream, params);
        if (!child) {
            fclose(stream);
            return NULL;
        }
        xchg_node *wrapper = malloc(sizeof(xchg_node));
        wrapper->type = XCHG_TYPE_DICT;
        wrapper->float_data = NULL;
        wrapper->byte_data = NULL;
        wrapper->float_value = 0.0f;
        wrapper->rows = 1;
        wrapper->columns = 1;
        wrapper->children = calloc(1, sizeof(void *));
        wrapper->keys = calloc(1, sizeof(void *));
        wrapper->children[0] = child;
        size_t name_size = strlen(name_buf) + 1;
        void *name_copy = malloc(name_size);
        memcpy(name_copy, name_buf, name_size);
        wrapper->keys[0] = name_copy;
        result = wrapper;
    }

    fclose(stream);
    return result;
}

int xchg_file_lock(FILE *stream)
{
    int i;
    BOOL ok;
    for (i = 0;; i++) {
        OVERLAPPED overlapped = {0};
        overlapped.Offset = 4;
        int fd = fileno(stream);
        void *handle = (void *)_get_osfhandle(fd);
        ok = LockFileEx(handle, LOCKFILE_EXCLUSIVE_LOCK | LOCKFILE_FAIL_IMMEDIATELY,
                        0, 4, 0, &overlapped);
        if (ok)
            break;
        Sleep(300);
        printf("wait_lock: IO waiting (%d)...\n", i);
    }
    return ok;
}
