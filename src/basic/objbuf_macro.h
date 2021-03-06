#ifndef OBJBUF_H_
#define OBJBUF_H_

#include <stdlib.h>

//
// objbuf_macro.h
// Define a buffer with functions to alloc, dealloc, resize, add, append, reset
//
// Example:
//
//   create_objbuf(charbuf,String,char)
//
// Creates:
//
//   typedef struct {
//     char *data;
//     size_t len, capacity;
//   } String;
//
//   void charbuf_alloc(String *buf, size_t capacity)
//   void charbuf_dealloc(String *buf)
//   void charbuf_ensure_capacity(String *buf, size_t capacity)
//   void charbuf_add(String *buf, char obj)
//   int charbuf_attempt_add(String *buf, char obj)
//   void charbuf_append(String *buf, char *obj, size_t n)
//   void charbuf_reset(String *buf)
//

// Round a number up to the nearest number that is a power of two
#ifndef roundup2pow
  #define roundup2pow(x) (1UL << (64 - leading_zeros(x)))
#endif

#define OBJBUF_INIT {.data = NULL, .len = 0, .capacity = 0}

#define create_objbuf(FUNC,buf_t,obj_t)                                        \
                                                                               \
typedef struct {                                                               \
  obj_t *data;                                                                 \
  size_t len, capacity;                                                        \
} buf_t;                                                                       \
                                                                               \
/* Define functions with unused attribute in case they're not used */          \
static inline void FUNC ## _alloc(buf_t *buf, size_t capacity)                 \
 __attribute__((unused));                                                      \
static inline void FUNC ## _dealloc(buf_t *buf)                                \
 __attribute__((unused));                                                      \
static inline void FUNC ## _ensure_capacity(buf_t *buf, size_t cap)            \
 __attribute__((unused));                                                      \
static inline size_t FUNC ## _add(buf_t *buf, obj_t obj)                       \
 __attribute__((unused));                                                      \
static inline int FUNC ## _attempt_add(buf_t *buf, obj_t obj)                  \
 __attribute__((unused));                                                      \
static inline size_t FUNC ## _append(buf_t *buf, obj_t *obj, size_t n)         \
 __attribute__((unused));                                                      \
static inline void FUNC ## _reset(buf_t *buf)                                  \
 __attribute__((unused));                                                      \
                                                                               \
                                                                               \
static inline void FUNC ## _alloc(buf_t *buf, size_t capacity) {               \
  buf->capacity = capacity;                                                    \
  buf->data = ctx_malloc(sizeof(obj_t) * buf->capacity);                       \
  buf->len = 0;                                                                \
}                                                                              \
                                                                               \
static inline void FUNC ## _dealloc(buf_t *buf) {                              \
  ctx_free(buf->data);                                                         \
  buf->data = NULL;                                                            \
  buf->len = buf->capacity = 0;                                                \
}                                                                              \
                                                                               \
static inline void FUNC ## _ensure_capacity(buf_t *buf, size_t cap) {          \
  if(cap > buf->capacity) {                                                    \
    buf->capacity = roundup2pow(cap);                                          \
    buf->data = ctx_realloc(buf->data, sizeof(obj_t) * buf->capacity);         \
  }                                                                            \
}                                                                              \
                                                                               \
/* Returns index of new object in buffer */                                    \
static inline size_t FUNC ## _add(buf_t *buf, obj_t obj) {                     \
  FUNC ## _ensure_capacity(buf, buf->len+1);                                   \
  size_t idx = buf->len;                                                       \
  memcpy(buf->data+idx, &obj, sizeof(obj));                                    \
  buf->len++;                                                                  \
  return idx;                                                                  \
}                                                                              \
                                                                               \
static inline int FUNC ## _attempt_add(buf_t *buf, obj_t obj) {                \
  if(buf->len >= buf->capacity) return 0;                                      \
  memcpy(buf->data+buf->len, &obj, sizeof(obj));                               \
  buf->len++; return 1;                                                        \
}                                                                              \
                                                                               \
/* Returns index of first new object in buffer */                              \
static inline size_t FUNC ## _append(buf_t *buf, obj_t *obj, size_t n) {       \
  FUNC ## _ensure_capacity(buf, buf->len+n);                                   \
  memcpy(buf->data+buf->len, obj, n*sizeof(obj_t));                            \
  size_t idx = buf->len;                                                       \
  buf->len += n;                                                               \
  return idx;                                                                  \
}                                                                              \
                                                                               \
static inline void FUNC ## _reset(buf_t *buf) { buf->len = 0; }                \
\


#endif /* OBJBUF_H_ */
