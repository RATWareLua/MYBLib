/* SPDX-License-Identifier: BSD-3-Clause */
/* NAIF text kernel reader. Format: NAIF "Kernel Required Reading" (text kernels).
 *
 * Storage: four flat arrays (variables, numbers, string descriptors, characters). A load
 * works on a copy of them and replaces the pool only on success, so a failed load leaves the
 * pool exactly as it was. */
#include <stdckdint.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "myb/io.h"
#include "platform/platform.h"

enum {
  NAME_MAX_LENGTH = 32,              /* NAIF limit for variable names */
  FILE_MAX_SIZE = 256 * 1024 * 1024, /* larger text kernels are rejected */
};

typedef struct kstring {
  size_t offset; /* into chars */
  size_t length;
} kstring;

typedef struct kvar {
  char name[NAME_MAX_LENGTH + 1];
  myb_kvar_type type;
  size_t first; /* index into numbers or strings */
  size_t count;
} kvar;

typedef struct storage {
  kvar* vars;
  size_t var_count;
  size_t var_capacity;
  double* numbers;
  size_t number_count;
  size_t number_capacity;
  kstring* strings;
  size_t string_count;
  size_t string_capacity;
  char* chars;
  size_t char_count;
  size_t char_capacity;
} storage;

struct myb_kpool {
  myb_allocator allocator;
  storage data;
};

/* ---- Storage ---------------------------------------------------------------------------- */

/* Makes room for needed elements, keeping the first used ones. */
static myb_status reserve(const myb_allocator* a, void** array, size_t used, size_t* capacity,
                          size_t needed, size_t element) {
  if (needed <= *capacity) {
    return MYB_OK;
  }
  size_t new_capacity = *capacity < 16 ? 16 : *capacity;
  /* bound: doubling, at most 64 iterations before ckd_mul overflows */
  while (new_capacity < needed) {
    if (ckd_mul(&new_capacity, new_capacity, (size_t)2)) {
      return MYB_E_OVERFLOW;
    }
  }
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(a, new_capacity, element, &block);
  if (status != MYB_OK) {
    return status;
  }
  if (used != 0) {
    memcpy(block, *array, used * element);
  }
  myb_internal_release(a, *array);
  *array = block;
  *capacity = new_capacity;
  return MYB_OK;
}

static void storage_free(const myb_allocator* a, storage* s) {
  myb_internal_release(a, s->vars);
  myb_internal_release(a, s->numbers);
  myb_internal_release(a, s->strings);
  myb_internal_release(a, s->chars);
  *s = (storage){};
}

static myb_status copy_array(const myb_allocator* a, const void* source, size_t count,
                             size_t element, void** target, size_t* capacity) {
  *target = NULL;
  *capacity = 0;
  if (count == 0) {
    return MYB_OK;
  }
  const myb_status status = myb_internal_alloc_array(a, count, element, target);
  if (status == MYB_OK) {
    memcpy(*target, source, count * element);
    *capacity = count;
  }
  return status;
}

static myb_status storage_clone(const myb_allocator* a, const storage* source, storage* target) {
  *target = (storage){.var_count = source->var_count,
                      .number_count = source->number_count,
                      .string_count = source->string_count,
                      .char_count = source->char_count};
  void* block = NULL;
  myb_status status =
      copy_array(a, source->vars, source->var_count, sizeof(kvar), &block, &target->var_capacity);
  target->vars = block;
  if (status == MYB_OK) {
    status = copy_array(a, source->numbers, source->number_count, sizeof(double), &block,
                        &target->number_capacity);
    target->numbers = block;
  }
  if (status == MYB_OK) {
    status = copy_array(a, source->strings, source->string_count, sizeof(kstring), &block,
                        &target->string_capacity);
    target->strings = block;
  }
  if (status == MYB_OK) {
    status = copy_array(a, source->chars, source->char_count, 1, &block, &target->char_capacity);
    target->chars = block;
  }
  if (status != MYB_OK) {
    storage_free(a, target);
  }
  return status;
}

static kvar* find_var(const storage* s, const char* name, size_t length) {
  if (length > NAME_MAX_LENGTH) {
    return NULL;
  }
  for (size_t i = 0; i < s->var_count; ++i) {
    if (strlen(s->vars[i].name) == length && memcmp(s->vars[i].name, name, length) == 0) {
      return &s->vars[i];
    }
  }
  return NULL;
}

/* ---- Parser state ------------------------------------------------------------------------- */

typedef struct parser {
  const myb_allocator* allocator;
  storage* s;
  const char* text;
  size_t size;
  size_t pos;
  long line;
  int in_data;
  myb_diag* diag;
} parser;

typedef enum token_kind {
  TOKEN_END,
  TOKEN_WORD,
  TOKEN_STRING,
  TOKEN_DATE,
  TOKEN_ASSIGN,
  TOKEN_APPEND,
  TOKEN_OPEN,
  TOKEN_CLOSE
} token_kind;

typedef struct token {
  token_kind kind;
  const char* start; /* WORD, DATE: the characters; STRING: inside the quotes, raw */
  size_t length;
  long line;
} token;

static myb_status fail(parser* p, long line, const char* what) {
  myb_internal_diag(p->diag, line, "line %ld: %s", line, what);
  return MYB_E_FORMAT;
}

static int is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\f' || c == '\v';
}

static int is_word_end(const parser* p, size_t i) {
  const char c = p->text[i];
  if (is_space(c) || c == '\n' || c == '=' || c == '(' || c == ')' || c == ',' || c == '\'') {
    return 1;
  }
  return c == '+' && i + 1 < p->size && p->text[i + 1] == '=';
}

/* Line [pos, end) without surrounding white space equals word. */
static int line_is(const parser* p, size_t pos, size_t end, const char* word) {
  /* bound: the line length */
  while (pos < end && is_space(p->text[pos])) {
    ++pos;
  }
  while (end > pos && is_space(p->text[end - 1])) {
    --end;
  }
  const size_t n = strlen(word);
  return end - pos == n && memcmp(p->text + pos, word, n) == 0;
}

static size_t line_end(const parser* p, size_t pos) {
  const char* newline = memchr(p->text + pos, '\n', p->size - pos);
  return newline != NULL ? (size_t)(newline - p->text) : p->size;
}

/* At the start of a line: handles \begindata / \begintext and skips commentary until data
 * or the end of the text. */
static void enter_line(parser* p) {
  /* bound: one iteration per line of the text */
  while (p->pos < p->size) {
    const size_t end = line_end(p, p->pos);
    if (line_is(p, p->pos, end, "\\begindata")) {
      p->in_data = 1;
    } else if (line_is(p, p->pos, end, "\\begintext")) {
      p->in_data = 0;
    } else if (p->in_data) {
      return;
    }
    p->pos = end < p->size ? end + 1 : end;
    p->line += 1;
  }
}

static myb_status next_token(parser* p, token* t) {
  /* bound: every iteration consumes at least one character or returns */
  while (p->pos < p->size) {
    const char c = p->text[p->pos];
    if (c == '\n') {
      p->pos += 1;
      p->line += 1;
      enter_line(p);
      continue;
    }
    if (is_space(c) || c == ',') {
      p->pos += 1;
      continue;
    }
    t->line = p->line;
    t->start = p->text + p->pos;
    t->length = 1;
    if (c == '=') {
      t->kind = TOKEN_ASSIGN;
      p->pos += 1;
      return MYB_OK;
    }
    if (c == '+' && p->pos + 1 < p->size && p->text[p->pos + 1] == '=') {
      t->kind = TOKEN_APPEND;
      t->length = 2;
      p->pos += 2;
      return MYB_OK;
    }
    if (c == '(' || c == ')') {
      t->kind = c == '(' ? TOKEN_OPEN : TOKEN_CLOSE;
      p->pos += 1;
      return MYB_OK;
    }
    if (c == '\'') {
      /* String on one line; '' stands for a quote. */
      size_t close = 0;
      int escaped = 0;
      for (size_t i = p->pos + 1; i < p->size && p->text[i] != '\n' && close == 0; ++i) {
        if (escaped) {
          escaped = 0;
        } else if (p->text[i] == '\'') {
          if (i + 1 < p->size && p->text[i + 1] == '\'') {
            escaped = 1;
          } else {
            close = i;
          }
        }
      }
      if (close == 0) {
        return fail(p, p->line, "unterminated string (strings end on their line)");
      }
      t->kind = TOKEN_STRING;
      t->start = p->text + p->pos + 1;
      t->length = close - p->pos - 1;
      p->pos = close + 1;
      return MYB_OK;
    }
    size_t i = p->pos;
    for (; i < p->size && !is_word_end(p, i); ++i) {
    }
    t->kind = c == '@' ? TOKEN_DATE : TOKEN_WORD;
    t->length = i - p->pos;
    p->pos = i;
    return MYB_OK;
  }
  t->kind = TOKEN_END;
  t->line = p->line;
  t->start = p->text + p->size;
  t->length = 0;
  return MYB_OK;
}

/* ---- Values --------------------------------------------------------------------------------- */

static myb_status push_number(parser* p, const token* t) {
  double value = 0.0;
  size_t used = 0;
  if (myb_parse_double(t->start, t->length, &value, &used) != MYB_OK || used != t->length) {
    myb_internal_diag(p->diag, t->line, "line %ld: '%.*s' is not a number", t->line,
                      (int)(t->length < 60 ? t->length : 60), t->start);
    return MYB_E_FORMAT;
  }
  storage* s = p->s;
  void* array = s->numbers;
  const myb_status status = reserve(p->allocator, &array, s->number_count, &s->number_capacity,
                                    s->number_count + 1, sizeof(double));
  s->numbers = array;
  if (status != MYB_OK) {
    return status;
  }
  s->numbers[s->number_count++] = value;
  return MYB_OK;
}

static myb_status push_string(parser* p, const token* t) {
  storage* s = p->s;
  void* array = s->chars;
  size_t needed;
  if (ckd_add(&needed, s->char_count, t->length)) {
    return MYB_E_OVERFLOW;
  }
  myb_status status = reserve(p->allocator, &array, s->char_count, &s->char_capacity, needed, 1);
  s->chars = array;
  if (status != MYB_OK) {
    return status;
  }
  const size_t offset = s->char_count;
  /* Copy, collapsing '' into ' inside quoted strings. */
  int skip = 0;
  for (size_t i = 0; i < t->length; ++i) {
    if (skip) {
      skip = 0;
      continue;
    }
    s->chars[s->char_count++] = t->start[i];
    skip = t->kind == TOKEN_STRING && t->start[i] == '\'';
  }
  array = s->strings;
  status = reserve(p->allocator, &array, s->string_count, &s->string_capacity, s->string_count + 1,
                   sizeof(kstring));
  s->strings = array;
  if (status != MYB_OK) {
    return status;
  }
  s->strings[s->string_count++] = (kstring){.offset = offset, .length = s->char_count - offset};
  return MYB_OK;
}

static myb_status push_value(parser* p, const token* t, myb_kvar_type* type) {
  const myb_kvar_type kind = t->kind == TOKEN_WORD ? MYB_KVAR_NUMBERS : MYB_KVAR_STRINGS;
  if (*type != 0 && *type != kind) {
    return fail(p, t->line, "numbers and strings mixed in one variable");
  }
  *type = kind;
  return kind == MYB_KVAR_NUMBERS ? push_number(p, t) : push_string(p, t);
}

/* Moves the values of v to the end of their array, so that appended values follow them. */
static myb_status move_to_end(parser* p, kvar* v) {
  storage* s = p->s;
  const int numeric = v->type == MYB_KVAR_NUMBERS;
  const size_t count_now = numeric ? s->number_count : s->string_count;
  if (v->first + v->count == count_now) {
    return MYB_OK;
  }
  void* array = numeric ? (void*)s->numbers : (void*)s->strings;
  size_t* capacity = numeric ? &s->number_capacity : &s->string_capacity;
  const size_t element = numeric ? sizeof(double) : sizeof(kstring);
  const myb_status status =
      reserve(p->allocator, &array, count_now, capacity, count_now + v->count, element);
  if (numeric) {
    s->numbers = array;
  } else {
    s->strings = array;
  }
  if (status != MYB_OK) {
    return status;
  }
  char* base = array;
  memcpy(base + count_now * element, base + v->first * element, v->count * element);
  v->first = count_now;
  if (numeric) {
    s->number_count += v->count;
  } else {
    s->string_count += v->count;
  }
  return MYB_OK;
}

/* NAME (= | +=) value | ( values ) */
static myb_status parse_assignment(parser* p, const token* name) {
  if (name->kind != TOKEN_WORD) {
    return fail(p, name->line, "expected a variable name");
  }
  if (name->length > NAME_MAX_LENGTH) {
    return fail(p, name->line, "variable name longer than 32 characters");
  }
  token op;
  myb_status status = next_token(p, &op);
  if (status != MYB_OK) {
    return status;
  }
  if (op.kind != TOKEN_ASSIGN && op.kind != TOKEN_APPEND) {
    return fail(p, op.line, "expected '=' or '+=' after the variable name");
  }
  storage* s = p->s;
  kvar* existing = find_var(s, name->start, name->length);
  myb_kvar_type type = 0;
  size_t first_number = s->number_count;
  size_t first_string = s->string_count;
  size_t kept = 0;
  if (op.kind == TOKEN_APPEND && existing != NULL) {
    status = move_to_end(p, existing);
    if (status != MYB_OK) {
      return status;
    }
    type = existing->type;
    kept = existing->count;
    first_number = existing->first;
    first_string = existing->first;
  }

  token t;
  status = next_token(p, &t);
  if (status != MYB_OK) {
    return status;
  }
  size_t added = 0;
  if (t.kind == TOKEN_OPEN) {
    const long open_line = t.line;
    /* bound: every value consumes a token; the list ends at ')' or the end of the text */
    for (;;) {
      status = next_token(p, &t);
      if (status != MYB_OK) {
        return status;
      }
      if (t.kind == TOKEN_CLOSE) {
        break;
      }
      if (t.kind == TOKEN_END) {
        return fail(p, open_line, "unterminated '(' list");
      }
      if (t.kind != TOKEN_WORD && t.kind != TOKEN_STRING && t.kind != TOKEN_DATE) {
        return fail(p, t.line, "expected a value");
      }
      status = push_value(p, &t, &type);
      if (status != MYB_OK) {
        return status;
      }
      ++added;
    }
  } else if (t.kind == TOKEN_WORD || t.kind == TOKEN_STRING || t.kind == TOKEN_DATE) {
    status = push_value(p, &t, &type);
    if (status != MYB_OK) {
      return status;
    }
    added = 1;
  } else {
    return fail(p, t.line, "expected a value or '('");
  }
  if (added == 0) {
    return fail(p, t.line, "empty value list");
  }

  if (existing == NULL) {
    void* array = s->vars;
    status = reserve(p->allocator, &array, s->var_count, &s->var_capacity, s->var_count + 1,
                     sizeof(kvar));
    s->vars = array;
    if (status != MYB_OK) {
      return status;
    }
    existing = &s->vars[s->var_count++];
    *existing = (kvar){};
    memcpy(existing->name, name->start, name->length);
    existing->name[name->length] = '\0';
  }
  existing->type = type;
  existing->first = type == MYB_KVAR_NUMBERS ? first_number : first_string;
  existing->count = kept + added;
  return MYB_OK;
}

static myb_status parse(parser* p) {
  enter_line(p);
  /* bound: every assignment consumes at least two tokens */
  for (;;) {
    token t;
    const myb_status status = next_token(p, &t);
    if (status != MYB_OK) {
      return status;
    }
    if (t.kind == TOKEN_END) {
      return MYB_OK;
    }
    const myb_status assigned = parse_assignment(p, &t);
    if (assigned != MYB_OK) {
      return assigned;
    }
  }
}

/* ---- Public API ------------------------------------------------------------------------------ */

myb_status myb_kpool_create(const myb_allocator* allocator, myb_kpool** pool) {
  MYB_CHECK(pool != NULL);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  const myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_kpool), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_kpool* created = block;
  *created = (myb_kpool){.allocator = a};
  *pool = created;
  return MYB_OK;
}

void myb_kpool_destroy(myb_kpool* pool) {
  if (pool == NULL) {
    return;
  }
  const myb_allocator a = pool->allocator;
  storage_free(&a, &pool->data);
  myb_internal_release(&a, pool);
}

myb_status myb_kpool_load_text(myb_kpool* pool, const char* text, size_t size, myb_diag* diag) {
  MYB_CHECK(pool != NULL);
  MYB_CHECK(text != NULL || size == 0);
  myb_internal_diag_clear(diag);
  storage work;
  myb_status status = storage_clone(&pool->allocator, &pool->data, &work);
  if (status != MYB_OK) {
    return status;
  }
  parser p = {.allocator = &pool->allocator,
              .s = &work,
              .text = text,
              .size = size,
              .pos = 0,
              .line = 1,
              .in_data = 0,
              .diag = diag};
  status = parse(&p);
  if (status != MYB_OK) {
    if (status != MYB_E_FORMAT) {
      myb_internal_diag(diag, p.line, "line %ld: %s", p.line, myb_status_name(status));
    }
    storage_free(&pool->allocator, &work);
    return status;
  }
  storage_free(&pool->allocator, &pool->data);
  pool->data = work;
  return MYB_OK;
}

myb_status myb_kpool_load_file(myb_kpool* pool, const char* path, myb_diag* diag) {
  MYB_CHECK(pool != NULL && path != NULL);
  myb_internal_diag_clear(diag);
  char* text = NULL;
  size_t size = 0;
  myb_status status = myb_platform_read_file(path, &pool->allocator, FILE_MAX_SIZE, &text, &size);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot read '%s' (%s)", path, myb_status_name(status));
    return status;
  }
  status = myb_kpool_load_text(pool, text, size, diag);
  myb_internal_release(&pool->allocator, text);
  return status;
}

static const kvar* lookup(const myb_kpool* pool, const char* name) {
  MYB_CHECK(pool != NULL && name != NULL);
  return find_var(&pool->data, name, strlen(name));
}

myb_status myb_kpool_describe(const myb_kpool* pool, const char* name, myb_kvar_type* type,
                              size_t* count) {
  const kvar* v = lookup(pool, name);
  if (v == NULL) {
    return MYB_E_NOT_FOUND;
  }
  if (type != NULL) {
    *type = v->type;
  }
  if (count != NULL) {
    *count = v->count;
  }
  return MYB_OK;
}

myb_status myb_kpool_numbers(const myb_kpool* pool, const char* name, const double** values,
                             size_t* count) {
  MYB_CHECK(values != NULL && count != NULL);
  const kvar* v = lookup(pool, name);
  if (v == NULL) {
    return MYB_E_NOT_FOUND;
  }
  if (v->type != MYB_KVAR_NUMBERS) {
    return MYB_E_FORMAT;
  }
  *values = pool->data.numbers + v->first;
  *count = v->count;
  return MYB_OK;
}

myb_status myb_kpool_string(const myb_kpool* pool, const char* name, size_t index,
                            const char** value, size_t* length) {
  MYB_CHECK(value != NULL && length != NULL);
  const kvar* v = lookup(pool, name);
  if (v == NULL) {
    return MYB_E_NOT_FOUND;
  }
  if (v->type != MYB_KVAR_STRINGS) {
    return MYB_E_FORMAT;
  }
  if (index >= v->count) {
    return MYB_E_RANGE;
  }
  const kstring str = pool->data.strings[v->first + index];
  *value = pool->data.chars + str.offset;
  *length = str.length;
  return MYB_OK;
}
