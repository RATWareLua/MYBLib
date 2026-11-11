/* SPDX-License-Identifier: BSD-3-Clause */
/* Orbit data messages of the CCSDS in the keyword = value notation: OEM, OPM, OMM
 * (CCSDS 502.0-B-2, 2009, and 502.0-B-3, 2023; chapters 3-5 and 7). */
#include <math.h>
#include <string.h>

#include "core/alloc.h"
#include "core/base.h"
#include "io/text.h"
#include "myb/ccsds.h"
#include "myb/io.h"
#include "platform/platform.h"

static const double metres_per_km = 1000.0;
static const double gm_unit = 1.0e9;                 /* m^3/s^2 per km^3/s^2 */
static const double degree = 0.017453292519943295;   /* rad */
static const size_t max_file = (size_t)1 << 31;

/* ---- Lines ----------------------------------------------------------------------------------- */

typedef enum line_kind {
  LINE_END = 0,  /* no more lines */
  LINE_PAIR = 1, /* KEY = value */
  LINE_MARK = 2, /* a keyword alone: META_START, ... */
  LINE_DATA = 3  /* anything else that is not empty or a comment */
} line_kind;

typedef struct line {
  line_kind kind;
  size_t key; /* offsets into the text */
  size_t key_length;
  size_t value;
  size_t value_length;
  long number;
  const char* alias; /* the keyword of a mark that the text spells otherwise (XML) */
  int raw;           /* XML: the value is a CDATA section, to be taken as written */
} line;

typedef struct cursor {
  const char* text;
  size_t size;
  size_t at;
  long line;
  int xml;    /* the XML notation */
  int flat;   /* XML: the sections are passed over and their keywords come one by one (the
                 parameter messages, which have no sections in the other notation) */
  int rooted; /* XML: the root element was read */
} cursor;

static int blank(char c) {
  return c == ' ' || c == '\t' || c == '\r';
}

static int keyword_char(char c) {
  return (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
}

/* ---- The XML notation ---------------------------------------------------------------------------
 * The same messages as elements ("XML Specification for Navigation Data Messages", CCSDS
 * 505.0-B): the root names the message and its version in the attributes id and version; a
 * keyword is an element with its value as text; the sections are elements in lower case.
 * The scanner gives them as the lines of the other notation — nothing is copied —, so that one
 * reader serves both. A value may be given as one CDATA section (taken as written: no entities
 * in it); a document type declaration is passed over with its internal subset, whose entities
 * are not known to the reader (&name; stays as written). */

static int xml_space(char c) {
  return blank(c) || c == '\n';
}

/* A message in XML starts with '<' (after a byte order mark and blanks). */
static int is_xml(const char* text, size_t size) {
  size_t at = 0;
  if (size >= 3 && (unsigned char)text[0] == 0xEF && (unsigned char)text[1] == 0xBB &&
      (unsigned char)text[2] == 0xBF) {
    at = 3;
  }
  /* bound: the size of the text */
  while (at < size && xml_space(text[at])) {
    ++at;
  }
  return at < size && text[at] == '<';
}

/* Past the next `end`, or to the end of the text. */
static void xml_skip(cursor* c, const char* end) {
  const size_t n = strlen(end);
  /* bound: the size of the text */
  while (c->at < c->size) {
    if (c->size - c->at >= n && memcmp(c->text + c->at, end, n) == 0) {
      c->at += n;
      return;
    }
    c->line += c->text[c->at] == '\n';
    ++c->at;
  }
}

static int xml_at(const cursor* c, size_t at, const char* what) {
  const size_t n = strlen(what);
  return at <= c->size && c->size - at >= n && memcmp(c->text + at, what, n) == 0;
}

/* Past a document type declaration: its '>' outside the quotes and outside the internal subset
 * [ ... ], where comments and processing instructions may hold any character. */
static void xml_skip_doctype(cursor* c) {
  int subset = 0;
  char quote = '\0';
  /* bound: the size of the text; every pass moves on */
  while (c->at < c->size) {
    const char ch = c->text[c->at];
    if (quote != '\0') {
      quote = ch == quote ? '\0' : quote;
    } else if (subset && xml_at(c, c->at, "<!--")) {
      xml_skip(c, "-->");
      continue;
    } else if (subset && xml_at(c, c->at, "<?")) {
      xml_skip(c, "?>");
      continue;
    } else if (ch == '"' || ch == '\'') {
      quote = ch;
    } else if (ch == '[') {
      subset = 1;
    } else if (ch == ']') {
      subset = 0;
    } else if (ch == '>' && !subset) {
      ++c->at;
      return;
    }
    c->line += ch == '\n';
    ++c->at;
  }
}

/* A value that is one CDATA section, [*from, *to) trimmed: narrowed to what the section holds.
 * 1 if it was one; -1 if a section is mixed with other characters (no value of a message is). */
static int xml_cdata(const cursor* c, size_t* from, size_t* to) {
  static const char open[] = "<![CDATA[";
  static const char close[] = "]]>";
  const char* s = c->text;
  size_t first = *to; /* of the first section in the value */
  /* bound: the characters of the value */
  for (size_t i = *from; i + sizeof open - 1 <= *to && first == *to; ++i) {
    if (memcmp(s + i, open, sizeof open - 1) == 0) {
      first = i;
    }
  }
  if (first == *to) {
    return 0;
  }
  if (first != *from || *to - *from < sizeof open - 1 + sizeof close - 1 ||
      memcmp(s + *to - (sizeof close - 1), close, sizeof close - 1) != 0) {
    return -1;
  }
  const size_t inner = *from + sizeof open - 1;
  const size_t inner_end = *to - (sizeof close - 1);
  /* bound: the characters of the section */
  for (size_t i = inner; i + sizeof close - 1 <= inner_end; ++i) {
    if (memcmp(s + i, close, sizeof close - 1) == 0) {
      return -1; /* it ended before: more follows it */
    }
  }
  *from = inner;
  *to = inner_end;
  return 1;
}

/* The value of the attribute `name` among the characters [from, to) of a tag. */
static int xml_attribute(const char* s, size_t from, size_t to, const char* name, size_t* value,
                         size_t* length) {
  const size_t n = strlen(name);
  /* bound: the characters of the tag */
  for (size_t i = from; i + n < to; ++i) {
    if (!xml_space(s[i]) || memcmp(s + i + 1, name, n) != 0) {
      continue;
    }
    size_t k = i + 1 + n;
    while (k < to && xml_space(s[k])) {
      ++k;
    }
    if (k >= to || s[k] != '=') {
      continue;
    }
    ++k;
    while (k < to && xml_space(s[k])) {
      ++k;
    }
    if (k >= to || (s[k] != '"' && s[k] != '\'')) {
      continue;
    }
    const char quote = s[k];
    size_t end = k + 1;
    while (end < to && s[end] != quote) {
      ++end;
    }
    if (end >= to) {
      return 0;
    }
    *value = k + 1;
    *length = end - (k + 1);
    return 1;
  }
  return 0;
}

/* The sections that the other notation marks, and the elements that hold a line of data. */
typedef struct xml_section {
  const char* element;
  const char* start; /* NULL: a line of data */
  const char* stop;
} xml_section;

static const xml_section xml_sections[] = {
    {"metadata", "META_START", "META_STOP"},
    {"covarianceMatrix", "COVARIANCE_START", "COVARIANCE_STOP"},
    {"stateVector", NULL, NULL},
    {"traj", "TRAJ_START", "TRAJ_STOP"},
    {"phys", "PHYS_START", "PHYS_STOP"},
    {"cov", "COV_START", "COV_STOP"},
    {"man", "MAN_START", "MAN_STOP"},
    {"pert", "PERT_START", "PERT_STOP"},
    {"od", "OD_START", "OD_STOP"},
    {"user", "USER_START", "USER_STOP"},
    {"trajLine", NULL, NULL},
    {"covLine", NULL, NULL},
    {"manLine", NULL, NULL}};

static void trim(const char* s, size_t* from, size_t* to) {
  /* bound: the characters between the two */
  while (*from < *to && xml_space(s[*from])) {
    ++*from;
  }
  while (*to > *from && xml_space(s[*to - 1])) {
    --*to;
  }
}

static void next_xml(cursor* c, line* out) {
  const char* s = c->text;
  /* bound: the size of the text; every pass moves on */
  while (c->at < c->size) {
    if (s[c->at] != '<') {
      c->line += s[c->at] == '\n';
      ++c->at;
      continue;
    }
    const size_t tag = c->at;
    if (c->size - tag >= 4 && memcmp(s + tag, "<!--", 4) == 0) {
      xml_skip(c, "-->");
      continue;
    }
    if (tag + 1 < c->size && s[tag + 1] == '?') {
      xml_skip(c, "?>");
      continue;
    }
    if (tag + 1 < c->size && s[tag + 1] == '!') {
      if (xml_at(c, tag, "<![CDATA[")) {
        xml_skip(c, "]]>"); /* character data between the elements: nothing reads it */
      } else if (xml_at(c, tag, "<!DOCTYPE")) {
        xml_skip_doctype(c);
      } else {
        xml_skip(c, ">");
      }
      continue;
    }
    const int closing = tag + 1 < c->size && s[tag + 1] == '/';
    const size_t name = tag + 1 + (closing ? 1u : 0u);
    size_t name_end = name;
    while (name_end < c->size && !xml_space(s[name_end]) && s[name_end] != '>' &&
           s[name_end] != '/') {
      ++name_end;
    }
    /* the end of the tag: a '>' outside the quotes of an attribute */
    size_t close = name_end;
    char quote = '\0';
    while (close < c->size && (quote != '\0' || s[close] != '>')) {
      if (quote != '\0') {
        quote = s[close] == quote ? '\0' : quote;
      } else if (s[close] == '"' || s[close] == '\'') {
        quote = s[close];
      }
      c->line += s[close] == '\n';
      ++close;
    }
    out->number = c->line + 1;
    out->alias = NULL;
    out->raw = 0;
    if (close >= c->size) { /* a tag that does not end: no keyword reads it */
      c->at = c->size;
      out->kind = LINE_DATA;
      out->key = tag;
      out->key_length = 0;
      out->value = tag;
      out->value_length = c->size - tag;
      return;
    }
    const int empty = s[close - 1] == '/';
    c->at = close + 1;
    const size_t n = name_end - name;
    if (!c->rooted) {
      /* the root: <oem id="CCSDS_OEM_VERS" version="2.0"> is the first keyword */
      c->rooted = 1;
      size_t id = 0;
      size_t id_length = 0;
      size_t version = 0;
      size_t version_length = 0;
      if (closing || !xml_attribute(s, name_end, close, "id", &id, &id_length) ||
          !xml_attribute(s, name_end, close, "version", &version, &version_length)) {
        out->kind = LINE_DATA;
        out->key = tag;
        out->key_length = 0;
        out->value = tag;
        out->value_length = close - tag;
        return;
      }
      out->kind = LINE_PAIR;
      out->key = id;
      out->key_length = id_length;
      out->value = version;
      out->value_length = version_length;
      return;
    }
    int keyword = n > 0;
    for (size_t i = 0; i < n; ++i) {
      keyword = keyword && keyword_char(s[name + i]);
    }
    if (keyword) {
      if (closing) {
        continue;
      }
      /* <KEYWORD>value</KEYWORD>: the value, and past the closing tag */
      size_t from = c->at;
      size_t to = from;
      int section = 0;
      if (empty) {
        from = to = close;
      } else {
        /* bound: the size of the text */
        while (to < c->size && (s[to] != '<' || xml_at(c, to, "<![CDATA["))) {
          if (s[to] == '<') { /* a section: to its end, whatever it holds */
            c->at = to;
            xml_skip(c, "]]>");
            to = c->at;
            continue;
          }
          c->line += s[to] == '\n';
          ++to;
        }
        c->at = to;
        xml_skip(c, ">");
        trim(s, &from, &to);
        section = xml_cdata(c, &from, &to);
      }
      trim(s, &from, &to);
      if (n == 7 && memcmp(s + name, "COMMENT", 7) == 0) {
        continue;
      }
      if (section < 0) { /* a section among other characters: no keyword reads it */
        out->kind = LINE_DATA;
        out->key = tag;
        out->key_length = 0;
        out->value = tag;
        out->value_length = close - tag;
        return;
      }
      out->raw = section;
      out->kind = LINE_PAIR;
      out->key = name;
      out->key_length = n;
      out->value = from;
      out->value_length = to - from;
      /* <USER_DEFINED parameter="X">: the keyword is the parameter */
      if (n == 12 && memcmp(s + name, "USER_DEFINED", 12) == 0) {
        (void)xml_attribute(s, name_end, close, "parameter", &out->key, &out->key_length);
      }
      return;
    }
    const xml_section* section = NULL;
    for (size_t k = 0; k < sizeof xml_sections / sizeof xml_sections[0]; ++k) {
      if (strlen(xml_sections[k].element) == n && memcmp(s + name, xml_sections[k].element, n) == 0) {
        section = &xml_sections[k];
      }
    }
    if (section == NULL || c->flat) {
      continue; /* header, body, segment, data and what later versions add */
    }
    if (section->start != NULL) {
      out->kind = LINE_MARK;
      out->alias = closing ? section->stop : section->start;
      out->key = name;
      out->key_length = n;
      out->value = close;
      out->value_length = 0;
      return;
    }
    if (closing || empty) {
      continue;
    }
    /* a line of data: what the element holds, up to its closing tag */
    size_t from = c->at;
    size_t to = from;
    /* bound: the size of the text */
    while (to < c->size) {
      if (s[to] == '<' && to + 1 < c->size && s[to + 1] == '/' && c->size - (to + 2) >= n &&
          memcmp(s + to + 2, s + name, n) == 0) {
        break;
      }
      c->line += s[to] == '\n';
      ++to;
    }
    c->at = to;
    xml_skip(c, ">");
    trim(s, &from, &to);
    /* the whole line as one section: what it holds (sections inside the elements of a line
     * are for next_word) */
    if (xml_cdata(c, &from, &to) > 0) {
      out->raw = 1;
      trim(s, &from, &to);
    }
    out->kind = LINE_DATA;
    out->key = from;
    out->key_length = 0;
    out->value = from;
    out->value_length = to - from;
    return;
  }
  out->kind = LINE_END;
}

/* &lt; &gt; &amp; &quot; &apos; in a value, in place: the new length. */
static size_t xml_decode(char* s, size_t n) {
  static const struct {
    const char* name;
    char c;
  } entities[5] = {{"&lt;", '<'}, {"&gt;", '>'}, {"&amp;", '&'}, {"&quot;", '"'}, {"&apos;", '\''}};
  size_t out = 0;
  /* bound: the characters of the value */
  for (size_t i = 0; i < n;) {
    size_t used = 0;
    for (int k = 0; k < 5 && s[i] == '&' && used == 0; ++k) {
      const size_t length = strlen(entities[k].name);
      if (n - i >= length && memcmp(s + i, entities[k].name, length) == 0) {
        s[out++] = entities[k].c;
        used = length;
      }
    }
    if (used == 0) {
      s[out++] = s[i];
      used = 1;
    }
    i += used;
  }
  return out;
}

/* The next line that carries something; comments and empty lines are skipped. */
static void next_line(cursor* c, line* out) {
  out->kind = LINE_END;
  out->alias = NULL;
  out->raw = 0;
  if (c->xml) {
    next_xml(c, out);
    return;
  }
  /* bound: the size of the text */
  while (c->at < c->size) {
    size_t start = c->at;
    size_t end = start;
    while (end < c->size && c->text[end] != '\n') {
      ++end;
    }
    c->at = end < c->size ? end + 1 : end;
    c->line++;
    while (start < end && blank(c->text[start])) {
      ++start;
    }
    while (end > start && blank(c->text[end - 1])) {
      --end;
    }
    if (start == end) {
      continue;
    }
    const char* s = c->text + start;
    const size_t n = end - start;
    if (n >= 7 && memcmp(s, "COMMENT", 7) == 0 && (n == 7 || blank(s[7]))) {
      continue;
    }
    out->number = c->line;
    size_t k = 0;
    while (k < n && keyword_char(s[k])) {
      ++k;
    }
    size_t after = k;
    while (after < n && blank(s[after])) {
      ++after;
    }
    out->key = start;
    out->key_length = k;
    if (k > 0 && after < n && s[after] == '=') {
      size_t v = after + 1;
      while (v < n && blank(s[v])) {
        ++v;
      }
      out->kind = LINE_PAIR;
      out->value = start + v;
      out->value_length = n - v;
    } else if (k == n) {
      out->kind = LINE_MARK;
      out->value = end;
      out->value_length = 0;
    } else {
      out->kind = LINE_DATA;
      out->key_length = 0;
      out->value = start;
      out->value_length = n;
    }
    return;
  }
}

static int key_is(const cursor* c, const line* l, const char* name) {
  if (l->alias != NULL) {
    return strcmp(l->alias, name) == 0;
  }
  const size_t n = strlen(name);
  return l->key_length == n && memcmp(c->text + l->key, name, n) == 0;
}

/* The next word of a line of data: blanks, line ends and the tags of the XML notation part
 * the words, and so do the two ends of a CDATA section. */
static int next_word(const char* s, size_t n, size_t* at, size_t* start, size_t* length) {
  size_t i = *at;
  /* bound: the characters of the line */
  while (i < n) {
    if (xml_space(s[i])) {
      ++i;
    } else if (n - i >= 9 && memcmp(s + i, "<![CDATA[", 9) == 0) {
      i += 9;
    } else if (n - i >= 3 && memcmp(s + i, "]]>", 3) == 0) {
      i += 3;
    } else if (s[i] == '<') {
      while (i < n && s[i] != '>') {
        ++i;
      }
      i += i < n;
    } else {
      break;
    }
  }
  size_t end = i;
  while (end < n && !xml_space(s[end]) && s[end] != '<' &&
         !(n - end >= 3 && memcmp(s + end, "]]>", 3) == 0)) {
    ++end;
  }
  *start = i;
  *length = end - i;
  *at = end;
  return end > i;
}

static myb_status fail(myb_diag* diag, long number, const char* what) {
  myb_internal_diag(diag, number, "%s", what);
  return MYB_E_FORMAT;
}

/* ---- Values ---------------------------------------------------------------------------------- */

/* A number, with an optional unit in brackets after it. */
static myb_status number_of(const char* s, size_t n, double* value) {
  if (n > 0 && s[n - 1] == ']') {
    size_t open = n;
    /* bound: the characters of the value */
    while (open > 0 && s[open - 1] != '[') {
      --open;
    }
    if (open == 0) {
      return MYB_E_FORMAT;
    }
    n = open - 1;
    while (n > 0 && blank(s[n - 1])) {
      --n;
    }
  }
  size_t used = 0;
  double v = 0.0;
  if (n == 0 || myb_parse_double(s, n, &v, &used) != MYB_OK || used != n || !isfinite(v)) {
    return MYB_E_FORMAT;
  }
  *value = v;
  return MYB_OK;
}

static myb_status integer_of(const char* s, size_t n, long* value) {
  double v = 0.0;
  if (number_of(s, n, &v) != MYB_OK || v != floor(v) || fabs(v) > 2.0e9) {
    return MYB_E_FORMAT;
  }
  *value = (long)v;
  return MYB_OK;
}

/* `count` decimal digits at s. */
static int digits_of(const char* s, size_t count, int* value) {
  int v = 0;
  for (size_t i = 0; i < count; ++i) {
    if (s[i] < '0' || s[i] > '9') {
      return 0;
    }
    v = 10 * v + (s[i] - '0');
  }
  *value = v;
  return 1;
}

static int leap_year(int64_t year) {
  return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

static int month_days(int64_t year, int month) {
  static const int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  return month == 2 && leap_year(year) ? 29 : days[month - 1];
}

/* YYYY-MM-DDThh:mm:ss[.d...][Z] or YYYY-DDDThh:mm:ss[.d...][Z]. */
static myb_status calendar_of(const char* s, size_t n, myb_calendar* out) {
  if (n > 0 && s[n - 1] == 'Z') {
    --n;
  }
  int year = 0;
  if (n < 17 || !digits_of(s, 4, &year) || s[4] != '-') {
    return MYB_E_FORMAT;
  }
  size_t at = 5;
  int month = 0;
  int day = 0;
  if (n > 8 && s[8] == 'T') { /* day of the year */
    int doy = 0;
    if (!digits_of(s + 5, 3, &doy) || doy < 1 || doy > (leap_year(year) ? 366 : 365)) {
      return MYB_E_FORMAT;
    }
    month = 1;
    /* bound: twelve months */
    while (month < 12 && doy > month_days(year, month)) {
      doy -= month_days(year, month);
      ++month;
    }
    day = doy;
    at = 8;
  } else {
    if (n < 19 || !digits_of(s + 5, 2, &month) || s[7] != '-' || !digits_of(s + 8, 2, &day) ||
        month < 1 || month > 12 || day < 1 || day > month_days(year, month)) {
      return MYB_E_FORMAT;
    }
    at = 10;
  }
  int hour = 0;
  int minute = 0;
  if (at + 9 > n || s[at] != 'T' || !digits_of(s + at + 1, 2, &hour) || s[at + 3] != ':' ||
      !digits_of(s + at + 4, 2, &minute) || s[at + 6] != ':' || hour > 23 || minute > 59) {
    return MYB_E_FORMAT;
  }
  const char* seconds = s + at + 7;
  const size_t left = n - (at + 7);
  int whole = 0;
  double second = 0.0;
  size_t used = 0;
  if (left < 2 || !digits_of(seconds, 2, &whole) || (left > 2 && seconds[2] != '.') ||
      myb_parse_double(seconds, left, &second, &used) != MYB_OK || used != left ||
      !(second >= 0.0) || !(second < 61.0)) {
    return MYB_E_FORMAT;
  }
  out->year = year;
  out->month = month;
  out->day = day;
  out->hour = hour;
  out->minute = minute;
  out->second = second;
  return MYB_OK;
}

static myb_status text_of(const char* s, size_t n, char out[MYB_ODM_TEXT]) {
  if (n == 0 || n >= MYB_ODM_TEXT) {
    return MYB_E_FORMAT;
  }
  memcpy(out, s, n);
  out[n] = '\0';
  return MYB_OK;
}

/* ---- Writing --------------------------------------------------------------------------------- */

static int valid_calendar(const myb_calendar* c) {
  return c->year >= 0 && c->year <= 9999 && c->month >= 1 && c->month <= 12 && c->day >= 1 &&
         c->day <= month_days(c->year, c->month) && c->hour >= 0 && c->hour <= 23 &&
         c->minute >= 0 && c->minute <= 59 && c->second >= 0.0 && c->second < 61.0;
}

static int valid_text(const char* s) {
  if (s == NULL || s[0] == '\0') {
    return 0;
  }
  for (size_t i = 0; s[i] != '\0'; ++i) {
    if (s[i] == '\n' || s[i] == '\r') {
      return 0;
    }
  }
  return 1;
}

/* YYYY-MM-DDThh:mm:ss.ddd: the seconds to 1e-15 s, trailing zeros beyond the third dropped. */
static void put_calendar(myb_internal_text* t, const myb_calendar* c) {
  myb_internal_text_int(t, c->year, 4, '0');
  myb_internal_text_char(t, '-');
  myb_internal_text_int(t, c->month, 2, '0');
  myb_internal_text_char(t, '-');
  myb_internal_text_int(t, c->day, 2, '0');
  myb_internal_text_char(t, 'T');
  myb_internal_text_int(t, c->hour, 2, '0');
  myb_internal_text_char(t, ':');
  myb_internal_text_int(t, c->minute, 2, '0');
  myb_internal_text_char(t, ':');
  /* the seconds: the shortest decimal that reads back to the same number, written without
   * an exponent, with three decimals at least */
  char shortest[MYB_DOUBLE_CHARS];
  const size_t n = myb_format_double(c->second, shortest);
  size_t i = 0;
  while (i < n && shortest[i] != '.' && shortest[i] != 'e') {
    ++i;
  }
  int decimals = 0;
  if (i < n && shortest[i] == '.') {
    ++i;
    while (i < n && shortest[i] >= '0' && shortest[i] <= '9') {
      ++decimals;
      ++i;
    }
  }
  if (i + 1 < n && shortest[i] == 'e' && shortest[i + 1] == '-') { /* below 1e-4 s */
    int exponent = 0;
    for (size_t k = i + 2; k < n; ++k) {
      exponent = 10 * exponent + (shortest[k] - '0');
    }
    decimals += exponent;
  }
  decimals = decimals < 3 ? 3 : decimals > 30 ? 30 : decimals;
  char seconds[48];
  size_t length = 0;
  const myb_status status = myb_format_fixed(c->second, decimals, seconds, sizeof seconds, &length);
  MYB_CHECK(status == MYB_OK);
  if (length >= 2 && seconds[1] == '.') {
    myb_internal_text_char(t, '0');
  }
  myb_internal_text_span(t, seconds, length);
}

/* value / scale, the number a message carries in its own unit (km for m, degrees for rad), as
 * the shortest decimal that the reader, multiplying by scale, turns back into value. The
 * quotient or one of its two neighbours does that for all but a few values in a hundred, which
 * come back one unit of the last place away: of those that do, the one with the fewest digits
 * is written, so that 7000 km read from a message is written as 7000, not 7000.000000000001. */
static void put_scaled(myb_internal_text* t, double value, double scale) {
  const double quotient = value / scale;
  const double candidates[3] = {quotient, nextafter(quotient, -(double)INFINITY),
                                nextafter(quotient, (double)INFINITY)};
  double best = quotient;
  size_t best_length = 0;
  for (int k = 0; k < 3; ++k) {
    if (candidates[k] * scale == value) {
      char text[MYB_DOUBLE_CHARS];
      const size_t n = myb_format_double(candidates[k], text);
      if (best_length == 0 || n < best_length) {
        best = candidates[k];
        best_length = n;
      }
    }
  }
  myb_internal_text_double(t, best);
}

/* What a writer writes into, and in which notation. */
typedef struct sink {
  myb_internal_text t;
  int xml;
} sink;

static void put_escaped(myb_internal_text* t, const char* s) {
  for (size_t i = 0; s[i] != '\0'; ++i) {
    if (s[i] == '&') {
      myb_internal_text_put(t, "&amp;");
    } else if (s[i] == '<') {
      myb_internal_text_put(t, "&lt;");
    } else if (s[i] == '>') {
      myb_internal_text_put(t, "&gt;");
    } else {
      myb_internal_text_char(t, s[i]);
    }
  }
}

static void put_key(sink* k, const char* key) {
  if (k->xml) {
    myb_internal_text_char(&k->t, '<');
    myb_internal_text_put(&k->t, key);
    myb_internal_text_char(&k->t, '>');
  } else {
    myb_internal_text_put(&k->t, key);
    myb_internal_text_put(&k->t, " = ");
  }
}

static void put_end(sink* k, const char* key) {
  if (k->xml) {
    myb_internal_text_put(&k->t, "</");
    myb_internal_text_put(&k->t, key);
    myb_internal_text_char(&k->t, '>');
  }
  myb_internal_text_char(&k->t, '\n');
}

static void put_text(sink* k, const char* key, const char* value) {
  put_key(k, key);
  if (k->xml) {
    put_escaped(&k->t, value);
  } else {
    myb_internal_text_put(&k->t, value);
  }
  put_end(k, key);
}

static void put_number(sink* k, const char* key, double value, double scale) {
  put_key(k, key);
  put_scaled(&k->t, value, scale);
  put_end(k, key);
}

static void put_integer(sink* k, const char* key, long value) {
  put_key(k, key);
  myb_internal_text_int(&k->t, value, 0, ' ');
  put_end(k, key);
}

static void put_time(sink* k, const char* key, const myb_calendar* c) {
  put_key(k, key);
  put_calendar(&k->t, c);
  put_end(k, key);
}

/* The start and the end of a section: an element in XML, the line `mark` (NULL: nothing) in
 * the other notation. */
static void section(sink* k, const char* element, const char* mark) {
  if (k->xml) {
    myb_internal_text_char(&k->t, '<');
    myb_internal_text_put(&k->t, element);
    myb_internal_text_put(&k->t, ">\n");
  } else if (mark != NULL) {
    myb_internal_text_put(&k->t, mark);
    myb_internal_text_char(&k->t, '\n');
  }
}

static void section_end(sink* k, const char* element, const char* mark) {
  if (k->xml) {
    myb_internal_text_put(&k->t, "</");
    myb_internal_text_put(&k->t, element);
    myb_internal_text_put(&k->t, ">\n");
  } else if (mark != NULL) {
    myb_internal_text_put(&k->t, mark);
    myb_internal_text_char(&k->t, '\n');
  }
}

/* The first line of a message: its kind and version. */
static void message(sink* k, const char* element, const char* key, const char* version) {
  if (!k->xml) {
    put_text(k, key, version);
    return;
  }
  myb_internal_text_put(&k->t, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<");
  myb_internal_text_put(&k->t, element);
  myb_internal_text_put(&k->t,
                        " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
                        " xsi:noNamespaceSchemaLocation="
                        "\"http://sanaregistry.org/r/ndmxml/ndmxml-1.0-master.xsd\" id=\"");
  myb_internal_text_put(&k->t, key);
  myb_internal_text_put(&k->t, "\" version=\"");
  myb_internal_text_put(&k->t, version);
  myb_internal_text_put(&k->t, "\">\n");
}

static myb_status message_end(sink* k, const char* element, size_t* length) {
  if (k->xml) {
    myb_internal_text_put(&k->t, "</");
    myb_internal_text_put(&k->t, element);
    myb_internal_text_put(&k->t, ">\n");
  }
  return myb_internal_text_finish(&k->t, length);
}

/* ---- Covariances ------------------------------------------------------------------------------
 * The lower triangle of the 6 x 6 matrix of a state, row by row: 21 numbers, in the file in
 * km^2, km^2/s and km^2/s^2. */

static const double covariance_unit = 1.0e6; /* m^2 per km^2 */

static const char* const covariance_keys[21] = {
    "CX_X",         "CY_X",         "CY_Y",         "CZ_X",         "CZ_Y",         "CZ_Z",
    "CX_DOT_X",     "CX_DOT_Y",     "CX_DOT_Z",     "CX_DOT_X_DOT", "CY_DOT_X",     "CY_DOT_Y",
    "CY_DOT_Z",     "CY_DOT_X_DOT", "CY_DOT_Y_DOT", "CZ_DOT_X",     "CZ_DOT_Y",     "CZ_DOT_Z",
    "CZ_DOT_X_DOT", "CZ_DOT_Y_DOT", "CZ_DOT_Z_DOT"};

/* The row and the column of the k-th number of the triangle. */
static void triangle(int k, int* row, int* column) {
  int i = 0;
  /* bound: six rows */
  while ((i + 1) * (i + 2) / 2 <= k) {
    ++i;
  }
  *row = i;
  *column = k - i * (i + 1) / 2;
}

static void covariance_set(double covariance[36], int k, double value) {
  int row = 0;
  int column = 0;
  triangle(k, &row, &column);
  covariance[6 * row + column] = value;
  covariance[6 * column + row] = value;
}

static int covariance_finite(const double covariance[36]) {
  for (int k = 0; k < 21; ++k) {
    int row = 0;
    int column = 0;
    triangle(k, &row, &column);
    if (!isfinite(covariance[6 * row + column])) {
      return 0;
    }
  }
  return 1;
}

/* The 21 keywords of a parameter message and of the XML notation. */
static void put_covariance(sink* k, const double covariance[36], double unit) {
  for (int i = 0; i < 21; ++i) {
    int row = 0;
    int column = 0;
    triangle(i, &row, &column);
    put_number(k, covariance_keys[i], covariance[6 * row + column], unit);
  }
}

/* ---- OEM ------------------------------------------------------------------------------------- */

typedef struct pair {
  const char* key;
  const char* value;
} pair;

typedef struct oem_segment {
  size_t pair_first;
  size_t pair_count;
  size_t state_first;
  size_t state_count;
  size_t covariance_first;
  size_t covariance_count;
} oem_segment;

struct myb_oem {
  myb_allocator allocator;
  char* text;
  pair* pairs; /* the header's first, then each segment's */
  size_t header_pairs;
  size_t pair_count;
  oem_segment* segments;
  size_t segment_count;
  myb_calendar* epochs;
  double* values; /* 9 per state: position, velocity, acceleration (SI) */
  unsigned char* accelerated;
  size_t state_count;
  myb_calendar* covariance_epochs;
  const char** covariance_frames; /* NULL where the message names none */
  double* covariances;            /* 36 each (SI) */
  size_t covariance_count;
};

void myb_oem_destroy(myb_oem* oem) {
  if (oem == NULL) {
    return;
  }
  const myb_allocator a = oem->allocator;
  myb_internal_release(&a, oem->text);
  myb_internal_release(&a, oem->pairs);
  myb_internal_release(&a, oem->segments);
  myb_internal_release(&a, oem->epochs);
  myb_internal_release(&a, oem->values);
  myb_internal_release(&a, oem->accelerated);
  myb_internal_release(&a, oem->covariance_epochs);
  myb_internal_release(&a, oem->covariance_frames);
  myb_internal_release(&a, oem->covariances);
  myb_internal_release(&a, oem);
}

/* One line of an ephemeris: the epoch and 6 or 9 numbers. */
static myb_status data_line(const char* s, size_t n, myb_calendar* epoch, double values[9],
                            int* count) {
  size_t at = 0;
  size_t start = 0;
  size_t length = 0;
  if (!next_word(s, n, &at, &start, &length) || calendar_of(s + start, length, epoch) != MYB_OK) {
    return MYB_E_FORMAT;
  }
  *count = 0;
  /* bound: the words of the line */
  while (next_word(s, n, &at, &start, &length)) {
    if (*count == 9 || number_of(s + start, length, &values[*count]) != MYB_OK) {
      return MYB_E_FORMAT;
    }
    ++*count;
  }
  return *count == 6 || *count == 9 ? MYB_OK : MYB_E_FORMAT;
}

typedef enum oem_mode { OEM_HEADER, OEM_META, OEM_DATA, OEM_COVARIANCE } oem_mode;

/* One pass over the message; with fill == 0 it only counts (oem->pair_count, ...). */
static myb_status oem_pass(myb_oem* oem, size_t size, int fill, myb_diag* diag) {
  cursor c = {.text = oem->text, .size = size, .xml = is_xml(oem->text, size)};
  oem_mode mode = OEM_HEADER;
  size_t pairs = 0;
  size_t segments = 0;
  size_t states = 0;
  size_t covariances = 0;
  int filled = 21; /* numbers of the covariance matrix being read; 21: none is open */
  line l;
  /* bound: the lines of the text */
  for (next_line(&c, &l); l.kind != LINE_END; next_line(&c, &l)) {
    if (l.kind == LINE_MARK && mode != OEM_COVARIANCE) {
      if (key_is(&c, &l, "META_START") && (mode == OEM_HEADER || mode == OEM_DATA)) {
        if (fill) {
          oem->segments[segments] = (oem_segment){
              .pair_first = pairs, .state_first = states, .covariance_first = covariances};
        }
        ++segments;
        mode = OEM_META;
      } else if (key_is(&c, &l, "META_STOP") && mode == OEM_META) {
        mode = OEM_DATA;
      } else if (key_is(&c, &l, "COVARIANCE_START") && mode == OEM_DATA) {
        mode = OEM_COVARIANCE;
      } else {
        return fail(diag, l.number, "OEM: a keyword out of place");
      }
      continue;
    }
    if (mode == OEM_COVARIANCE) {
      if (l.kind == LINE_MARK && key_is(&c, &l, "COVARIANCE_STOP")) {
        if (filled != 21) {
          return fail(diag, l.number, "OEM: a covariance matrix does not have 21 numbers");
        }
        mode = OEM_DATA;
        continue;
      }
      if (l.kind == LINE_PAIR && key_is(&c, &l, "EPOCH")) {
        if (filled != 21) {
          return fail(diag, l.number, "OEM: a covariance matrix does not have 21 numbers");
        }
        myb_calendar epoch;
        if (calendar_of(oem->text + l.value, l.value_length, &epoch) != MYB_OK) {
          return fail(diag, l.number, "OEM: not a date of the form YYYY-MM-DDThh:mm:ss");
        }
        if (fill) {
          oem->covariance_epochs[covariances] = epoch;
          oem->covariance_frames[covariances] = NULL;
          oem->segments[segments - 1].covariance_count++;
        }
        ++covariances;
        filled = 0;
        continue;
      }
      if (filled == 21) {
        return fail(diag, l.number, "OEM: a covariance matrix starts with its EPOCH");
      }
      if (l.kind == LINE_PAIR && key_is(&c, &l, "COV_REF_FRAME")) {
        if (fill) {
          oem->text[l.value + l.value_length] = '\0';
          oem->covariance_frames[covariances - 1] = oem->text + l.value;
        }
        continue;
      }
      /* numbers: the rows of the triangle, or one keyword each in XML; a row of one number
       * that looks like a keyword ("1") is a number all the same */
      const size_t from = l.kind == LINE_MARK ? l.key : l.value;
      const size_t extent = l.kind == LINE_MARK ? l.key_length : l.value_length;
      size_t at = 0;
      size_t start = 0;
      size_t length = 0;
      /* bound: the words of the line */
      while (next_word(oem->text + from, extent, &at, &start, &length)) {
        double v = 0.0;
        if (filled == 21 || number_of(oem->text + from + start, length, &v) != MYB_OK) {
          return fail(diag, l.number, "OEM: not the numbers of a covariance matrix");
        }
        if (fill) {
          covariance_set(oem->covariances + 36 * (covariances - 1), filled, v * covariance_unit);
        }
        ++filled;
      }
      continue;
    }
    if (l.kind == LINE_PAIR) {
      if (mode != OEM_HEADER && mode != OEM_META) {
        return fail(diag, l.number, "OEM: a keyword among the ephemeris lines");
      }
      if (pairs == 0 && !key_is(&c, &l, "CCSDS_OEM_VERS")) {
        return fail(diag, l.number, "OEM: the message does not start with CCSDS_OEM_VERS");
      }
      if (fill) {
        const size_t value_length =
            c.xml && !l.raw ? xml_decode(oem->text + l.value, l.value_length) : l.value_length;
        oem->text[l.key + l.key_length] = '\0';
        oem->text[l.value + value_length] = '\0';
        oem->pairs[pairs] = (pair){oem->text + l.key, oem->text + l.value};
        if (mode == OEM_META) {
          oem->segments[segments - 1].pair_count++;
        }
      }
      ++pairs;
      if (mode == OEM_HEADER) {
        oem->header_pairs = pairs;
      }
      continue;
    }
    if (mode != OEM_DATA) {
      return fail(diag, l.number, "OEM: an ephemeris line outside the data of a segment");
    }
    if (fill) {
      double v[9] = {0.0};
      int count = 0;
      if (data_line(oem->text + l.value, l.value_length, &oem->epochs[states], v, &count) !=
          MYB_OK) {
        return fail(diag, l.number, "OEM: not an epoch and 6 or 9 numbers");
      }
      for (int i = 0; i < 9; ++i) {
        oem->values[9 * states + (size_t)i] = v[i] * metres_per_km;
      }
      oem->accelerated[states] = count == 9;
      oem->segments[segments - 1].state_count++;
    }
    ++states;
  }
  if (mode == OEM_META || mode == OEM_COVARIANCE) {
    return fail(diag, c.line, "OEM: a section is not closed");
  }
  if (segments == 0) {
    return fail(diag, c.line, "OEM: no META_START");
  }
  oem->pair_count = pairs;
  oem->segment_count = segments;
  oem->state_count = states;
  oem->covariance_count = covariances;
  return MYB_OK;
}

myb_status myb_oem_load(const char* text, size_t size, const myb_allocator* allocator,
                        myb_oem** oem, myb_diag* diag) {
  MYB_CHECK(text != NULL && oem != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_oem), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_oem* o = block;
  memset(o, 0, sizeof *o);
  o->allocator = a;
  status = myb_internal_alloc_array(&a, size + 1, 1, &block);
  if (status != MYB_OK) {
    goto fail;
  }
  o->text = block;
  memcpy(o->text, text, size);
  o->text[size] = '\0';
  status = oem_pass(o, size, 0, diag);
  if (status != MYB_OK) {
    goto fail;
  }
  status = myb_internal_alloc_array(&a, o->pair_count, sizeof(pair), &block);
  if (status == MYB_OK) {
    o->pairs = block;
    status = myb_internal_alloc_array(&a, o->segment_count, sizeof(oem_segment), &block);
  }
  if (status == MYB_OK) {
    o->segments = block;
    status = myb_internal_alloc_array(&a, o->state_count, sizeof(myb_calendar), &block);
  }
  if (status == MYB_OK) {
    o->epochs = block;
    status = myb_internal_alloc_array(&a, 9 * o->state_count, sizeof(double), &block);
  }
  if (status == MYB_OK) {
    o->values = block;
    status = myb_internal_alloc_array(&a, o->state_count, 1, &block);
  }
  if (status == MYB_OK) {
    o->accelerated = block;
    status = myb_internal_alloc_array(&a, o->covariance_count, sizeof(myb_calendar), &block);
  }
  if (status == MYB_OK) {
    o->covariance_epochs = block;
    status = myb_internal_alloc_array(&a, o->covariance_count, sizeof(const char*), &block);
  }
  if (status == MYB_OK) {
    o->covariance_frames = block;
    status = myb_internal_alloc_array(&a, 36 * o->covariance_count, sizeof(double), &block);
  }
  if (status != MYB_OK) {
    goto fail;
  }
  o->covariances = block;
  status = oem_pass(o, size, 1, diag);
  if (status != MYB_OK) {
    goto fail;
  }
  /* what every segment must say about itself */
  static const char* const mandatory[7] = {"OBJECT_NAME", "OBJECT_ID",  "CENTER_NAME", "REF_FRAME",
                                           "TIME_SYSTEM", "START_TIME", "STOP_TIME"};
  for (size_t s = 0; s < o->segment_count; ++s) {
    for (int k = 0; k < 7; ++k) {
      if (myb_oem_meta(o, s, mandatory[k]) == NULL) {
        myb_internal_diag(diag, 0, "OEM: segment %zu has no %s", s + 1, mandatory[k]);
        status = MYB_E_FORMAT;
        goto fail;
      }
    }
  }
  *oem = o;
  return MYB_OK;
fail:
  myb_oem_destroy(o);
  return status;
}

myb_status myb_oem_load_file(const char* path, const myb_allocator* allocator, myb_oem** oem,
                             myb_diag* diag) {
  MYB_CHECK(path != NULL && oem != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  char* data = NULL;
  size_t size = 0;
  const myb_status status = myb_platform_read_file(path, &a, max_file, &data, &size);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot read '%s'", path);
    return status;
  }
  const myb_status loaded = myb_oem_load(data, size, allocator, oem, diag);
  myb_internal_release(&a, data);
  return loaded;
}

static const char* find(const pair* pairs, size_t count, const char* key) {
  for (size_t k = 0; k < count; ++k) {
    if (strcmp(pairs[k].key, key) == 0) {
      return pairs[k].value;
    }
  }
  return NULL;
}

const char* myb_oem_header(const myb_oem* oem, const char* key) {
  MYB_CHECK(oem != NULL && key != NULL);
  return find(oem->pairs, oem->header_pairs, key);
}

size_t myb_oem_segment_count(const myb_oem* oem) {
  MYB_CHECK(oem != NULL);
  return oem->segment_count;
}

const char* myb_oem_meta(const myb_oem* oem, size_t segment, const char* key) {
  MYB_CHECK(oem != NULL && key != NULL && segment < oem->segment_count);
  const oem_segment* s = &oem->segments[segment];
  return find(oem->pairs + s->pair_first, s->pair_count, key);
}

size_t myb_oem_state_count(const myb_oem* oem, size_t segment) {
  MYB_CHECK(oem != NULL && segment < oem->segment_count);
  return oem->segments[segment].state_count;
}

int myb_oem_state(const myb_oem* oem, size_t segment, size_t index, myb_calendar* epoch,
                  double state[6], double acceleration[3]) {
  MYB_CHECK(oem != NULL && segment < oem->segment_count &&
            index < oem->segments[segment].state_count);
  const size_t k = oem->segments[segment].state_first + index;
  if (epoch != NULL) {
    *epoch = oem->epochs[k];
  }
  if (state != NULL) {
    memcpy(state, oem->values + 9 * k, 6 * sizeof(double));
  }
  if (acceleration != NULL) {
    memcpy(acceleration, oem->values + 9 * k + 6, 3 * sizeof(double));
  }
  return oem->accelerated[k] ? 3 : 1;
}

size_t myb_oem_covariance_count(const myb_oem* oem, size_t segment) {
  MYB_CHECK(oem != NULL && segment < oem->segment_count);
  return oem->segments[segment].covariance_count;
}

void myb_oem_covariance(const myb_oem* oem, size_t segment, size_t index, myb_calendar* epoch,
                        const char** ref_frame, double covariance[36]) {
  MYB_CHECK(oem != NULL && segment < oem->segment_count &&
            index < oem->segments[segment].covariance_count);
  const size_t k = oem->segments[segment].covariance_first + index;
  if (epoch != NULL) {
    *epoch = oem->covariance_epochs[k];
  }
  if (ref_frame != NULL) {
    *ref_frame = oem->covariance_frames[k];
  }
  if (covariance != NULL) {
    memcpy(covariance, oem->covariances + 36 * k, 36 * sizeof(double));
  }
}

static const char* const state_keys[6] = {"X", "Y", "Z", "X_DOT", "Y_DOT", "Z_DOT"};

myb_status myb_oem_write_as(const myb_oem_data* data, myb_odm_notation notation, char* buffer,
                            size_t size, size_t* length) {
  MYB_CHECK(data != NULL && (buffer != NULL || size == 0));
  if (data->segment_count == 0 || data->segments == NULL || !valid_text(data->originator) ||
      !valid_calendar(&data->creation_date) ||
      (notation != MYB_ODM_KVN && notation != MYB_ODM_XML)) {
    return MYB_E_DOMAIN;
  }
  for (size_t s = 0; s < data->segment_count; ++s) {
    const myb_oem_segment_data* seg = &data->segments[s];
    if (seg->count == 0 || seg->epochs == NULL || seg->states == NULL ||
        !valid_text(seg->object_name) || !valid_text(seg->object_id) ||
        !valid_text(seg->center_name) || !valid_text(seg->ref_frame) ||
        !valid_text(seg->time_system) ||
        (seg->interpolation != NULL && !valid_text(seg->interpolation)) ||
        (seg->covariance_count > 0 &&
         (seg->covariance_epochs == NULL || seg->covariances == NULL)) ||
        (seg->covariance_frame != NULL && !valid_text(seg->covariance_frame))) {
      return MYB_E_DOMAIN;
    }
    for (size_t k = 0; k < seg->count; ++k) {
      if (!valid_calendar(&seg->epochs[k])) {
        return MYB_E_DOMAIN;
      }
      for (int i = 0; i < 6; ++i) {
        if (!isfinite(seg->states[6 * k + (size_t)i])) {
          return MYB_E_DOMAIN;
        }
      }
    }
    for (size_t k = 0; k < seg->covariance_count; ++k) {
      if (!valid_calendar(&seg->covariance_epochs[k]) ||
          !covariance_finite(seg->covariances + 36 * k)) {
        return MYB_E_DOMAIN;
      }
    }
  }
  sink out = {.t = {.out = buffer, .size = size}, .xml = notation == MYB_ODM_XML};
  sink* k = &out;
  message(k, "oem", "CCSDS_OEM_VERS", "2.0");
  section(k, "header", NULL);
  put_time(k, "CREATION_DATE", &data->creation_date);
  put_text(k, "ORIGINATOR", data->originator);
  section_end(k, "header", NULL);
  section(k, "body", NULL);
  for (size_t s = 0; s < data->segment_count; ++s) {
    const myb_oem_segment_data* seg = &data->segments[s];
    section(k, "segment", NULL);
    if (!k->xml) {
      myb_internal_text_char(&k->t, '\n');
    }
    section(k, "metadata", "META_START");
    put_text(k, "OBJECT_NAME", seg->object_name);
    put_text(k, "OBJECT_ID", seg->object_id);
    put_text(k, "CENTER_NAME", seg->center_name);
    put_text(k, "REF_FRAME", seg->ref_frame);
    put_text(k, "TIME_SYSTEM", seg->time_system);
    put_time(k, "START_TIME", &seg->epochs[0]);
    put_time(k, "STOP_TIME", &seg->epochs[seg->count - 1]);
    if (seg->interpolation != NULL) {
      put_text(k, "INTERPOLATION", seg->interpolation);
      put_integer(k, "INTERPOLATION_DEGREE", seg->interpolation_degree);
    }
    section_end(k, "metadata", "META_STOP");
    if (!k->xml) {
      myb_internal_text_char(&k->t, '\n');
    }
    section(k, "data", NULL);
    for (size_t e = 0; e < seg->count; ++e) {
      if (k->xml) {
        section(k, "stateVector", NULL);
        put_time(k, "EPOCH", &seg->epochs[e]);
        for (int i = 0; i < 6; ++i) {
          put_number(k, state_keys[i], seg->states[6 * e + (size_t)i], metres_per_km);
        }
        section_end(k, "stateVector", NULL);
        continue;
      }
      put_calendar(&k->t, &seg->epochs[e]);
      for (int i = 0; i < 6; ++i) {
        myb_internal_text_char(&k->t, ' ');
        put_scaled(&k->t, seg->states[6 * e + (size_t)i], metres_per_km);
      }
      myb_internal_text_char(&k->t, '\n');
    }
    if (seg->covariance_count > 0 && !k->xml) {
      myb_internal_text_put(&k->t, "\nCOVARIANCE_START\n");
    }
    for (size_t e = 0; e < seg->covariance_count; ++e) {
      const double* covariance = seg->covariances + 36 * e;
      section(k, "covarianceMatrix", NULL);
      put_time(k, "EPOCH", &seg->covariance_epochs[e]);
      if (seg->covariance_frame != NULL) {
        put_text(k, "COV_REF_FRAME", seg->covariance_frame);
      }
      if (k->xml) {
        put_covariance(k, covariance, covariance_unit);
      } else {
        /* the triangle, a row a line */
        for (int row = 0; row < 6; ++row) {
          for (int column = 0; column <= row; ++column) {
            if (column > 0) {
              myb_internal_text_char(&k->t, ' ');
            }
            put_scaled(&k->t, covariance[6 * row + column], covariance_unit);
          }
          myb_internal_text_char(&k->t, '\n');
        }
      }
      section_end(k, "covarianceMatrix", NULL);
    }
    if (seg->covariance_count > 0 && !k->xml) {
      myb_internal_text_put(&k->t, "COVARIANCE_STOP\n");
    }
    section_end(k, "data", NULL);
    section_end(k, "segment", NULL);
  }
  section_end(k, "body", NULL);
  return message_end(k, "oem", length);
}

myb_status myb_oem_write(const myb_oem_data* data, char* buffer, size_t size, size_t* length) {
  return myb_oem_write_as(data, MYB_ODM_KVN, buffer, size, length);
}

/* ---- OPM ------------------------------------------------------------------------------------- */

/* A number of a parameter message, scaled to SI. */
static myb_status scaled(const cursor* c, const line* l, double scale, double* out,
                         myb_diag* diag) {
  double v = 0.0;
  if (number_of(c->text + l->value, l->value_length, &v) != MYB_OK) {
    return fail(diag, l->number, "not a number");
  }
  *out = v * scale;
  return MYB_OK;
}

static myb_status dated(const cursor* c, const line* l, myb_calendar* out, myb_diag* diag) {
  if (calendar_of(c->text + l->value, l->value_length, out) != MYB_OK) {
    return fail(diag, l->number, "not a date of the form YYYY-MM-DDThh:mm:ss");
  }
  return MYB_OK;
}

static myb_status named(const cursor* c, const line* l, char out[MYB_ODM_TEXT], myb_diag* diag) {
  if (text_of(c->text + l->value, l->value_length, out) != MYB_OK) {
    return fail(diag, l->number, "an empty value or one longer than 63 characters");
  }
  if (c->xml && !l->raw) {
    out[xml_decode(out, l->value_length)] = '\0';
  }
  return MYB_OK;
}

/* A keyword of the covariance of a parameter message: the frame or one of the 21 numbers
 * (*seen gathers which). Returns 0 if the keyword is neither. */
static int covariance_keyword(const cursor* c, const line* l, double unit,
                              char frame[MYB_ODM_TEXT], double covariance[36], unsigned* seen,
                              myb_status* status, myb_diag* diag) {
  if (key_is(c, l, "COV_REF_FRAME")) {
    *status = named(c, l, frame, diag);
    return 1;
  }
  if (l->key_length < 4 || c->text[l->key] != 'C') {
    return 0;
  }
  for (int k = 0; k < 21; ++k) {
    if (key_is(c, l, covariance_keys[k])) {
      double v = 0.0;
      *status = scaled(c, l, unit, &v, diag);
      covariance_set(covariance, k, v);
      *seen |= 1u << (unsigned)k;
      return 1;
    }
  }
  return 0;
}

myb_status myb_opm_parse(const char* text, size_t size, myb_opm* opm, myb_diag* diag) {
  MYB_CHECK(text != NULL && opm != NULL);
  myb_internal_diag_clear(diag);
  enum {
    S_VERSION = 1u << 0, S_CREATION = 1u << 1, S_ORIGINATOR = 1u << 2, S_NAME = 1u << 3,
    S_ID = 1u << 4, S_CENTER = 1u << 5, S_FRAME = 1u << 6, S_TIME = 1u << 7, S_EPOCH = 1u << 8,
    S_STATE = 1u << 9, /* six bits from here */
    S_ELEMENTS = 1u << 15, /* seven bits from here */
    S_MANDATORY = (1u << 15) - 1u
  };
  static const char* const missing[15] = {
      "CCSDS_OPM_VERS", "CREATION_DATE", "ORIGINATOR", "OBJECT_NAME", "OBJECT_ID",
      "CENTER_NAME",    "REF_FRAME",     "TIME_SYSTEM", "EPOCH",      "X",
      "Y",              "Z",             "X_DOT",       "Y_DOT",      "Z_DOT"};
  myb_opm* const out = opm;
  memset(out, 0, sizeof *out);
#define o (*out)
  unsigned seen = 0;
  unsigned covariance_seen = 0;
  cursor c = {.text = text, .size = size, .xml = is_xml(text, size), .flat = 1};
  line l;
  myb_status status = MYB_OK;
  myb_opm_maneuver* man = NULL;
  /* bound: the lines of the text */
  for (next_line(&c, &l); l.kind != LINE_END && status == MYB_OK; next_line(&c, &l)) {
    if (l.kind != LINE_PAIR) {
      return fail(diag, l.number, "OPM: not a keyword = value line");
    }
    if (seen == 0 && !key_is(&c, &l, "CCSDS_OPM_VERS")) {
      return fail(diag, l.number, "OPM: the message does not start with CCSDS_OPM_VERS");
    }
    if (key_is(&c, &l, "CCSDS_OPM_VERS")) {
      seen |= S_VERSION;
    } else if (key_is(&c, &l, "CREATION_DATE")) {
      status = dated(&c, &l, &o.creation_date, diag);
      seen |= S_CREATION;
    } else if (key_is(&c, &l, "ORIGINATOR")) {
      status = named(&c, &l, o.originator, diag);
      seen |= S_ORIGINATOR;
    } else if (key_is(&c, &l, "OBJECT_NAME")) {
      status = named(&c, &l, o.object_name, diag);
      seen |= S_NAME;
    } else if (key_is(&c, &l, "OBJECT_ID")) {
      status = named(&c, &l, o.object_id, diag);
      seen |= S_ID;
    } else if (key_is(&c, &l, "CENTER_NAME")) {
      status = named(&c, &l, o.center_name, diag);
      seen |= S_CENTER;
    } else if (key_is(&c, &l, "REF_FRAME")) {
      status = named(&c, &l, o.ref_frame, diag);
      seen |= S_FRAME;
    } else if (key_is(&c, &l, "TIME_SYSTEM")) {
      status = named(&c, &l, o.time_system, diag);
      seen |= S_TIME;
    } else if (key_is(&c, &l, "EPOCH")) {
      status = dated(&c, &l, &o.epoch, diag);
      seen |= S_EPOCH;
    } else if (key_is(&c, &l, "SEMI_MAJOR_AXIS")) {
      status = scaled(&c, &l, metres_per_km, &o.semi_major_axis, diag);
      seen |= S_ELEMENTS;
    } else if (key_is(&c, &l, "ECCENTRICITY")) {
      status = scaled(&c, &l, 1.0, &o.eccentricity, diag);
      seen |= S_ELEMENTS << 1;
    } else if (key_is(&c, &l, "INCLINATION")) {
      status = scaled(&c, &l, degree, &o.inclination, diag);
      seen |= S_ELEMENTS << 2;
    } else if (key_is(&c, &l, "RA_OF_ASC_NODE")) {
      status = scaled(&c, &l, degree, &o.ra_of_asc_node, diag);
      seen |= S_ELEMENTS << 3;
    } else if (key_is(&c, &l, "ARG_OF_PERICENTER")) {
      status = scaled(&c, &l, degree, &o.arg_of_pericenter, diag);
      seen |= S_ELEMENTS << 4;
    } else if (key_is(&c, &l, "TRUE_ANOMALY") || key_is(&c, &l, "MEAN_ANOMALY")) {
      status = scaled(&c, &l, degree, &o.anomaly, diag);
      o.anomaly_is_mean = key_is(&c, &l, "MEAN_ANOMALY");
      seen |= S_ELEMENTS << 5;
    } else if (key_is(&c, &l, "GM")) {
      status = scaled(&c, &l, gm_unit, &o.gm, diag);
      seen |= S_ELEMENTS << 6;
    } else if (key_is(&c, &l, "MASS")) {
      status = scaled(&c, &l, 1.0, &o.mass, diag);
      o.has_spacecraft = 1;
    } else if (key_is(&c, &l, "SOLAR_RAD_AREA")) {
      status = scaled(&c, &l, 1.0, &o.solar_rad_area, diag);
      o.has_spacecraft = 1;
    } else if (key_is(&c, &l, "SOLAR_RAD_COEFF")) {
      status = scaled(&c, &l, 1.0, &o.solar_rad_coeff, diag);
      o.has_spacecraft = 1;
    } else if (key_is(&c, &l, "DRAG_AREA")) {
      status = scaled(&c, &l, 1.0, &o.drag_area, diag);
      o.has_spacecraft = 1;
    } else if (key_is(&c, &l, "DRAG_COEFF")) {
      status = scaled(&c, &l, 1.0, &o.drag_coeff, diag);
      o.has_spacecraft = 1;
    } else if (key_is(&c, &l, "MAN_EPOCH_IGNITION")) {
      if (o.maneuver_count == MYB_OPM_MAX_MANEUVERS) {
        return fail(diag, l.number, "OPM: too many manoeuvres");
      }
      man = &o.maneuvers[o.maneuver_count++];
      status = dated(&c, &l, &man->epoch_ignition, diag);
    } else if (l.key_length > 4 && memcmp(text + l.key, "MAN_", 4) == 0) {
      if (man == NULL) {
        return fail(diag, l.number, "OPM: a manoeuvre must start with MAN_EPOCH_IGNITION");
      }
      if (key_is(&c, &l, "MAN_DURATION")) {
        status = scaled(&c, &l, 1.0, &man->duration, diag);
      } else if (key_is(&c, &l, "MAN_DELTA_MASS")) {
        status = scaled(&c, &l, 1.0, &man->delta_mass, diag);
      } else if (key_is(&c, &l, "MAN_REF_FRAME")) {
        status = named(&c, &l, man->ref_frame, diag);
      } else if (key_is(&c, &l, "MAN_DV_1")) {
        status = scaled(&c, &l, metres_per_km, &man->dv[0], diag);
      } else if (key_is(&c, &l, "MAN_DV_2")) {
        status = scaled(&c, &l, metres_per_km, &man->dv[1], diag);
      } else if (key_is(&c, &l, "MAN_DV_3")) {
        status = scaled(&c, &l, metres_per_km, &man->dv[2], diag);
      }
    } else if (!covariance_keyword(&c, &l, covariance_unit, o.cov_ref_frame, o.covariance,
                                   &covariance_seen, &status, diag)) {
      for (int i = 0; i < 6; ++i) {
        if (key_is(&c, &l, state_keys[i])) {
          status = scaled(&c, &l, metres_per_km, &o.state[i], diag);
          seen |= (unsigned)S_STATE << (unsigned)i;
        }
      }
      /* anything else (user-defined keywords, later versions) is passed over */
    }
  }
  if (status != MYB_OK) {
    return status;
  }
  for (int k = 0; k < 15; ++k) {
    if (!(seen & (1u << k))) {
      myb_internal_diag(diag, 0, "OPM: no %s", missing[k]);
      return MYB_E_FORMAT;
    }
  }
  const unsigned elements = (seen >> 15) & 0x7Fu;
  if (elements != 0 && elements != 0x7Fu) {
    return fail(diag, 0, "OPM: the osculating elements are incomplete");
  }
  o.has_elements = elements != 0;
  if (covariance_seen != 0 && covariance_seen != 0x1FFFFFu) {
    return fail(diag, 0, "OPM: the covariance matrix is incomplete");
  }
  o.has_covariance = covariance_seen != 0;
#undef o
  return MYB_OK;
}

myb_status myb_opm_write_as(const myb_opm* opm, myb_odm_notation notation, char* buffer,
                            size_t size, size_t* length) {
  MYB_CHECK(opm != NULL && (buffer != NULL || size == 0));
  if (!valid_text(opm->originator) || !valid_text(opm->object_name) ||
      !valid_text(opm->object_id) || !valid_text(opm->center_name) ||
      !valid_text(opm->ref_frame) || !valid_text(opm->time_system) ||
      !valid_calendar(&opm->creation_date) || !valid_calendar(&opm->epoch) ||
      opm->maneuver_count > MYB_OPM_MAX_MANEUVERS ||
      (opm->has_covariance && !covariance_finite(opm->covariance)) ||
      (notation != MYB_ODM_KVN && notation != MYB_ODM_XML)) {
    return MYB_E_DOMAIN;
  }
  for (int i = 0; i < 6; ++i) {
    if (!isfinite(opm->state[i])) {
      return MYB_E_DOMAIN;
    }
  }
  for (size_t m = 0; m < opm->maneuver_count; ++m) {
    if (!valid_calendar(&opm->maneuvers[m].epoch_ignition) ||
        !valid_text(opm->maneuvers[m].ref_frame)) {
      return MYB_E_DOMAIN;
    }
  }
  sink out = {.t = {.out = buffer, .size = size}, .xml = notation == MYB_ODM_XML};
  sink* k = &out;
  message(k, "opm", "CCSDS_OPM_VERS", "2.0");
  section(k, "header", NULL);
  put_time(k, "CREATION_DATE", &opm->creation_date);
  put_text(k, "ORIGINATOR", opm->originator);
  section_end(k, "header", NULL);
  section(k, "body", NULL);
  section(k, "segment", NULL);
  section(k, "metadata", NULL);
  put_text(k, "OBJECT_NAME", opm->object_name);
  put_text(k, "OBJECT_ID", opm->object_id);
  put_text(k, "CENTER_NAME", opm->center_name);
  put_text(k, "REF_FRAME", opm->ref_frame);
  put_text(k, "TIME_SYSTEM", opm->time_system);
  section_end(k, "metadata", NULL);
  section(k, "data", NULL);
  section(k, "stateVector", NULL);
  put_time(k, "EPOCH", &opm->epoch);
  for (int i = 0; i < 6; ++i) {
    put_number(k, state_keys[i], opm->state[i], metres_per_km);
  }
  section_end(k, "stateVector", NULL);
  if (opm->has_elements) {
    section(k, "keplerianElements", NULL);
    put_number(k, "SEMI_MAJOR_AXIS", opm->semi_major_axis, metres_per_km);
    put_number(k, "ECCENTRICITY", opm->eccentricity, 1.0);
    put_number(k, "INCLINATION", opm->inclination, degree);
    put_number(k, "RA_OF_ASC_NODE", opm->ra_of_asc_node, degree);
    put_number(k, "ARG_OF_PERICENTER", opm->arg_of_pericenter, degree);
    put_number(k, opm->anomaly_is_mean ? "MEAN_ANOMALY" : "TRUE_ANOMALY", opm->anomaly, degree);
    put_number(k, "GM", opm->gm, gm_unit);
    section_end(k, "keplerianElements", NULL);
  }
  if (opm->has_spacecraft) {
    section(k, "spacecraftParameters", NULL);
    put_number(k, "MASS", opm->mass, 1.0);
    put_number(k, "SOLAR_RAD_AREA", opm->solar_rad_area, 1.0);
    put_number(k, "SOLAR_RAD_COEFF", opm->solar_rad_coeff, 1.0);
    put_number(k, "DRAG_AREA", opm->drag_area, 1.0);
    put_number(k, "DRAG_COEFF", opm->drag_coeff, 1.0);
    section_end(k, "spacecraftParameters", NULL);
  }
  if (opm->has_covariance) {
    section(k, "covarianceMatrix", NULL);
    if (opm->cov_ref_frame[0] != '\0') {
      put_text(k, "COV_REF_FRAME", opm->cov_ref_frame);
    }
    put_covariance(k, opm->covariance, covariance_unit);
    section_end(k, "covarianceMatrix", NULL);
  }
  for (size_t m = 0; m < opm->maneuver_count; ++m) {
    const myb_opm_maneuver* man = &opm->maneuvers[m];
    section(k, "maneuverParameters", NULL);
    put_time(k, "MAN_EPOCH_IGNITION", &man->epoch_ignition);
    put_number(k, "MAN_DURATION", man->duration, 1.0);
    put_number(k, "MAN_DELTA_MASS", man->delta_mass, 1.0);
    put_text(k, "MAN_REF_FRAME", man->ref_frame);
    put_number(k, "MAN_DV_1", man->dv[0], metres_per_km);
    put_number(k, "MAN_DV_2", man->dv[1], metres_per_km);
    put_number(k, "MAN_DV_3", man->dv[2], metres_per_km);
    section_end(k, "maneuverParameters", NULL);
  }
  section_end(k, "data", NULL);
  section_end(k, "segment", NULL);
  section_end(k, "body", NULL);
  return message_end(k, "opm", length);
}

myb_status myb_opm_write(const myb_opm* opm, char* buffer, size_t size, size_t* length) {
  return myb_opm_write_as(opm, MYB_ODM_KVN, buffer, size, length);
}

/* ---- OMM ------------------------------------------------------------------------------------- */

myb_status myb_omm_parse(const char* text, size_t size, myb_omm* omm, myb_diag* diag) {
  MYB_CHECK(text != NULL && omm != NULL);
  myb_internal_diag_clear(diag);
  enum {
    M_VERSION = 1u << 0, M_CREATION = 1u << 1, M_ORIGINATOR = 1u << 2, M_NAME = 1u << 3,
    M_ID = 1u << 4, M_CENTER = 1u << 5, M_FRAME = 1u << 6, M_TIME = 1u << 7, M_THEORY = 1u << 8,
    M_EPOCH = 1u << 9, M_ECC = 1u << 10, M_INC = 1u << 11, M_NODE = 1u << 12, M_ARG = 1u << 13,
    M_ANOMALY = 1u << 14, M_MOTION = 1u << 15, M_AXIS = 1u << 16
  };
  static const char* const missing[15] = {
      "CCSDS_OMM_VERS", "CREATION_DATE", "ORIGINATOR",          "OBJECT_NAME",  "OBJECT_ID",
      "CENTER_NAME",    "REF_FRAME",     "TIME_SYSTEM",         "MEAN_ELEMENT_THEORY",
      "EPOCH",          "ECCENTRICITY",  "INCLINATION",         "RA_OF_ASC_NODE",
      "ARG_OF_PERICENTER", "MEAN_ANOMALY"};
  myb_omm* const out = omm;
  memset(out, 0, sizeof *out);
#define o (*out)
  unsigned seen = 0;
  unsigned covariance_seen = 0;
  cursor c = {.text = text, .size = size, .xml = is_xml(text, size), .flat = 1};
  line l;
  myb_status status = MYB_OK;
  /* bound: the lines of the text */
  for (next_line(&c, &l); l.kind != LINE_END && status == MYB_OK; next_line(&c, &l)) {
    if (l.kind != LINE_PAIR) {
      return fail(diag, l.number, "OMM: not a keyword = value line");
    }
    if (seen == 0 && !key_is(&c, &l, "CCSDS_OMM_VERS")) {
      return fail(diag, l.number, "OMM: the message does not start with CCSDS_OMM_VERS");
    }
    if (key_is(&c, &l, "CCSDS_OMM_VERS")) {
      seen |= M_VERSION;
    } else if (key_is(&c, &l, "CREATION_DATE")) {
      status = dated(&c, &l, &o.creation_date, diag);
      seen |= M_CREATION;
    } else if (key_is(&c, &l, "ORIGINATOR")) {
      status = named(&c, &l, o.originator, diag);
      seen |= M_ORIGINATOR;
    } else if (key_is(&c, &l, "OBJECT_NAME")) {
      status = named(&c, &l, o.object_name, diag);
      seen |= M_NAME;
    } else if (key_is(&c, &l, "OBJECT_ID")) {
      status = named(&c, &l, o.object_id, diag);
      seen |= M_ID;
    } else if (key_is(&c, &l, "CENTER_NAME")) {
      status = named(&c, &l, o.center_name, diag);
      seen |= M_CENTER;
    } else if (key_is(&c, &l, "REF_FRAME")) {
      status = named(&c, &l, o.ref_frame, diag);
      seen |= M_FRAME;
    } else if (key_is(&c, &l, "TIME_SYSTEM")) {
      status = named(&c, &l, o.time_system, diag);
      seen |= M_TIME;
    } else if (key_is(&c, &l, "MEAN_ELEMENT_THEORY")) {
      status = named(&c, &l, o.mean_element_theory, diag);
      seen |= M_THEORY;
    } else if (key_is(&c, &l, "EPOCH")) {
      status = dated(&c, &l, &o.epoch, diag);
      seen |= M_EPOCH;
    } else if (key_is(&c, &l, "MEAN_MOTION")) {
      status = scaled(&c, &l, 1.0, &o.mean_motion, diag);
      o.has_mean_motion = 1;
      seen |= M_MOTION;
    } else if (key_is(&c, &l, "SEMI_MAJOR_AXIS")) {
      status = scaled(&c, &l, 1.0, &o.semi_major_axis, diag);
      seen |= M_AXIS;
    } else if (key_is(&c, &l, "ECCENTRICITY")) {
      status = scaled(&c, &l, 1.0, &o.eccentricity, diag);
      seen |= M_ECC;
    } else if (key_is(&c, &l, "INCLINATION")) {
      status = scaled(&c, &l, 1.0, &o.inclination, diag);
      seen |= M_INC;
    } else if (key_is(&c, &l, "RA_OF_ASC_NODE")) {
      status = scaled(&c, &l, 1.0, &o.ra_of_asc_node, diag);
      seen |= M_NODE;
    } else if (key_is(&c, &l, "ARG_OF_PERICENTER")) {
      status = scaled(&c, &l, 1.0, &o.arg_of_pericenter, diag);
      seen |= M_ARG;
    } else if (key_is(&c, &l, "MEAN_ANOMALY")) {
      status = scaled(&c, &l, 1.0, &o.mean_anomaly, diag);
      seen |= M_ANOMALY;
    } else if (key_is(&c, &l, "GM")) {
      status = scaled(&c, &l, 1.0, &o.gm, diag);
      o.has_gm = 1;
    } else if (key_is(&c, &l, "EPHEMERIS_TYPE")) {
      long v = 0;
      status = integer_of(text + l.value, l.value_length, &v) == MYB_OK
                   ? MYB_OK
                   : fail(diag, l.number, "not an integer");
      o.ephemeris_type = (int)v;
      o.has_tle = 1;
    } else if (key_is(&c, &l, "CLASSIFICATION_TYPE")) {
      if (l.value_length != 1) {
        return fail(diag, l.number, "OMM: the classification is one letter");
      }
      o.classification_type = text[l.value];
      o.has_tle = 1;
    } else if (key_is(&c, &l, "NORAD_CAT_ID") || key_is(&c, &l, "ELEMENT_SET_NO") ||
               key_is(&c, &l, "REV_AT_EPOCH")) {
      long v = 0;
      status = integer_of(text + l.value, l.value_length, &v) == MYB_OK
                   ? MYB_OK
                   : fail(diag, l.number, "not an integer");
      if (key_is(&c, &l, "NORAD_CAT_ID")) {
        o.norad_cat_id = v;
      } else if (key_is(&c, &l, "ELEMENT_SET_NO")) {
        o.element_set_no = v;
      } else {
        o.rev_at_epoch = v;
      }
      o.has_tle = 1;
    } else if (key_is(&c, &l, "BSTAR")) {
      status = scaled(&c, &l, 1.0, &o.bstar, diag);
      o.has_tle = 1;
    } else if (key_is(&c, &l, "MEAN_MOTION_DOT")) {
      status = scaled(&c, &l, 1.0, &o.mean_motion_dot, diag);
      o.has_tle = 1;
    } else if (key_is(&c, &l, "MEAN_MOTION_DDOT")) {
      status = scaled(&c, &l, 1.0, &o.mean_motion_ddot, diag);
      o.has_tle = 1;
    } else {
      (void)covariance_keyword(&c, &l, 1.0, o.cov_ref_frame, o.covariance, &covariance_seen,
                               &status, diag);
    }
    /* anything else (spacecraft parameters, user-defined keywords) is passed over */
  }
  if (status != MYB_OK) {
    return status;
  }
  for (int k = 0; k < 15; ++k) {
    if (!(seen & (1u << k))) {
      myb_internal_diag(diag, 0, "OMM: no %s", missing[k]);
      return MYB_E_FORMAT;
    }
  }
  if (!(seen & (M_MOTION | M_AXIS))) {
    return fail(diag, 0, "OMM: neither MEAN_MOTION nor SEMI_MAJOR_AXIS");
  }
  if (o.has_tle && o.classification_type == '\0') {
    o.classification_type = 'U';
  }
  if (covariance_seen != 0 && covariance_seen != 0x1FFFFFu) {
    return fail(diag, 0, "OMM: the covariance matrix is incomplete");
  }
  o.has_covariance = covariance_seen != 0;
#undef o
  return MYB_OK;
}

myb_status myb_omm_write_as(const myb_omm* omm, myb_odm_notation notation, char* buffer,
                            size_t size, size_t* length) {
  MYB_CHECK(omm != NULL && (buffer != NULL || size == 0));
  if (!valid_text(omm->originator) || !valid_text(omm->object_name) ||
      !valid_text(omm->object_id) || !valid_text(omm->center_name) ||
      !valid_text(omm->ref_frame) || !valid_text(omm->time_system) ||
      !valid_text(omm->mean_element_theory) || !valid_calendar(&omm->creation_date) ||
      !valid_calendar(&omm->epoch) || !isfinite(omm->eccentricity) ||
      !isfinite(omm->inclination) || !isfinite(omm->ra_of_asc_node) ||
      !isfinite(omm->arg_of_pericenter) || !isfinite(omm->mean_anomaly) ||
      !isfinite(omm->has_mean_motion ? omm->mean_motion : omm->semi_major_axis) ||
      (omm->has_covariance && !covariance_finite(omm->covariance)) ||
      (notation != MYB_ODM_KVN && notation != MYB_ODM_XML)) {
    return MYB_E_DOMAIN;
  }
  sink out = {.t = {.out = buffer, .size = size}, .xml = notation == MYB_ODM_XML};
  sink* k = &out;
  message(k, "omm", "CCSDS_OMM_VERS", "2.0");
  section(k, "header", NULL);
  put_time(k, "CREATION_DATE", &omm->creation_date);
  put_text(k, "ORIGINATOR", omm->originator);
  section_end(k, "header", NULL);
  section(k, "body", NULL);
  section(k, "segment", NULL);
  section(k, "metadata", NULL);
  put_text(k, "OBJECT_NAME", omm->object_name);
  put_text(k, "OBJECT_ID", omm->object_id);
  put_text(k, "CENTER_NAME", omm->center_name);
  put_text(k, "REF_FRAME", omm->ref_frame);
  put_text(k, "TIME_SYSTEM", omm->time_system);
  put_text(k, "MEAN_ELEMENT_THEORY", omm->mean_element_theory);
  section_end(k, "metadata", NULL);
  section(k, "data", NULL);
  section(k, "meanElements", NULL);
  put_time(k, "EPOCH", &omm->epoch);
  if (omm->has_mean_motion) {
    put_number(k, "MEAN_MOTION", omm->mean_motion, 1.0);
  } else {
    put_number(k, "SEMI_MAJOR_AXIS", omm->semi_major_axis, 1.0);
  }
  put_number(k, "ECCENTRICITY", omm->eccentricity, 1.0);
  put_number(k, "INCLINATION", omm->inclination, 1.0);
  put_number(k, "RA_OF_ASC_NODE", omm->ra_of_asc_node, 1.0);
  put_number(k, "ARG_OF_PERICENTER", omm->arg_of_pericenter, 1.0);
  put_number(k, "MEAN_ANOMALY", omm->mean_anomaly, 1.0);
  if (omm->has_gm) {
    put_number(k, "GM", omm->gm, 1.0);
  }
  section_end(k, "meanElements", NULL);
  if (omm->has_tle) {
    const char classification[2] = {omm->classification_type, '\0'};
    section(k, "tleParameters", NULL);
    put_integer(k, "EPHEMERIS_TYPE", omm->ephemeris_type);
    put_text(k, "CLASSIFICATION_TYPE", classification);
    put_integer(k, "NORAD_CAT_ID", omm->norad_cat_id);
    put_integer(k, "ELEMENT_SET_NO", omm->element_set_no);
    put_integer(k, "REV_AT_EPOCH", omm->rev_at_epoch);
    put_number(k, "BSTAR", omm->bstar, 1.0);
    put_number(k, "MEAN_MOTION_DOT", omm->mean_motion_dot, 1.0);
    put_number(k, "MEAN_MOTION_DDOT", omm->mean_motion_ddot, 1.0);
    section_end(k, "tleParameters", NULL);
  }
  if (omm->has_covariance) {
    section(k, "covarianceMatrix", NULL);
    if (omm->cov_ref_frame[0] != '\0') {
      put_text(k, "COV_REF_FRAME", omm->cov_ref_frame);
    }
    put_covariance(k, omm->covariance, 1.0);
    section_end(k, "covarianceMatrix", NULL);
  }
  section_end(k, "data", NULL);
  section_end(k, "segment", NULL);
  section_end(k, "body", NULL);
  return message_end(k, "omm", length);
}

myb_status myb_omm_write(const myb_omm* omm, char* buffer, size_t size, size_t* length) {
  return myb_omm_write_as(omm, MYB_ODM_KVN, buffer, size, length);
}

/* ---- OCM ------------------------------------------------------------------------------------- */

enum { OCM_KINDS = 7 };

static const char* const ocm_start[OCM_KINDS] = {"TRAJ_START", "PHYS_START", "COV_START",
                                                 "MAN_START",  "PERT_START", "OD_START",
                                                 "USER_START"};
static const char* const ocm_stop[OCM_KINDS] = {"TRAJ_STOP", "PHYS_STOP", "COV_STOP", "MAN_STOP",
                                                "PERT_STOP", "OD_STOP",   "USER_STOP"};
static const char* const ocm_element[OCM_KINDS] = {"traj", "phys", "cov", "man",
                                                   "pert", "od",   "user"};
/* the element of a line of data; NULL: the kind has no lines */
static const char* const ocm_line_element[OCM_KINDS] = {"trajLine", NULL, "covLine", "manLine",
                                                        NULL,       NULL, NULL};
static const char user_prefix[] = "USER_DEFINED_";

typedef struct ocm_block {
  int kind;
  size_t pair_first;
  size_t pair_count;
  size_t line_first;
  size_t line_count;
} ocm_block;

struct myb_ocm {
  myb_allocator allocator;
  char* text;
  pair* pairs; /* the header's, the metadata's, then each block's */
  size_t header_pairs;
  size_t meta_pairs;
  size_t pair_count;
  ocm_block* blocks;
  size_t block_count;
  const char** lines;
  size_t line_count;
};

void myb_ocm_destroy(myb_ocm* ocm) {
  if (ocm == NULL) {
    return;
  }
  const myb_allocator a = ocm->allocator;
  myb_internal_release(&a, ocm->text);
  myb_internal_release(&a, ocm->pairs);
  myb_internal_release(&a, ocm->blocks);
  myb_internal_release(&a, (void*)ocm->lines);
  myb_internal_release(&a, ocm);
}

enum { OCM_HEADER = -3, OCM_META = -2, OCM_BODY = -1 }; /* 0 ..: inside a block of that kind */

static myb_status ocm_pass(myb_ocm* ocm, size_t size, int fill, myb_diag* diag) {
  cursor c = {.text = ocm->text, .size = size, .xml = is_xml(ocm->text, size)};
  int mode = OCM_HEADER;
  size_t pairs = 0;
  size_t blocks = 0;
  size_t lines = 0;
  line l;
  /* bound: the lines of the text */
  for (next_line(&c, &l); l.kind != LINE_END; next_line(&c, &l)) {
    if (l.kind == LINE_MARK) {
      int handled = 0;
      if (key_is(&c, &l, "META_START") && mode == OCM_HEADER) {
        mode = OCM_META;
        handled = 1;
      } else if (key_is(&c, &l, "META_STOP") && mode == OCM_META) {
        mode = OCM_BODY;
        handled = 1;
      }
      for (int kind = 0; kind < OCM_KINDS && !handled; ++kind) {
        if (mode == OCM_BODY && key_is(&c, &l, ocm_start[kind])) {
          if (fill) {
            ocm->blocks[blocks] =
                (ocm_block){.kind = kind, .pair_first = pairs, .line_first = lines};
          }
          ++blocks;
          mode = kind;
          handled = 1;
        } else if (mode == kind && key_is(&c, &l, ocm_stop[kind])) {
          mode = OCM_BODY;
          handled = 1;
        }
      }
      if (handled) {
        continue;
      }
      if (mode < 0 || ocm_line_element[mode] == NULL || l.alias != NULL) {
        return fail(diag, l.number, "OCM: a keyword out of place");
      }
      /* a line of data of one word that looks like a keyword */
      l.kind = LINE_DATA;
      l.value = l.key;
      l.value_length = l.key_length;
    }
    if (l.kind == LINE_PAIR) {
      if (mode == OCM_BODY) {
        return fail(diag, l.number, "OCM: a keyword outside the sections");
      }
      if (pairs == 0 && !key_is(&c, &l, "CCSDS_OCM_VERS")) {
        return fail(diag, l.number, "OCM: the message does not start with CCSDS_OCM_VERS");
      }
      if (fill) {
        const size_t value_length =
            c.xml && !l.raw ? xml_decode(ocm->text + l.value, l.value_length) : l.value_length;
        ocm->text[l.key + l.key_length] = '\0';
        ocm->text[l.value + value_length] = '\0';
        ocm->pairs[pairs] = (pair){ocm->text + l.key, ocm->text + l.value};
        if (mode >= 0) {
          ocm->blocks[blocks - 1].pair_count++;
        }
      }
      ++pairs;
      if (mode == OCM_HEADER) {
        ocm->header_pairs = pairs;
      } else if (mode == OCM_META) {
        ocm->meta_pairs = pairs - ocm->header_pairs;
      }
      continue;
    }
    if (mode < 0 || ocm_line_element[mode] == NULL) {
      return fail(diag, l.number, "OCM: a line of data where there is none");
    }
    if (fill) {
      const size_t length =
          c.xml && !l.raw ? xml_decode(ocm->text + l.value, l.value_length) : l.value_length;
      ocm->text[l.value + length] = '\0';
      ocm->lines[lines] = ocm->text + l.value;
      ocm->blocks[blocks - 1].line_count++;
    }
    ++lines;
  }
  if (mode != OCM_BODY) {
    return fail(diag, c.line, mode == OCM_HEADER ? "OCM: no META_START"
                                                 : "OCM: a section is not closed");
  }
  ocm->pair_count = pairs;
  ocm->block_count = blocks;
  ocm->line_count = lines;
  return MYB_OK;
}

myb_status myb_ocm_load(const char* text, size_t size, const myb_allocator* allocator,
                        myb_ocm** ocm, myb_diag* diag) {
  MYB_CHECK(text != NULL && ocm != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  void* block = NULL;
  myb_status status = myb_internal_alloc_array(&a, 1, sizeof(myb_ocm), &block);
  if (status != MYB_OK) {
    return status;
  }
  myb_ocm* o = block;
  memset(o, 0, sizeof *o);
  o->allocator = a;
  status = myb_internal_alloc_array(&a, size + 1, 1, &block);
  if (status != MYB_OK) {
    goto fail;
  }
  o->text = block;
  memcpy(o->text, text, size);
  o->text[size] = '\0';
  status = ocm_pass(o, size, 0, diag);
  if (status != MYB_OK) {
    goto fail;
  }
  status = myb_internal_alloc_array(&a, o->pair_count, sizeof(pair), &block);
  if (status == MYB_OK) {
    o->pairs = block;
    status = myb_internal_alloc_array(&a, o->block_count, sizeof(ocm_block), &block);
  }
  if (status == MYB_OK) {
    o->blocks = block;
    status = myb_internal_alloc_array(&a, o->line_count, sizeof(const char*), &block);
  }
  if (status != MYB_OK) {
    goto fail;
  }
  o->lines = block;
  o->header_pairs = 0;
  o->meta_pairs = 0;
  status = ocm_pass(o, size, 1, diag);
  if (status != MYB_OK) {
    goto fail;
  }
  static const char* const mandatory[2] = {"TIME_SYSTEM", "EPOCH_TZERO"};
  for (int k = 0; k < 2; ++k) {
    if (myb_ocm_meta(o, mandatory[k]) == NULL) {
      myb_internal_diag(diag, 0, "OCM: the metadata have no %s", mandatory[k]);
      status = MYB_E_FORMAT;
      goto fail;
    }
  }
  *ocm = o;
  return MYB_OK;
fail:
  myb_ocm_destroy(o);
  return status;
}

myb_status myb_ocm_load_file(const char* path, const myb_allocator* allocator, myb_ocm** ocm,
                             myb_diag* diag) {
  MYB_CHECK(path != NULL && ocm != NULL);
  myb_internal_diag_clear(diag);
  const myb_allocator a = myb_internal_allocator(allocator);
  char* data = NULL;
  size_t size = 0;
  const myb_status status = myb_platform_read_file(path, &a, max_file, &data, &size);
  if (status != MYB_OK) {
    myb_internal_diag(diag, 0, "cannot read '%s'", path);
    return status;
  }
  const myb_status loaded = myb_ocm_load(data, size, allocator, ocm, diag);
  myb_internal_release(&a, data);
  return loaded;
}

/* As find; a user-defined keyword is kept without its prefix when the message is in XML. */
static const char* ocm_find(const pair* pairs, size_t count, const char* key) {
  const size_t prefix = sizeof user_prefix - 1;
  const int user = strncmp(key, user_prefix, prefix) == 0;
  for (size_t k = 0; k < count; ++k) {
    if (strcmp(pairs[k].key, key) == 0 || (user && strcmp(pairs[k].key, key + prefix) == 0)) {
      return pairs[k].value;
    }
  }
  return NULL;
}

const char* myb_ocm_header(const myb_ocm* ocm, const char* key) {
  MYB_CHECK(ocm != NULL && key != NULL);
  return find(ocm->pairs, ocm->header_pairs, key);
}

const char* myb_ocm_meta(const myb_ocm* ocm, const char* key) {
  MYB_CHECK(ocm != NULL && key != NULL);
  return find(ocm->pairs + ocm->header_pairs, ocm->meta_pairs, key);
}

size_t myb_ocm_block_count(const myb_ocm* ocm, myb_ocm_kind kind) {
  MYB_CHECK(ocm != NULL && (int)kind >= 0 && (int)kind < OCM_KINDS);
  size_t count = 0;
  for (size_t k = 0; k < ocm->block_count; ++k) {
    count += ocm->blocks[k].kind == (int)kind;
  }
  return count;
}

/* Block `index` of a kind; the caller's index is within the count. */
static const ocm_block* ocm_block_of(const myb_ocm* ocm, myb_ocm_kind kind, size_t index) {
  MYB_CHECK(ocm != NULL && (int)kind >= 0 && (int)kind < OCM_KINDS);
  size_t seen = 0;
  for (size_t k = 0; k < ocm->block_count; ++k) {
    if (ocm->blocks[k].kind == (int)kind && seen++ == index) {
      return &ocm->blocks[k];
    }
  }
  MYB_CHECK(0 && "no such block");
  return NULL;
}

const char* myb_ocm_value(const myb_ocm* ocm, myb_ocm_kind kind, size_t index, const char* key) {
  MYB_CHECK(key != NULL);
  const ocm_block* b = ocm_block_of(ocm, kind, index);
  return ocm_find(ocm->pairs + b->pair_first, b->pair_count, key);
}

size_t myb_ocm_line_count(const myb_ocm* ocm, myb_ocm_kind kind, size_t index) {
  return ocm_block_of(ocm, kind, index)->line_count;
}

const char* myb_ocm_line_text(const myb_ocm* ocm, myb_ocm_kind kind, size_t index, size_t line) {
  const ocm_block* b = ocm_block_of(ocm, kind, index);
  MYB_CHECK(line < b->line_count);
  return ocm->lines[b->line_first + line];
}

myb_status myb_ocm_line(const myb_ocm* ocm, myb_ocm_kind kind, size_t index, size_t line,
                        int* absolute, myb_calendar* epoch, double* offset, double* values,
                        size_t max, size_t* count) {
  MYB_CHECK(values != NULL || max == 0);
  const char* s = myb_ocm_line_text(ocm, kind, index, line);
  const size_t n = strlen(s);
  size_t at = 0;
  size_t start = 0;
  size_t length = 0;
  if (!next_word(s, n, &at, &start, &length)) {
    return MYB_E_FORMAT;
  }
  myb_calendar label = {0};
  double seconds = 0.0;
  const int dated_line = calendar_of(s + start, length, &label) == MYB_OK;
  if (!dated_line && number_of(s + start, length, &seconds) != MYB_OK) {
    return MYB_E_FORMAT;
  }
  size_t found = 0;
  /* bound: the words of the line */
  while (next_word(s, n, &at, &start, &length)) {
    double v = 0.0;
    if (number_of(s + start, length, &v) != MYB_OK) {
      return MYB_E_FORMAT;
    }
    if (found < max) {
      values[found] = v;
    }
    ++found;
  }
  if (absolute != NULL) {
    *absolute = dated_line;
  }
  if (epoch != NULL) {
    *epoch = label;
  }
  if (offset != NULL) {
    *offset = seconds;
  }
  if (count != NULL) {
    *count = found;
  }
  return MYB_OK;
}

myb_status myb_ocm_state(const myb_ocm* ocm, size_t index, size_t line, int* absolute,
                         myb_calendar* epoch, double* offset, double state[9], int* count) {
  MYB_CHECK(state != NULL);
  const char* type = myb_ocm_value(ocm, MYB_OCM_TRAJ, index, "TRAJ_TYPE");
  const int expected = type == NULL                 ? 0
                       : strcmp(type, "CARTP") == 0   ? 3
                       : strcmp(type, "CARTPV") == 0  ? 6
                       : strcmp(type, "CARTPVA") == 0 ? 9
                                                      : 0;
  if (expected == 0) {
    return MYB_E_DOMAIN;
  }
  double v[9] = {0.0};
  size_t found = 0;
  const myb_status status =
      myb_ocm_line(ocm, MYB_OCM_TRAJ, index, line, absolute, epoch, offset, v, 9, &found);
  if (status != MYB_OK) {
    return status;
  }
  if (found != (size_t)expected) {
    return MYB_E_FORMAT;
  }
  for (int i = 0; i < 9; ++i) {
    state[i] = v[i] * metres_per_km;
  }
  if (count != NULL) {
    *count = expected;
  }
  return MYB_OK;
}

static int valid_keyword(const char* s) {
  if (s == NULL || s[0] == '\0') {
    return 0;
  }
  for (size_t i = 0; s[i] != '\0'; ++i) {
    if (!keyword_char(s[i])) {
      return 0;
    }
  }
  return 1;
}

static int valid_pairs(const myb_ocm_pair* pairs, size_t count) {
  if (count > 0 && pairs == NULL) {
    return 0;
  }
  for (size_t k = 0; k < count; ++k) {
    if (!valid_keyword(pairs[k].key) || !valid_text(pairs[k].value)) {
      return 0;
    }
  }
  return 1;
}

static void put_pairs(sink* k, const myb_ocm_pair* pairs, size_t count) {
  const size_t prefix = sizeof user_prefix - 1;
  for (size_t i = 0; i < count; ++i) {
    if (k->xml && strncmp(pairs[i].key, user_prefix, prefix) == 0 &&
        pairs[i].key[prefix] != '\0') {
      myb_internal_text_put(&k->t, "<USER_DEFINED parameter=\"");
      myb_internal_text_put(&k->t, pairs[i].key + prefix);
      myb_internal_text_put(&k->t, "\">");
      put_escaped(&k->t, pairs[i].value);
      myb_internal_text_put(&k->t, "</USER_DEFINED>\n");
    } else {
      put_text(k, pairs[i].key, pairs[i].value);
    }
  }
}

myb_status myb_ocm_write(const myb_ocm_data* data, myb_odm_notation notation, char* buffer,
                         size_t size, size_t* length) {
  MYB_CHECK(data != NULL && (buffer != NULL || size == 0));
  if (!valid_text(data->originator) || !valid_calendar(&data->creation_date) ||
      !valid_pairs(data->meta, data->meta_count) ||
      (data->block_count > 0 && data->blocks == NULL) ||
      (notation != MYB_ODM_KVN && notation != MYB_ODM_XML)) {
    return MYB_E_DOMAIN;
  }
  int time_system = 0;
  int tzero = 0;
  for (size_t i = 0; i < data->meta_count; ++i) {
    time_system |= strcmp(data->meta[i].key, "TIME_SYSTEM") == 0;
    tzero |= strcmp(data->meta[i].key, "EPOCH_TZERO") == 0;
  }
  if (!time_system || !tzero) {
    return MYB_E_DOMAIN;
  }
  for (size_t b = 0; b < data->block_count; ++b) {
    const myb_ocm_block_data* block = &data->blocks[b];
    if ((int)block->kind < 0 || (int)block->kind >= OCM_KINDS ||
        !valid_pairs(block->pairs, block->pair_count) ||
        (block->line_count > 0 &&
         (block->lines == NULL || ocm_line_element[block->kind] == NULL))) {
      return MYB_E_DOMAIN;
    }
    for (size_t i = 0; i < block->line_count; ++i) {
      if (!valid_text(block->lines[i])) {
        return MYB_E_DOMAIN;
      }
    }
  }
  sink out = {.t = {.out = buffer, .size = size}, .xml = notation == MYB_ODM_XML};
  sink* k = &out;
  message(k, "ocm", "CCSDS_OCM_VERS", "3.0");
  section(k, "header", NULL);
  put_time(k, "CREATION_DATE", &data->creation_date);
  put_text(k, "ORIGINATOR", data->originator);
  section_end(k, "header", NULL);
  section(k, "body", NULL);
  section(k, "segment", NULL);
  section(k, "metadata", "META_START");
  put_pairs(k, data->meta, data->meta_count);
  section_end(k, "metadata", "META_STOP");
  section(k, "data", NULL);
  for (size_t b = 0; b < data->block_count; ++b) {
    const myb_ocm_block_data* block = &data->blocks[b];
    section(k, ocm_element[block->kind], ocm_start[block->kind]);
    put_pairs(k, block->pairs, block->pair_count);
    for (size_t i = 0; i < block->line_count; ++i) {
      if (k->xml) {
        myb_internal_text_char(&k->t, '<');
        myb_internal_text_put(&k->t, ocm_line_element[block->kind]);
        myb_internal_text_char(&k->t, '>');
        put_escaped(&k->t, block->lines[i]);
        myb_internal_text_put(&k->t, "</");
        myb_internal_text_put(&k->t, ocm_line_element[block->kind]);
        myb_internal_text_put(&k->t, ">\n");
      } else {
        myb_internal_text_put(&k->t, block->lines[i]);
        myb_internal_text_char(&k->t, '\n');
      }
    }
    section_end(k, ocm_element[block->kind], ocm_stop[block->kind]);
  }
  section_end(k, "data", NULL);
  section_end(k, "segment", NULL);
  section_end(k, "body", NULL);
  return message_end(k, "ocm", length);
}

myb_status myb_ocm_format_line(const myb_calendar* epoch, double offset, const double* values,
                               size_t count, char* buffer, size_t size, size_t* length) {
  MYB_CHECK((values != NULL || count == 0) && (buffer != NULL || size == 0));
  if ((epoch != NULL && !valid_calendar(epoch)) || (epoch == NULL && !isfinite(offset))) {
    return MYB_E_DOMAIN;
  }
  for (size_t i = 0; i < count; ++i) {
    if (!isfinite(values[i])) {
      return MYB_E_DOMAIN;
    }
  }
  myb_internal_text t = {.out = buffer, .size = size};
  if (epoch != NULL) {
    put_calendar(&t, epoch);
  } else {
    myb_internal_text_double(&t, offset);
  }
  for (size_t i = 0; i < count; ++i) {
    myb_internal_text_char(&t, ' ');
    myb_internal_text_double(&t, values[i]);
  }
  return myb_internal_text_finish(&t, length);
}
