#ifndef TEST_H
#define TEST_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lazer.h"

/* automake test exit status */
#define TEST_RC_PASS 0
#define TEST_RC_SKIP 77
#define TEST_RC_FAIL 1
#define TEST_RC_ERROR 99

#define TEST_PASS() exit (TEST_RC_PASS)

#define TEST_SKIP_IF(cond)                                                    \
  do                                                                          \
    {                                                                         \
      if ((cond))                                                             \
        {                                                                     \
          fprintf (stderr, "[SKIP] (%s:%d).\n", __FILE__, __LINE__);          \
          exit (TEST_RC_SKIP);                                                \
        }                                                                     \
    }                                                                         \
  while (0)

#define TEST_EXPECT(cond)                                                     \
  do                                                                          \
    {                                                                         \
      if (!(cond))                                                            \
        {                                                                     \
          fprintf (stderr, "[FAIL] (%s:%d).\n", __FILE__, __LINE__);          \
          exit (TEST_RC_FAIL);                                                \
        }                                                                     \
    }                                                                         \
  while (0)

#define TEST_ASSERT(cond)                                                     \
  do                                                                          \
    {                                                                         \
      if (!(cond))                                                            \
        {                                                                     \
          fprintf (stderr, "[ERROR] (%s:%d).\n", __FILE__, __LINE__);         \
          exit (TEST_RC_ERROR);                                               \
        }                                                                     \
    }                                                                         \
  while (0)

void test_hexstr2buf (uint8_t *, size_t *, const char *);
void test_buf2hexstr (char *hexstr, const uint8_t *, size_t);

#endif
