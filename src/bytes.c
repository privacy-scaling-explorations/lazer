#include "lazer.h"

#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static inline int ishexdigit (const int);
static inline uint8_t hexdigit2bits (int);
static int bits2hexdigit (uint8_t);

void
bytes_urandom (uint8_t *bytes, const size_t len)
{
  int rv;

  ASSERT_ERR (len <= 256);

  rv = getentropy (bytes, len);
  ERR (rv != 0, "getentropy failed (size %llu).\n", (unsigned long long)len);
}

size_t
bytes_out_str (FILE *stream, const uint8_t *bytes, size_t len)
{
  const uint8_t *in = bytes;
  size_t i;
  int c;

  for (i = 0; i < len; i++)
    {
      c = bits2hexdigit (in[i] >> 4);
      if (UNLIKELY (fputc (c, stream) != c))
        return 0;

      c = bits2hexdigit (in[i] & 0xf);
      if (UNLIKELY (fputc (c, stream) != c))
        return 0;
    }
  return 2 * len;
}

size_t
bytes_inp_str (uint8_t *bytes, size_t len, FILE *stream)
{
  size_t i = 0;
  int c;

  do
    {
      c = fgetc (stream);
    }
  while (isspace (c));
  ungetc (c, stream);

  for (; i < 2 * len; i += 2)
    {
      c = fgetc (stream);
      if (UNLIKELY (!ishexdigit (c)))
        goto err;

      bytes[i / 2] = hexdigit2bits (c) << 4;

      c = fgetc (stream);
      if (UNLIKELY (!ishexdigit (c)))
        goto err;

      bytes[i / 2] += hexdigit2bits (c);
    }

  return 2 * len;
err:
  ungetc (c, stream);
  return 0;
}

size_t
bytes_out_raw (FILE *stream, const uint8_t *bytes, size_t len)
{
  return fwrite (bytes, 1, len, stream);
}

size_t
bytes_inp_raw (uint8_t *bytes, size_t len, FILE *stream)
{
  return fread (bytes, 1, len, stream);
}

void
bytes_clear (uint8_t *bytes, const size_t nbytes)
{
  explicit_bzero (bytes, nbytes);
}

static inline int
ishexdigit (const int d)
{
  return ((d >= '0' && d <= '9') || (d >= 'A' && d <= 'F')
          || (d >= 'a' && d <= 'f'));
}

static inline uint8_t
hexdigit2bits (int d)
{
  const int noff = '0' - 0;
  const int uoff = 'A' - 10;
  const int loff = 'a' - 10;

  return (d >= 'a' ? d - loff : (d >= 'A' ? d - uoff : d - noff));
}

static inline int
bits2hexdigit (uint8_t b)
{
  const char noff = '0' - 0;
  const char loff = 'a' - 10;

  return (b >= 10 ? b + loff : b + noff);
}
