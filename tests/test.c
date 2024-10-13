#include "test.h"

static int ishexdigit (const char);
static unsigned char hexdigit2byte (char);
static char byte2hexdigit (unsigned char);

void
test_hexstr2buf (uint8_t *buf, size_t *buflen, const char *hexstr)
{
  const char *ptr;
  size_t len, i;

  TEST_ASSERT (buf != NULL && buflen != NULL && hexstr != NULL
               && strlen (hexstr) % 2 == 0);

  *buflen = 0;
  ptr = hexstr;

  /* Skip possible leading '0x'. */
  if (strlen (ptr) > 2 && ptr[0] == '0' && ptr[1] == 'x')
    ptr += 2;
  len = strlen (ptr);

  for (i = 0; i + 1 < len; i += 2)
    {
      TEST_ASSERT (ishexdigit (ptr[i]) && ishexdigit (ptr[i + 1]));

      buf[i / 2] = hexdigit2byte (ptr[i]) << 4;
      buf[i / 2] += hexdigit2byte (ptr[i + 1]);
      (*buflen)++;
    }
}

void
test_buf2hexstr (char *hexstr, const uint8_t *buf, size_t buflen)
{
  size_t i;

  TEST_ASSERT (hexstr != NULL && buf != NULL);

  for (i = 0; i < buflen; i++)
    {
      hexstr[2 * i] = byte2hexdigit (buf[i] >> 4);
      hexstr[2 * i + 1] = byte2hexdigit (buf[i] & 0xf);
    }
  hexstr[2 * i] = '\0';
}

static int
ishexdigit (const char d)
{
  return ((d >= '0' && d <= '9') || (d >= 'A' && d <= 'F')
          || (d >= 'a' && d <= 'f'));
}

static unsigned char
hexdigit2byte (char d)
{
  const char noff = '0' - 0;
  const char uoff = 'A' - 10;
  const char loff = 'a' - 10;

  return (d >= 'a' ? d - loff : (d >= 'A' ? d - uoff : d - noff));
}

static char
byte2hexdigit (unsigned char b)
{
  const char noff = '0' - 0;
  const char loff = 'a' - 10;

  TEST_ASSERT ((b & 0xf0) == 0);

  return (b >= 10 ? b + loff : b + noff);
}
