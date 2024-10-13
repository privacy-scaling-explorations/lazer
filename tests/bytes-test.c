#include "test.h"

#define BUFLEN 256

static void zero_bufs (void);

const uint8_t raw[] = { 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0 };
const char str[] = "123456789abcdef0";

uint8_t r1[BUFLEN], r2[BUFLEN];
char s1[BUFLEN], s2[BUFLEN];

int
main (void)
{
  size_t len;
  FILE *f;

  lazer_init();

  bytes_urandom (r1, BUFLEN);
  bytes_clear (r1, BUFLEN);

  /* raw to str */
  zero_bufs ();

  test_buf2hexstr (s1, raw, sizeof (raw));

  f = fmemopen (s2, sizeof (s2), "w");
  TEST_ASSERT (f != NULL);
  bytes_out_str (f, raw, sizeof (raw));
  fclose (f);

  TEST_EXPECT (strcmp (s1, str) == 0);
  TEST_EXPECT (strcmp (s2, str) == 0);

  /* str to raw */
  zero_bufs ();

  test_hexstr2buf (r1, &len, str);
  TEST_EXPECT (len == sizeof (raw));

  f = fmemopen ((void *)str, sizeof (str), "r");
  TEST_ASSERT (f != NULL);
  bytes_inp_str (r2, sizeof (raw), f);
  fclose (f);

  TEST_EXPECT (memcmp (r1, raw, sizeof (raw)) == 0);
  TEST_EXPECT (memcmp (r2, raw, sizeof (raw)) == 0);

  /* raw to raw write */
  zero_bufs ();

  f = fmemopen (r1, sizeof (r1), "w");
  TEST_ASSERT (f != NULL);
  bytes_out_raw (f, raw, sizeof (raw));
  fclose (f);

  TEST_EXPECT (memcmp (r1, raw, sizeof (raw)) == 0);

  /* raw to raw read */
  zero_bufs ();

  f = fmemopen ((void *)raw, sizeof (raw), "r");
  TEST_ASSERT (f != NULL);
  bytes_inp_raw (r1, sizeof (raw), f);
  fclose (f);

  TEST_EXPECT (memcmp (r1, raw, sizeof (raw)) == 0);

  TEST_PASS ();
}

static void
zero_bufs (void)
{
  memset (r1, 0, sizeof (r1));
  memset (r2, 0, sizeof (r2));
  memset (s1, 0, sizeof (s1));
  memset (s2, 0, sizeof (s2));
}
