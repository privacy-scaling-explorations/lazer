#include "test.h"

/* output bytes 480-511 for zero-length message */
const char test1[]
    = "43E41B45A653F2A5C4492C1ADD544512DDA2529833462B71A41A45BE97290B6F";

/* output bytes 480-511 for message of 200 0xA3 bytes */
const char test2[]
    = "44C9FB359FD56AC0A9A75A743CFF6862F17D7259AB075216C0699511643B6439";

int
main (void)
{
  uint8_t a[32], b[32], buf[20];
  shake128_state_t state;
  size_t len;
  int i;

  lazer_init();

  /* test 1 */

  test_hexstr2buf (a, &len, test1);
  TEST_ASSERT (len == 32);

  shake128_init (state);
  for (i = 0; i < 512; i += 32) /* discard bytes 0-479 */
    shake128_squeeze (state, b, 32);
  shake128_clear (state);

  TEST_EXPECT (memcmp (a, b, 32) == 0);

  /* test 2 */

  test_hexstr2buf (a, &len, test2);
  TEST_ASSERT (len == 32);

  memset (buf, 0xA3, 20);

  shake128_init (state);
  for (i = 0; i < 200; i += 20)
    shake128_absorb (state, buf, 20);
  for (i = 0; i < 512; i += 32) /* discard bytes 0-479 */
    shake128_squeeze (state, b, 32);
  shake128_clear (state);

  TEST_EXPECT (memcmp (a, b, 32) == 0);

  TEST_PASS ();
}
