#include "test.h"

#if RNG == RNG_SHAKE128RNG
const char seedstr[]
    = "22c1fd7342356b0a1a0ef75e7346c2df8a76148407f7f1132e47ed9d59ae4147";
const char domstr[] = "a6211d51ed664050";
const char streamstr[] = "f99123f410a594dda2238d0007ec8d01";

int
main (void)
{
  uint8_t seed[32], out[16], stream[16];
  uint64_t dom = 0;
  rng_state_t state;
  size_t len;

  lazer_init();

  test_hexstr2buf (seed, &len, seedstr);
  TEST_ASSERT (len == 32);

  test_hexstr2buf ((uint8_t *)&dom, &len, domstr);
  TEST_ASSERT (len == 8);
  dom = le64toh (dom);

  test_hexstr2buf (stream, &len, streamstr);
  TEST_ASSERT (len == 16);

  rng_init (state, seed, dom);

  rng_urandom (state, out, 3);
  rng_urandom (state, out + 3, 10);
  rng_urandom (state, out + 3 + 10, 16 - 3 - 10);
  TEST_EXPECT (memcmp (out, stream, 16) == 0);

  rng_clear (state);
  TEST_PASS ();
}

#elif RNG == RNG_AES256CTR
const char seedstr[]
    = "FF7A617CE69148E4F1726E2F43581DE2AA62D9F805532EDFF1EED687FB54153D";
const char domstr[] = "001CC5B751A51D70";
const char streamstr[] = "913cd4d68a9feed715e3bd37489e266f8a3c490cefe47e14bbde"
                         "6ade9317f9619c99e38a";

int
main (void)
{
  uint8_t seed[32], out[36], stream[36];
  uint64_t dom = 0;
  rng_state_t state;
  size_t len;

  lazer_init();

  test_hexstr2buf (seed, &len, seedstr);
  TEST_ASSERT (len == 32);

  test_hexstr2buf ((uint8_t *)&dom, &len, domstr);
  TEST_ASSERT (len == 8);
  dom = le64toh (dom);

  test_hexstr2buf (stream, &len, streamstr);
  TEST_ASSERT (len == 36);

  rng_init (state, seed, dom);

  rng_urandom (state, out, 10);
  rng_urandom (state, out + 10, 10);
  rng_urandom (state, out + 20, 10);
  rng_urandom (state, out + 30, 6);

  TEST_EXPECT (memcmp (out, stream, sizeof (stream)) == 0);

  rng_clear (state);
  TEST_PASS ();
}

#else
#error "Invalid rng option."
#endif
