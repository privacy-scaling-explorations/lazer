#include "test.h"

#define LIMB0(x) ((x)->limbs[0])
#define LIMB1(x) ((x)->limbs[1])
#define LIMB2(x) ((x)->limbs[2])
#define LIMB3(x) ((x)->limbs[3])
#define IS_NEG(x) ((x)->neg)

#define IT 400000

static void encode (void);
static void encode_uniform (void);
static void encode_gaussian (void);
static void encode_ghint (void);

int
main (void)
{
  lazer_init();

  encode ();
  encode_uniform ();
  encode_gaussian ();
  encode_ghint ();

  TEST_PASS ();
}

static void
encode_uniform (void)
{
  unsigned int nbits1, nbits2, nbits3, nbits, i;
  CODER_STATE_T (state);
  uint8_t buf[100000];
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  INTVEC_T (v1, 4, 1);
  INTVEC_T (v2, 4, 1);
  INTVEC_T (v3, 4, 1);
  int rc;
  INT_T (mod, 1);
  unsigned int mbits;

  int_set_i64 (mod, 17);
  mbits = 5;

  bytes_urandom (seed, sizeof (seed));

  for (i = 0; i < IT; i++)
    {
      intvec_urandom (v1, mod, mbits, seed, dom);
      dom++;

      coder_enc_begin (state, buf);

      coder_enc_urandom (state, v1, mod, mbits);
      nbits1 = coder_get_offset (state);
      coder_enc_urandom (state, v1, mod, mbits);
      nbits2 = coder_get_offset (state);

      coder_enc_end (state);
      nbits3 = coder_get_offset (state);

      coder_dec_begin (state, buf);

      coder_dec_urandom (state, v2, mod, mbits);
      nbits = coder_get_offset (state);
      TEST_EXPECT (nbits == nbits1);
      TEST_EXPECT (intvec_eq (v1, v2));
      coder_dec_urandom (state, v3, mod, mbits);
      nbits = coder_get_offset (state);
      TEST_EXPECT (nbits == nbits2);
      TEST_EXPECT (intvec_eq (v1, v3));
      rc = coder_dec_end (state);
      TEST_EXPECT (rc == 1);
      nbits3 = coder_get_offset (state);
      TEST_EXPECT (nbits3 == nbits3);
    }
}

static void
encode_gaussian (void)
{
  const unsigned int log2o = 10;
  unsigned int nbits1, nbits2, nbits3, nbits, i;
  CODER_STATE_T (state);
  uint8_t buf[100000];
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  INTVEC_T (v1, 4, 1);
  INTVEC_T (v2, 4, 1);
  INTVEC_T (v3, 4, 1);
  int rc;

  bytes_urandom (seed, sizeof (seed));

  for (i = 0; i < IT; i++)
    {
      intvec_grandom (v1, log2o, seed, dom);
      dom++;

      coder_enc_begin (state, buf);

      coder_enc_grandom (state, v1, log2o);
      nbits1 = coder_get_offset (state);
      coder_enc_grandom (state, v1, log2o);
      nbits2 = coder_get_offset (state);

      coder_enc_end (state);
      nbits3 = coder_get_offset (state);

      coder_dec_begin (state, buf);

      coder_dec_grandom (state, v2, log2o);
      nbits = coder_get_offset (state);
      TEST_EXPECT (nbits == nbits1);
      TEST_EXPECT (intvec_eq (v1, v2));
      coder_dec_grandom (state, v3, log2o);
      nbits = coder_get_offset (state);
      TEST_EXPECT (nbits == nbits2);
      TEST_EXPECT (intvec_eq (v1, v3));
      rc = coder_dec_end (state);
      TEST_EXPECT (rc == 1);
      nbits3 = coder_get_offset (state);
      TEST_EXPECT (nbits3 == nbits3);
    }
}

static void
encode_ghint (void)
{
  unsigned int nbits1, nbits2, nbits3, nbits, i;
  CODER_STATE_T (state);
  uint8_t buf[100000];
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  INTVEC_T (v1, 4, 1);
  INTVEC_T (v2, 4, 1);
  INTVEC_T (v3, 4, 1);
  int rc;
  INT_T (mod, 1);
  unsigned int mbits;

  int_set_i64 (mod, 17);
  mbits = 5;

  bytes_urandom (seed, sizeof (seed));

  for (i = 0; i < IT; i++)
    {
      intvec_urandom (v1, mod, mbits, seed, dom);
      dom++;

      coder_enc_begin (state, buf);

      coder_enc_ghint (state, v1);
      nbits1 = coder_get_offset (state);
      coder_enc_ghint (state, v1);
      nbits2 = coder_get_offset (state);

      coder_enc_end (state);
      nbits3 = coder_get_offset (state);

      coder_dec_begin (state, buf);

      coder_dec_ghint (state, v2);
      nbits = coder_get_offset (state);
      TEST_EXPECT (nbits == nbits1);
      TEST_EXPECT (intvec_eq (v1, v2));
      coder_dec_ghint (state, v3);
      nbits = coder_get_offset (state);
      TEST_EXPECT (nbits == nbits2);
      TEST_EXPECT (intvec_eq (v1, v3));
      rc = coder_dec_end (state);
      TEST_EXPECT (rc == 1);
      nbits3 = coder_get_offset (state);
      TEST_EXPECT (nbits3 == nbits3);
    }
}

static void
encode (void)
{
  CODER_STATE_T (state);
  const char encv1_1[]
      = "badda91889d2fe3776eb7f7d9dffe94bed37924d947a8102e900";
  const char encv1_2[] = "badda91889d2fe3776eb7f7d9dffe94bed37924d947a8102e974"
                         "bb533112a5fd6fecd6fffa3affd397da6f249b28f50205d201";
  const char encv1_3[] = "badda91889d2fe3776eb7f7d9dffe94bed37924d947a8102e974"
                         "bb533112a5fd6fecd6fffa3affd397da6f249b28f50205d2e976"
                         "a762244afbdfd8adfff575fea72fb5df483651ea050aa403";
  const char encv1_4[]
      = "badda91889d2fe3776eb7f7d9dffe94bed37924d947a8102e974bb533112a5fd6fecd"
        "6fffa3affd397da6f249b28f50205d2e976a762244afbdfd8adfff575fea72fb5df48"
        "3651ea050aa4c306000e0002";
  const char encv1_5[]
      = "badda91889d2fe3776eb7f7d9dffe94bed37924d947a8102e974bb533112a5fd6fecd"
        "6fffa3affd397da6f249b28f50205d2e976a762244afbdfd8adfff575fea72fb5df48"
        "3651ea050aa4c306000e00620300070001";
  const char encv1_6[]
      = "badda91889d2fe3776eb7f7d9dffe94bed37924d947a8102e974bb533112a5fd6fecd"
        "6fffa3affd397da6f249b28f50205d2e976a762244afbdfd8adfff575fea72fb5df48"
        "3651ea050aa4c306000e00620300070003";
  const int64_t ghint[] = { 0, 1, -1, 9, -8 };
  size_t nbits, len;
  uint8_t buf[100000], buf2[100000], bytes3[3], tmp[3];
  int rc;

  INT_T (zero, 3);
  INT_T (m, 3);

  /* range 6 * 2^64 + 2^64 - 2 */
  IS_NEG (m) = 0;
  LIMB0 (m) = (~(limb_t)0) - 1;
  LIMB1 (m) = 6;
  LIMB2 (m) = 0;

  LIMB0 (m) += 1; /* modulus = range + 1 */

  int_set_i64 (zero, 0);

  memset (buf, 0xff, sizeof (buf));
  bytes_urandom (bytes3, sizeof (bytes3));

  INT_T (v1_0, 3);
  IS_NEG (v1_0) = 0;
  LIMB0 (v1_0) = 4034893802436681146ULL;
  LIMB1 (v1_0) = 6ULL;
  LIMB2 (v1_0) = 0;
  INT_T (v1_1, 3);
  IS_NEG (v1_1) = 0;
  LIMB0 (v1_1) = 12212988080355802478ULL;
  LIMB1 (v1_1) = 5ULL;
  LIMB2 (v1_1) = 0;
  INT_T (v1_2, 3);
  IS_NEG (v1_2) = 0;
  LIMB0 (v1_2) = 11820266675930286303ULL;
  LIMB1 (v1_2) = 3ULL;
  LIMB2 (v1_2) = 0;

  INTVEC_T (v1, 3, 3);
  INTVEC_T (v2, 3, 3);
  INTVEC_T (v3, 3, 3);
  INTVEC_T (v4, 3, 3);
  INTVEC_T (v5, 5, 1);
  INTVEC_T (v5_, 5, 1);

  intvec_set_elem (v1, 0, v1_0);
  intvec_set_elem (v1, 1, v1_1);
  intvec_set_elem (v1, 2, v1_2);
  intvec_set_i64 (v5, ghint);

  coder_enc_begin (state, buf);

  coder_enc_bytes (state, bytes3, sizeof (bytes3));

  coder_enc_urandom (state, v1, m, 67);
  nbits = coder_get_offset (state);

  TEST_EXPECT (nbits == 3 * 8 + 3 * 67);
  test_hexstr2buf (buf2, &len, encv1_1);
  TEST_EXPECT (memcmp (buf + 3, buf2, len) == 0);

  coder_enc_urandom (state, v1, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (nbits == 3 * 8 + 6 * 67);
  test_hexstr2buf (buf2, &len, encv1_2);
  TEST_EXPECT (memcmp (buf + 3, buf2, len) == 0);

  coder_enc_urandom (state, v1, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (nbits == 3 * 8 + 9 * 67);
  test_hexstr2buf (buf2, &len, encv1_3);
  TEST_EXPECT (memcmp (buf + 3, buf2, len) == 0);

  coder_enc_ghint (state, v5);
  nbits = coder_get_offset (state);
  TEST_EXPECT (nbits == 3 * 8 + 9 * 67 + 3 * 2 + 3 + 14 + 3 + 13);
  test_hexstr2buf (buf2, &len, encv1_4);
  TEST_EXPECT (memcmp (buf + 3, buf2, len) == 0);

  coder_enc_ghint (state, v5);
  nbits = coder_get_offset (state);
  TEST_EXPECT (nbits
               == 3 * 8 + 9 * 67 + 3 * 2 + 3 + 14 + 3 + 13 + 3 * 2 + 3 + 14 + 3
                      + 13);
  test_hexstr2buf (buf2, &len, encv1_5);
  TEST_EXPECT (memcmp (buf + 3, buf2, len) == 0);

  coder_enc_end (state);
  nbits = coder_get_offset (state);
  TEST_EXPECT (nbits == 89 * 8);
  test_hexstr2buf (buf2, &len, encv1_6);
  TEST_EXPECT (memcmp (buf + 3, buf2, len) == 0);


  test_hexstr2buf (buf2, &len, encv1_3);
  coder_dec_begin (state, buf2);

  rc = coder_dec_urandom (state, v2, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 3 * 67);
  TEST_EXPECT (intvec_eq (v1, v2) == 1);

  rc = coder_dec_urandom (state, v3, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 6 * 67);
  TEST_EXPECT (intvec_eq (v3, v1) == 1);

  rc = coder_dec_urandom (state, v4, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 9 * 67);
  TEST_EXPECT (intvec_eq (v4, v1) == 1);

  rc = coder_dec_end (state);
  TEST_EXPECT (rc == 0);

  coder_dec_begin (state, buf);

  rc = coder_dec_bytes (state, tmp, sizeof (bytes3));
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 3 * 8);
  TEST_EXPECT (memcmp (tmp, bytes3, sizeof (bytes3)) == 0);

  rc = coder_dec_urandom (state, v2, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 3 * 8 + 3 * 67);
  TEST_EXPECT (intvec_eq (v1, v2) == 1);

  rc = coder_dec_urandom (state, v3, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 3 * 8 + 6 * 67);
  TEST_EXPECT (intvec_eq (v3, v1) == 1);

  rc = coder_dec_urandom (state, v4, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 3 * 8 + 9 * 67);
  TEST_EXPECT (intvec_eq (v4, v1) == 1);

  coder_dec_ghint (state, v5_);
  TEST_EXPECT (intvec_eq (v5_, v5) == 1);

  coder_dec_ghint (state, v5_);
  TEST_EXPECT (intvec_eq (v5_, v5) == 1);

  rc = coder_dec_end (state);
  TEST_EXPECT (rc == 1);
  nbits = coder_get_offset (state);
  TEST_EXPECT (nbits == 89 * 8);

  test_hexstr2buf (buf2, &len, encv1_3);
  coder_dec_begin (state, buf2);

  rc = coder_dec_urandom (state, v2, m, 67);
  nbits = coder_get_offset (state);
  TEST_EXPECT (rc == 0);
  TEST_EXPECT (nbits == 3 * 67);
  TEST_EXPECT (intvec_eq (v1, v2) == 1);

  rc = coder_dec_end (state);
  TEST_EXPECT (rc == 0);
}
