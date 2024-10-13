#include "abdlop-params1.h"
#include "test.h"

static void test (abdlop_params_srcptr params);

int
main (void)
{
  lazer_init();

  test (params1);

  TEST_PASS ();
}

static void
test (abdlop_params_srcptr params)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  const unsigned int deg = params->ring->d;
  const unsigned int nlimbs = params->ring->q->nlimbs;
  INTVEC_T (r, deg, nlimbs);
  INTVEC_T (r0, deg, nlimbs);
  INTVEC_T (r1, deg, nlimbs);
  INTVEC_T (r1prime, deg, nlimbs);
  INTVEC_T (z, deg, nlimbs);
  INTVEC_T (y, deg, nlimbs);
  INTVEC_T (yprime, deg, nlimbs);
  INT_T (qminus1, nlimbs);
  INT_T (negqminus1, nlimbs);
  INT_T (gammaby2, nlimbs);
  INT_T (neggammaby2, nlimbs);
  unsigned int i;

  int_set (qminus1, params->dcompress->qminus1);
  int_neg (negqminus1, qminus1);
  int_set (gammaby2, params->dcompress->gammaby2);
  int_neg (neggammaby2, gammaby2);

  for (i = 0; i < 4000; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      dom = 1;
      intvec_urandom_bnd (r, negqminus1, qminus1, seed, dom);
      dom = 2;
      intvec_urandom_bnd (z, neggammaby2, gammaby2, seed, dom);

      dcompress_make_ghint (y, z, r, params->dcompress);
      dcompress_use_ghint (r1prime, y, r, params->dcompress);

      intvec_add (r, r, z);
      intvec_mod (r, r, params->ring->q);
      dcompress_decompose (r1, r0, r, params->dcompress);

      TEST_EXPECT (intvec_eq (r1prime, r1) == 1);
    }
}
