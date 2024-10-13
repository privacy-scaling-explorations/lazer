#include "abdlop-params1.h"
#include "lazer.h"
#include "lnp-tbox-params1.h"
#include "test.h"

static void test0 (lnp_tbox_params_srcptr params);
static void test1 (abdlop_params_srcptr params);
static void test2 (abdlop_params_srcptr params);
static void tracemap (void);
static void isoring (void);

int
main (void)
{
  lazer_init ();

  test0 (tbox_params1);
  test1 (params1);
  test2 (params1);
  tracemap ();
  isoring ();
  TEST_PASS ();
}

static void
test0 (lnp_tbox_params_srcptr params)
{
  polyring_srcptr Rq = params->tbox->ring;
  int_srcptr q = Rq->q;
  unsigned int d = Rq->d;
  unsigned int i;
  INT_T (q2, q->nlimbs);
  INT_T (negq2, q->nlimbs);
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  INT_T (one, q->nlimbs);
  INT_T (tmp, q->nlimbs);
  uint8_t seed[32] = { 0 };
  const int64_t coeffs[]
      = { 549755813943,  -15, 549755813944,  -14, 549755813945,  -13,
          549755813946,  -12, 549755813947,  -11, 549755813948,  -10,
          549755813949,  -9,  549755813950,  -8,  549755813951,  -7,
          549755813952,  -6,  549755813953,  -5,  549755813954,  -4,
          549755813955,  -3,  549755813956,  -2,  549755813957,  -1,
          549755813958,  0,   -549755813958, 1,   -549755813957, 2,
          -549755813956, 3,   -549755813955, 4,   -549755813954, 5,
          -549755813953, 6,   -549755813952, 7,   -549755813951, 8,
          -549755813950, 9,   -549755813949, 10,  -549755813948, 11,
          -549755813947, 12,  -549755813946, 13,  -549755813945, 14,
          -549755813944, 15,  -549755813943, 16 };

  int_set_one (one);

  int_set (q2, q);
  int_sub (q2, q2, one);
  int_rshift (q2, q2, 1); // (q-1)/2
  int_neg (negq2, q2);

  poly_set_zero (a);
  poly_set (b, a);
  poly_tocrt (a);
  poly_fromcrt (a);
  TEST_EXPECT (poly_eq (a, b));

  poly_set_one (a);
  poly_set (b, a);
  poly_tocrt (a);
  poly_fromcrt (a);
  TEST_EXPECT (poly_eq (a, b));

  for (i = 0; i < d; i++)
    poly_set_coeff (a, i, q2);
  poly_set (b, a);
  poly_tocrt (a);
  poly_fromcrt (a);
  TEST_EXPECT (poly_eq (a, b));

  poly_neg_self (a);
  poly_set (b, a);
  poly_tocrt (a);
  poly_fromcrt (a);
  TEST_EXPECT (poly_eq (a, b));

  for (i = 0; i < 100; i++)
    {
      poly_urandom_bnd (a, negq2, q2, seed, 0);
      poly_set (b, a);
      poly_fromcrt (a);
      TEST_EXPECT (poly_eq (a, b));
    }

  poly_set_zero (a);
  poly_set_coeff (a, 0, q2);
  poly_set (b, a);
  poly_mul (a, a, b);
  poly_fromcrt (a);
  int_set_i64 (tmp, -274877906979);
  poly_set_zero (b);
  poly_set_coeff (b, 0, tmp);
  TEST_EXPECT (poly_eq (a, b));

  poly_set_zero (a);
  poly_set_coeff (a, 0, negq2);
  poly_set (b, a);
  poly_mul (a, a, b);
  poly_fromcrt (a);
  int_set_i64 (tmp, -274877906979);
  poly_set_zero (b);
  poly_set_coeff (b, 0, tmp);
  TEST_EXPECT (poly_eq (a, b));

  poly_set_zero (a);
  poly_set_coeff (a, 0, q2);
  poly_neg (b, a);
  poly_mul (a, a, b);
  poly_fromcrt (a);
  int_set_i64 (tmp, 274877906979);
  poly_set_zero (b);
  poly_set_coeff (b, 0, tmp);
  TEST_EXPECT (poly_eq (a, b));

  for (i = 0; i < d; i++)
    poly_set_coeff (a, i, q2);
  poly_set (b, a);
  poly_mul (a, a, b);
  poly_fromcrt (a);
  poly_set_coeffvec_i64 (b, coeffs);
  //poly_dump (a);
  //poly_dump (b);
  TEST_EXPECT (poly_eq (a, b));

  TEST_PASS ();
}

static void
test1 (abdlop_params_srcptr params)
{
  polyring_srcptr Rq = params->ring;
  unsigned int deg = polyring_get_deg (Rq);
  int_srcptr mod = polyring_get_mod (Rq);
  POLY_T (a, Rq);
  POLY_T (ar, Rq);
  POLY_T (al, Rq);
  INTVEC_T (v, deg, mod->nlimbs);
  INTVEC_T (vr, deg, mod->nlimbs);
  INT_T (s, 1);
  unsigned int i;
  long c[deg], cr[deg], cl[deg];

  TEST_ASSERT (deg >= 16);

  int_set_i64 (s, -2);

  memset (c, 0, sizeof (c));
  memset (cr, 0, sizeof (cr));
  memset (cl, 0, sizeof (cl));
  c[0] = -1;
  c[1] = 5;
  c[deg - 2] = -99;
  c[deg - 1] = 0;
  cl[2] = -1;
  cl[3] = 5;
  cl[0] = 99;
  cl[1] = 0;
  cr[deg - 2] = 1;
  cr[deg - 1] = -5;
  cr[deg - 4] = -99;
  cr[deg - 3] = 0;

  intvec_set_i64 (v, c);
  intvec_set_i64 (vr, cr);
  poly_set_coeffvec2 (a, v);
  poly_set_coeffvec (ar, vr);
  poly_set_coeffvec_i64 (al, cl);

  poly_rrot (a, a, 2);
  TEST_EXPECT (poly_eq (a, ar) == 1);

  poly_lrot (a, a, 4);
  TEST_EXPECT (poly_eq (a, al) == 1);

  poly_rrot (a, a, 2);
  poly_scale (a, s, a);
  TEST_EXPECT (int_get_i64 (intvec_get_elem (poly_get_coeffvec (a), 0)) == 2);
  TEST_EXPECT (int_get_i64 (intvec_get_elem (poly_get_coeffvec (a), 1))
               == -10);
  for (i = 2; i < deg - 2; i++)
    TEST_EXPECT (int_get_i64 (intvec_get_elem (poly_get_coeffvec (a), i))
                 == 0);
  TEST_EXPECT (int_get_i64 (intvec_get_elem (poly_get_coeffvec (a), deg - 2))
               == 198);
  TEST_EXPECT (int_get_i64 (intvec_get_elem (poly_get_coeffvec (a), deg - 1))
               == 0);
}

static void
test2 (abdlop_params_srcptr params)
{
  polyring_srcptr Rq = params->ring;
  int_srcptr mod = polyring_get_mod (Rq);
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (a, Rq);
  POLY_T (b, Rq);

  int_set_i64 (lo, -1000);
  int_set_i64 (hi, 1000);

  seed[0] = 1;
  poly_urandom_bnd (a, lo, hi, seed, dom);
  poly_set (b, a);
  poly_dump (a);
  poly_dump (b);

  poly_tocrt (a);
  poly_fromcrt (a);
  poly_dump (a);

  TEST_EXPECT (poly_eq (a, b));
}

static void
tracemap (void)
{
  int64_t coeffs[] = { 1, -1, -8, 8 };
  int64_t coeffs_tm[] = { 1, 4, 0, -4 };
  INT_T (q, 1);
  INT_T (inv2, 1);
  polyring_t Rq = { { q, 4, 0, 0, NULL, 0, NULL, NULL, inv2 } };

  int_set_i64 (q, 17);
  int_set_i64 (inv2, -8);

  POLY_T (a, Rq);
  POLY_T (a2, Rq);
  POLY_T (r, Rq);

  poly_set_coeffvec_i64 (a, coeffs_tm);
  poly_set_coeffvec_i64 (r, coeffs);
  poly_dump (r);
  poly_tracemap (r, r);
  poly_dump (r);
  poly_dump (a);
  TEST_EXPECT (poly_eq (r, a) == 1);

  poly_set_coeffvec_i64 (a, coeffs_tm);
  poly_set_coeffvec_i64 (a2, coeffs);
  poly_dump (a2);
  poly_tracemap (r, a2);
  poly_dump (r);
  poly_dump (a);
  TEST_EXPECT (poly_eq (r, a) == 1);
}

static void
isoring (void)
{
  const uint8_t seed[32] = { 0 };
  const unsigned int dprime = 256;
  const unsigned int d = 64;
  const unsigned int k = dprime / d;

  INT_T (q, 7);    /* dummy */
  INT_T (inv2, 1); /* dummy */
  polyring_t Rprime = { { q, dprime, 0, 0, NULL, 0, NULL, NULL, inv2 } };
  polyring_t R = { { q, d, 0, 0, NULL, 0, NULL, NULL, inv2 } };
  POLYVEC_T (ivec, R, k);
  POLYVEC_T (ovec, R, k);
  POLY_T (ia, Rprime);
  POLY_T (oa, Rprime);

  poly_brandom (ia, 1, seed, 0);
  poly_toisoring (ovec, ia);
  poly_fromisoring (oa, ovec);
  TEST_EXPECT (poly_eq (ia, oa) == 1);

  polyvec_brandom (ivec, 1, seed, 1);
  poly_fromisoring (oa, ivec);
  poly_toisoring (ovec, oa);
  TEST_EXPECT (polyvec_eq (ivec, ovec) == 1);
}
