#include "abdlop-params1.h"
#include "test.h"

static void misc (void);
static void isoring (void);

int
main (void)
{
  lazer_init ();

  misc ();
  isoring ();
  TEST_PASS ();
}

static void
misc (void)
{
  INT_T (z, 1);
  INT_T (z2, 2);
  INT_T (q, 1);
  INT_T (nq, 1);
  INT_T (inv2, 1);
  int_set_i64 (q, 17);
  int_set_i64 (nq, -17);
  int_set_i64 (inv2, -8);

  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  polyring_t rq = { { q, 4, 0, 0, NULL, 0, NULL, NULL, inv2 } };
  const int64_t ci[]
      = { -2, 0, 3, 8, 10, -3, 1, 7, 0, -5, 5, 6, 3, 2, -1, -1 };
  const int32_t ci32[]
      = { -2, 0, 3, 8, 10, -3, 1, 7, 0, -5, 5, 6, 3, 2, -1, -1 };
  const int64_t cisubv[] = { 10, -3, 1, 7, 3, 2, -1, -1 };
  int64_t i[16];
  int32_t i32[16];
  POLYVEC_T (tmp, rq, 4);
  POLYVEC_T (v, rq, 4);
  POLYVEC_T (v2, rq, 4);
  POLYVEC_T (v3, rq, 4);
  POLYVEC_T (vsub, rq, 2);
  POLY_T (cp0, rq);
  POLY_T (cp1, rq);
  poly_ptr pp;
  polyvec_t subv;

  poly_set_coeffvec_i64 (cp0, &ci[0]);
  poly_set_coeffvec_i64 (cp1, &ci[4]);

  polyvec_set_elem (v, 1, cp1);
  pp = polyvec_get_elem (v, 1);
  TEST_EXPECT (poly_eq (pp, cp1) == 1);
  TEST_EXPECT (poly_eq (pp, cp0) == 0);

  polyvec_set_coeffvec_i64 (v, ci);
  polyvec_get_coeffvec_i64 (i, v);
  TEST_EXPECT (memcmp (ci, i, sizeof (ci)) == 0);

  polyvec_set_coeffvec_i32 (v, ci32);
  polyvec_get_coeffvec_i32 (i32, v);
  TEST_EXPECT (memcmp (ci32, i32, sizeof (ci32)) == 0);

  polyvec_set (v2, v);
  TEST_EXPECT (polyvec_eq (v2, v) == 1);

  polyvec_set_coeffvec_i64 (vsub, cisubv);
  polyvec_get_subvec (subv, v, 1, 2, 2);
  polyvec_dump (subv);
  polyvec_dump (vsub);
  TEST_EXPECT (polyvec_eq (subv, vsub) == 1);

  polyvec_rshift (v2, v, 1);
  polyvec_set_coeffvec_i32 (
      tmp, (int32_t[]){ -1, 0, 1, 4, 5, -1, 0, 3, 0, -2, 2, 3, 1, 1, 0, 0 });
  TEST_EXPECT (polyvec_eq (v2, tmp) == 1);

  polyvec_lshift (v2, v2, 1);
  polyvec_set_coeffvec_i32 (
      tmp, (int32_t[]){ -2, 0, 2, 8, 10, -2, 0, 6, 0, -4, 4, 6, 2, 2, 0, 0 });
  TEST_EXPECT (polyvec_eq (v2, tmp) == 1);

  polyvec_rrot (v2, v2, 1);
  polyvec_set_coeffvec_i32 (
      tmp, (int32_t[]){ 0, 2, 8, 2, -2, 0, 6, -10, -4, 4, 6, 0, 2, 0, 0, -2 });
  TEST_EXPECT (polyvec_eq (v2, tmp) == 1);

  polyvec_lrot (v2, v, 1);
  polyvec_set_coeffvec_i32 (tmp, (int32_t[]){ -8, -2, 0, 3, -7, 10, -3, 1, -6,
                                              0, -5, 5, 1, 3, 2, -1 });
  TEST_EXPECT (polyvec_eq (v2, tmp) == 1);

  polyvec_add (v2, v, v, 0);
  polyvec_scale (v2, q, v2);
  polyvec_sub (v2, v, v2, 0);
  polyvec_mod (v2, v2);
  polyvec_redc (v2, v2);
  polyvec_set_coeffvec_i64 (tmp, (int64_t[]){ -2, 0, 3, 8, -7, -3, 1, 7, 0, -5,
                                              5, 6, 3, 2, -1, -1 });
  TEST_EXPECT (polyvec_eq (v2, tmp) == 1);

  polyvec_urandom (v2, q, 5, seed, dom);
  polyvec_urandom_bnd (v2, nq, q, seed, dom);

  polyvec_urandom_autostable (v2, 2, 3, seed, dom);
  polyvec_auto (v3, v2);
  TEST_EXPECT (polyvec_eq (v2, v3) == 1);

  polyvec_set (v2, v); /* v is unstable under auto */
  polyvec_auto_self (v2);
  TEST_EXPECT (polyvec_eq (v2, v) == 0);
  polyvec_auto_self (v2);
  TEST_EXPECT (polyvec_eq (v2, v) == 1);

  polyvec_redc (v3, v);
  polyvec_auto (v2, v3);
  polyvec_add (v2, v, v2, 0);
  polyvec_scale (v2, inv2, v2);
  poly_tracemap (cp0, polyvec_get_elem (v3, 1));
  TEST_EXPECT (poly_eq (cp0, polyvec_get_elem (v2, 1)) == 1);

  polyvec_linf (z, v);
  TEST_EXPECT (int_get_i64 (z) == 10);

  polyvec_l2sqr (z2, v);
  TEST_EXPECT (int_get_i64 (z2) == 337);
}

static void
isoring (void)
{
  const unsigned int nelems = 10;
  const uint8_t seed[32] = { 0 };
  const unsigned int dprime = 256;
  const unsigned int d = 64;
  const unsigned int k = dprime / d;

  INT_T (q, 7);    /* dummy */
  INT_T (inv2, 1); /* dummy */
  polyring_t Rprime = { { q, dprime, 0, 0, NULL, 0, NULL, NULL, inv2 } };
  polyring_t R = { { q, d, 0, 0, NULL, 0, NULL, NULL, inv2 } };
  POLYVEC_T (ivec, R, k * nelems);
  POLYVEC_T (ovec, R, k * nelems);
  POLYVEC_T (ia, Rprime, nelems);
  POLYVEC_T (oa, Rprime, nelems);

  polyvec_brandom (ia, 1, seed, 0);
  polyvec_toisoring (ovec, ia);
  polyvec_fromisoring (oa, ovec);
  TEST_EXPECT (polyvec_eq (ia, oa) == 1);

  polyvec_brandom (ivec, 1, seed, 1);
  polyvec_fromisoring (oa, ivec);
  polyvec_toisoring (ovec, oa);
  TEST_EXPECT (polyvec_eq (ivec, ovec) == 1);
}
