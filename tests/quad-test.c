#include "test.h"

static void lin (void);
static void lin2 (void);
static void quad (void);

int
main (void)
{
  lazer_init ();

  lin ();
  lin2 ();
  // quad ();
  TEST_PASS ();
}

static void
lin (void)
{
  const int64_t s0coeffs[]
      = { 1, 6, 4, 2, 5, 0, 6, 3, 0, 5, 3, 1, 1, 1, 4, 4, 1, 0, 6, 0, 0, 3,
          1, 3, 4, 1, 5, 3, 5, 2, 3, 6, 5, 6, 2, 4, 1, 4, 1, 3, 1, 3, 0, 6,
          4, 2, 6, 6, 4, 0, 6, 0, 3, 6, 5, 0, 5, 4, 2, 5, 1, 1, 2, 0, 0, 4,
          0, 1, 1, 4, 1, 1, 1, 0, 2, 0, 4, 5, 0, 6, 2, 1, 3, 3, 3, 5, 2, 6,
          1, 6, 1, 4, 0, 2, 4, 0, 4, 1, 5, 5, 4, 0, 5, 1, 1, 1, 2, 6, 4, 3,
          4, 2, 2, 3, 4, 5, 4, 3, 2, 3, 1, 1, 5, 5, 1, 3, 4, 3 };
  const int64_t s1coeffs[]
      = { 6, 3, 4, 5, 2, 3, 5, 2, 3, 6, 6, 3, 5, 2, 1, 4, 2, 6, 1, 0, 2, 2,
          1, 6, 2, 6, 4, 5, 5, 0, 5, 3, 2, 5, 1, 1, 0, 5, 1, 5, 6, 3, 1, 3,
          6, 4, 2, 4, 2, 5, 0, 3, 3, 1, 0, 0, 3, 3, 4, 5, 4, 3, 2, 2, 2, 2,
          0, 3, 2, 4, 6, 0, 6, 3, 6, 0, 6, 0, 2, 6, 1, 2, 6, 6, 5, 0, 6, 2,
          4, 4, 5, 2, 6, 5, 5, 6, 6, 5, 6, 4, 4, 1, 3, 2, 6, 2, 4, 5, 0, 1,
          2, 6, 3, 5, 1, 6, 4, 0, 3, 4, 5, 4, 4, 3, 5, 5, 4, 3 };
  const int64_t r10coeffs[]
      = { 1, 4, 0, 6, 1, 4, 5, 2, 5, 5, 6, 2, 3, 5, 4, 0, 5, 4, 4, 5, 3, 6,
          2, 1, 2, 4, 4, 2, 4, 2, 1, 3, 5, 6, 1, 5, 2, 3, 2, 0, 6, 3, 5, 1,
          4, 0, 5, 1, 4, 2, 4, 5, 6, 5, 0, 4, 1, 1, 3, 2, 1, 5, 1, 4, 5, 2,
          2, 4, 1, 3, 3, 5, 0, 3, 3, 0, 5, 2, 2, 4, 1, 0, 2, 4, 2, 1, 1, 4,
          6, 5, 1, 2, 3, 2, 1, 1, 0, 4, 2, 1, 4, 4, 0, 1, 4, 3, 1, 5, 4, 5,
          6, 4, 5, 0, 1, 6, 2, 0, 6, 3, 3, 0, 1, 5, 4, 5, 1, 4 };
  const int64_t r11coeffs[]
      = { 3, 1, 6, 0, 5, 0, 4, 6, 3, 6, 3, 5, 6, 6, 5, 3, 5, 6, 0, 0, 2, 4,
          4, 3, 2, 0, 2, 6, 4, 4, 6, 2, 5, 1, 5, 1, 2, 4, 0, 0, 3, 5, 4, 4,
          6, 1, 4, 0, 3, 1, 0, 4, 2, 2, 1, 3, 0, 1, 2, 4, 4, 1, 6, 2, 4, 4,
          6, 0, 2, 5, 4, 1, 1, 3, 3, 0, 3, 5, 3, 0, 6, 1, 3, 2, 2, 0, 4, 4,
          2, 0, 1, 0, 3, 2, 4, 4, 2, 2, 0, 5, 1, 5, 1, 5, 1, 1, 4, 2, 4, 3,
          5, 5, 1, 3, 0, 5, 3, 6, 1, 1, 5, 3, 5, 1, 2, 5, 1, 3 };
  const int64_t r0coeffs[]
      = { 5, 1, 1, 2, 1, 3, 0, 0, 3, 5, 0, 0, 1, 0, 5, 3, 2, 4, 2, 3, 0, 0,
          6, 5, 5, 4, 2, 2, 2, 0, 1, 1, 3, 3, 3, 1, 1, 5, 5, 2, 3, 4, 4, 0,
          1, 4, 5, 4, 5, 0, 4, 1, 0, 2, 5, 3, 4, 2, 2, 6, 4, 4, 0, 2, 4, 2,
          2, 3, 3, 3, 1, 3, 4, 3, 6, 0, 1, 5, 2, 1, 0, 2, 0, 3, 0, 1, 1, 6,
          6, 6, 3, 2, 6, 1, 1, 3, 2, 6, 0, 3, 5, 3, 5, 5, 2, 6, 6, 0, 0, 0,
          5, 6, 4, 4, 0, 2, 1, 2, 4, 5, 1, 5, 5, 5, 2, 2, 0, 1 };
  const unsigned int dprime = 128;
  const unsigned int d = 64;
  const unsigned int k = dprime / d;
  const unsigned int n = 2;
  unsigned int i;
  poly_ptr poly;

  INT_T (q, 1);
  int_set_i64 (q, 7);
  INT_T (inv2, 1);
  int_set_i64 (inv2, -3);

  INT_T (Pprime, 1);
  int_set_i64 (Pprime, 1125899906840833);
  INT_T (P, 1);
  int_set_i64 (P, 1125899906840833);
  INT_T (Ppmodq_, 1);
  int_set_i64 (Ppmodq_, 1);
  int_srcptr Ppmodq[] = { Ppmodq_ };

  polyring_t Rprime
      = { { q, dprime, 3, 7, moduli_d128, 1, Pprime, Ppmodq, inv2 } };
  polyring_t R = { { q, d, 3, 6, moduli_d64, 1, P, Ppmodq, inv2 } };

  POLYVEC_T (sprime, Rprime, n);
  POLYVEC_T (r1prime, Rprime, n);
  POLY_T (r0prime, Rprime);
  POLYVEC_T (s, R, k * n);
  POLYVEC_T (r10, R, k * n);
  POLYVEC_T (r11, R, k * n);
  POLY_T (r00, R);
  POLY_T (r01, R);
  POLY_T (zero, R);
  polyvec_ptr r1[n];
  poly_ptr r0[n];

  poly_set_zero (zero);

  r1[0] = r10;
  r1[1] = r11;
  r0[0] = r00;
  r0[1] = r01;

  /* sprime */
  poly = polyvec_get_elem (sprime, 0);
  poly_set_coeffvec_i64 (poly, s0coeffs);
  poly = polyvec_get_elem (sprime, 1);
  poly_set_coeffvec_i64 (poly, s1coeffs);

  /* r1prime */
  poly = polyvec_get_elem (r1prime, 0);
  poly_set_coeffvec_i64 (poly, r10coeffs);
  poly = polyvec_get_elem (r1prime, 1);
  poly_set_coeffvec_i64 (poly, r11coeffs);

  /* r0 */
  poly_set_coeffvec_i64 (r0prime, r0coeffs);

  polyvec_toisoring (s, sprime);
  quad_toisoring (NULL, r1, r0, NULL, r1prime, r0prime);

  for (i = 0; i < k; i++)
    {
      poly_adddot (r0[i], r1[i], s, 0);
      TEST_EXPECT (poly_eq (r0[i], zero) == 1);
    }
}

static void
lin2 (void)
{
  const int64_t s0coeffs[]
      = { 1, 6, 4, 2, 5, 0, 6, 3, 0, 5, 3, 1, 1, 1, 4, 4, 1, 0, 6, 0, 0, 3,
          1, 3, 4, 1, 5, 3, 5, 2, 3, 6, 5, 6, 2, 4, 1, 4, 1, 3, 1, 3, 0, 6,
          4, 2, 6, 6, 4, 0, 6, 0, 3, 6, 5, 0, 5, 4, 2, 5, 1, 1, 2, 0, 0, 4,
          0, 1, 1, 4, 1, 1, 1, 0, 2, 0, 4, 5, 0, 6, 2, 1, 3, 3, 3, 5, 2, 6,
          1, 6, 1, 4, 0, 2, 4, 0, 4, 1, 5, 5, 4, 0, 5, 1, 1, 1, 2, 6, 4, 3,
          4, 2, 2, 3, 4, 5, 4, 3, 2, 3, 1, 1, 5, 5, 1, 3, 4, 3 };
  const int64_t s1coeffs[]
      = { 6, 3, 4, 5, 2, 3, 5, 2, 3, 6, 6, 3, 5, 2, 1, 4, 2, 6, 1, 0, 2, 2,
          1, 6, 2, 6, 4, 5, 5, 0, 5, 3, 2, 5, 1, 1, 0, 5, 1, 5, 6, 3, 1, 3,
          6, 4, 2, 4, 2, 5, 0, 3, 3, 1, 0, 0, 3, 3, 4, 5, 4, 3, 2, 2, 2, 2,
          0, 3, 2, 4, 6, 0, 6, 3, 6, 0, 6, 0, 2, 6, 1, 2, 6, 6, 5, 0, 6, 2,
          4, 4, 5, 2, 6, 5, 5, 6, 6, 5, 6, 4, 4, 1, 3, 2, 6, 2, 4, 5, 0, 1,
          2, 6, 3, 5, 1, 6, 4, 0, 3, 4, 5, 4, 4, 3, 5, 5, 4, 3 };
  const int64_t r10coeffs[]
      = { 1, 4, 0, 6, 1, 4, 5, 2, 5, 5, 6, 2, 3, 5, 4, 0, 5, 4, 4, 5, 3, 6,
          2, 1, 2, 4, 4, 2, 4, 2, 1, 3, 5, 6, 1, 5, 2, 3, 2, 0, 6, 3, 5, 1,
          4, 0, 5, 1, 4, 2, 4, 5, 6, 5, 0, 4, 1, 1, 3, 2, 1, 5, 1, 4, 5, 2,
          2, 4, 1, 3, 3, 5, 0, 3, 3, 0, 5, 2, 2, 4, 1, 0, 2, 4, 2, 1, 1, 4,
          6, 5, 1, 2, 3, 2, 1, 1, 0, 4, 2, 1, 4, 4, 0, 1, 4, 3, 1, 5, 4, 5,
          6, 4, 5, 0, 1, 6, 2, 0, 6, 3, 3, 0, 1, 5, 4, 5, 1, 4 };
  const int64_t r11coeffs[]
      = { 3, 1, 6, 0, 5, 0, 4, 6, 3, 6, 3, 5, 6, 6, 5, 3, 5, 6, 0, 0, 2, 4,
          4, 3, 2, 0, 2, 6, 4, 4, 6, 2, 5, 1, 5, 1, 2, 4, 0, 0, 3, 5, 4, 4,
          6, 1, 4, 0, 3, 1, 0, 4, 2, 2, 1, 3, 0, 1, 2, 4, 4, 1, 6, 2, 4, 4,
          6, 0, 2, 5, 4, 1, 1, 3, 3, 0, 3, 5, 3, 0, 6, 1, 3, 2, 2, 0, 4, 4,
          2, 0, 1, 0, 3, 2, 4, 4, 2, 2, 0, 5, 1, 5, 1, 5, 1, 1, 4, 2, 4, 3,
          5, 5, 1, 3, 0, 5, 3, 6, 1, 1, 5, 3, 5, 1, 2, 5, 1, 3 };
  const int64_t r0coeffs[]
      = { 5, 1, 1, 2, 1, 3, 0, 0, 3, 5, 0, 0, 1, 0, 5, 3, 2, 4, 2, 3, 0, 0,
          6, 5, 5, 4, 2, 2, 2, 0, 1, 1, 3, 3, 3, 1, 1, 5, 5, 2, 3, 4, 4, 0,
          1, 4, 5, 4, 5, 0, 4, 1, 0, 2, 5, 3, 4, 2, 2, 6, 4, 4, 0, 2, 4, 2,
          2, 3, 3, 3, 1, 3, 4, 3, 6, 0, 1, 5, 2, 1, 0, 2, 0, 3, 0, 1, 1, 6,
          6, 6, 3, 2, 6, 1, 1, 3, 2, 6, 0, 3, 5, 3, 5, 5, 2, 6, 6, 0, 0, 0,
          5, 6, 4, 4, 0, 2, 1, 2, 4, 5, 1, 5, 5, 5, 2, 2, 0, 1 };
  const unsigned int dprime = 128;
  const unsigned int d = 64;
  const unsigned int k = dprime / d;
  const unsigned int n = 2;
  unsigned int i;
  polyvec_t row;
  poly_ptr poly;

  INT_T (q, 1);
  int_set_i64 (q, 7);
  INT_T (inv2, 1);
  int_set_i64 (inv2, -3);

  INT_T (Pprime, 1);
  int_set_i64 (Pprime, 1125899906840833);
  INT_T (P, 1);
  int_set_i64 (P, 1125899906840833);
  INT_T (Ppmodq_, 1);
  int_set_i64 (Ppmodq_, 1);
  int_srcptr Ppmodq[] = { Ppmodq_ };

  polyring_t Rprime
      = { { q, dprime, 3, 7, moduli_d128, 1, Pprime, Ppmodq, inv2 } };
  polyring_t R = { { q, d, 3, 6, moduli_d64, 1, P, Ppmodq, inv2 } };

  POLYVEC_T (sprime, Rprime, n);
  POLYVEC_T (r1prime, Rprime, n);
  POLY_T (r0prime, Rprime);
  POLYVEC_T (s, R, k * n);
  POLYVEC_T (r10, R, k * n);
  POLYVEC_T (r11, R, k * n);
  POLY_T (r00, R);
  POLY_T (r01, R);
  POLYVEC_T (zero, R, k * 2);
  POLYMAT_T (Aprime, Rprime, 2, n);
  POLYVEC_T (tprime, Rprime, 2);
  POLYMAT_T (A, R, k * 2, k * n);
  POLYVEC_T (t, R, k * 2);

  polyvec_set_zero (zero);

  /* sprime */
  poly = polyvec_get_elem (sprime, 0);
  poly_set_coeffvec_i64 (poly, s0coeffs);
  poly = polyvec_get_elem (sprime, 1);
  poly_set_coeffvec_i64 (poly, s1coeffs);

  /* r1prime */
  poly = polyvec_get_elem (r1prime, 0);
  poly_set_coeffvec_i64 (poly, r10coeffs);
  poly = polyvec_get_elem (r1prime, 1);
  poly_set_coeffvec_i64 (poly, r11coeffs);

  /* r0 */
  poly_set_coeffvec_i64 (r0prime, r0coeffs);

  for (i = 0; i < 2; i++)
    {
      polymat_get_row (row, Aprime, i);
      polyvec_set (row, r1prime);
      poly = polyvec_get_elem (tprime, i);
      poly_set (poly, r0prime);
    }

  polyvec_toisoring (s, sprime);
  lin_toisoring (A, t, Aprime, tprime);

  polyvec_addmul (t, A, s, 0);
  TEST_EXPECT (polyvec_eq (t, zero) == 1);
}

static void
quad (void)
{
  const int64_t s0coeffs[]
      = { 1, 6, 4, 2, 5, 0, 6, 3, 0, 5, 3, 1, 1, 1, 4, 4, 1, 0, 6, 0, 0, 3,
          1, 3, 4, 1, 5, 3, 5, 2, 3, 6, 5, 6, 2, 4, 1, 4, 1, 3, 1, 3, 0, 6,
          4, 2, 6, 6, 4, 0, 6, 0, 3, 6, 5, 0, 5, 4, 2, 5, 1, 1, 2, 0, 0, 4,
          0, 1, 1, 4, 1, 1, 1, 0, 2, 0, 4, 5, 0, 6, 2, 1, 3, 3, 3, 5, 2, 6,
          1, 6, 1, 4, 0, 2, 4, 0, 4, 1, 5, 5, 4, 0, 5, 1, 1, 1, 2, 6, 4, 3,
          4, 2, 2, 3, 4, 5, 4, 3, 2, 3, 1, 1, 5, 5, 1, 3, 4, 3 };
  const int64_t s1coeffs[]
      = { 6, 3, 4, 5, 2, 3, 5, 2, 3, 6, 6, 3, 5, 2, 1, 4, 2, 6, 1, 0, 2, 2,
          1, 6, 2, 6, 4, 5, 5, 0, 5, 3, 2, 5, 1, 1, 0, 5, 1, 5, 6, 3, 1, 3,
          6, 4, 2, 4, 2, 5, 0, 3, 3, 1, 0, 0, 3, 3, 4, 5, 4, 3, 2, 2, 2, 2,
          0, 3, 2, 4, 6, 0, 6, 3, 6, 0, 6, 0, 2, 6, 1, 2, 6, 6, 5, 0, 6, 2,
          4, 4, 5, 2, 6, 5, 5, 6, 6, 5, 6, 4, 4, 1, 3, 2, 6, 2, 4, 5, 0, 1,
          2, 6, 3, 5, 1, 6, 4, 0, 3, 4, 5, 4, 4, 3, 5, 5, 4, 3 };
  const int64_t R200coeffs[]
      = { 5, 0, 6, 5, 4, 2, 4, 3, 2, 1, 1, 5, 5, 6, 4, 1, 6, 0, 0, 4, 0, 0,
          5, 0, 1, 0, 0, 1, 0, 5, 2, 3, 1, 5, 0, 3, 0, 5, 1, 5, 5, 2, 6, 3,
          5, 2, 1, 2, 4, 4, 6, 0, 1, 3, 1, 5, 6, 4, 4, 4, 6, 4, 6, 4, 6, 1,
          1, 4, 0, 1, 0, 0, 1, 0, 0, 6, 2, 5, 2, 3, 6, 0, 6, 0, 6, 4, 5, 6,
          1, 3, 0, 4, 2, 3, 5, 5, 3, 0, 2, 0, 6, 0, 3, 0, 3, 1, 2, 1, 0, 6,
          0, 4, 0, 3, 1, 0, 6, 3, 6, 2, 0, 4, 2, 2, 6, 2, 0, 6 };
  const int64_t R201coeffs[]
      = { 3, 0, 0, 4, 6, 4, 4, 4, 6, 4, 0, 1, 3, 6, 0, 6, 3, 3, 3, 6, 1, 6,
          1, 4, 5, 0, 6, 2, 5, 2, 5, 4, 3, 0, 4, 6, 3, 5, 0, 3, 3, 1, 1, 0,
          0, 5, 0, 6, 5, 0, 2, 1, 1, 3, 3, 4, 5, 4, 2, 2, 6, 0, 3, 2, 1, 6,
          6, 4, 2, 0, 6, 5, 5, 3, 5, 6, 4, 2, 2, 0, 2, 0, 5, 4, 4, 2, 3, 4,
          0, 6, 6, 5, 1, 4, 6, 3, 4, 1, 2, 1, 2, 5, 3, 4, 4, 1, 1, 1, 2, 6,
          1, 4, 2, 3, 4, 1, 3, 6, 2, 2, 5, 5, 4, 1, 6, 5, 2, 2 };
  const int64_t R210coeffs[]
      = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  const int64_t R211coeffs[]
      = { 6, 1, 5, 1, 4, 2, 4, 3, 6, 1, 3, 0, 2, 4, 1, 2, 0, 4, 2, 4, 2, 6,
          0, 4, 3, 5, 0, 6, 2, 5, 5, 4, 3, 4, 5, 4, 3, 1, 2, 2, 1, 2, 3, 0,
          2, 6, 3, 1, 5, 0, 6, 5, 4, 0, 4, 3, 4, 0, 5, 6, 2, 1, 3, 1, 3, 4,
          0, 1, 0, 2, 5, 2, 3, 5, 2, 0, 3, 1, 1, 5, 3, 4, 1, 0, 5, 1, 1, 0,
          4, 1, 5, 4, 6, 4, 6, 1, 6, 0, 2, 0, 2, 4, 3, 2, 6, 2, 3, 4, 0, 2,
          5, 6, 0, 0, 3, 4, 6, 4, 4, 4, 3, 5, 4, 1, 4, 4, 4, 2 };
  const int64_t r10coeffs[]
      = { 1, 4, 0, 6, 1, 4, 5, 2, 5, 5, 6, 2, 3, 5, 4, 0, 5, 4, 4, 5, 3, 6,
          2, 1, 2, 4, 4, 2, 4, 2, 1, 3, 5, 6, 1, 5, 2, 3, 2, 0, 6, 3, 5, 1,
          4, 0, 5, 1, 4, 2, 4, 5, 6, 5, 0, 4, 1, 1, 3, 2, 1, 5, 1, 4, 5, 2,
          2, 4, 1, 3, 3, 5, 0, 3, 3, 0, 5, 2, 2, 4, 1, 0, 2, 4, 2, 1, 1, 4,
          6, 5, 1, 2, 3, 2, 1, 1, 0, 4, 2, 1, 4, 4, 0, 1, 4, 3, 1, 5, 4, 5,
          6, 4, 5, 0, 1, 6, 2, 0, 6, 3, 3, 0, 1, 5, 4, 5, 1, 4 };
  const int64_t r11coeffs[]
      = { 3, 1, 6, 0, 5, 0, 4, 6, 3, 6, 3, 5, 6, 6, 5, 3, 5, 6, 0, 0, 2, 4,
          4, 3, 2, 0, 2, 6, 4, 4, 6, 2, 5, 1, 5, 1, 2, 4, 0, 0, 3, 5, 4, 4,
          6, 1, 4, 0, 3, 1, 0, 4, 2, 2, 1, 3, 0, 1, 2, 4, 4, 1, 6, 2, 4, 4,
          6, 0, 2, 5, 4, 1, 1, 3, 3, 0, 3, 5, 3, 0, 6, 1, 3, 2, 2, 0, 4, 4,
          2, 0, 1, 0, 3, 2, 4, 4, 2, 2, 0, 5, 1, 5, 1, 5, 1, 1, 4, 2, 4, 3,
          5, 5, 1, 3, 0, 5, 3, 6, 1, 1, 5, 3, 5, 1, 2, 5, 1, 3 };
  const int64_t r0coeffs[]
      = { 2, 2, 2, 5, 3, 1, 4, 4, 3, 1, 6, 3, 3, 3, 2, 2, 0, 6, 1, 0, 6, 4,
          6, 6, 3, 5, 0, 0, 4, 6, 5, 6, 2, 3, 4, 0, 0, 1, 4, 3, 2, 3, 1, 6,
          6, 5, 2, 4, 1, 1, 1, 5, 4, 0, 5, 4, 0, 6, 2, 4, 2, 0, 2, 6, 5, 6,
          4, 5, 1, 5, 5, 6, 5, 1, 0, 1, 0, 1, 3, 5, 1, 2, 5, 3, 5, 1, 6, 5,
          2, 0, 6, 6, 4, 3, 0, 4, 5, 0, 0, 2, 3, 6, 1, 6, 0, 2, 1, 3, 1, 5,
          2, 5, 0, 3, 2, 2, 1, 4, 2, 4, 0, 1, 0, 2, 0, 3, 5, 3 };
  const unsigned int dprime = 128;
  const unsigned int d = 64;
  const unsigned int k = dprime / d;
  const unsigned int n = 2;
  unsigned int i;
  poly_ptr poly;

  INT_T (q, 1);
  int_set_i64 (q, 7);
  INT_T (inv2, 1);
  int_set_i64 (inv2, -3);

  INT_T (Pprime, 1);
  int_set_i64 (Pprime, 1125899906840833);
  INT_T (P, 1);
  int_set_i64 (P, 1125899906840833);
  INT_T (Ppmodq_, 1);
  int_set_i64 (Ppmodq_, 1);
  int_srcptr Ppmodq[] = { Ppmodq_ };

  polyring_t Rprime
      = { { q, dprime, 3, 7, moduli_d128, 1, Pprime, Ppmodq, inv2 } };
  polyring_t R = { { q, d, 3, 6, moduli_d64, 1, P, Ppmodq, inv2 } };

  POLYVEC_T (sprime, Rprime, n);
  POLYMAT_T (R2prime, Rprime, n, n);
  POLYVEC_T (r1prime, Rprime, n);
  POLY_T (r0prime, Rprime);
  POLYVEC_T (s, R, k * n);
  POLYMAT_T (R20, R, k * n, k * n);
  POLYMAT_T (R21, R, k * n, k * n);
  POLYVEC_T (r10, R, k * n);
  POLYVEC_T (r11, R, k * n);
  POLY_T (r00, R);
  POLY_T (r01, R);
  POLY_T (zero, R);
  polymat_ptr R2[n];
  polyvec_ptr r1[n];
  poly_ptr r0[n];

  poly_set_zero (zero);

  R2[0] = R20;
  R2[1] = R21;
  r1[0] = r10;
  r1[1] = r11;
  r0[0] = r00;
  r0[1] = r01;

  /* sprime */
  poly = polyvec_get_elem (sprime, 0);
  poly_set_coeffvec_i64 (poly, s0coeffs);
  poly = polyvec_get_elem (sprime, 1);
  poly_set_coeffvec_i64 (poly, s1coeffs);

  /* R2prime */
  poly = polymat_get_elem (R2prime, 0, 0);
  poly_set_coeffvec_i64 (poly, R200coeffs);
  poly = polymat_get_elem (R2prime, 0, 1);
  poly_set_coeffvec_i64 (poly, R201coeffs);
  poly = polymat_get_elem (R2prime, 1, 0);
  poly_set_coeffvec_i64 (poly, R210coeffs);
  poly = polymat_get_elem (R2prime, 1, 1);
  poly_set_coeffvec_i64 (poly, R211coeffs);

  /* r1prime */
  poly = polyvec_get_elem (r1prime, 0);
  poly_set_coeffvec_i64 (poly, r10coeffs);
  poly = polyvec_get_elem (r1prime, 1);
  poly_set_coeffvec_i64 (poly, r11coeffs);

  /* r0 */
  poly_set_coeffvec_i64 (r0prime, r0coeffs);

  polyvec_toisoring (s, sprime);
  quad_toisoring (R2, r1, r0, R2prime, r1prime, r0prime);

  for (i = 0; i < k; i++)
    {
      printf ("s\n");
      polyvec_dump (s);
      printf ("\nr0\n");
      poly_dump (r0[i]);
      printf ("\nr1\n");
      polyvec_dump (r1[i]);
      printf ("\nR2\n");
      polymat_dump (R2[i]);

      poly_adddot (r0[i], r1[i], s, 0);
      poly_fromcrt (r0[i]);
      polyvec_mul (r1[i], R2[i], s);
      polyvec_fromcrt (r1[i]);
      poly_adddot (r0[i], s, r1[i], 0);
      TEST_EXPECT (poly_eq (r0[i], zero) == 1);
    }
}
