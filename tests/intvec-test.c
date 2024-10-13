#include "test.h"

#define LIMB0(x) ((x)->limbs[0])
#define LIMB1(x) ((x)->limbs[1])
#define LIMB2(x) ((x)->limbs[2])
#define LIMB3(x) ((x)->limbs[3])
#define IS_NEG(x) ((x)->neg)

static void arith (void);
static void arith2 (void);
static void arith3 (void);
static void autom (void);
static void mulsgn (void);

int
main (void)
{
  lazer_init();

  arith ();
  arith2 ();
  arith3 ();
  autom ();
  mulsgn ();

  TEST_PASS ();
}

static void
arith (void)
{
  int_srcptr srcptri;
  int_ptr ptri;
  intvec_t subv, c, d;

  INT_T (ni, 2);
  INT_T (pi, 1);
  INTVEC_T (a, 4, 2);
  INTVEC_T (b, 4, 2);

  IS_NEG (ni) = 1;
  LIMB0 (ni) = 8;
  LIMB1 (ni) = 7;

  IS_NEG (pi) = 0;
  LIMB0 (pi) = 9;

  intvec_alloc (c, 3, 2);
  intvec_alloc (d, 3, 2);
  intvec_set_ones (c);
  intvec_set_ones (d);
  TEST_EXPECT (intvec_eq (c, d) == 1);
  intvec_free(c);
  intvec_free(d);

  intvec_set_elem (a, 0, ni);
  intvec_set_elem (a, 1, pi);
  intvec_set_elem (a, 2, ni);
  intvec_set_elem (a, 3, pi);
  TEST_EXPECT (intvec_eq (a, a) == 1);

  intvec_set (b, a);
  TEST_EXPECT (intvec_eq (a, b) == 1);

  intvec_set_elem (b, 0, ni);
  intvec_set_elem (b, 1, pi);
  intvec_set_elem (b, 2, pi);
  intvec_set_elem (b, 3, ni);
  TEST_EXPECT (intvec_eq (a, b) == 0);

  ptri = intvec_get_elem (a, 2);
  TEST_EXPECT (int_eq (ptri, ni) == 1);

  srcptri = intvec_get_elem_src (b, 3);
  TEST_EXPECT (int_eq (srcptri, ni) == 1);

  intvec_dump (a);
  intvec_dump (b);

  intvec_set_i64 (a, (long[]){ -1, 2, 0, 0 });
  ptri = intvec_get_elem (a, 0);
  TEST_EXPECT (int_get_i64 (ptri) == -1);
  ptri = intvec_get_elem (a, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 2);
  ptri = intvec_get_elem (a, 2);
  TEST_EXPECT (int_get_i64 (ptri) == 0);
  ptri = intvec_get_elem (a, 3);
  TEST_EXPECT (int_get_i64 (ptri) == 0);
  intvec_dump (a);

  intvec_get_subvec (subv, a, 1, 2, 1);
  intvec_dump (subv);
  intvec_get_subvec (subv, a, 0, 2, 2);
  intvec_dump (subv);
}

static void
arith2 (void)
{
  INTVEC_T (r, 2, 4);
  INTVEC_T (v, 2, 2);
  INT_T (s, 2);

  IS_NEG (s) = 1;
  LIMB0 (s) = 2;
  LIMB1 (s) = 1;

  intvec_set_i64 (v, (const long[]){ -2, 3 });

  intvec_scale (r, s, v);
  TEST_EXPECT (IS_NEG (intvec_get_elem_src (r, 0)) == 0);
  TEST_EXPECT (LIMB0 (intvec_get_elem_src (r, 0)) == 4);
  TEST_EXPECT (LIMB1 (intvec_get_elem_src (r, 0)) == 2);
  TEST_EXPECT (IS_NEG (intvec_get_elem_src (r, 1)) == 1);
  TEST_EXPECT (LIMB0 (intvec_get_elem_src (r, 1)) == 6);
  TEST_EXPECT (LIMB1 (intvec_get_elem_src (r, 1)) == 3);

  intvec_rrot (v, v, 1);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 0)) == 3);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 1)) == 2);

  intvec_rrot (v, v, 1);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 0)) == 2);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 1)) == -3);

  intvec_lrot (v, v, 1);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 0)) == 3);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 1)) == 2);

  intvec_lrot (v, v, 1);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 0)) == -2);
  TEST_EXPECT (int_get_i64 (intvec_get_elem_src (v, 1)) == 3);
}

static void
arith3 (void)
{
  const int64_t iv[] = { -12, 100 };
  const int64_t im[] = { -2, 3, 0, -10 };
  INTVEC_T (r, 2, 2);
  INTMAT_T (m, 2, 2, 1);
  INTVEC_T (v, 2, 1);

  intvec_set_i64 (v, iv);
  intmat_set_i64 (m, im);

  intvec_mul_matvec (r, m, v);
  TEST_EXPECT (intvec_get_elem_i64 (r, 0) == 324);
  TEST_EXPECT (intvec_get_elem_i64 (r, 1) == -1000);
}

static void
autom (void)
{
  const int64_t iv1[] = { -2, 3, 1, -1, 4, 5, -4, 2 };
  const int64_t iv1_auto[] = { -2, -2, 4, -5, -4, 1, -1, -3 };
  const int64_t iv2[] = { 0, 3, 1, 2, -1, -5, -4, 0 };
  const int64_t iv2_auto[] = { 0, 0, 4, 5, 1, -2, -1, -3 };
  const int64_t iv3[] = { -1, -2, 3, 1, 0, -1, -3, 2 };
  const int64_t iv3_auto[] = { -1, -2, 3, 1, 0, -1, -3, 2 };
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  int64_t ri[8];
  int64_t omega;
  INT_T (omega_, 1);
  INTVEC_T (v, 8, 1);
  INTVEC_T (r, 8, 1);
  int hit_bnd;

  intvec_set_i64 (v, iv1);
  intvec_dump (v);
  intvec_auto (r, v);
  intvec_dump (r);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv1_auto, sizeof (ri)) == 0);

  intvec_set_i64 (v, iv2);
  intvec_auto (r, v);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv2_auto, sizeof (ri)) == 0);

  intvec_set_i64 (v, iv3);
  intvec_auto (r, v);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv3_auto, sizeof (ri)) == 0);

  /* in-place */
  intvec_set_i64 (r, iv1);
  intvec_dump (r);
  intvec_auto_self (r);
  intvec_dump (r);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv1_auto, sizeof (ri)) == 0);

  intvec_set_i64 (r, iv2);
  intvec_auto_self (r);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv2_auto, sizeof (ri)) == 0);

  intvec_set_i64 (r, iv3);
  intvec_auto_self (r);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv3_auto, sizeof (ri)) == 0);

  /* create challenge */
  omega = 2;
  int_set_i64 (omega_, omega);

  hit_bnd = 0;
  dom = 0;
  for (i = 0; i < 10; i++)
    {
      intvec_urandom_autostable (r, omega, 3, seed, dom);
      dom++;

      intvec_auto (v, r);
      TEST_EXPECT (intvec_eq (v, r) == 1);

      TEST_EXPECT (intvec_le (r, omega_) == 1);
      if (intvec_lt (r, omega_) != 1)
        hit_bnd = 1;

      intvec_mul_sgn_self (r, -1);
      TEST_EXPECT (intvec_le (r, omega_) == 1);
      if (intvec_lt (r, omega_) != 1)
        hit_bnd = 1;
    }
  TEST_EXPECT (hit_bnd == 1);

  omega = 8;
  int_set_i64 (omega_, omega);

  hit_bnd = 0;
  dom = 0;
  for (i = 0; i < 10; i++)
    {
      intvec_urandom_autostable (r, omega, 5, seed, dom);
      dom++;

      intvec_auto (v, r);
      TEST_EXPECT (intvec_eq (v, r) == 1);

      TEST_EXPECT (intvec_le (r, omega_) == 1);
      if (intvec_lt (r, omega_) != 1)
        hit_bnd = 1;

      intvec_mul_sgn_self (r, -1);
      TEST_EXPECT (intvec_le (r, omega_) == 1);
      if (intvec_lt (r, omega_) != 1)
        hit_bnd = 1;
    }
  TEST_EXPECT (hit_bnd == 1);
}

static void
mulsgn (void)
{
  const int64_t iv1[] = { -2, 3, 1, -1, 4, 5, -4, 2 };
  const int64_t iv1_mul1[] = { -2, 3, 1, -1, 4, 5, -4, 2 };
  const int64_t iv2[] = { 0, 3, 1, 2, -1, -5, -4, 0 };
  const int64_t iv2_muln1[] = { 0, -3, -1, -2, 1, 5, 4, 0 };
  int64_t ri[8];
  INTVEC_T (r, 8, 1);

  /* in-place */
  intvec_set_i64 (r, iv1);
  intvec_dump (r);
  intvec_mul_sgn_self (r, 1);
  intvec_dump (r);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv1_mul1, sizeof (ri)) == 0);

  intvec_set_i64 (r, iv2);
  intvec_mul_sgn_self (r, -1);
  intvec_get_i64 (ri, r);
  TEST_EXPECT (memcmp (ri, iv2_muln1, sizeof (ri)) == 0);
}
