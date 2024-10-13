#include "test.h"

#define LIMB0(x) ((x)->limbs[0])
#define LIMB1(x) ((x)->limbs[1])
#define LIMB2(x) ((x)->limbs[2])
#define LIMB3(x) ((x)->limbs[3])
#define IS_NEG(x) ((x)->neg)

static void arith (void);

int
main (void)
{
  lazer_init();

  arith ();

  TEST_PASS ();
}

static void
arith (void)
{
  const int64_t i64vec_src[] = { -1, 2, 3, 4, 4, 0 };
  int64_t i64vec[6];
  int_srcptr srcptri;
  int_ptr ptri;
  intvec_t subv;
  intmat_t subm, c, d;
  INTMAT_T (a, 2, 3, 2);
  INTMAT_T (b, 2, 3, 2);

  INT_T (ni, 2);
  INT_T (pi, 2);

  IS_NEG (ni) = 1;
  LIMB0 (ni) = 8;
  LIMB1 (ni) = 7;

  IS_NEG (pi) = 0;
  LIMB0 (pi) = 9;

  intmat_alloc (c, 2, 2, 2);
  intmat_alloc (d, 2, 2, 2);
  intmat_set_one (c);
  intmat_set_one (d);
  TEST_EXPECT (intmat_eq (c, d) == 1);
  intmat_free (c);
  intmat_free (d);

  intmat_set_elem (a, 0, 0, ni);
  intmat_set_elem (a, 0, 1, ni);
  intmat_set_elem (a, 0, 2, pi);
  intmat_set_elem (a, 1, 0, ni);
  intmat_set_elem (a, 1, 1, pi);
  intmat_set_elem (a, 1, 2, pi);
  TEST_EXPECT (intmat_eq (a, a) == 1);

  intmat_set (b, a);
  TEST_EXPECT (intmat_eq (a, b) == 1);

  intmat_set_elem (b, 0, 0, pi);
  intmat_set_elem (b, 1, 1, ni);
  TEST_EXPECT (intmat_eq (a, b) == 0);

  ptri = intmat_get_elem (a, 0, 2);
  TEST_EXPECT (int_eq (ptri, pi) == 1);

  srcptri = intmat_get_elem_src (b, 1, 2);
  TEST_EXPECT (int_eq (srcptri, pi) == 1);

  intmat_dump (a);
  intmat_dump (b);

  intmat_set_i64 (a, i64vec_src);
  intmat_dump (a);
  ptri = intmat_get_elem (a, 0, 0);
  TEST_EXPECT (int_get_i64 (ptri) == -1);
  ptri = intmat_get_elem (a, 1, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  memset (i64vec, 0, sizeof (i64vec));
  intmat_get_i64 (i64vec, a);
  TEST_EXPECT (memcmp (i64vec, i64vec_src, sizeof (i64vec)) == 0);

  intmat_get_submat (subm, a, 1, 1, 1, 2, 1, 1);
  intmat_dump (subm);
  ptri = intmat_get_elem (subm, 0, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 4);
  ptri = intmat_get_elem (subm, 0, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 0);

  intmat_get_submat (subm, a, 0, 1, 2, 1, 1, 2);
  intmat_dump (subm);
  ptri = intmat_get_elem (subm, 0, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 2);
  ptri = intmat_get_elem (subm, 1, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  intmat_get_row (subv, a, 1);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 3);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 4);
  ptri = intvec_get_elem (subv, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 4);
  ptri = intvec_get_elem (subv, 2);
  TEST_EXPECT (int_get_i64 (ptri) == 0);

  intmat_get_col (subv, a, 1);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 2);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 2);
  ptri = intvec_get_elem (subv, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  intmat_get_diag (subv, a, 0);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 2);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == -1);
  ptri = intvec_get_elem (subv, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  intmat_get_diag (subv, a, 1);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 2);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 2);
  ptri = intvec_get_elem (subv, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 0);

  intmat_get_diag (subv, a, -1);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 1);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  intmat_get_antidiag (subv, a, 0);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 2);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 3);
  ptri = intvec_get_elem (subv, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  intmat_get_antidiag (subv, a, 1);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 2);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 2);
  ptri = intvec_get_elem (subv, 1);
  TEST_EXPECT (int_get_i64 (ptri) == 4);

  intmat_get_antidiag (subv, a, -1);
  intvec_dump (subv);
  TEST_EXPECT (intvec_get_nelems (subv) == 1);
  ptri = intvec_get_elem (subv, 0);
  TEST_EXPECT (int_get_i64 (ptri) == 0);
}
