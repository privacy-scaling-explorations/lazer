#include "test.h"

int
main (void)
{
  lazer_init();

  spolymat_t A, B, C;
  INT_T (q, 1);
  polyring_t Rq = { { q, 64, 0, 0, NULL, 0, NULL, NULL, NULL } };
  poly_ptr poly;
  POLY_T (zero, Rq);
  POLY_T (one, Rq);

  int_set_i64 (q, 7);
  poly_set_zero (zero);
  poly_set_one (one);

  spolymat_alloc (A, Rq, 20, 20, 200);
  spolymat_alloc (B, Rq, 20, 20, 20);
  spolymat_alloc (C, Rq, 20, 20, 8);

  poly = spolymat_insert_elem (A, 5, 11);
  poly_set_zero (poly);
  poly = spolymat_insert_elem (A, 5, 1);
  poly_set_zero (poly);
  poly = spolymat_insert_elem (A, 1, 1);
  poly_set_zero (poly);
  poly = spolymat_insert_elem (A, 0, 4);
  poly_set_zero (poly);
  poly = spolymat_insert_elem (A, 2, 2);
  poly_set_zero (poly);
  poly = spolymat_insert_elem (A, 2, 4);
  poly_set_zero (poly);
  poly = spolymat_insert_elem (A, 2, 3);
  poly_set_one (poly);

  poly = spolymat_insert_elem (B, 5, 11);
  poly_set_one (poly);
  poly = spolymat_insert_elem (B, 11, 0);
  poly_set_one (poly);

  poly = spolymat_get_elem (A, 0);
  TEST_EXPECT (poly_eq (poly, zero));

  spolymat_sort (A);
  spolymat_sort (B);

  poly = spolymat_get_elem (A, 3);
  TEST_EXPECT (poly_eq (poly, one));

  spolymat_add (C, A, B, 0);

  poly = spolymat_get_elem (C, 6);
  TEST_EXPECT (poly_eq (poly, one));
  poly = spolymat_get_elem (C, 0);
  TEST_EXPECT (poly_eq (poly, zero));
  poly = spolymat_get_elem (C, 7);
  TEST_EXPECT (poly_eq (poly, one));

  spolymat_free (A);
  spolymat_free (B);
  spolymat_free (C);
  TEST_PASS ();
}
