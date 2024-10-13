#include "test.h"

#define LIMB0(x) ((x)->limbs[0])
#define LIMB1(x) ((x)->limbs[1])
#define LIMB2(x) ((x)->limbs[2])
#define LIMB3(x) ((x)->limbs[3])
#define IS_NEG(x) ((x)->neg)

static void cmp (void);
static void arith (void);
static void import_export (void);

int
main (void)
{
  lazer_init();

  import_export ();
  cmp ();
  arith ();

  TEST_PASS ();
}

static void
cmp (void)
{
  int x;
  INT_T (a, 2);
  INT_T (b, 2);
  INT_T (c, 2);
  INT_T (negz, 2);

  int_set_i64 (a, -1);
  int_set_i64 (b, 0);
  int_set_i64 (c, 1);

  int_set_i64 (negz, -1);
  LIMB0 (negz) = 0; /* negative zero */

  x = int_eqzero (a);
  TEST_EXPECT (x == 0);
  x = int_eqzero (b);
  TEST_EXPECT (x == 1);
  x = int_eqzero (c);
  TEST_EXPECT (x == 0);
  x = int_eqzero (negz);
  TEST_EXPECT (x == 1);

  x = int_eq (a, b);
  TEST_EXPECT (x == 0);
  x = int_eq (a, c);
  TEST_EXPECT (x == 0);
  x = int_eq (negz, b);
  TEST_EXPECT (x == 1);
  x = int_eq (a, a);
  TEST_EXPECT (x == 1);
  x = int_eq (c, c);
  TEST_EXPECT (x == 1);

  x = int_lt (a, b);
  TEST_EXPECT (x == 1);
  x = int_lt (b, c);
  TEST_EXPECT (x == 1);
  x = int_lt (a, c);
  TEST_EXPECT (x == 1);
  x = int_lt (negz, b);
  TEST_EXPECT (x == 0);
  x = int_lt (b, a);
  TEST_EXPECT (x == 0);
  x = int_lt (c, b);
  TEST_EXPECT (x == 0);
  x = int_lt (a, a);
  TEST_EXPECT (x == 0);

  x = int_le (a, b);
  TEST_EXPECT (x == 1);
  x = int_le (b, c);
  TEST_EXPECT (x == 1);
  x = int_le (a, c);
  TEST_EXPECT (x == 1);
  x = int_le (negz, b);
  TEST_EXPECT (x == 1);
  x = int_le (b, a);
  TEST_EXPECT (x == 0);
  x = int_le (c, b);
  TEST_EXPECT (x == 0);
  x = int_le (a, a);
  TEST_EXPECT (x == 1);

  x = int_gt (a, b);
  TEST_EXPECT (x == 0);
  x = int_gt (b, c);
  TEST_EXPECT (x == 0);
  x = int_gt (a, c);
  TEST_EXPECT (x == 0);
  x = int_gt (negz, b);
  TEST_EXPECT (x == 0);
  x = int_gt (b, a);
  TEST_EXPECT (x == 1);
  x = int_gt (c, b);
  TEST_EXPECT (x == 1);
  x = int_gt (a, a);
  TEST_EXPECT (x == 0);

  x = int_ge (a, b);
  TEST_EXPECT (x == 0);
  x = int_ge (b, c);
  TEST_EXPECT (x == 0);
  x = int_ge (a, c);
  TEST_EXPECT (x == 0);
  x = int_ge (negz, b);
  TEST_EXPECT (x == 1);
  x = int_ge (b, a);
  TEST_EXPECT (x == 1);
  x = int_ge (c, b);
  TEST_EXPECT (x == 1);
  x = int_ge (a, a);
  TEST_EXPECT (x == 1);

  x = int_abseq (a, b);
  TEST_EXPECT (x == 0);
  x = int_abseq (a, c);
  TEST_EXPECT (x == 1);
  x = int_abseq (negz, b);
  TEST_EXPECT (x == 1);

  x = int_abslt (a, b);
  TEST_EXPECT (x == 0);
  x = int_abslt (b, a);
  TEST_EXPECT (x == 1);
  x = int_abslt (a, c);
  TEST_EXPECT (x == 0);
  x = int_abslt (negz, b);
  TEST_EXPECT (x == 0);

  x = int_absle (a, b);
  TEST_EXPECT (x == 0);
  x = int_absle (b, a);
  TEST_EXPECT (x == 1);
  x = int_absle (a, c);
  TEST_EXPECT (x == 1);
  x = int_absle (b, negz);
  TEST_EXPECT (x == 1);

  x = int_absgt (a, b);
  TEST_EXPECT (x == 1);
  x = int_absgt (b, a);
  TEST_EXPECT (x == 0);
  x = int_absgt (a, c);
  TEST_EXPECT (x == 0);
  x = int_absgt (negz, b);
  TEST_EXPECT (x == 0);

  x = int_absge (a, b);
  TEST_EXPECT (x == 1);
  x = int_absge (b, a);
  TEST_EXPECT (x == 0);
  x = int_absge (a, c);
  TEST_EXPECT (x == 1);
  x = int_absge (negz, b);
  TEST_EXPECT (x == 1);
}

static void
arith (void)
{
  int x;
  INT_T (a, 2);
  INT_T (b, 2);
  INT_T (c, 2);
  INT_T (d, 4);
  INT_T (e, 3);
  INT_T (f, 1);
  INT_T (m, 1);
  INT_T (negz, 2);
  INT_T (negz2, 4);
  int_t g, h;

  int_alloc (g, 2);
  int_alloc (h, 2);
  int_set_i64 (g, -3);
  int_set_i64 (h, -1);
  int_add (g, h, g);
  TEST_EXPECT (IS_NEG (g) == 1);
  TEST_EXPECT (LIMB0 (g) == 4);
  TEST_EXPECT (LIMB1 (g) == 0);
  int_free (g);
  int_free (h);

  int_set_i64 (negz, -1);
  LIMB0 (negz) = 0; /* negative zero */
  int_set_i64 (negz2, -1);
  LIMB0 (negz2) = 0; /* negative zero */

  int_set_i64 (a, -5);
  int_set (b, a);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 5);

  int_set_i64 (a, -1);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 1);
  int_neg (a, a);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);
  int_neg (b, a);
  TEST_EXPECT (IS_NEG (b) == 1);
  TEST_EXPECT (LIMB0 (b) == 1);

  int_set_i64 (a, -1);
  int_abs (a, a);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (d, -1);
  LIMB0 (d) = 4;
  LIMB1 (d) = 3;
  LIMB2 (d) = 2;
  LIMB3 (d) = 1;
  int_rshift (d, d, 1);
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == (1ULL << 63) + 2);
  TEST_EXPECT (LIMB1 (d) == 1);
  TEST_EXPECT (LIMB2 (d) == (1ULL << 63) + 1);
  TEST_EXPECT (LIMB3 (d) == 0);

  int_set_i64 (d, 1);
  LIMB0 (d) = 1;
  LIMB1 (d) = 2;
  LIMB2 (d) = 3;
  LIMB3 (d) = 4;
  int_rshift (e, d, 65);
  TEST_EXPECT (IS_NEG (e) == 0);
  TEST_EXPECT (LIMB0 (e) == (1ULL << 63) + 1);
  TEST_EXPECT (LIMB1 (e) == 1);
  TEST_EXPECT (LIMB2 (e) == 2);

  int_set_i64 (a, 1);
  LIMB0 (a) = 3;
  LIMB1 (a) = 5;
  LIMB0 (d) = 1;
  LIMB1 (d) = 1;
  LIMB2 (d) = 1;
  LIMB3 (d) = 1;
  int_rshift (d, a, 2);
  TEST_EXPECT (IS_NEG (d) == 0);
  TEST_EXPECT (LIMB0 (d) == (1ULL << 62));
  TEST_EXPECT (LIMB1 (d) == 1);
  TEST_EXPECT (LIMB2 (d) == 0);
  TEST_EXPECT (LIMB3 (d) == 0);

  int_set_i64 (d, -1);
  LIMB0 (d) = ~0;
  LIMB1 (d) = ~0;
  LIMB2 (d) = ~0;
  LIMB3 (d) = ~0;
  int_rshift (d, d, 1000);
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == 0);
  TEST_EXPECT (LIMB1 (d) == 0);
  TEST_EXPECT (LIMB2 (d) == 0);
  TEST_EXPECT (LIMB3 (d) == 0);

  int_set_i64 (a, -1);
  LIMB0 (a) = ~0;
  LIMB1 (a) = 5;
  int_lshift (d, a, 2);
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == ((~0UL) & (~0x3UL)));
  TEST_EXPECT (LIMB1 (d) == (5 << 2) + 3);
  TEST_EXPECT (LIMB2 (d) == 0);
  TEST_EXPECT (LIMB3 (d) == 0);

  int_set_i64 (a, -1);
  LIMB0 (a) = ~0;
  LIMB1 (a) = 5;
  int_lshift (d, a, 66);
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == 0);
  TEST_EXPECT (LIMB1 (d) == ((~0UL) & (~0x3UL)));
  TEST_EXPECT (LIMB2 (d) == (5 << 2) + 3);
  TEST_EXPECT (LIMB3 (d) == 0);

  int_set_i64 (a, 1);
  LIMB0 (a) = ~0;
  LIMB1 (a) = 5;
  int_lshift (a, a, 66);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);
  TEST_EXPECT (LIMB1 (a) == ((~0UL) & (~0x3UL)));

  int_set_i64 (a, 7);
  int_set_i64 (b, 3);
  int_set_i64 (c, 11);
  int_mul (d, a, b);
  TEST_EXPECT (IS_NEG (d) == 0);
  TEST_EXPECT (LIMB0 (d) == 21);
  int_mod (a, d, c);
  TEST_EXPECT (IS_NEG (d) == 0);
  TEST_EXPECT (LIMB0 (d) == 21);

  int_set_i64 (a, 4);
  int_set_i64 (b, -5);
  int_set_i64 (c, 11);
  int_set_i64 (f, 11);
  int_mul (d, a, b);
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == 20);
  int_mod (m, d, f);
  TEST_EXPECT (IS_NEG (m) == 1);
  TEST_EXPECT (LIMB0 (m) == 9);

  int_set_i64 (a, -4);
  int_mul (d, a, a);
  TEST_EXPECT (IS_NEG (d) == 0);
  TEST_EXPECT (LIMB0 (d) == 16);

  int_set_i64 (a, -4);
  int_set_i64 (b, 5);
  int_mul (d, a, b);
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == 20);

  int_set_i64 (a, -3);
  int_mul (d, a, a);
  TEST_EXPECT (IS_NEG (d) == 0);
  TEST_EXPECT (LIMB0 (d) == 9);

  int_set_i64 (a, 0);
  int_set_i64 (b, -1);
  int_mul (d, a, b); /* negative zero */
  TEST_EXPECT (IS_NEG (d) == 1);
  TEST_EXPECT (LIMB0 (d) == 0);

  int_set_i64 (a, 0);
  int_set_i64 (b, -1);
  int_mul (d, a, b); /* negative zero */
  x = int_sgn (d);
  TEST_EXPECT (x == 1);

  int_set_i64 (a, 0);
  x = int_sgn (d);
  TEST_EXPECT (x == 1);

  int_set_i64 (a, 1);
  x = int_sgn (d);
  TEST_EXPECT (x == 1);

  int_set_i64 (a, -1);
  x = int_sgn (d);
  TEST_EXPECT (x == 1);

  int_set_i64 (b, 7);
  int_set_i64 (c, 17);
  int_add (a, b, c);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 24);

  int_set_i64 (b, -17);
  int_set_i64 (c, -7);
  int_add (c, b, c);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 24);

  int_set_i64 (b, -7);
  int_set_i64 (c, 17);
  int_add (b, b, c);
  TEST_EXPECT (IS_NEG (b) == 0);
  TEST_EXPECT (LIMB0 (b) == 10);

  int_set_i64 (b, -17);
  int_set_i64 (c, 7);
  int_add (c, b, c);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 10);

  int_set_i64 (b, 17);
  int_set_i64 (c, -7);
  int_add (b, b, c);
  TEST_EXPECT (IS_NEG (b) == 0);
  TEST_EXPECT (LIMB0 (b) == 10);

  int_set_i64 (b, 7);
  int_set_i64 (c, -17);
  int_add (c, b, c);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 10);

  int_set_i64 (c, -17);
  int_add (c, c, c);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 34);

  int_set_i64 (c, 17);
  int_add (c, c, c);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 34);

  int_set_i64 (b, 7);
  int_set_i64 (c, 17);
  int_sub (a, b, c);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 10);

  int_set_i64 (b, 17);
  int_set_i64 (c, 7);
  int_sub (a, b, c);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 10);

  int_set_i64 (b, -7);
  int_set_i64 (c, -17);
  int_sub (a, b, c);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 10);
  
  int_set_i64 (b, -17);
  int_set_i64 (c, -7);
  int_sub (a, b, c);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 10);

  int_set_i64 (b, -7);
  int_set_i64 (c, 17);
  int_sub (b, b, c);
  TEST_EXPECT (IS_NEG (b) == 1);
  TEST_EXPECT (LIMB0 (b) == 24);

  int_set_i64 (b, -17);
  int_set_i64 (c, 7);
  int_sub (c, b, c);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 24);

  int_set_i64 (b, 17);
  int_set_i64 (c, -7);
  int_sub (b, b, c);
  TEST_EXPECT (IS_NEG (b) == 0);
  TEST_EXPECT (LIMB0 (b) == 24);

  int_set_i64 (b, 7);
  int_set_i64 (c, -17);
  int_sub (c, b, c);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 24);

  //XXXint_set_i64 (c, -17);
  //int_sub (c, c, c);
  //TEST_EXPECT (IS_NEG (c) == 0); /* pos. or neg. 0 ? XXX */
  //TEST_EXPECT (LIMB0 (c) == 0);

  int_set_i64 (c, 17);
  int_sub (c, c, c);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 0);

  int_redc (a, negz, b);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, 0);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, 16);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (a, 9);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 8);

  int_set_i64 (a, 8);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 8);

  int_set_i64 (b, 17);
  int_redp (a, negz, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, -8);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 9);

  int_set_i64 (a, -1);
  int_set_i64 (b, 17);
  int_redp (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 16);

  int_set_i64 (a, 0);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, 0);
  int_set_i64 (b, 17);
  int_neg (a, a);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, 1);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (a, 16);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 16);

  int_set_i64 (a, 0);
  LIMB0 (a) = 1;
  LIMB1 (a) = 1; /* a = 2^64 + 1 */
  int_set_i64 (c, 31);
  int_invmod (b, a, c);
  TEST_EXPECT (IS_NEG (b) == 0);
  TEST_EXPECT (LIMB0 (b) == 11);

  int_set_i64 (a, -1);
  LIMB0 (a) = 1;
  LIMB1 (a) = 1; /* a = -2^64 + 1 */
  int_set_i64 (c, 31);
  int_invmod (b, a, c);
  TEST_EXPECT (IS_NEG (b) == 1);
  TEST_EXPECT (LIMB0 (b) == 11);

  int_set_i64 (c, 31);
  int_invmod (b, negz, c);
  TEST_EXPECT (IS_NEG (b) == 1);
  TEST_EXPECT (LIMB0 (b) == 0);

  int_set_i64 (c, 31);
  int_mod (a, negz2, c);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, -16);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (a, -9);
  int_set_i64 (b, 17);
  int_redc (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 8);

  int_set_i64 (a, -8);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (LIMB0 (a) == 8);

  int_set_i64 (a, -1);
  int_set_i64 (b, 17);
  int_redc (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 1);

  int_set_i64 (b, 17);
  int_redc (c, negz, b);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 0);

  int_set_i64 (a, 0);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, 1);
  int_set_i64 (b, 17);
  int_redc (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (a, 8);
  int_set_i64 (b, 17);
  int_redc (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 8);

  int_set_i64 (a, 9);
  int_set_i64 (b, 17);
  int_redc (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 8);

  int_set_i64 (a, 16);
  int_set_i64 (b, 17);
  int_redc (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 1);
  TEST_EXPECT (LIMB0 (c) == 1);

  int_set_i64 (a, -16);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (a, -9);
  int_set_i64 (b, 17);
  int_redp (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 8);

  int_set_i64 (a, -8);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 9);

  int_set_i64 (a, -1);
  int_set_i64 (b, 17);
  int_redp (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 16);

  int_set_i64 (b, 17);
  int_redp (c, negz, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 0);

  int_set_i64 (a, 0);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 0);

  int_set_i64 (a, 1);
  int_set_i64 (b, 17);
  int_redp (a, a, b);
  TEST_EXPECT (IS_NEG (a) == 0);
  TEST_EXPECT (LIMB0 (a) == 1);

  int_set_i64 (a, 8);
  int_set_i64 (b, 17);
  int_redp (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 8);

  int_set_i64 (a, 9);
  int_set_i64 (b, 17);
  int_redp (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 9);

  int_set_i64 (a, 16);
  int_set_i64 (b, 17);
  int_redp (c, a, b);
  TEST_EXPECT (IS_NEG (c) == 0);
  TEST_EXPECT (LIMB0 (c) == 16);

  int_set_i64 (e, -1);
  LIMB0 (e) = 1;
  LIMB1 (e) = 4;
  LIMB2 (e) = 0;
  int_set_i64 (c, 1);
  LIMB0 (c) = 0;
  LIMB1 (c) = 1;
  int_div (a, b, e, c);
  TEST_EXPECT (IS_NEG (a) == 1);
  TEST_EXPECT (IS_NEG (b) == 1);
  TEST_EXPECT (LIMB0 (a) == 4);
  TEST_EXPECT (LIMB0 (b) == 1);
}

static void
import_export (void)
{
  uint8_t buf[1000];
  size_t len;
  INT_T (a, 2);
  INT_T (b, 2);

  int_set_i64 (a, -11);

  int_export (buf, &len, a);
  int_import (b, buf, len);
  TEST_EXPECT (IS_NEG (b) == 0);
  TEST_EXPECT (LIMB0 (b) == 11);
}
