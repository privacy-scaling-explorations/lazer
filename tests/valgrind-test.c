#include "abdlop-params1.h"
#include "lazer.h"
#include "test.h"
#include <valgrind/memcheck.h>

#if ASSERT == ASSERT_ENABLED
int
main (void)
{
  lazer_init();

  TEST_SKIP_IF (VALGRIND == VALGRIND_DISABLED);
  TEST_SKIP_IF (ASSERT == ASSERT_ENABLED);
}
#else

polyring_srcptr Rq = params1_ring;
int_srcptr mod = params1_q;
uint8_t seed[32];
uint32_t dom = 0;
unsigned int i;

static void test_poly_tocrt (void);
static void test_poly_fromcrt (void);
static void test_poly_add (void);
static void test_poly_sub (void);
static void test_poly_mul (void);
static void test_poly_scale (void);
static void test_poly_rrot (void);
static void test_poly_lrot (void);
static void test_poly_rshift (void);
static void test_poly_lshift (void);
static void test_poly_mod (void);
static void test_poly_redc (void);
static void test_poly_eq (void);

static inline void
_undef_int (int_t a)
{
  bytes_urandom (seed, sizeof (seed));
  int_grandom (a, 10, seed, dom);
  VALGRIND_MAKE_MEM_UNDEFINED (a->limbs, a->nlimbs * sizeof (a->limbs[0]));
  VALGRIND_MAKE_MEM_UNDEFINED (&(a->neg), sizeof (a->neg));
}

static inline void
_undef_coeffs (poly_t a)
{
  bytes_urandom (seed, sizeof (seed));
  poly_grandom (a, 10, seed, dom);
  VALGRIND_MAKE_MEM_UNDEFINED (a->coeffs->bytes,
                               a->coeffs->nlimbs * a->coeffs->nelems
                                   * sizeof (a->coeffs->elems[0].limbs[0]));
}

static inline void
_undef_crtrep (poly_t a)
{
  bytes_urandom (seed, sizeof (seed));
  poly_grandom (a, 10, seed, dom);
  poly_tocrt (a);
  VALGRIND_MAKE_MEM_UNDEFINED (a->crtrep, a->ring->nmoduli * a->ring->d
                                              * sizeof (crtcoeff_t));
}

int
main (void)
{
  test_poly_tocrt ();
  test_poly_fromcrt ();
  test_poly_add ();
  test_poly_sub ();
  test_poly_mul ();
  test_poly_scale ();
  test_poly_rrot ();
  test_poly_lrot ();
  test_poly_rshift ();
  test_poly_lshift ();
  test_poly_mod ();
  test_poly_redc ();
  test_poly_eq ();
  TEST_PASS ();
}

static void
test_poly_tocrt (void)
{
  POLY_T (a, Rq);

  _undef_coeffs (a);
  poly_tocrt (a);
}

static void
test_poly_fromcrt (void)
{
  POLY_T (a, Rq);

  _undef_crtrep (a);
  poly_fromcrt (a);
}

static void
test_poly_add (void)
{
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  _undef_coeffs (b);
  poly_add (r, a, b, 0);
  poly_add (r, a, b, 1);

  _undef_crtrep (a);
  _undef_crtrep (b);
  poly_add (r, a, b, 0);
  poly_add (r, a, b, 0);
}

static void
test_poly_sub (void)
{
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  _undef_coeffs (b);
  poly_sub (r, a, b, 0);
  poly_sub (r, a, b, 1);

  _undef_crtrep (a);
  _undef_crtrep (b);
  poly_sub (r, a, b, 0);
  poly_sub (r, a, b, 0);
}

static void
test_poly_mul (void)
{
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  POLY_T (r, Rq);

  _undef_crtrep (a);
  _undef_crtrep (b);
  poly_mul (r, a, b);
}

static void
test_poly_scale (void)
{
  INT_T (a, Rq->q->nlimbs);
  POLY_T (b, Rq);
  POLY_T (r, Rq);

  _undef_int (a);
  _undef_coeffs (b);
  poly_scale (r, a, b);
}

static void
test_poly_rrot (void)
{
  POLY_T (a, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  poly_rrot (r, a, 1);
}

static void
test_poly_lrot (void)
{
  POLY_T (a, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  poly_lrot (r, a, 1);
}

static void
test_poly_rshift (void)
{
  POLY_T (a, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  poly_rshift (r, a, 1);
}

static void
test_poly_lshift (void)
{
  POLY_T (a, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  poly_lshift (r, a, 1);
}

static void
test_poly_mod (void)
{
  POLY_T (a, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  poly_mod (r, a);
}

static void
test_poly_redc (void)
{
  POLY_T (a, Rq);
  POLY_T (r, Rq);

  _undef_coeffs (a);
  poly_redc (r, a);
}

static void
test_poly_eq (void)
{
  POLY_T (a, Rq);
  POLY_T (b, Rq);

  _undef_coeffs (a);
  _undef_coeffs (b);
  poly_eq (a, b);
  poly_eq (a, a);
}

#endif
