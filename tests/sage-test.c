#include "abdlop-params1.h"
#include "lazer.h"
#include "test.h"

polyring_srcptr Rq = params1_ring;
int_srcptr mod = params1_q;

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

// XXX
// d=4
// q = 17
// Rq =
// QuotientRing(Integers(q)['x'],Integers(q)['x'].ideal(x^d+1),'x')
// vector(Rq,((2),(1,1)))
// matrix(Rq,[[(2),(3),(2,2)],[(1,1,1),(1,0,1),(2)]]

int
main (void)
{
  lazer_init();

  fprintf (stdout, "import sys\n");
  fprintf (stdout, "\n");

  fprintf (stdout, "q = ");
  int_out_str (stdout, 10, Rq->q);
  fprintf (stdout, "\n");

  fprintf (stdout, "d = %u\n", Rq->d);

  fprintf (stdout, "Rq = QuotientRing(Integers(q)['x'], "
                   "Integers(q)['x'].ideal(x^d+1),'x')\n");
  fprintf (stdout, "\n");

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

  fprintf (stdout, "sys.exit(int(%d))\n", TEST_RC_PASS);
  TEST_PASS ();
}

static void
test_poly_add (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  int crt;

  fprintf (stdout, "print ('poly_add')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  for (i = 0; i < 50; i++)
    {
      bytes_urandom (seed, sizeof (seed));

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      dom = 2;
      poly_urandom_bnd (b, lo, hi, seed, dom);
      fprintf (stdout, "b = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (b));
      fprintf (stdout, ")\n");

      if (seed[0] & 0x1)
        poly_tocrt (a);
      if (seed[0] & 0x2)
        poly_tocrt (b);
      crt = !!(seed[0] & 0x4);

      poly_add (r, a, b, crt);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != a+b:\n");
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_sub (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  int crt;

  fprintf (stdout, "print ('poly_sub')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  for (i = 0; i < 50; i++)
    {
      bytes_urandom (seed, sizeof (seed));

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      dom = 2;
      poly_urandom_bnd (b, lo, hi, seed, dom);
      fprintf (stdout, "b = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (b));
      fprintf (stdout, ")\n");

      if (seed[0] & 0x1)
        poly_tocrt (a);
      if (seed[0] & 0x2)
        poly_tocrt (b);
      crt = !!(seed[0] & 0x4);

      poly_sub (r, a, b, crt);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != a-b:\n");
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_mul (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);
  POLY_T (b, Rq);

  fprintf (stdout, "print ('poly_mul')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      dom = 2;
      poly_urandom_bnd (b, lo, hi, seed, dom);
      fprintf (stdout, "b = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (b));
      fprintf (stdout, ")\n");

      poly_tocrt (a);
      poly_tocrt (b);
      poly_mul (r, a, b);
      poly_fromcrt (r);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != a*b:\n");
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_scale (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  INT_T (s, Rq->q->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_scale')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q / 2 - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q / 2 - 1));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));

      dom = 1;
      int_urandom_bnd (s, lo, hi, seed, dom);
      fprintf (stdout, "s = ");
      int_out_str (stdout, 10, s);
      fprintf (stdout, "\n");

      dom = 2;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_scale (r, s, a);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != s*a:\n");
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_rshift (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i, shift;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  INT_T (s, Rq->q->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_rshift')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 18)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 18));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      shift = seed[1] % 5;

      int_set_i64 (s, 1 << shift);

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      poly_scale (a, s, a);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_rshift (r, a, shift);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != a/%d:\n", 1 << shift);
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_lshift (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i, shift;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_lshift')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 18)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 18));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      shift = seed[1] % 5;

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_lshift (r, a, shift);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != %d*a:\n", 1 << shift);
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_rrot (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i, rot;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_rrot')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      rot = seed[1] % 5;

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_rrot (r, a, rot);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != a/(Rq(x)^%u):\n", rot);
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_lrot (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i, rot;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_lrot')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      rot = seed[1] % 5;

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_lrot (r, a, rot);
      fprintf (stdout, "r = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");

      fprintf (stdout, "if r != a*Rq(x)^%u:\n", rot);
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_mod (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_mod')\n");

  int_set_i64 (lo, INT64_MIN);
  int_set_i64 (hi, INT64_MAX);

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_mod (r, a);
      fprintf (stdout, "vr = (");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");
      fprintf (stdout, "r = Rq(vr)\n");

      fprintf (stdout, "if r != a or max(max(vr),-min(vr)) >= q:\n");
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}

static void
test_poly_redc (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom;
  unsigned int i;
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (r, Rq);
  POLY_T (a, Rq);

  fprintf (stdout, "print ('poly_redc')\n");

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 2)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 2));

  for (i = 0; i < 20; i++)
    {
      bytes_urandom (seed, sizeof (seed));

      dom = 1;
      poly_urandom_bnd (a, lo, hi, seed, dom);
      poly_add (a, a, a, 0);
      fprintf (stdout, "a = Rq(");
      intvec_out_str (stdout, 10, poly_get_coeffvec (a));
      fprintf (stdout, ")\n");

      poly_redc (r, a);
      fprintf (stdout, "vr = (");
      intvec_out_str (stdout, 10, poly_get_coeffvec (r));
      fprintf (stdout, ")\n");
      fprintf (stdout, "r = Rq(vr)\n");

      fprintf (stdout, "if r != a or max(max(vr),-min(vr)) > (q-1)/2:\n");
      fprintf (stdout, "\tsys.exit(int(%d))\n", TEST_RC_FAIL);
      fprintf (stdout, "\n");
    }
}
