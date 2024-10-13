#include "abdlop-params1.h"
#include "abdlop-params2.h"
#include "abdlop-params3.h"
#include "abdlop-params4.h"
#include "lazer.h"
#include "test.h"

static void test_abdlop (uint8_t seed[32], const abdlop_params_t params);

int
main (void)
{
  unsigned int i;
  uint8_t seed[32];

  lazer_init();

  for (i = 0; i < 2; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_abdlop (seed, params1);
    }
  for (i = 0; i < 2; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_abdlop (seed, params2);
    }
  for (i = 0; i < 2; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_abdlop (seed, params3);
    }
  for (i = 0; i < 2; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_abdlop (seed, params4);
    }

  TEST_PASS ();
}

static void
test_abdlop (uint8_t seed[32], const abdlop_params_t params)
{
  uint8_t hashp[32] = { 0 };
  uint8_t hashv[32] = { 0 };
  uint8_t hashcomm[32] = { 0 };
  polyring_srcptr Rq = params->ring;
  INT_T (lo, Rq->q->nlimbs);
  INT_T (hi, Rq->q->nlimbs);
  uint8_t buf[2];
  uint32_t dom;
  unsigned int i;
  int b;
  polymat_t A1, A2prime, Bprime, A1err, A2primeerr;
  polyvec_t z1err, z21err, herr, tA1err, s1, s2, m, tA1, tA2, tB, z1, z21, h;
  poly_t c;
  const unsigned int l = params->l + params->lext;

  poly_alloc (c, Rq);
  polyvec_alloc (z1err, Rq, params->m1);
  polyvec_alloc (z21err, Rq, params->m2 - params->kmsis);
  polyvec_alloc (herr, Rq, params->kmsis);
  polyvec_alloc (tA1err, Rq, params->kmsis);
  polyvec_alloc (s1, Rq, params->m1);
  polyvec_alloc (s2, Rq, params->m2);
  if (l > 0)
    polyvec_alloc (m, Rq, params->l + params->lext);
  polyvec_alloc (tA1, Rq, params->kmsis);
  polyvec_alloc (tA2, Rq, params->kmsis);
  if (l > 0)
    polyvec_alloc (tB, Rq, params->l + params->lext);
  polyvec_alloc (z1, Rq, params->m1);
  polyvec_alloc (z21, Rq, params->m2 - params->kmsis);
  polyvec_alloc (h, Rq, params->kmsis);
  polymat_alloc (A1, Rq, params->kmsis, params->m1);
  polymat_alloc (A2prime, Rq, params->kmsis, params->m2 - params->kmsis);
  if (l > 0)
    polymat_alloc (Bprime, Rq, params->l + params->lext,
                   params->m2 - params->kmsis);
  polymat_alloc (A1err, Rq, params->kmsis, params->m1);
  polymat_alloc (A2primeerr, Rq, params->kmsis, params->m2 - params->kmsis);

  dom = 0;
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s1, lo, hi, seed, dom++);
  polyvec_urandom_bnd (s2, lo, hi, seed, dom++);
  if (l > 0)
    polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

  /* generate public parameters */

  abdlop_keygen (A1, A2prime, Bprime, seed, params);

  /* generate proof */

  memcpy (hashp, seed, 32);
  abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, params);
  abdlop_hashcomm (hashp, tA1, tB, params);
  abdlop_prove (hashp, c, z1, z21, h, tA2, s1, s2, A1, A2prime, seed, params);

  /* expect successful verification */

  memcpy (hashv, seed, 32);
  abdlop_hashcomm (hashv, tA1, tB, params);
  b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
  TEST_EXPECT (b == 1);
  TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);

  memcpy (hashcomm, seed, 32);
  abdlop_hashcomm (hashcomm, tA1, tB, params);

  for (i = 0; i < 10; i++)
    {
      /* expect verification failures */

      bytes_urandom (buf, sizeof (buf));
      memcpy (hashv, hashcomm, 32);
      hashv[buf[0] % 32] ^= (1 << (buf[1] % 8));
      b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z1err, 1, seed, dom++);
      polyvec_add (z1err, z1err, z1, 0);
      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1err, z21, h, tA1, A1, A2prime, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z21err, 1, seed, dom++);
      polyvec_add (z21err, z21err, z21, 0);
      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1, z21err, h, tA1, A1, A2prime, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (herr, 1, seed, dom++);
      polyvec_add (herr, herr, h, 0);
      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1, z21, herr, tA1, A1, A2prime, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (tA1err, 1, seed, dom++); /* sometimes fails XXX */
      polyvec_add (tA1err, tA1err, tA1, 0);
      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1, z21, h, tA1err, A1, A2prime, params);
      TEST_EXPECT (b == 0);

      polymat_brandom (A1err, 1, seed, dom++);
      polymat_add (A1err, A1err, A1, 0);
      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1err, A2prime, params);
      TEST_EXPECT (b == 0);

      polymat_brandom (A2primeerr, 1, seed, dom++);
      polymat_add (A2primeerr, A2primeerr, A2prime, 0);
      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2primeerr, params);
      TEST_EXPECT (b == 0);

      /* expect successful verification */

      memcpy (hashv, hashcomm, 32);
      b = abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
      TEST_EXPECT (b == 1);
      TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
    }

  poly_free (c);
  polyvec_free (z1err);
  polyvec_free (z21err);
  polyvec_free (herr);
  polyvec_free (tA1err);
  polyvec_free (s1);
  polyvec_free (s2);
  if (l > 0)
    polyvec_free (m);
  polyvec_free (tA1);
  polyvec_free (tA2);
  if (l > 0)
    polyvec_free (tB);
  polyvec_free (z1);
  polyvec_free (z21);
  polyvec_free (h);
  polymat_free (A1);
  polymat_free (A2prime);
  if (l > 0)
    polymat_free (Bprime);
  polymat_free (A1err);
  polymat_free (A2primeerr);
}
