#include "lazer.h"
#include "lnp-quad-params1.h"
#include "lnp-quad-params2.h"
#include "lnp-quad-params3.h"
#include "test.h"
#include <mpfr.h>

#define N 3 /* number of quadratic equations */

static void test_lnp_quad_many (uint8_t seed[32],
                                const abdlop_params_t params);

int
main (void)
{
  unsigned int i;
  uint8_t seed[32] = { 0 };

  lazer_init ();

  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_lnp_quad_many (seed, params1);
    }
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_lnp_quad_many (seed, params2);
    }
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_lnp_quad_many (seed, params3);
    }

  mpfr_free_cache ();
  TEST_PASS ();
}

static void
test_lnp_quad_many (uint8_t seed[32], const abdlop_params_t params)
{
  uint8_t hashp[32] = { 0 };
  uint8_t hashv[32] = { 0 };
  polyring_srcptr Rq = params->ring;
  INT_T (lo, Rq->q->nlimbs);
  INT_T (hi, Rq->q->nlimbs);
  int b;
  uint8_t buf[2];
  uint32_t dom;
  unsigned int i, j, k;
  spolymat_t R2i[N];
  spolyvec_t r1i[N];
  poly_t r0i[N];
  spolyvec_t r1err, r1err_;
  spolymat_t R2err, R2err_;
  polymat_t A1, A2prime, Bprime, Bprimeerr, A1err, A2primeerr;
  polyvec_t s1, s2, m, tA1, tA2, tB, tBerr, z1, z21, h, s, tmp, z1err, z21err,
      herr, tA1err;
  poly_t r0err, c;
  polyvec_t asub, asub_auto, bsub, bsub_auto, subv;
  spolymat_ptr R2[N];
  spolyvec_ptr r1[N];
  poly_ptr r0[N];
  const unsigned int nelems = 2 * (params->m1 + params->l);

  poly_alloc (c, Rq);
  poly_alloc (r0err, Rq);
  polyvec_alloc (s1, Rq, params->m1);
  polyvec_alloc (s2, Rq, params->m2);
  polyvec_alloc (m, Rq, params->l + params->lext);
  polyvec_alloc (tA1, Rq, params->kmsis);
  polyvec_alloc (tA2, Rq, params->kmsis);
  polyvec_alloc (tB, Rq, params->l + params->lext);
  polyvec_alloc (tBerr, Rq, params->l + params->lext);
  polyvec_alloc (z1, Rq, params->m1);
  polyvec_alloc (z21, Rq, params->m2 - params->kmsis);
  polyvec_alloc (h, Rq, params->kmsis);
  polyvec_alloc (s, Rq, 2 * (params->m1 + params->l));
  polyvec_alloc (tmp, Rq, 2 * (params->m1 + params->l));
  polyvec_alloc (z1err, Rq, params->m1);
  polyvec_alloc (z21err, Rq, params->m2 - params->kmsis);
  polyvec_alloc (herr, Rq, params->kmsis);
  polyvec_alloc (tA1err, Rq, params->kmsis);
  polymat_alloc (A1, Rq, params->kmsis, params->m1);
  polymat_alloc (A2prime, Rq, params->kmsis, params->m2 - params->kmsis);
  polymat_alloc (Bprime, Rq, params->l + params->lext,
                 params->m2 - params->kmsis);
  polymat_alloc (Bprimeerr, Rq, params->l + params->lext,
                 params->m2 - params->kmsis);
  polymat_alloc (A1err, Rq, params->kmsis, params->m1);
  polymat_alloc (A2primeerr, Rq, params->kmsis, params->m2 - params->kmsis);
  spolymat_alloc (R2err, Rq, nelems, nelems,
                  (nelems * nelems - nelems) / 2 + nelems);
  spolymat_alloc (R2err_, Rq, nelems, nelems,
                  (nelems * nelems - nelems) / 2 + nelems);
  spolyvec_alloc (r1err, Rq, nelems, nelems);
  spolyvec_alloc (r1err_, Rq, nelems, nelems);
  for (i = 0; i < N; i++)
    {
      spolymat_alloc (R2i[i], Rq, nelems, nelems,
                      (nelems * nelems - nelems) / 2 + nelems);
      R2[i] = R2i[i];
      spolyvec_alloc (r1i[i], Rq, nelems, nelems);
      r1[i] = r1i[i];
      poly_alloc (r0i[i], Rq);
      r0[i] = r0i[i];

      for (j = 0; j < nelems; j++)
        {
          for (k = j; k < nelems; k++)
            spolymat_insert_elem (R2i[i], j, k);

          spolyvec_insert_elem (r1i[i], j);
        }
      spolyvec_sort (r1i[i]);
      spolymat_sort (R2i[i]);
    }
  for (j = 0; j < nelems; j++)
    {
      for (k = j; k < nelems; k++)
        spolymat_insert_elem (R2err, j, k);

      spolyvec_insert_elem (r1err, j);
    }
  spolyvec_sort (r1err);
  spolymat_sort (R2err);

  dom = 0;
  int_set_i64 (lo, -3);
  int_set_i64 (hi, 3);
  polyvec_urandom_bnd (s1, lo, hi, seed, dom++);
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s2, lo, hi, seed, dom++);
  polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

  /* s = (<s1>,<m>) */

  polyvec_get_subvec (asub, s, 0, params->m1, 2);
  polyvec_get_subvec (asub_auto, s, 1, params->m1, 2);
  polyvec_set (asub, s1);
  polyvec_auto (asub_auto, s1);

  if (params->l > 0)
    {
      polyvec_get_subvec (bsub, s, params->m1 * 2, params->l, 2);
      polyvec_get_subvec (bsub_auto, s, params->m1 * 2 + 1, params->l, 2);
      polyvec_get_subvec (subv, m, 0, params->l, 1);
      polyvec_set (bsub, subv);
      polyvec_auto (bsub_auto, subv);
    }

  /* generate quadratic equations (in s) randomly */
  for (i = 0; i < N; i++)
    {
      spolymat_urandom (R2[i], Rq->q, Rq->log2q, seed, dom++);
      spolyvec_urandom (r1[i], Rq->q, Rq->log2q, seed, dom++);

      polyvec_dot2 (r0[i], r1[i], s);
      polyvec_mulsparse (tmp, R2[i], s);
      polyvec_fromcrt (tmp);
      poly_adddot (r0[i], s, tmp, 0);
      poly_neg_self (r0[i]);
      poly_fromcrt (r0[i]);
    }

  /* generate public parameters */

  abdlop_keygen (A1, A2prime, Bprime, seed, params);

  /* generate proof */

  memset (hashp, 0xff, 32);
  abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, params);
  lnp_quad_many_prove (hashp, tB, c, z1, z21, h, s1, m, s2, tA2, A1, A2prime,
                       Bprime, R2, r1, N, seed, params);

  /* expect successful verification */

  memset (hashv, 0xff, 32);
  b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime, Bprime,
                            R2, r1, r0, N, params);
  TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
  TEST_EXPECT (b == 1);

  for (i = 0; i < 1; i++)
    {
      /* expect verification failures */

      bytes_urandom (buf, sizeof (buf));
      memset (hashv, 0xff, 32);
      hashv[buf[0] % 32] ^= (1 << (buf[1] % 8));
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z1err, 1, seed, dom++);
      polyvec_add (z1err, z1err, z1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1err, z21, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z21err, 1, seed, dom++);
      polyvec_add (z21err, z21err, z21, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21err, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (herr, 1, seed, dom++);
      polyvec_add (herr, herr, h, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, herr, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (tA1err, 1, seed, dom++);
      polyvec_add (tA1err, tA1err, tA1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1err, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      if (params->l > 0)
        {
          polyvec_brandom (tBerr, 1, seed, dom++);
          polyvec_add (tBerr, tBerr, tB, 0);
          memset (hashv, 0xff, 32);
          b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tBerr, A1,
                                    A2prime, Bprime, R2, r1, r0, N, params);
          TEST_EXPECT (b == 0);
        }

      polymat_brandom (A1err, 1, seed, dom++);
      polymat_add (A1err, A1err, A1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1err, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      polymat_brandom (A2primeerr, 1, seed, dom++);
      polymat_add (A2primeerr, A2primeerr, A2prime, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2primeerr,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 0);

      if (params->l > 0)
        {
          polymat_brandom (Bprimeerr, 1, seed, dom++);
          polymat_add (Bprimeerr, Bprimeerr, Bprime, 0);
          memset (hashv, 0xff, 32);
          b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime,
                                    Bprimeerr, R2, r1, r0, N, params);
          TEST_EXPECT (b == 0);
        }

      spolymat_brandom (R2err, 1, seed, dom++);
      spolymat_add (R2err_, R2[0], R2err, 0);
      R2[0] = R2err_;
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      R2[0] = R2i[0];
      TEST_EXPECT (b == 0);

      spolyvec_brandom (r1err, 1, seed, dom++);
      spolyvec_add (r1err_, r1[1], r1err, 0);
      r1[1] = r1err_;
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      r1[1] = r1i[1];
      TEST_EXPECT (b == 0);

      poly_brandom (r0err, 1, seed, dom++);
      poly_add (r0[2], r0[2], r0err, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      poly_sub (r0[2], r0[2], r0err, 0);
      TEST_EXPECT (b == 0);

      /* expect successful verification */

      memset (hashv, 0xff, 32);
      b = lnp_quad_many_verify (hashv, c, z1, z21, h, tA1, tB, A1, A2prime,
                                Bprime, R2, r1, r0, N, params);
      TEST_EXPECT (b == 1);
      TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
    }

  poly_free (c);
  poly_free (r0err);
  polyvec_free (s1);
  polyvec_free (s2);
  polyvec_free (m);
  polyvec_free (tA1);
  polyvec_free (tA2);
  polyvec_free (tB);
  polyvec_free (tBerr);
  polyvec_free (z1);
  polyvec_free (z21);
  polyvec_free (h);
  polyvec_free (s);
  polyvec_free (tmp);
  polyvec_free (z1err);
  polyvec_free (z21err);
  polyvec_free (herr);
  polyvec_free (tA1err);
  polymat_free (A1);
  polymat_free (A2prime);
  polymat_free (Bprime);
  polymat_free (Bprimeerr);
  polymat_free (A1err);
  polymat_free (A2primeerr);
  spolymat_free (R2err);
  spolymat_free (R2err_);
  spolyvec_free (r1err);
  spolyvec_free (r1err_);
  for (i = 0; i < N; i++)
    {
      spolymat_free (R2i[i]);
      spolyvec_free (r1i[i]);
      poly_free (r0i[i]);
    }
}
