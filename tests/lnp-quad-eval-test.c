#include "lazer.h"
#include "lnp-quad-eval-params1.h"
#include "lnp-quad-eval-params2.h"
#include "test.h"
#include <mpfr.h>

#define N 3 /* number of quadratic equations */
#define M 3 /* number of quadratic eval equations */

static void test_lnp_quad_eval (uint8_t seed[32],
                                const lnp_quad_eval_params_t params);

int
main (void)
{
  unsigned int i;
  uint8_t seed[32] = { 0 };

  lazer_init ();

  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_lnp_quad_eval (seed, params1);
    }
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_lnp_quad_eval (seed, params2);
    }

  mpfr_free_cache ();
  TEST_PASS ();
}

static void
test_lnp_quad_eval (uint8_t seed[32], const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr abdlop = params->quad_eval;
  uint8_t hashp[32] = { 0 };
  uint8_t hashv[32] = { 0 };
  polyring_srcptr Rq = abdlop->ring;
  const unsigned int lambda = params->lambda;
  const unsigned int N_ = lambda / 2;
  INT_T (lo, Rq->q->nlimbs);
  INT_T (hi, Rq->q->nlimbs);
  int b;
  uint8_t buf[2];
  uint32_t dom;
  unsigned int i, j, k;
  spolymat_t R2i[N + lambda / 2], Rprime2i[M];
  spolyvec_t r1i[N + lambda / 2], rprime1i[M];
  poly_t r0i[N + lambda / 2], rprime0i[M];
  polyvec_t asub, asub_auto, bsub, bsub_auto, subv;
  int_ptr coeff;
  polymat_t A1err, A2primeerr, A1, A2prime, Bprime, Bprimeerr;
  polyvec_t s1, s2, m, tA1, tA2, tB, tBerr, z1, z21, hint, h, s, tmp, z1err,
      z21err, hinterr, tA1err, herr;
  poly_t r0err, rprime0err, c, cerr;
  spolymat_ptr R2[N + lambda / 2], Rprime2[M];
  spolyvec_ptr r1[N + lambda / 2], rprime1[M];
  poly_ptr r0[N + lambda / 2], rprime0[M];
  poly_ptr poly;
  spolyvec_t r1err, r1err_, rprime1err, rprime1err_;
  spolymat_t R2err, Rprime2err, R2err_, Rprime2err_;
  const unsigned int n = 2 * (abdlop->m1 + abdlop->l) + params->lambda;
  const unsigned int np = 2 * (abdlop->m1 + abdlop->l);

  dom = 0;

  poly_alloc (r0err, Rq);
  poly_alloc (rprime0err, Rq);
  poly_alloc (c, Rq);
  poly_alloc (cerr, Rq);
  polyvec_alloc (s1, Rq, abdlop->m1);
  polyvec_alloc (s2, Rq, abdlop->m2);
  polyvec_alloc (m, Rq, abdlop->l + params->lambda / 2 + 1);
  polyvec_alloc (tA1, Rq, abdlop->kmsis);
  polyvec_alloc (tA2, Rq, abdlop->kmsis);
  polyvec_alloc (tB, Rq, abdlop->l + abdlop->lext);
  polyvec_alloc (tBerr, Rq, abdlop->l + abdlop->lext);
  polyvec_alloc (z1, Rq, abdlop->m1);
  polyvec_alloc (z21, Rq, abdlop->m2 - abdlop->kmsis);
  polyvec_alloc (hint, Rq, abdlop->kmsis);
  polyvec_alloc (h, Rq, params->lambda / 2);
  polyvec_alloc (s, Rq, 2 * (abdlop->m1 + abdlop->l));
  polyvec_alloc (tmp, Rq, 2 * (abdlop->m1 + abdlop->l));
  spolyvec_alloc (r1err, Rq, n, n);
  spolyvec_alloc (r1err_, Rq, n, n);
  spolyvec_alloc (rprime1err, Rq, np, np);
  spolyvec_alloc (rprime1err_, Rq, np, np);
  polyvec_alloc (herr, Rq, params->lambda / 2);
  polyvec_alloc (z1err, Rq, abdlop->m1);
  polyvec_alloc (z21err, Rq, abdlop->m2 - abdlop->kmsis);
  polyvec_alloc (hinterr, Rq, abdlop->kmsis);
  polyvec_alloc (tA1err, Rq, abdlop->kmsis);
  polymat_alloc (A1err, Rq, abdlop->kmsis, abdlop->m1);
  polymat_alloc (A2primeerr, Rq, abdlop->kmsis, abdlop->m2 - abdlop->kmsis);
  spolymat_alloc (R2err, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (Rprime2err, Rq, np, np, (np * np - np) / 2 + np);
  spolymat_alloc (R2err_, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (Rprime2err_, Rq, np, np, (np * np - np) / 2 + np);
  polymat_alloc (A1, Rq, abdlop->kmsis, abdlop->m1);
  polymat_alloc (A2prime, Rq, abdlop->kmsis, abdlop->m2 - abdlop->kmsis);
  polymat_alloc (Bprime, Rq, abdlop->l + abdlop->lext,
                 abdlop->m2 - abdlop->kmsis);
  polymat_alloc (Bprimeerr, Rq, abdlop->l + abdlop->lext,
                 abdlop->m2 - abdlop->kmsis);
  for (i = 0; i < N + lambda / 2; i++)
    {
      spolymat_alloc (R2i[i], Rq, n, n, (np * np - np) / 2 + np);
      R2[i] = R2i[i];
      spolyvec_alloc (r1i[i], Rq, n, n);
      r1[i] = r1i[i];
      poly_alloc (r0i[i], Rq);
      r0[i] = r0i[i];

      for (j = 0; j < np; j++)
        {
          for (k = j; k < np; k++)
            {
              poly = spolymat_insert_elem (R2i[i], j, k);
              poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
            }
          poly = spolyvec_insert_elem (r1i[i], j);
          poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
        }
      spolyvec_sort (r1i[i]);
      spolymat_sort (R2i[i]);
    }
  for (i = 0; i < M; i++)
    {
      spolymat_alloc (Rprime2i[i], Rq, np, np, (np * np - np) / 2 + np);
      Rprime2[i] = Rprime2i[i];
      spolyvec_alloc (rprime1i[i], Rq, np, np);
      rprime1[i] = rprime1i[i];
      poly_alloc (rprime0i[i], Rq);
      rprime0[i] = rprime0i[i];

      for (j = 0; j < np; j++)
        {
          for (k = j; k < np; k++)
            {
              poly = spolymat_insert_elem (Rprime2i[i], j, k);
              poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
            }
          poly = spolyvec_insert_elem (rprime1i[i], j);
          poly_urandom (poly, Rq->q, Rq->log2q, seed, dom++);
        }
      spolyvec_sort (rprime1i[i]);
      spolymat_sort (Rprime2i[i]);
    }
  for (j = 0; j < np; j++)
    {
      for (k = j; k < np; k++)
        {
          spolymat_insert_elem (R2err, j, k);
          spolymat_insert_elem (Rprime2err, j, k);
        }
      spolyvec_insert_elem (r1err, j);
      spolyvec_insert_elem (rprime1err, j);
    }
  spolyvec_sort (r1err);
  spolyvec_sort (rprime1err);
  spolymat_sort (R2err);
  spolymat_sort (Rprime2err);

  int_set_i64 (lo, -3);
  int_set_i64 (hi, 3);
  polyvec_urandom_bnd (s1, lo, hi, seed, dom++);
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s2, lo, hi, seed, dom++);
  polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

  /* s = (<s1>,<m>) */

  polyvec_get_subvec (asub, s, 0, abdlop->m1, 2);
  polyvec_get_subvec (asub_auto, s, 1, abdlop->m1, 2);
  polyvec_set (asub, s1);
  polyvec_auto (asub_auto, s1);
  if (abdlop->l > 0)
    {
      polyvec_get_subvec (bsub, s, abdlop->m1 * 2, abdlop->l, 2);
      polyvec_get_subvec (bsub_auto, s, abdlop->m1 * 2 + 1, abdlop->l, 2);
      polyvec_get_subvec (subv, m, 0, abdlop->l, 1);
      polyvec_set (bsub, subv);
      polyvec_auto (bsub_auto, subv);
    }

  /* generate quadratic equations (in s) randomly */

  for (i = N_; i < N_ + N; i++)
    {
      /* R2, r1 already randomized */

      polyvec_dot2 (r0[i], r1[i], s);
      polyvec_mulsparse (tmp, R2i[i], s);
      polyvec_fromcrt (tmp);
      poly_adddot (r0[i], s, tmp, 0);
      poly_neg_self (r0[i]);
      poly_fromcrt (r0[i]);
    }

  for (i = 0; i < M; i++)
    {
      /* R2' already randomized */
      spolyvec_urandom (rprime1[i], Rq->q, Rq->log2q, seed, dom++);

      polyvec_dot2 (rprime0[i], rprime1[i], s);
      polyvec_mulsparse (tmp, Rprime2[i], s);
      polyvec_fromcrt (tmp);
      poly_adddot (rprime0[i], s, tmp, 0);
      poly_neg_self (rprime0[i]);
      poly_fromcrt (rprime0[i]);

      /* only constant coeff needs to be zero */
      poly_brandom (r0err, 1, seed, dom++);
      coeff = poly_get_coeff (r0err, 0);
      int_set_i64 (coeff, 0);

      poly_add (rprime0[i], rprime0[i], r0err, 0);
    }

  /* generate public parameters */

  abdlop_keygen (A1, A2prime, Bprime, seed, abdlop);

  /* generate proof */

  memset (hashp, 0xff, 32);
  abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, abdlop);
  lnp_quad_eval_prove (hashp, tB, h, c, z1, z21, hint, s1, m, s2, tA2, A1,
                       A2prime, Bprime, R2, r1, N, Rprime2, rprime1, rprime0,
                       M, seed, params);

  /* expect successful verification */

  memset (hashv, 0xff, 32);
  b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1, A2prime,
                            Bprime, R2, r1, r0, N, Rprime2, rprime1, rprime0,
                            M, params);
  TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
  TEST_EXPECT (b == 1);

  for (i = 0; i < 1; i++)
    {
      /* expect verification failures */

      bytes_urandom (buf, sizeof (buf));
      memset (hashv, 0xff, 32);
      hashv[buf[0] % 32] ^= (1 << (buf[1] % 8));
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (herr, 1, seed, dom++);
      polyvec_add (herr, herr, h, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, herr, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      poly_brandom (cerr, 1, seed, dom++);
      poly_add (cerr, cerr, c, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, cerr, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z1err, 1, seed, dom++);
      polyvec_add (z1err, z1err, z1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1err, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z21err, 1, seed, dom++);
      polyvec_add (z21err, z21err, z21, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21err, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (hinterr, 1, seed, dom++);
      polyvec_add (hinterr, hinterr, hint, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hinterr, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (tA1err, 1, seed, dom++);
      polyvec_add (tA1err, tA1err, tA1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1err, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      if (abdlop->l > 0)
        {
          polyvec_brandom (tBerr, 1, seed, dom++);
          polyvec_add (tBerr, tBerr, tB, 0);
          memset (hashv, 0xff, 32);
          b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tBerr, A1,
                                    A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                    rprime1, rprime0, M, params);
          TEST_EXPECT (b == 0);
        }

      polymat_brandom (A1err, 1, seed, dom++);
      polymat_add (A1err, A1err, A1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1err,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polymat_brandom (A2primeerr, 1, seed, dom++);
      polymat_add (A2primeerr, A2primeerr, A2prime, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2primeerr, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      polymat_brandom (Bprimeerr, 1, seed, dom++);
      polymat_add (Bprimeerr, Bprimeerr, Bprime, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprimeerr, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 0);

      spolymat_brandom (R2err, 1, seed, dom++);
      spolymat_add (R2err_, R2[0], R2err, 0);
      R2[N_] = R2err_;
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      R2[N_] = R2i[N_];
      TEST_EXPECT (b == 0);

      spolyvec_brandom (r1err, 1, seed, dom++);
      spolyvec_add (r1err_, r1[1], r1err, 0);
      r1[N_ + 1] = r1err_;
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      r1[N_ + 1] = r1i[N_ + 1];
      TEST_EXPECT (b == 0);

      poly_brandom (r0err, 1, seed, dom++);
      poly_add (r0[N_ + 2], r0[N_ + 2], r0err, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      poly_sub (r0[N_ + 2], r0[N_ + 2], r0err, 0);
      TEST_EXPECT (b == 0);

      spolymat_brandom (Rprime2err, 1, seed, dom++);
      spolymat_add (Rprime2err_, Rprime2[1], Rprime2err, 0);
      Rprime2[1] = Rprime2err_;
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      Rprime2[1] = Rprime2i[1];
      TEST_EXPECT (b == 0);

      spolyvec_brandom (rprime1err, 1, seed, dom++);
      spolyvec_add (rprime1err_, rprime1[1], rprime1err, 0);
      rprime1[1] = rprime1err_;
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      rprime1[1] = rprime1i[1];
      TEST_EXPECT (b == 0);

      poly_brandom (rprime0err, 1, seed, dom++);
      poly_add (rprime0[2], rprime0[2], rprime0err, 0);
      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      poly_sub (rprime0[2], rprime0[2], rprime0err, 0);
      TEST_EXPECT (b == 0);

      /* expect successful verification */

      memset (hashv, 0xff, 32);
      b = lnp_quad_eval_verify (hashv, h, c, z1, z21, hint, tA1, tB, A1,
                                A2prime, Bprime, R2, r1, r0, N, Rprime2,
                                rprime1, rprime0, M, params);
      TEST_EXPECT (b == 1);
      TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
    }

  poly_free (r0err);
  poly_free (rprime0err);
  poly_free (c);
  poly_free (cerr);
  polyvec_free (s1);
  polyvec_free (s2);
  polyvec_free (m);
  polyvec_free (tA1);
  polyvec_free (tA2);
  polyvec_free (tB);
  polyvec_free (tBerr);
  polyvec_free (z1);
  polyvec_free (z21);
  polyvec_free (hint);
  polyvec_free (h);
  polyvec_free (s);
  polyvec_free (tmp);
  spolyvec_free (r1err);
  spolyvec_free (r1err_);
  spolyvec_free (rprime1err);
  spolyvec_free (rprime1err_);
  polyvec_free (herr);
  polyvec_free (z1err);
  polyvec_free (z21err);
  polyvec_free (hinterr);
  polyvec_free (tA1err);
  polymat_free (A1err);
  polymat_free (A2primeerr);
  spolymat_free (R2err);
  spolymat_free (R2err_);
  spolymat_free (Rprime2err);
  spolymat_free (Rprime2err_);
  polymat_free (A1);
  polymat_free (A2prime);
  polymat_free (Bprime);
  polymat_free (Bprimeerr);
  for (i = 0; i < N + lambda / 2; i++)
    {
      spolymat_free (R2i[i]);
      spolyvec_free (r1i[i]);
      poly_free (r0i[i]);
    }
  for (i = 0; i < M; i++)
    {
      spolymat_free (Rprime2i[i]);
      spolyvec_free (rprime1i[i]);
      poly_free (rprime0i[i]);
    }
}
