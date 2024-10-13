#include "lazer.h"
#include "stopwatch.h"

#ifdef XXX
/*
 * r = U^T*auto(a) = U*auto(a)
 * for each dim 2 subvec:
 * (a,b) -> auto((b,a))
 */
static void
_shuffleautovec (polyvec_t r, polyvec_t a)
{
  poly_ptr rp, ap;
  unsigned int i;

  ASSERT_ERR (r != a);
  ASSERT_ERR (polyvec_get_ring (r) == polyvec_get_ring (a));
  ASSERT_ERR (polyvec_get_nelems (r) == polyvec_get_nelems (a));
  ASSERT_ERR (polyvec_get_nelems (r) % 2 == 0);

  for (i = 0; i < polyvec_get_nelems (r); i++)
    {
      rp = polyvec_get_elem (r, i);
      ap = polyvec_get_elem (a, FLOOR (i, 2) * 2 + (1 - i % 2));
      poly_auto (rp, ap);
    }
}
#endif

/*
 * r = U^T*auto(a) = U*auto(a)
 * for each dim 2 subvec:
 * (a,b) -> auto((b,a))
 */
static void
_shuffleautovecsparse (spolyvec_t r)
{
  poly_ptr rp;
  unsigned int i, elem;

  _SVEC_FOREACH_ELEM (r, i)
  {
    rp = spolyvec_get_elem (r, i);
    elem = spolyvec_get_elem_ (r, i);

    poly_auto_self (rp);
    spolyvec_set_elem_ (r, i, elem % 2 == 0 ? elem + 1 : elem - 1);
  }

  r->sorted = 0; // XXX simpler sort possible
  spolyvec_sort (r);
}

/*
 *
 * r = U^T*auto(a)*U = U*auto(a)*U
 * for each 2x2 submat on or above the main diagonal:
 * [[a,b],[c,d]] -> auto([[d,c],[b,a]])
 * r != a
 */
static void
_shuffleauto2x2submatssparse (spolymat_t a)
{
  poly_ptr ap;
  unsigned int i, arow, acol;

  ASSERT_ERR (spolymat_get_nrows (a) % 2 == 0);
  ASSERT_ERR (spolymat_get_ncols (a) % 2 == 0);
  ASSERT_ERR (spolymat_is_upperdiag (a));

  _SMAT_FOREACH_ELEM (a, i)
  {
    ap = spolymat_get_elem (a, i);
    arow = spolymat_get_row (a, i);
    acol = spolymat_get_col (a, i);

    if (arow % 2 == 0 && acol % 2 == 0)
      {
        spolymat_set_row (a, i, arow + 1);
        spolymat_set_col (a, i, acol + 1);
        poly_auto_self (ap);
      }
    else if (arow % 2 == 1 && acol % 2 == 1)
      {
        spolymat_set_row (a, i, arow - 1);
        spolymat_set_col (a, i, acol - 1);
        poly_auto_self (ap);
      }
    else if (arow % 2 == 1 && acol % 2 == 0)
      {
        spolymat_set_row (a, i, arow - 1);
        spolymat_set_col (a, i, acol + 1);
        poly_auto_self (ap);
      }
    else
      {
        /*
         * arow % 2 == 0 && acol % 2 == 1
         * This element's automorphism may land in the subdiagonal -1
         * if the 2x2 submat is on the main diagonal.
         * Check for this case and keep the matrix upper diagonal.
         */
        if (arow + 1 > acol - 1)
          {
            poly_auto_self (ap);
          }
        else
          {
            spolymat_set_row (a, i, arow + 1);
            spolymat_set_col (a, i, acol - 1);
            poly_auto_self (ap);
          }
      }
  }
  a->sorted = 0;
  spolymat_sort (a);
  ASSERT_ERR (spolymat_is_upperdiag (a));
}

static void
_schwartz_zippel_int (spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
                      spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                      poly_ptr rprime0i[], unsigned int M, polyvec_t h,
                      const intmat_t v, const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int lambda = params->lambda;
  const unsigned int N_ = lambda / 2;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u1, u2, u3;
  spolymat_t t0, t1, t2, t3;
  polyvec_t tmp;
  poly_ptr poly, hi;
  poly_t tpoly;
  unsigned int i, j;

  poly_alloc (tpoly, Rq);
  polyvec_alloc (tmp, Rq, 2 * (m1 + l));

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u1, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolyvec_alloc (u3, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t1, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t3, Rq, n, n, (n * n - n) / 2 + n);

  /* compute R2i, r1i, r0i for lambda/2 additional equations */
  STOPWATCH_START (stopwatch_lnp_quad_eval_schwartz_zippel_quad,
                   "lnp_quad_eval_prove_schwartz_zippel_quad");
  for (i = 0; i < N_; i++)
    {
      /* R2i */

      spolymat_set_empty (t0);
      for (j = 0; j < M; j++)
        {
          if (Rprime2i[j] == NULL)
            continue;

          spolymat_fromcrt (Rprime2i[j]);

          spolymat_set (t1, Rprime2i[j]);
          _shuffleauto2x2submatssparse (t1);
          spolymat_add (t2, Rprime2i[j], t1, 0);

          spolymat_scale (t3, intmat_get_elem (v, 2 * i, j), t2);
          spolymat_add (t1, t0, t3, 0);
          spolymat_set (t0, t1);

          spolymat_lrot (t1, t2, d / 2);

          spolymat_scale (t3, intmat_get_elem (v, 2 * i + 1, j), t1);
          spolymat_add (t1, t0, t3, 0);
          spolymat_set (t0, t1);
        }
      spolymat_scale (t1, Rq->inv2, t0);
      spolymat_mod (R2i[i], t1);
    }
  STOPWATCH_STOP (stopwatch_lnp_quad_eval_schwartz_zippel_quad);

  STOPWATCH_START (stopwatch_lnp_quad_eval_schwartz_zippel_lin,
                   "lnp_quad_eval_prove_schwartz_zippel_lin");
  for (i = 0; i < N_; i++)
    {
      /* r1i */

#ifdef XXX
      polyvec_get_subvec (subv, r1i[i], 0, 2 * (m1 + l), 1);
      polyvec_set_zero (subv);
      for (j = 0; j < M; j++)
        {
          if (rprime1i[j] == NULL)
            continue;

          spolyvec_fromcrt (rprime1i[j]);

          _shuffleautovec (tmp, rprime1i[j]);
          polyvec_add (tmp, rprime1i[j], tmp, 0);

          polyvec_addscale (subv, intmat_get_elem (v, 2 * i, j), tmp, 0);
          polyvec_lrot (tmp, tmp, d / 2);
          polyvec_addscale (subv, intmat_get_elem (v, 2 * i + 1, j), tmp, 0);
        }
      polyvec_scale (subv, Rq->inv2, subv);
      polyvec_mod (subv, subv);
#endif
      spolyvec_set_empty (u0);
      for (j = 0; j < M; j++)
        {
          if (rprime1i[j] == NULL)
            continue;

          spolyvec_fromcrt (rprime1i[j]);

          spolyvec_set (u1, rprime1i[j]);
          _shuffleautovecsparse (u1);
          spolyvec_add (u2, rprime1i[j], u1, 0);

          spolyvec_scale (u3, intmat_get_elem (v, 2 * i, j), u2);
          spolyvec_add (u1, u0, u3, 0);
          spolyvec_set (u0, u1);

          spolyvec_lrot (u1, u2, d / 2);

          spolyvec_scale (u3, intmat_get_elem (v, 2 * i + 1, j), u1);
          spolyvec_add (u1, u0, u3, 0);
          spolyvec_set (u0, u1);
        }
      spolyvec_scale (u1, Rq->inv2, u0);
      spolyvec_mod (r1i[i], u1);

#ifdef XXX
      polyvec_get_subvec (subv, r1i[i], 2 * (m1 + l), lambda, 1);
      for (j = 0; j < lambda; j++)
        {
          poly = polyvec_get_elem (subv, j);
          if (j == 2 * i)
            poly_set_one (poly);
          else
            poly_set_zero (poly);
        }
#endif
      for (j = 0; j < lambda; j++)
        {
          poly = spolyvec_insert_elem (r1i[i], 2 * (m1 + l) + j);
          if (j == 2 * i)
            poly_set_one (poly);
          else
            poly_set_zero (poly);
        }
      r1i[i]->sorted = 1; /* above for loop appends */
    }
  STOPWATCH_STOP (stopwatch_lnp_quad_eval_schwartz_zippel_lin);

  STOPWATCH_START (stopwatch_lnp_quad_eval_schwartz_zippel_const,
                   "lnp_quad_eval_prove_schwartz_zippel_const");
  if (r0i != NULL)
    {
      for (i = 0; i < N_; i++)
        {
          /* r0i = -hi + sum_(J in [0,M-1]) (tmp2[i,j]*Tr(rprime0j */

          hi = polyvec_get_elem (h, i);

          poly_set_zero (r0i[i]);
          for (j = 0; j < M; j++)
            {
              if (rprime0i[j] == NULL)
                continue;

              poly_fromcrt (rprime0i[j]);

              poly_auto (tpoly, rprime0i[j]);
              poly_add (tpoly, tpoly, rprime0i[j], 1);

              poly_addscale (r0i[i], intmat_get_elem (v, 2 * i, j), tpoly, 0);
              poly_lrot (tpoly, tpoly, d / 2);
              poly_addscale (r0i[i], intmat_get_elem (v, 2 * i + 1, j), tpoly,
                             0);
            }
          poly_scale (r0i[i], Rq->inv2, r0i[i]);
          poly_mod (r0i[i], r0i[i]);

          poly_sub (r0i[i], r0i[i], hi, 0);
        }
    }
  STOPWATCH_STOP (stopwatch_lnp_quad_eval_schwartz_zippel_const);

#ifdef XXX
  /* zero R1i scratch space */
  for (i = N_; i < N_ + N; i++)
    {
      if (r1i[i] == NULL)
        continue;

      polyvec_get_subvec (subv, r1i[i], 2 * (m1 + l), lambda, 1);
      polyvec_set_zero (subv);
    }
#endif

  spolyvec_free (u0);
  spolyvec_free (u1);
  spolyvec_free (u2);
  spolyvec_free (u3);
  spolymat_free (t0);
  spolymat_free (t1);
  spolymat_free (t2);
  spolymat_free (t3);
  poly_free (tpoly);
  polyvec_free (tmp);
}

static void
_lnp_quad_eval_prove (uint8_t hash[32], polyvec_t tB, polyvec_t h,
                      polyvec_t s1, polyvec_t m, polyvec_t s2,
                      polymat_t Bprime, spolymat_ptr R2i[], spolyvec_ptr r1i[],
                      spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                      poly_ptr rprime0i[], unsigned int M,
                      const uint8_t seed_quad_eval[32],
                      const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  const unsigned int lambda = params->lambda;
  polyring_srcptr Rq = quad_eval->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int m1 = quad_eval->m1;
  const unsigned int m2 = quad_eval->m2;
  const unsigned int l = quad_eval->l;
  const unsigned int kmsis = quad_eval->kmsis;
  const unsigned int n = 2 * (m1 + l);
  INTMAT_T (v, lambda, M, Rq->q->nlimbs);
  polyvec_t asub, bsub, asub_auto, bsub_auto, tg, subv, s2_, s, tr, rottr, tmp;
  spolyvec_t rprime1crt;
  spolymat_t Rprime2crt;
  shake128_state_t hstate;
  coder_state_t cstate;
  polymat_t Bextprime;
  poly_ptr poly, poly2, hi;
  intvec_ptr coeffs;
  intvec_t isubv;
  unsigned int i, j;
  size_t outlen;
  /* buff for encoding of tg */
  uint8_t out[CEIL (log2q * d * lambda / 2, 8) + 1];

  polyvec_alloc (s, Rq, n);
  polyvec_alloc (tr, Rq, M);
  polyvec_alloc (rottr, Rq, M);
  polyvec_alloc (tmp, Rq, n);
  spolyvec_alloc (rprime1crt, Rq, n, n);
  spolymat_alloc (Rprime2crt, Rq, n, n, (n * n - n) / 2 + n);

  /* tB = (tB_,tg,t) */

  polyvec_get_subvec (tg, tB, l, lambda / 2, 1);

  /* Bprime = (Bprime_,Bext,bext) */

  polymat_get_submat (Bextprime, Bprime, l, 0, lambda / 2, m2 - kmsis, 1, 1);

  /* s = (<s1>,<m>) */

  polyvec_get_subvec (asub, s, 0, m1, 2);
  polyvec_get_subvec (asub_auto, s, 1, m1, 2);
  polyvec_set (asub, s1);
  polyvec_auto (asub_auto, s1);

  if (l > 0)
    {
      polyvec_get_subvec (bsub, s, m1 * 2, l, 2);
      polyvec_get_subvec (bsub_auto, s, m1 * 2 + 1, l, 2);
      polyvec_get_subvec (subv, m, 0, l, 1);
      polyvec_set (bsub, subv);
      polyvec_auto (bsub_auto, subv);
    }

#if ASSERT == ASSERT_ENABLED
  for (j = 0; j < M; j++)
    {
      poly = polyvec_get_elem (tr, j);
      if (!(rprime0i[j] == NULL))
        poly_set (poly, rprime0i[j]);
      else
        poly_set_zero (poly);

      if (!(rprime1i[j] == NULL))
        poly_adddot2 (poly, rprime1i[j], s, 0);

      if (!(Rprime2i[j] == NULL))
        {
          polyvec_mulsparse (tmp, Rprime2i[j], s);
          polyvec_fromcrt (tmp);
          poly_adddot (poly, s, tmp, 0);
        }

      //XXX printf ("M=%u/%u\n", j,M); //XXX
      ASSERT_ERR (int_eqzero (poly_get_coeff (poly, 0)) == 1);
    }
#endif

  /* generate uniformly random h=g with zero constant coefficient. */

  for (i = 0; i < lambda / 2; i++)
    {
      poly = polyvec_get_elem (h, i);
      coeffs = poly_get_coeffvec (poly);
      intvec_get_subvec (isubv, coeffs, 1, d - 1, 1);

      intvec_set_elem_i64 (coeffs, 0, 0);
      intvec_urandom (isubv, q, log2q, seed_quad_eval, i);
    }

  /* append g to message m */
  polyvec_get_subvec (subv, m, l, lambda / 2, 1);
  polyvec_set (subv, h);

  /* tg = Bexptprime*s2 + g */

  polyvec_set (tg, h);
  polyvec_get_subvec (s2_, s2, 0, m2 - kmsis, 1);
  polyvec_addmul (tg, Bextprime, s2_, 0); // tg correct

  /* encode and hash tg */

  polyvec_mod (tg, tg);
  polyvec_redp (tg, tg);

  coder_enc_begin (cstate, out);
  coder_enc_urandom3 (cstate, tg, q, log2q);
  coder_enc_end (cstate);

  outlen = coder_get_offset (cstate);
  ASSERT_ERR (outlen % 8 == 0);
  ASSERT_ERR (outlen / 8 <= CEIL (log2q * d * lambda / 2, 8) + 1);
  outlen >>= 3; /* nbits to nbytes */

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, outlen);
  shake128_squeeze (hstate, hash, 32);
  shake128_clear (hstate);

  /* generate challenges v */

  intmat_urandom (v, q, log2q, hash, 0);

  STOPWATCH_START (stopwatch_lnp_quad_eval_prove_compute_h,
                   "lnp_quad_eval_prove_compute_h");

  /* compute h */

  for (j = 0; j < M; j++)
    {
      poly = polyvec_get_elem (tr, j);
      if (!(rprime0i[j] == NULL))
        poly_set (poly, rprime0i[j]);
      else
        poly_set_zero (poly);

      if (!(rprime1i[j] == NULL))
        {
          /* keep rprime1 in coeff-rep, by copying it */
          spolyvec_set (rprime1crt, rprime1i[j]);

          poly_adddot2 (poly, rprime1crt, s, 0);
        }

      if (!(Rprime2i[j] == NULL))
        {
          spolymat_set (Rprime2crt, Rprime2i[j]);

          polyvec_mulsparse (tmp, Rprime2crt, s);
          polyvec_fromcrt (tmp);
          poly_adddot (poly, s, tmp, 0);
        }

      poly_tracemap (poly, poly);
    }
  polyvec_lrot (rottr, tr, d / 2);

  for (i = 0; i < lambda / 2; i++)
    {
      hi = polyvec_get_elem (h, i); /* = gi */

      for (j = 0; j < M; j++)
        {
          poly = polyvec_get_elem (tr, j);
          poly2 = polyvec_get_elem (rottr, j);
          poly_addscale (hi, intmat_get_elem (v, 2 * i, j), poly, 0);
          poly_addscale (hi, intmat_get_elem (v, 2 * i + 1, j), poly2, 0);
        }
    }

  STOPWATCH_STOP (stopwatch_lnp_quad_eval_prove_compute_h);

  _schwartz_zippel_int (R2i, r1i, NULL, Rprime2i, rprime1i, rprime0i, M, h, v,
                        params);

  polyvec_free (s);
  polyvec_free (tr);
  polyvec_free (rottr);
  polyvec_free (tmp);
  spolyvec_free (rprime1crt);
  spolymat_free (Rprime2crt);
}

/*
 * hash hash of tA1, tB
 * tB = (tB-,tg,t)
 *
 * scratch space:
 * R2i, r1i, r0i are dim N + lambda/2 arrays.
 * R2i is 2*(m1+l)+lambda x 2*(m1+l)+lambda
 */
void
lnp_quad_eval_prove (uint8_t hash[32], polyvec_t tB, polyvec_t h, poly_t c,
                     polyvec_t z1, polyvec_t z21, polyvec_t hint, polyvec_t s1,
                     polyvec_t m, polyvec_t s2, polyvec_t tA2, polymat_t A1,
                     polymat_t A2prime, polymat_t Bprime, spolymat_ptr R2i[],
                     spolyvec_ptr r1i[], unsigned int N,
                     spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                     poly_ptr rprime0i[], unsigned int M,
                     const uint8_t seed[32],
                     const lnp_quad_eval_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  abdlop_params_srcptr quad_eval = params->quad_eval;
  const unsigned int lext = params->quad_eval->lext;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int m2 = quad_eval->m2;
  const unsigned int l = quad_eval->l;
  const unsigned int kmsis = quad_eval->kmsis;
  unsigned int i;
#endif
  abdlop_params_srcptr quad_many = params->quad_many;
  const unsigned int lambda = params->lambda;
  const unsigned int N_ = lambda / 2;
  shake128_state_t hstate;
  uint8_t expseed[64];
  const uint8_t *seed_quad_eval = expseed;
  const uint8_t *seed_quad_many = expseed + 32;

  STOPWATCH_START (stopwatch_lnp_quad_eval_prove, "lnp_quad_eval_prove");

  ASSERT_ERR (M > 0); /* use quad_many if no eval eq is needed. */
  ASSERT_ERR (lambda % 2 == 0);
  ASSERT_ERR (lext == lambda / 2 + 1);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polyvec_get_nelems (h) == lambda / 2);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_nelems (hint) == kmsis);
  ASSERT_ERR (polyvec_get_nelems (s1) == m1);
  ASSERT_ERR (polyvec_get_nelems (s2) == m2);
  ASSERT_ERR (polyvec_get_nelems (m) == l + lext);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (polyvec_get_nelems (tA2) == kmsis);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (polymat_get_ncols (Bprime) == m2 - kmsis);
#if ASSERT == ASSERT_ENABLED
  for (i = N_; i < N_ + N; i++)
    {
      ASSERT_ERR (i >= N || R2i[i] == NULL || spolymat_is_upperdiag (R2i[i]));
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_nrows (R2i[i]) == 2 * (m1 + l) + lambda);
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_ncols (R2i[i]) == 2 * (m1 + l) + lambda);
      ASSERT_ERR (r1i[i] == NULL
                  || r1i[i]->nelems_max == 2 * (m1 + l) + lambda);
    }
  for (i = 0; i < M; i++)
    {
      ASSERT_ERR (Rprime2i[i] == NULL || spolymat_is_upperdiag (Rprime2i[i]));
      ASSERT_ERR (Rprime2i[i] == NULL
                  || spolymat_get_nrows (Rprime2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (Rprime2i[i] == NULL
                  || spolymat_get_ncols (Rprime2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (rprime1i[i] == NULL
                  || rprime1i[i]->nelems_max == 2 * (m1 + l));
    }
#endif

  /*
   * Expand input seed into two seeds: one for quad_eval
   * and one for the sub-protocol quad_many.
   */
  shake128_init (hstate);
  shake128_absorb (hstate, seed, 32);
  shake128_squeeze (hstate, expseed, sizeof (expseed));
  shake128_clear (hstate);

  _lnp_quad_eval_prove (hash, tB, h, s1, m, s2, Bprime, R2i, r1i, Rprime2i,
                        rprime1i, rprime0i, M, seed_quad_eval, params);
  lnp_quad_many_prove (hash, tB, c, z1, z21, hint, s1, m, s2, tA2, A1, A2prime,
                       Bprime, R2i, r1i, N_ + N, seed_quad_many, quad_many);
  STOPWATCH_STOP (stopwatch_lnp_quad_eval_prove);
}

static int
_lnp_quad_eval_verify (uint8_t hash[32], polyvec_t h, polyvec_t tB,
                       spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
                       spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                       poly_ptr rprime0i[], unsigned int M,
                       const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  const unsigned int lambda = params->lambda;
  polyring_srcptr Rq = quad_eval->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  polyvec_t tg, s;
  INTMAT_T (v, lambda, M, Rq->q->nlimbs);
  unsigned int i;
  poly_ptr poly;
  int_ptr coeff;
  shake128_state_t hstate;
  coder_state_t cstate;
  size_t outlen;
  /* buff for encoding of tg */
  uint8_t out[CEIL (log2q * d * lambda / 2, 8) + 1];
  int b = 0;

  polyvec_alloc (s, Rq, 2 * (m1 + l));

  /* check if h's constant coeffs are zero. */
  for (i = 0; i < lambda / 2; i++)
    {
      poly = polyvec_get_elem (h, i);

      coeff = poly_get_coeff (poly, 0);
      if (int_eqzero (coeff) != 1)
        goto ret;
    }

  /* tB = (tB_,tg,t) */
  polyvec_get_subvec (tg, tB, l, lambda / 2, 1);

  /* encode and hash tg */

  polyvec_mod (tg, tg);
  polyvec_redp (tg, tg);

  coder_enc_begin (cstate, out);
  coder_enc_urandom3 (cstate, tg, q, log2q);
  coder_enc_end (cstate);

  outlen = coder_get_offset (cstate);
  ASSERT_ERR (outlen % 8 == 0);
  ASSERT_ERR (outlen / 8 <= CEIL (log2q * d * lambda / 2, 8) + 1);
  outlen >>= 3; /* nbits to nbytes */

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, outlen);
  shake128_squeeze (hstate, hash, 32);
  shake128_clear (hstate);

  /* generate challenges v */

  intmat_urandom (v, q, log2q, hash, 0);

  _schwartz_zippel_int (R2i, r1i, r0i, Rprime2i, rprime1i, rprime0i, M, h, v,
                        params);

  b = 1;
ret:
  polyvec_free (s);
  return b;
}

int
lnp_quad_eval_verify (uint8_t hash[32], polyvec_t h, poly_t c, polyvec_t z1,
                      polyvec_t z21, polyvec_t hint, polyvec_t tA1,
                      polyvec_t tB, polymat_t A1, polymat_t A2prime,
                      polymat_t Bprime, spolymat_ptr R2i[], spolyvec_ptr r1i[],
                      poly_ptr r0i[], unsigned int N, spolymat_ptr Rprime2i[],
                      spolyvec_ptr rprime1i[], poly_ptr rprime0i[],
                      unsigned int M, const lnp_quad_eval_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  abdlop_params_srcptr quad_eval = params->quad_eval;
  const unsigned int m2 = params->quad_eval->m2;
  const unsigned int kmsis = params->quad_eval->kmsis;
  const unsigned int lext = params->quad_eval->lext;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  unsigned int i;
#endif
  abdlop_params_srcptr quad_many = params->quad_many;
  const unsigned int lambda = params->lambda;
  const unsigned int N_ = lambda / 2;
  int b;

  STOPWATCH_START (stopwatch_lnp_quad_eval_verify, "lnp_quad_eval_verify");

  ASSERT_ERR (M > 0); /* for M=0 one may want to work with quad_many */
  ASSERT_ERR (lambda % 2 == 0);
  ASSERT_ERR (lext == lambda / 2 + 1);
  ASSERT_ERR (polyvec_get_nelems (h) == lambda / 2);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_nelems (hint) == kmsis);
  ASSERT_ERR (polyvec_get_nelems (tA1) == kmsis);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_nelems (h) == lambda / 2);
#if ASSERT == ASSERT_ENABLED
  for (i = 0; i < N + lambda / 2; i++)
    {
      ASSERT_ERR (i >= N || R2i[i] == NULL || spolymat_is_upperdiag (R2i[i]));
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_nrows (R2i[i]) == 2 * (m1 + l) + lambda);
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_ncols (R2i[i]) == 2 * (m1 + l) + lambda);
      ASSERT_ERR (r1i[i] == NULL
                  || r1i[i]->nelems_max == 2 * (m1 + l) + lambda);
    }
  for (i = 0; i < M; i++)
    {
      ASSERT_ERR (Rprime2i[i] == NULL || spolymat_is_upperdiag (Rprime2i[i]));
      ASSERT_ERR (Rprime2i[i] == NULL
                  || spolymat_get_nrows (Rprime2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (Rprime2i[i] == NULL
                  || spolymat_get_ncols (Rprime2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (rprime1i[i] == NULL
                  || rprime1i[i]->nelems_max == 2 * (m1 + l));
    }
#endif

  b = _lnp_quad_eval_verify (hash, h, tB, R2i, r1i, r0i, Rprime2i, rprime1i,
                             rprime0i, M, params);
  if (b != 1)
    goto ret;

  b = lnp_quad_many_verify (hash, c, z1, z21, hint, tA1, tB, A1, A2prime,
                            Bprime, R2i, r1i, r0i, N + N_, quad_many);
  if (b != 1)
    goto ret;

  b = 1;
ret:
  STOPWATCH_STOP (stopwatch_lnp_quad_eval_verify);
  return b;
}
