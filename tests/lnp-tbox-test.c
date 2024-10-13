#include "lazer.h"
#include "lnp-tbox-params1.h"
#include "test.h"
#include <mpfr.h>

#define N 1 /* number of quadratic equations */
#define M 4 /* number of quadratic eval equations */

static void test_lnp_tbox (uint8_t seed[32], const lnp_tbox_params_t params);

int
main (void)
{
  unsigned int i;
  uint8_t seed[32] = { 0 };
  seed[0] = 2;

  lazer_init();

  for (i = 0; i < 1; i++)
    {
      test_lnp_tbox (seed, tbox_params1);
    }

  mpfr_free_cache ();
  TEST_PASS ();
}

/* R2 != R2_ */
static void
_scatter_smat (spolymat_ptr R2, spolymat_ptr R2_, unsigned int m1,
               unsigned int Z, unsigned int l)
{
  const unsigned int nelems = R2_->nelems;
  unsigned int i, row, col;
  poly_ptr poly, poly2;

  ASSERT_ERR (R2->nelems_max >= R2_->nelems_max);
  ASSERT_ERR (spolymat_is_upperdiag (R2_));

  (void)l; /* unused */

  for (i = 0; i < nelems; i++)
    {
      poly = spolymat_get_elem (R2_, i);
      row = spolymat_get_row (R2_, i);
      col = spolymat_get_col (R2_, i);

      ASSERT_ERR (row < 2 * (m1 + l));
      ASSERT_ERR (col < 2 * (m1 + l));
      ASSERT_ERR (col >= row);

      if (col >= 2 * m1)
        col += 2 * Z;
      if (row >= 2 * m1)
        row += 2 * Z;

      poly2 = spolymat_insert_elem (R2, row, col);
      poly_set (poly2, poly);
    }
  R2->sorted = 0;
  spolymat_sort (R2);
  ASSERT_ERR (spolymat_is_upperdiag (R2));
}

/* r1, r1_ may not overlap */
static void
_scatter_vec (spolyvec_ptr r1, spolyvec_ptr r1_, unsigned int m1,
              unsigned int Z)
{
  const unsigned int nelems = r1_->nelems;
  unsigned int i, elem;
  poly_ptr poly, poly2;

  ASSERT_ERR (r1->nelems_max >= r1_->nelems_max);

  for (i = 0; i < nelems; i++)
    {
      poly = spolyvec_get_elem (r1_, i);
      elem = spolyvec_get_elem_ (r1_, i);

      if (elem >= 2 * m1)
        elem += 2 * Z;

      poly2 = spolyvec_insert_elem (r1, elem);
      poly_set (poly2, poly);
    }
  r1->sorted = 1;
}

static void
test_lnp_tbox (uint8_t seed[32], const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  abdlop_params_srcptr quad = params->quad_eval->quad_many;
  abdlop_params_srcptr quad_eval = params->quad_eval->quad_eval;
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int m1 = tbox->m1 - params->Z;
  const unsigned int Z = params->Z;
  const unsigned int l = tbox->l;
  uint8_t hashp[32] = { 0 };
  uint8_t hashv[32] = { 0 };
  INT_T (lo, Rq->q->nlimbs);
  INT_T (hi, Rq->q->nlimbs);
  uint32_t dom;
  unsigned int i, j, k, nelems, nrows, ncols;
  spolymat_t R2i[N];
  spolyvec_t r1i[N];
  poly_t r0i[N];
  spolymat_t R2primei[M];
  spolyvec_t r1primei[M];
  poly_t r0primei[M];
  polymat_t Esi[Z], Emi[Z];
  polyvec_t vi[Z], vi_[Z];
  spolymat_ptr R2[N];
  spolyvec_ptr r1[N];
  poly_ptr r0[N];
  spolymat_ptr R2prime[M];
  spolyvec_ptr r1prime[M];
  poly_ptr r0prime[M];
  polymat_ptr Es[Z], Em[Z];
  polyvec_ptr v[Z];
  polymat_t Ps, Pm;
  polyvec_t f, f_, tBerr, z1err, z21err, hinterr, tA1err;
  spolyvec_t r1err, r1err_, r1primeerr, r1primeerr_;
  polymat_t Ds, Dm;
  polyvec_t u, u_;
  polymat_t Bprime, A2prime, A1, Bprimeerr, A1err, A2primeerr, Es0err, Em0err,
      Pserr, Pmerr, Dserr, Dmerr, errmat;
  polyvec_t z4, z3, hint, z21, z1, h, tB, tA2, tA1, m, m_, s2, s1, s1_, s, tmp,
      z3err, z4err, herr, verr, ferr, uerr;
  poly_t c, r0_, errpoly;
  polyvec_t asub, asub_auto, bsub, bsub_auto, upsilon, errvec;
  spolymat_t R2err, R2primeerr, R2err_, R2primeerr_;
  int_ptr coeff;
  spolymat_t R2_;
  spolyvec_t r1_;
  uint8_t buf[2];
  int b;
  const unsigned int n = 2 * (tbox->m1 + quad->l);
  const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
  const unsigned int n_ = 2 * (m1 + l);

  dom = 0;

  poly_alloc (c, Rq);
  poly_alloc (r0_, Rq);
  polyvec_alloc (tmp, Rq, 2 * (m1 + l));
  polyvec_alloc (s, Rq, 2 * (m1 + l));
  polyvec_alloc (f, Rq, params->nbin);
  polyvec_alloc (f_, Rq, params->nbin);
  polyvec_alloc (u, Rq, params->nprime);
  polyvec_alloc (u_, Rq, params->nprime);
  polyvec_alloc (s1, Rq, tbox->m1);
  polyvec_alloc (s2, Rq, tbox->m2);
  polyvec_alloc (m, Rq, tbox->l + tbox->lext);
  polyvec_alloc (tA1, Rq, tbox->kmsis);
  polyvec_alloc (tA2, Rq, tbox->kmsis);
  polyvec_alloc (tB, Rq, tbox->l + tbox->lext);
  polyvec_alloc (h, Rq, quade->lambda / 2);
  polyvec_alloc (z1, Rq, tbox->m1);
  polyvec_alloc (z21, Rq, tbox->m2 - tbox->kmsis);
  polyvec_alloc (hint, Rq, tbox->kmsis);
  polyvec_alloc (z3, Rq, 256 / d);
  polyvec_alloc (z4, Rq, 256 / d);
  spolyvec_alloc (r1_, Rq, n_, n_);
  spolyvec_alloc (r1err, Rq, n, n);
  spolyvec_alloc (r1err_, Rq, n, n);
  spolyvec_alloc (r1primeerr, Rq, np, np);
  spolyvec_alloc (r1primeerr_, Rq, np, np);
  polymat_alloc (Ds, Rq, params->nprime, m1);
  polymat_alloc (Dm, Rq, params->nprime, tbox->l);
  polymat_alloc (Ps, Rq, params->nbin, m1);
  polymat_alloc (Pm, Rq, params->nbin, tbox->l);
  polymat_alloc (A1, Rq, tbox->kmsis, tbox->m1);
  polymat_alloc (A2prime, Rq, tbox->kmsis, tbox->m2 - tbox->kmsis);
  polymat_alloc (Bprime, Rq, tbox->l + tbox->lext, tbox->m2 - tbox->kmsis);
  spolymat_alloc (R2_, Rq, n_, n_, (n_ * n_ - n_) / 2 + n_);

  spolymat_alloc (R2err, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (R2primeerr, Rq, np, np, (np * np - np) / 2 + np);
  spolymat_alloc (R2err_, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (R2primeerr_, Rq, np, np, (np * np - np) / 2 + np);

  /* error poly */
  poly_alloc (errpoly, Rq);

  /* error vector large enough to hold smaller error subvecs */
  nelems = 0;
  nelems = MAX (nelems, quade->lambda / 2);
  nelems = MAX (nelems, 256 / d);
  nelems = MAX (nelems, tbox->l + tbox->lext);
  nelems = MAX (nelems, tbox->m1);
  nelems = MAX (nelems, tbox->m2 - tbox->kmsis);
  nelems = MAX (nelems, tbox->kmsis);
  nelems = MAX (nelems, params->nex);
  nelems = MAX (nelems, params->nprime);
  polyvec_alloc (errvec, Rq, nelems);

  polyvec_get_subvec (herr, errvec, 0, quade->lambda / 2, 1);
  polyvec_get_subvec (z3err, errvec, 0, 256 / d, 1);
  polyvec_get_subvec (z4err, errvec, 0, 256 / d, 1);
  polyvec_get_subvec (tBerr, errvec, 0, tbox->l + tbox->lext, 1);
  polyvec_get_subvec (z1err, errvec, 0, tbox->m1, 1);
  polyvec_get_subvec (z21err, errvec, 0, tbox->m2 - tbox->kmsis, 1);
  polyvec_get_subvec (hinterr, errvec, 0, tbox->kmsis, 1);
  polyvec_get_subvec (tA1err, errvec, 0, tbox->kmsis, 1);
  polyvec_get_subvec (verr, errvec, 0, params->n[0], 1);
  polyvec_get_subvec (uerr, errvec, 0, params->nprime, 1);
  polyvec_get_subvec (ferr, errvec, 0, params->nbin, 1);

  /* error matrix large enough to hold smaller error submats */
  nrows = 0;
  nrows = MAX (nrows, tbox->kmsis);
  nrows = MAX (nrows, tbox->l + tbox->lext);
  nrows = MAX (nrows, params->nex);
  nrows = MAX (nrows, params->nprime);
  ncols = 0;
  ncols = MAX (ncols, tbox->m1);
  ncols = MAX (ncols, tbox->l);
  ncols = MAX (ncols, tbox->m2 - tbox->kmsis);
  polymat_alloc (errmat, Rq, nrows, ncols);

  polymat_get_submat (A1err, errmat, 0, 0, tbox->kmsis, tbox->m1, 1, 1);
  polymat_get_submat (A2primeerr, errmat, 0, 0, tbox->kmsis,
                      tbox->m2 - tbox->kmsis, 1, 1);
  polymat_get_submat (Bprimeerr, errmat, 0, 0, tbox->l + tbox->lext,
                      tbox->m2 - tbox->kmsis, 1, 1);
  polymat_get_submat (Es0err, errmat, 0, 0, params->n[0], m1, 1, 1);
  polymat_get_submat (Em0err, errmat, 0, 0, params->n[0], l, 1, 1);
  polymat_get_submat (Pserr, errmat, 0, 0, params->nbin, m1, 1, 1);
  polymat_get_submat (Pmerr, errmat, 0, 0, params->nbin, l, 1, 1);
  polymat_get_submat (Dserr, errmat, 0, 0, params->nprime, m1, 1, 1);
  polymat_get_submat (Dmerr, errmat, 0, 0, params->nprime, l, 1, 1);

  for (j = 0; j < np; j++)
    {
      for (k = j; k < np; k++)
        {
          spolymat_insert_elem (R2err, j, k);
          spolymat_insert_elem (R2primeerr, j, k);
        }
      spolyvec_insert_elem (r1err, j);
      spolyvec_insert_elem (r1primeerr, j);
    }
  spolyvec_sort (r1err);
  spolyvec_sort (r1primeerr);
  spolymat_sort (R2err);
  spolymat_sort (R2primeerr);

  for (i = 0; i < N; i++)
    {
      spolymat_alloc (R2i[i], Rq, n, n, (n * n - n) / 2 + n);
      R2[i] = R2i[i];
      spolyvec_alloc (r1i[i], Rq, 2 * (tbox->m1 + quad->l),
                      2 * (tbox->m1 + quad->l));
      r1[i] = r1i[i];
      poly_alloc (r0i[i], Rq);
      r0[i] = r0i[i];
    }

  for (i = 0; i < M; i++)
    {
      spolymat_alloc (R2primei[i], Rq, np, np, (np * np - np) / 2 + np);
      R2prime[i] = R2primei[i];
      spolyvec_alloc (r1primei[i], Rq, 2 * (tbox->m1 + quad_eval->l),
                      2 * (tbox->m1 + quad_eval->l));
      r1prime[i] = r1primei[i];
      poly_alloc (r0primei[i], Rq);
      r0prime[i] = r0primei[i];
    }

  for (i = 0; i < Z; i++)
    {
      polymat_alloc (Esi[i], Rq, params->n[i], m1);
      Es[i] = Esi[i];
      polymat_alloc (Emi[i], Rq, params->n[i], tbox->l);
      Em[i] = Emi[i];
      polyvec_alloc (vi[i], Rq, params->n[i]);
      v[i] = vi[i];
      polyvec_alloc (vi_[i], Rq, params->n[i]);
    }

  int_set_i64 (lo, -3);
  int_set_i64 (hi, 3);
  polyvec_urandom_bnd (s1, lo, hi, seed, dom++);

  polyvec_get_subvec (upsilon, s1, m1, params->Z, 1);

  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s2, lo, hi, seed, dom++);
  polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

  polyvec_get_subvec (s1_, s1, 0, m1, 1);
  polyvec_get_subvec (m_, m, 0, l, 1);

  /* s = (<s1>,<m>) */

  polyvec_get_subvec (asub, s, 0, m1, 2);
  polyvec_get_subvec (asub_auto, s, 1, m1, 2);
  polyvec_set (asub, s1_);
  polyvec_auto (asub_auto, s1_);
  if (l > 0)
    {
      polyvec_get_subvec (bsub, s, m1 * 2, l, 2);
      polyvec_get_subvec (bsub_auto, s, m1 * 2 + 1, l, 2);
      polyvec_set (bsub, m_);
      polyvec_auto (bsub_auto, m_);
    }

  /* generate quadratic equations (in s) randomly */

  for (i = 0; i < N; i++)
    {
      spolymat_set_empty (R2_);
      for (k = 0; k < n_; k++)
        {
          for (j = k; j < n_; j++)
            spolymat_insert_elem (R2_, k, j);
        }
      spolyvec_set_empty (r1_);
      for (j = 0; j < n_; j++)
        {
          spolyvec_insert_elem (r1_, j);
        }

      spolymat_urandom (R2_, Rq->q, Rq->log2q, seed, dom++);
      spolyvec_urandom (r1_, Rq->q, Rq->log2q, seed, dom++);

      polyvec_dot2 (r0[i], r1_, s);
      polyvec_mulsparse (tmp, R2_, s);
      polyvec_fromcrt (tmp);
      poly_adddot (r0[i], s, tmp, 0);
      poly_neg_self (r0[i]);

      spolymat_fromcrt (R2_);
      spolyvec_fromcrt (r1_);
      poly_fromcrt (r0[i]);
      ASSERT_ERR (spolymat_is_upperdiag (R2_));

      _scatter_vec (r1[i], r1_, m1, Z);
      _scatter_smat (R2i[i], R2_, m1, Z, l);

      ASSERT_ERR (spolymat_is_upperdiag (R2i[i]));
    }

  for (i = 0; i < M; i++)
    {
      spolymat_set_empty (R2_);
      for (k = 0; k < n_; k++)
        {
          for (j = k; j < n_; j++)
            spolymat_insert_elem (R2_, k, j);
        }
      spolyvec_set_empty (r1_);
      for (j = 0; j < n_; j++)
        {
          spolyvec_insert_elem (r1_, j);
        }

      spolymat_urandom (R2_, Rq->q, Rq->log2q, seed, dom++);
      spolyvec_urandom (r1_, Rq->q, Rq->log2q, seed, dom++);

      polyvec_dot2 (r0prime[i], r1_, s);
      polyvec_mulsparse (tmp, R2_, s);
      polyvec_fromcrt (tmp);
      poly_adddot (r0prime[i], s, tmp, 0);
      poly_neg_self (r0prime[i]);

      spolymat_fromcrt (R2_);
      spolyvec_fromcrt (r1_);
      poly_fromcrt (r0prime[i]);
      ASSERT_ERR (spolymat_is_upperdiag (R2_));

      _scatter_vec (r1prime[i], r1_, m1, Z);
      _scatter_smat (R2primei[i], R2_, m1, Z, l);

      /* only constant coeff needs to be zero */
      poly_brandom (errpoly, 1, seed, dom++);
      coeff = poly_get_coeff (errpoly, 0);
      int_set_i64 (coeff, 0);

      poly_add (r0prime[i], r0prime[i], errpoly, 0);
      ASSERT_ERR (spolymat_is_upperdiag (R2primei[i]));
    }

  /* generate equations (in s1,m) for binary vector randomly */
  int_set_i64 (lo, 0);
  int_set_i64 (hi, 1);

  polymat_urandom (Ps, Rq->q, Rq->log2q, seed, 10000);
  polymat_urandom (Pm, Rq->q, Rq->log2q, seed, 10001);
  polyvec_urandom_bnd (f_, lo, hi, seed, 10002);
  polyvec_mul (f, Ps, s1_);
  polyvec_addmul (f, Pm, m_, 0);
  polyvec_neg_self (f);
  polyvec_add (f, f, f_, 0);

  /* generate equations (in s1,m) for exact 2-norm randomly */
  for (i = 0; i < Z; i++)
    {
      INT_T (l2sqr, 2 * Rq->q->nlimbs);
      INT_T (l2sqr2, 2 * Rq->q->nlimbs);

#if 0
      polymat_set_zero (Es[i]);
      polymat_set_zero (Em[i]);
      polyvec_set_zero (v[i]);
      polyvec_set_zero (vi_[i]);

      //Es[i] = NULL;
      //Em[i] = NULL;
      //v[i] = NULL;
#else
      polymat_urandom (Es[i], Rq->q, Rq->log2q, seed, 20000);
      polymat_urandom (Em[i], Rq->q, Rq->log2q, seed, 20001);
      polyvec_brandom (vi_[i], 1, seed, 20002);
      polyvec_mul (v[i], Es[i], s1_);
      polyvec_addmul (v[i], Em[i], m_, 0);
      polyvec_neg_self (v[i]);
      polyvec_add (v[i], v[i], vi_[i], 0);
#endif

      int_set (l2sqr2, params->l2Bsqr[i]);
      polyvec_l2sqr (l2sqr, vi_[i]);
      int_sub (l2sqr2, l2sqr2, l2sqr);
      int_binexp (polyvec_get_elem (upsilon, i), NULL, l2sqr2);
    }

  polymat_urandom (Ds, Rq->q, Rq->log2q, seed, 30000);
  polymat_urandom (Dm, Rq->q, Rq->log2q, seed, 30001);
  polyvec_brandom (u_, 4, seed, 30002);
  polyvec_mul (u, Ds, s1_);
  polyvec_addmul (u, Dm, m_, 0);
  polyvec_neg_self (u);
  polyvec_add (u, u, u_, 0);

  /* generate public parameters */

  abdlop_keygen (A1, A2prime, Bprime, seed, tbox);

  /* generate proof */

  memset (hashp, 0xff, 32);
  abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, tbox);

  lnp_tbox_prove (hashp, tB, h, c, z1, z21, hint, z3, z4, s1, m, s2, tA2, A1,
                  A2prime, Bprime, R2, r1, N, R2prime, r1prime, r0prime, M, Es,
                  Em, v, Ps, Pm, f, Ds, Dm, u, seed, params);

  /* expect successful verification */
  memset (hashv, 0xff, 32);
  b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                       A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                       r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u, params);
  TEST_EXPECT (b == 1);
  TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);

  for (i = 0; i < 1; i++)
    {
      /* expect verification failures */
      bytes_urandom (buf, sizeof (buf));
      memset (hashv, 0xff, 32);
      hashv[buf[0] % 32] ^= (1 << (buf[1] % 8));
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (herr, 1, seed, dom++);
      polyvec_add (herr, herr, h, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, herr, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      poly_brandom (errpoly, 1, seed, dom++);
      poly_add (errpoly, errpoly, c, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, errpoly, z1, z21, hint, z3, z4, tA1, tB,
                           A1, A2prime, Bprime, R2, r1, r0, N, R2prime,
                           r1prime, r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm,
                           u, params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z1err, 1, seed, dom++);
      polyvec_add (z1err, z1err, z1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1err, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z21err, 1, seed, dom++);
      polyvec_add (z21err, z21err, z21, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21err, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (hinterr, 1, seed, dom++);
      polyvec_add (hinterr, hinterr, hint, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hinterr, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z3err, 1, seed, dom++);
      polyvec_add (z3err, z3err, z3, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3err, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (z4err, 1, seed, dom++);
      polyvec_add (z4err, z4err, z3, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4err, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (tA1err, 1, seed, dom++);
      polyvec_add (tA1err, tA1err, tA1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1err, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      if (tbox->l > 0)
        {
          polyvec_brandom (tBerr, 1, seed, dom++);
          polyvec_add (tBerr, tBerr, tB, 0);
          memset (hashv, 0xff, 32);
          b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tBerr,
                               A1, A2prime, Bprime, R2, r1, r0, N, R2prime,
                               r1prime, r0prime, M, Es, Em, v, Ps, Pm, f, Ds,
                               Dm, u, params);
          TEST_EXPECT (b == 0);
        }

      polymat_brandom (A1err, 1, seed, dom++);
      polymat_add (A1err, A1err, A1, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1err,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polymat_brandom (A2primeerr, 1, seed, dom++);
      polymat_add (A2primeerr, A2primeerr, A2prime, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2primeerr, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polymat_brandom (Bprimeerr, 1, seed, dom++);
      polymat_add (Bprimeerr, Bprimeerr, Bprime, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprimeerr, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      spolymat_brandom (R2err, 1, seed, dom++);
      spolymat_add (R2err_, R2[0], R2err, 0);
      R2[0] = R2err_;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      R2[0] = R2i[0];
      TEST_EXPECT (b == 0);

      spolyvec_brandom (r1err, 1, seed, dom++);
      spolyvec_add (r1err_, r1err, r1[0], 0);
      r1[0] = r1err_;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      r1[0] = r1i[0];
      TEST_EXPECT (b == 0);

      poly_brandom (errpoly, 1, seed, dom++);
      poly_add (errpoly, errpoly, r0[0], 0);
      r0[0] = errpoly;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      r0[0] = r0i[0];
      TEST_EXPECT (b == 0);

      spolymat_brandom (R2primeerr, 1, seed, dom++);
      spolymat_add (R2primeerr_, R2prime[0], R2primeerr, 0);
      R2prime[0] = R2primeerr_;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      R2prime[0] = R2primei[0];
      TEST_EXPECT (b == 0);

      spolyvec_brandom (r1primeerr, 1, seed, dom++);
      spolyvec_add (r1primeerr_, r1primeerr, r1prime[0], 0);
      r1prime[0] = r1primeerr_;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      r1prime[0] = r1primei[0];
      TEST_EXPECT (b == 0);

      poly_brandom (errpoly, 1, seed, dom++);
      poly_add (errpoly, errpoly, r0prime[0], 0);
      r0prime[0] = errpoly;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      r0prime[0] = r0primei[0];
      TEST_EXPECT (b == 0);

      polymat_brandom (Es0err, 1, seed, dom++);
      polymat_add (Es0err, Es0err, Es[0], 0);
      Es[0] = Es0err;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      Es[0] = Esi[0];
      TEST_EXPECT (b == 0);

      polymat_brandom (Em0err, 1, seed, dom++);
      polymat_add (Em0err, Em0err, Es[0], 0);
      Em[0] = Em0err;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      Em[0] = Emi[0];
      TEST_EXPECT (b == 0);

      polyvec_brandom (verr, 1, seed, dom++);
      polyvec_add (verr, verr, v[0], 0);
      v[0] = verr;
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      v[0] = vi[0];
      TEST_EXPECT (b == 0);

      polymat_brandom (Pserr, 1, seed, dom++);
      polymat_add (Pserr, Pserr, Ps, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Pserr, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polymat_brandom (Pmerr, 1, seed, dom++);
      polymat_add (Pmerr, Pmerr, Pm, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pmerr, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (ferr, 1, seed, dom++);
      polyvec_add (ferr, ferr, f, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, ferr, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polymat_brandom (Dserr, 1, seed, dom++);
      polymat_add (Dserr, Dserr, Ds, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Dserr, Dm, u,
                           params);
      TEST_EXPECT (b == 0);

      polymat_brandom (Dmerr, 1, seed, dom++);
      polymat_add (Dmerr, Dmerr, Dm, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dmerr, u,
                           params);
      TEST_EXPECT (b == 0);

      polyvec_brandom (uerr, 1, seed, dom++);
      polyvec_add (uerr, uerr, u, 0);
      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, uerr,
                           params);
      TEST_EXPECT (b == 0);

      /* expect successful verification */

      memset (hashv, 0xff, 32);
      b = lnp_tbox_verify (hashv, h, c, z1, z21, hint, z3, z4, tA1, tB, A1,
                           A2prime, Bprime, R2, r1, r0, N, R2prime, r1prime,
                           r0prime, M, Es, Em, v, Ps, Pm, f, Ds, Dm, u,
                           params);
      TEST_EXPECT (b == 1);
      TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);
    }

  poly_free (c);
  poly_free (r0_);
  poly_free (errpoly);
  polyvec_free (tmp);
  polyvec_free (s);
  polyvec_free (f);
  polyvec_free (f_);
  polyvec_free (u);
  polyvec_free (u_);
  polyvec_free (s1);
  polyvec_free (s2);
  polyvec_free (m);
  polyvec_free (tA1);
  polyvec_free (tA2);
  polyvec_free (tB);
  polyvec_free (h);
  polyvec_free (z1);
  polyvec_free (z21);
  polyvec_free (hint);
  polyvec_free (z3);
  polyvec_free (z4);
  spolyvec_free (r1_);
  spolyvec_free (r1err);
  spolyvec_free (r1err_);
  spolyvec_free (r1primeerr);
  spolyvec_free (r1primeerr_);
  polyvec_free (errvec);
  polymat_free (Ds);
  polymat_free (Dm);
  polymat_free (Ps);
  polymat_free (Pm);
  polymat_free (A1);
  polymat_free (A2prime);
  polymat_free (Bprime);
  polymat_free (errmat);
  spolymat_free (R2_);
  spolymat_free (R2err);
  spolymat_free (R2primeerr);
  spolymat_free (R2err_);
  spolymat_free (R2primeerr_);
  for (i = 0; i < N; i++)
    {
      spolymat_free (R2[i]);
      spolyvec_free (r1[i]);
      poly_free (r0[i]);
    }
  for (i = 0; i < M; i++)
    {
      spolymat_free (R2prime[i]);
      spolyvec_free (r1prime[i]);
      poly_free (r0prime[i]);
    }
  for (i = 0; i < Z; i++)
    {
      polymat_free (Esi[i]);
      polymat_free (Emi[i]);
      polyvec_free (vi[i]);
      polyvec_free (vi_[i]);
    }
}
