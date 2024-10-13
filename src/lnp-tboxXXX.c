#include "lnp-tbox.h"
#include "lazer.h"
#include "stopwatch.h"

/* expand i-th row of R from cseed and i */
static inline void
_expand_R_i (int8_t *Ri, unsigned int ncols, unsigned int i,
             const uint8_t cseed[32])
{
  _brandom (Ri, ncols, 1, cseed, i);
}

/* expand i-th row of Rprime from cseed and 256 + i */
static inline void
_expand_Rprime_i (int8_t *Rprimei, unsigned int ncols, unsigned int i,
                  const uint8_t cseed[32])
{
  _brandom (Rprimei, ncols, 1, cseed, 256 + i);
}

/* swap row and col iff row > col */
static inline void
_diag (unsigned int *row, unsigned int *col, unsigned int r, unsigned int c)
{
  if (r > c)
    {
      *row = c;
      *col = r;
    }
  else
    {
      *row = r;
      *col = c;
    }
}

/* beta3^2 - 1 == 0 */
static void
_quadeq_beta3 (spolymat_ptr R2[], poly_ptr r0[],
               const lnp_tbox_params_t params)
{
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int lambda = quade->lambda;
  const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  poly_ptr poly;
  int_ptr coeff;
  spolymat_ptr R2i = R2[QUADEQ_BETA3_OFF];

  /* R2 */
  spolymat_set_empty (R2i);

  poly = spolymat_insert_elem (R2i, ibeta, ibeta);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_set (coeff, params->inv4);

  poly = spolymat_insert_elem (R2i, ibeta, ibeta + 1);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_set (coeff, Rq->inv2);

  poly = spolymat_insert_elem (R2i, ibeta + 1, ibeta + 1);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_set (coeff, params->inv4);

  if (r0 != NULL)
    {
      poly_ptr r0i = r0[QUADEQ_BETA3_OFF];

      /* r0 */
      poly_set_zero (r0i);

      coeff = poly_get_coeff (r0i, 0);
      int_set_i64 (coeff, -1);
    }

  // spolymat_sort (R2i); // not required, elems are inserted in order
  R2i->sorted = 1;
  ASSERT_ERR (spolymat_is_upperdiag (R2i));
}

/* beta4^2 - 1 == 0 */
static void
_quadeq_beta4 (spolymat_ptr R2[], poly_ptr r0[],
               const lnp_tbox_params_t params)
{
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int lambda = quade->lambda;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  poly_ptr poly;
  int_ptr coeff;
  spolymat_ptr R2i = R2[QUADEQ_BETA4_OFF];

  /* R2 */
  spolymat_set_empty (R2i);

  poly = spolymat_insert_elem (R2i, ibeta, ibeta);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_neg (coeff, params->inv4);

  poly = spolymat_insert_elem (R2i, ibeta, ibeta + 1);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_set (coeff, Rq->inv2);

  poly = spolymat_insert_elem (R2i, ibeta + 1, ibeta + 1);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_neg (coeff, params->inv4);

  if (r0 != NULL)
    {
      poly_ptr r0i = r0[QUADEQ_BETA4_OFF];

      /* r0 */
      poly_set_zero (r0i);

      coeff = poly_get_coeff (r0i, 0);
      int_set_i64 (coeff, -1);
    }

  // spolymat_sort (R2i); // not required, elems are inserted in order
  R2i->sorted = 1;
  ASSERT_ERR (spolymat_is_upperdiag (R2i));
}

/* const.coeff of X^i * beta3 == 0 for i in [1,d-1] */
static void
_evaleq_beta3 (spolyvec_ptr r1prime[], const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  unsigned int i, j;
  spolyvec_ptr r1primei;
  poly_ptr poly;
  int_ptr coeff;

  for (j = 0; j < d - 1; j++)
    {
      r1primei = r1prime[EVALEQ_BETA3_OFF + j];
      i = j + 1;

      /* r1' */
      spolyvec_set_empty (r1primei);

      poly = spolyvec_insert_elem (r1primei, ibeta);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, i);
      int_set (coeff, Rq->inv2);
      poly = spolyvec_insert_elem (r1primei, ibeta + 1);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, i);
      int_set (coeff, Rq->inv2);

      r1primei->sorted = 1;
    }
}

/* const.coeff of X^i * beta4 == 0 for i in [1,d-1] */
static void
_evaleq_beta4 (spolyvec_ptr r1prime[], const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  spolyvec_ptr r1primei;
  unsigned int i, j;
  poly_ptr poly;
  int_ptr coeff;

  for (j = 0; j < d - 1; j++)
    {
      r1primei = r1prime[EVALEQ_BETA4_OFF + j];
      i = j + 1;

      /* r1' */
      spolyvec_set_empty (r1primei);

      poly = spolyvec_insert_elem (r1primei, ibeta);
      poly_set_zero (poly);
      if (i < d / 2)
        {
          coeff = poly_get_coeff (poly, i + d / 2);
          int_neg (coeff, Rq->inv2);
        }
      else
        {
          coeff = poly_get_coeff (poly, i + d / 2 - d);
          int_set (coeff, Rq->inv2);
        }
      poly = spolyvec_insert_elem (r1primei, ibeta + 1);
      poly_set_zero (poly);
      if (i < d / 2)
        {
          coeff = poly_get_coeff (poly, i + d / 2);
          int_set (coeff, Rq->inv2);
        }
      else
        {
          coeff = poly_get_coeff (poly, i + d / 2 - d);
          int_neg (coeff, Rq->inv2);
        }

      r1primei->sorted = 1;
    }
}

/* z3 = y3 + beta3*R*s3 is well-formed. */
static void
_evaleq_z3 (spolymat_ptr R2prime[], spolyvec_ptr r1prime[], poly_ptr r0prime[],
            polymat_ptr Es[], polymat_ptr Em[], polyvec_ptr v[], polymat_t Ps,
            polymat_t Pm, polyvec_t f, const intvec_t z3coeffs,
            const uint8_t cseed[32], const lnp_tbox_params_t params,
            polyvec_t tmp_polyvec, UNUSED unsigned int j_)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  const unsigned int nbin = params->nbin;
  const unsigned int *ni = params->n;
  unsigned int coloff, rowoff, quot, rem, i, j, k, row, col;
  INT_T (a_, int_get_nlimbs (q));
  INT_T (a, 2 * int_get_nlimbs (q));
  INT_T (at, 2 * int_get_nlimbs (q));
  INTVEC_T (row_, nex * d, int_get_nlimbs (q));
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
  polyvec_t ri, subv, subv2;
  poly_ptr poly, poly2;
  int_ptr coeff;
  intvec_t subrow;
  int8_t Ri[nex * d];

  ASSERT_ERR (j_ < 256);

  polyvec_get_subvec (ri, tmp_polyvec, 0, nex, 1);

  for (j = 0; j < 256; j++)
    {
      R2primei = R2prime[EVALEQ_Z3_OFF + j];
      r1primei = r1prime[EVALEQ_Z3_OFF + j];
      r0primei = r0prime[EVALEQ_Z3_OFF + j];

      /* r0': neg. j-th coeff of z3 */
      poly_set_zero (r0primei);

      coeff = poly_get_coeff (r0primei, 0);
      int_set (coeff, intvec_get_elem (z3coeffs, j));
      int_neg_self (coeff);

      /*
       * r1': a*beta3 + o(e_j)*y3,
       *      a = <j-th row of R,(f,v,0)>
       */
      spolyvec_set_empty (r1primei);

      _expand_R_i (Ri, nex * d, j, cseed);
      intvec_set_i8 (row_, Ri);

      int_set_zero (a);
      if (f != NULL)
        {
          coloff = 0;

          for (i = 0; i < nbin; i++)
            {
              poly = polyvec_get_elem (f, i);
              intvec_get_subvec (subrow, row_, coloff, d, 1);
              coloff += d;
              intvec_dot (at, subrow, poly_get_coeffvec (poly));
              int_add (a, a, at);
            }
        }
      if (v != NULL)
        {
          coloff = nbin * d;

          for (i = 0; i < Z; i++)
            {
              if (v[i] != NULL)
                {
                  for (k = 0; k < ni[i]; k++)
                    {
                      poly = polyvec_get_elem (v[i], k);
                      intvec_get_subvec (subrow, row_, coloff, d, 1);
                      coloff += d;
                      intvec_dot (at, subrow, poly_get_coeffvec (poly));
                      int_add (a, a, at);
                    }
                }
            }
        }
      if (f != NULL || v != NULL)
        {
          int_mod (a_, a, q);
          int_mul (a, Rq->inv2, a_);
          int_mod (a_, a, q);
          int_redc (a_, a_, q);
        }

      poly = spolyvec_insert_elem (r1primei, ibeta);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, 0);
      int_set (coeff, a_);
      poly = spolyvec_insert_elem (r1primei, ibeta + 1);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, 0);
      int_set (coeff, a_);

      quot = j / d;
      rem = j - d * quot;

      poly = spolyvec_insert_elem (r1primei, 2 * (m1 + Z + l + quot));
      poly_set_zero (poly);
      if (rem == 0)
        {
          coeff = poly_get_coeff (poly, 0);
          int_set_i64 (coeff, 1);
        }
      else
        {
          coeff = poly_get_coeff (poly, d - rem);
          int_set_i64 (coeff, -1);
        }

      /* R2' */
      spolymat_set_empty (R2primei);

      polyvec_set_coeffvec (ri, row_);
      polyvec_auto_self (ri);
      polyvec_scale (ri, Rq->inv2, ri);

      for (k = 0; k < 2 * (m1 + Z + l); k += 2)
        {
          /* diagonalize */
          if (ibeta <= k)
            {
              row = ibeta;
              col = k;
            }
          else
            {
              row = k;
              col = ibeta;
            }
          poly = spolymat_insert_elem (R2primei, row, col);
          poly_set_zero (poly);

          if (ibeta + 1 <= k)
            {
              row = ibeta + 1;
              col = k;
            }
          else
            {
              row = k;
              col = ibeta + 1;
            }
          poly2 = spolymat_insert_elem (R2primei, row, col);

          if (k < 2 * m1)
            {
              rowoff = 0;

              if (Ps != NULL)
                {
                  polyvec_get_subvec (subv2, ri, rowoff, nbin, 1);

                  polymat_get_col (subv, Ps, k / 2);
                  poly_adddot (poly, subv2, subv, 0);
                }
              rowoff += nbin;

              if (Es != NULL)
                {
                  for (i = 0; i < Z; i++)
                    {
                      if (Es[i] != NULL)
                        {
                          polyvec_get_subvec (subv2, ri, rowoff, ni[i], 1);
                          polymat_get_col (subv, Es[i], k / 2);
                          poly_adddot (poly, subv2, subv, 0);
                        }
                      rowoff += ni[i];
                    }
                }
            }
          else if (k < 2 * (m1 + Z))
            {
              rowoff = (nex - Z) + (k - 2 * m1) / 2;
              poly_set (poly, polyvec_get_elem (ri, rowoff));
            }
          else
            {
              rowoff = 0;

              if (Pm != NULL)
                {
                  polyvec_get_subvec (subv2, ri, rowoff, nbin, 1);
                  polymat_get_col (subv, Pm, (k - 2 * (m1 + Z)) / 2);
                  poly_adddot (poly, subv2, subv, 0);
                }
              rowoff += nbin;

              if (Em != NULL)
                {
                  for (i = 0; i < Z; i++)
                    {
                      if (Em[i] != NULL)
                        {
                          polyvec_get_subvec (subv2, ri, rowoff, ni[i], 1);
                          polymat_get_col (subv, Em[i],
                                           (k - 2 * (m1 + Z)) / 2);
                          poly_adddot (poly, subv2, subv, 0);
                        }
                      rowoff += ni[i];
                    }
                }
            }

          poly_set (poly2, poly);
        }
      spolyvec_sort (r1primei); // XXX already compute sorted
      spolymat_fromcrt (R2primei);
      spolymat_sort (R2primei); // XXX already compute sorted
      ASSERT_ERR (spolymat_is_upperdiag (R2primei));
    }
}

/* z4 = y4 + beta4*Rprime*s4 is well-formed. */
static void
_evaleq_z4 (spolymat_ptr R2prime[], spolyvec_ptr r1prime[], poly_ptr r0prime[],
            polymat_t Ds, polymat_t Dm, polyvec_t u, const intvec_t z4coeffs,
            const uint8_t cseed[32], const lnp_tbox_params_t params,
            polyvec_t tmp_polyvec, UNUSED unsigned int j_)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  unsigned int coloff, rowoff, quot, rem, i, j, k, row, col;
  INT_T (a_, int_get_nlimbs (q));
  INT_T (a, 2 * int_get_nlimbs (q));
  INT_T (at, 2 * int_get_nlimbs (q));
  INTVEC_T (row_, nprime * d, int_get_nlimbs (q));
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
  polyvec_t ri_, subv, subv2;
  poly_ptr poly, poly2;
  int_ptr coeff;
  intvec_t subrow;
  int8_t Rprimei[nprime * d];

  ASSERT_ERR (j_ < 256);

  polyvec_get_subvec (ri_, tmp_polyvec, 0, nprime, 1);

  for (j = 0; j < 256; j++)
    {
      R2primei = R2prime[EVALEQ_Z4_OFF + j];
      r1primei = r1prime[EVALEQ_Z4_OFF + j];
      r0primei = r0prime[EVALEQ_Z4_OFF + j];

      /* r0': neg. j-th coeff of z3 */
      poly_set_zero (r0primei);

      coeff = poly_get_coeff (r0primei, 0);
      int_set (coeff, intvec_get_elem (z4coeffs, j));
      int_neg_self (coeff);

      /*
       * r1': a*beta4 + o(e_j)*y4,
       *      a = <j-th row of R',u>
       */
      spolyvec_set_empty (r1primei);

      _expand_Rprime_i (Rprimei, nprime * d, j, cseed);
      intvec_set_i8 (row_, Rprimei);

      int_set_zero (a);

      if (u != NULL)
        {
          coloff = 0;
          for (i = 0; i < nprime; i++)
            {
              poly = polyvec_get_elem (u, i);
              intvec_get_subvec (subrow, row_, coloff, d, 1);
              coloff += d;
              intvec_dot (at, subrow, poly_get_coeffvec (poly));
              int_add (a, a, at);
            }
          int_mod (a_, a, q);
          int_mul (a, Rq->inv2, a_);
          int_mod (a_, a, q);
          int_redc (a_, a_, q);
        }

      poly = spolyvec_insert_elem (r1primei, ibeta);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, d / 2);
      int_set (coeff, a_);
      int_neg_self (coeff);
      poly = spolyvec_insert_elem (r1primei, ibeta + 1);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, d / 2);
      int_set (coeff, a_);

      quot = j / d;
      rem = j - d * quot;

      poly
          = spolyvec_insert_elem (r1primei, 2 * (m1 + Z + l + 256 / d + quot));
      poly_set_zero (poly);
      if (rem == 0)
        {
          coeff = poly_get_coeff (poly, 0);
          int_set_i64 (coeff, 1);
        }
      else
        {
          coeff = poly_get_coeff (poly, d - rem);
          int_set_i64 (coeff, -1);
        }

      /* R2' */
      spolymat_set_empty (R2primei);

      polyvec_set_coeffvec (ri_, row_);
      polyvec_auto_self (ri_);
      polyvec_scale (ri_, Rq->inv2, ri_);
      polyvec_neg_self (ri_);
      polyvec_lrot (ri_, ri_, d / 2);

      for (k = 0; k < 2 * (m1 + Z + l); k += 2)
        {
          /* diagonalize */
          if (ibeta <= k)
            {
              row = ibeta;
              col = k;
            }
          else
            {
              row = k;
              col = ibeta;
            }
          poly = spolymat_insert_elem (R2primei, row, col);
          poly_set_zero (poly);

          if (ibeta + 1 <= k)
            {
              row = ibeta + 1;
              col = k;
            }
          else
            {
              row = k;
              col = ibeta + 1;
            }
          poly2 = spolymat_insert_elem (R2primei, row, col);
          poly_set_zero (poly2); // XXX

          if (k < 2 * m1)
            {
              if (Ds != NULL)
                {
                  rowoff = 0;
                  polyvec_get_subvec (subv2, ri_, rowoff, nprime, 1);
                  rowoff += nprime;
                  polymat_get_col (subv, Ds, k / 2);
                  poly_adddot (poly, subv2, subv, 0);
                }
            }
          else if (k < 2 * (m1 + Z))
            {
              continue;
            }
          else
            {
              if (Dm != NULL)
                {
                  rowoff = 0;
                  polyvec_get_subvec (subv2, ri_, rowoff, nprime, 1);
                  rowoff += nprime;
                  polymat_get_col (subv, Dm, (k - 2 * (m1 + Z)) / 2);
                  poly_adddot (poly, subv2, subv, 0);
                }
            }

          poly_set (poly2, poly);
          poly_neg_self (poly2);
        }
      spolyvec_sort (r1primei); // XXX already compute sorted
      spolymat_fromcrt (R2primei);
      spolymat_sort (R2primei); // XXX already compute sorted
      ASSERT_ERR (spolymat_is_upperdiag (R2primei));
    }
}

/* <upsilon_j,upsilon_j - 1> = 0 */
static void
_evaleq_upsilon (spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                 const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  unsigned int i, j, iupsilonj;
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr poly;
  int_ptr coeff;

  for (j = 0; j < Z; j++)
    {
      R2primei = R2prime[EVALEQ_UPSILON_OFF + j];
      r1primei = r1prime[EVALEQ_UPSILON_OFF + j];
      iupsilonj = (m1 + j) * 2;

      /* R2' */
      spolymat_set_empty (R2primei);

      poly = spolymat_insert_elem (R2primei, iupsilonj, iupsilonj + 1);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, 0);
      int_set_one (coeff);

      /* r1' */
      spolyvec_set_empty (r1primei);

      poly = spolyvec_insert_elem (r1primei, iupsilonj + 1);
      for (i = 0; i < d; i++)
        {
          coeff = poly_get_coeff (poly, i);
          int_set_one (coeff);
          int_neg_self (coeff);
        }

      // spolymat_sort (R2primei); // already sorted
      r1primei->sorted = 1;
      R2primei->sorted = 1;
      ASSERT_ERR (spolymat_is_upperdiag (R2primei));
    }
}

/* <Ps*s1+Pm*m+f,Ps*s1+Pm*m+f-1> = 0 */
static void
_evaleq_bin (spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
             poly_ptr r0prime[], polymat_t Ps, polymat_t Pm, polyvec_t f,
             const lnp_tbox_params_t params, polyvec_t tmp_polyvec,
             polyvec_t tmp_polyvec2)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nbin = params->nbin;
  INT_T (a_, int_get_nlimbs (q));
  INT_T (a, 2 * int_get_nlimbs (q));
  INT_T (at, 2 * int_get_nlimbs (q));
  polyvec_t subv;
  polyvec_ptr f_;
  intvec_ptr coeffs, coeffs2;
  unsigned int i, j, k, row, col;
  int_ptr coeff;
  poly_ptr poly, poly2;
  spolymat_ptr R2primei = R2prime[EVALEQ_BIN_OFF];
  spolyvec_ptr r1primei = r1prime[EVALEQ_BIN_OFF];
  poly_ptr r0primei = r0prime[EVALEQ_BIN_OFF];
  polyvec_t tmp_polyvec3;

  polyvec_alloc (tmp_polyvec3, Rq, MAX (m1, l));

  f_ = tmp_polyvec;

  /* r0' */
  poly_set_zero (r0primei);

  for (j = 0; j < nbin; j++)
    {
      poly = polyvec_get_elem (f_, j);
      for (i = 0; i < d; i++)
        {
          coeff = poly_get_coeff (poly, i);
          int_set_one (coeff);
          int_neg_self (coeff);
        }
    }

  if (f != NULL)
    {
      polyvec_add (f_, f_, f, 0);
      polyvec_redc (f_, f_); /* f - 1 XXX correct */

      int_set_zero (a);
      for (j = 0; j < nbin; j++)
        {
          poly = polyvec_get_elem (f_, j);
          poly2 = polyvec_get_elem (f, j);
          coeffs = poly_get_coeffvec (poly);
          coeffs2 = poly_get_coeffvec (poly2);

          intvec_dot (at, coeffs, coeffs2);
          int_add (a, a, at);
        }
      int_mod (a_, a, q);
      int_redc (a_, a_, q);
      coeff = poly_get_coeff (r0primei, 0);
      int_set (coeff, a_); /* <f,f-1> XXX correct */
    }

  /* R2', r1' */
  spolymat_set_empty (R2primei);
  spolyvec_set_empty (r1primei);

  if (Ps != NULL)
    {
      for (k = 0; k < m1; k++)
        {
          polymat_get_col (subv, Ps, k);
          polyvec_auto (tmp_polyvec2, subv);

          /* o(s1)^T*o(Ps)^T*Ps*s1 */
          for (j = 0; j < m1; j++)
            {
              polymat_get_col (subv, Ps, j);
              _diag (&row, &col, 2 * k + 1, 2 * j);
              poly = spolymat_insert_elem (R2primei, row, col);
              polyvec_dot (poly, tmp_polyvec2, subv);
            }

          /* o(s1)^T*o(Ps)^T*Pm*m */
          if (Pm != NULL)
            {
              for (j = 0; j < l; j++)
                {
                  polymat_get_col (subv, Pm, j);
                  _diag (&row, &col, 2 * k + 1, 2 * (m1 + Z + j));
                  poly = spolymat_insert_elem (R2primei, row, col);
                  polyvec_dot (poly, tmp_polyvec2, subv);
                }
            }

          /* o(s1)^T*o(Ps)^T*(f-1) */
          poly = spolyvec_insert_elem (r1primei, 2 * k + 1);
          polyvec_dot (poly, tmp_polyvec2, f_);
        }
    }
  if (Pm != NULL)
    {
      for (k = 0; k < l; k++)
        {
          polymat_get_col (subv, Pm, k);
          polyvec_auto (tmp_polyvec2, subv);

          /* o(m)^T*o(Pm)^T*Ps*s1 */
          if (Ps != NULL)
            {
              for (j = 0; j < m1; j++)
                {
                  polymat_get_col (subv, Ps, j);
                  _diag (&row, &col, 2 * (m1 + Z + k) + 1, 2 * j);
                  poly = spolymat_insert_elem (R2primei, row, col);
                  polyvec_dot (poly, tmp_polyvec2, subv);
                }
            }

          /* o(m)^T*o(Pm)^T*Pm*m */
          for (j = 0; j < l; j++)
            {
              polymat_get_col (subv, Pm, j);
              _diag (&row, &col, 2 * (m1 + Z + k) + 1, 2 * (m1 + Z + j));
              poly = spolymat_insert_elem (R2primei, row, col);
              polyvec_dot (poly, tmp_polyvec2, subv);
            }

          /* o(m)^T*o(Pm)^T*(f-1) */
          poly = spolyvec_insert_elem (r1primei, 2 * (m1 + Z + k) + 1);
          polyvec_dot (poly, tmp_polyvec2, f_);
        }
    }

  if (f != NULL)
    {
      polyvec_auto (f_, f);

      /* o(f)^T*Ps*s1 */
      if (Ps != NULL)
        {
          polyvec_get_subvec (subv, tmp_polyvec3, 0, m1, 1);
          polyvec_mul2 (subv, f_, Ps);

          _VEC_FOREACH_ELEM (subv, i)
          {
            poly = polyvec_get_elem (tmp_polyvec3, i);
            poly2 = spolyvec_insert_elem (r1primei, 0 + 2 * i);
            poly_set (poly2, poly);
          }
        }

      /* o(f)^T*Pm*m */
      if (Pm != NULL)
        {
          polyvec_get_subvec (subv, tmp_polyvec3, 0, l, 1);
          polyvec_mul2 (subv, f_, Pm);

          _VEC_FOREACH_ELEM (subv, i)
          {
            poly = polyvec_get_elem (tmp_polyvec3, i);
            poly2 = spolyvec_insert_elem (r1primei, 2 * (m1 + Z) + 2 * i);
            poly_set (poly2, poly);
          }
        }
    }

  spolyvec_fromcrt (r1primei);
  spolyvec_sort (r1primei); // compute already sorted ?
  spolymat_fromcrt (R2primei);
  spolymat_sort (R2primei); // compute already sorted ?
  ASSERT_ERR (spolymat_is_upperdiag (R2primei));

  polyvec_free (tmp_polyvec3);
}

/* <bin(Bi^2),upsiloni> + <Es*s1+Em*m+v, Es*s1+Em*m+v> - Bi^2= 0 */
static void
_evaleq_l2 (spolymat_ptr R2prime[], spolyvec_ptr r1prime[], poly_ptr r0prime[],
            polymat_ptr Es[], polymat_ptr Em[], polyvec_ptr v[],
            const lnp_tbox_params_t params, polyvec_t tmp_polyvec)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int *ni = params->n;
  INT_T (a_, int_get_nlimbs (q));
  INT_T (a, 2 * int_get_nlimbs (q));
  INT_T (at, 2 * int_get_nlimbs (q));
  unsigned int i, j, k, row, col;
  polyvec_t subv2, subv, tmp_polyvec3;
  intvec_ptr coeffs;
  int_ptr coeff;
  poly_ptr poly, poly2;
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;

  polyvec_alloc (tmp_polyvec3, Rq, MAX (m1, l));

  for (i = 0; i < Z; i++)
    {
      R2primei = R2prime[EVALEQ_L2_OFF + i];
      r1primei = r1prime[EVALEQ_L2_OFF + i];
      r0primei = r0prime[EVALEQ_L2_OFF + i];

      polyvec_get_subvec (subv2, tmp_polyvec, 0, ni[i], 1);

      /* r0' */
      poly_set_zero (r0primei);

      if (v != NULL)
        {
          int_set_zero (a);
          for (j = 0; j < ni[i]; j++)
            {
              if (v[i] == NULL)
                continue;

              poly = polyvec_get_elem (v[i], j);
              coeffs = poly_get_coeffvec (poly);

              intvec_dot (at, coeffs, coeffs);
              int_add (a, a, at);
            }

          int_mod (a_, a, q);
          int_redc (a_, a_, q);
        }
      else
        {
          int_set_zero (a_);
        }

      int_sub (a_, a_, params->l2Bsqr[i]);
      int_redc (a_, a_, q);
      coeff = poly_get_coeff (r0primei, 0);
      int_set (coeff, a_); /* <v,v> - Bi^2 */

      /* R2', r1' */
      spolymat_set_empty (R2primei);
      spolyvec_set_empty (r1primei);

      if (Es != NULL && Es[i] != NULL)
        {
          for (k = 0; k < m1; k++)
            {
              polymat_get_col (subv, Es[i], k);
              polyvec_auto (subv2, subv);

              /* o(s1)^T*o(Es)^T*Es*s1 */
              for (j = 0; j < m1; j++)
                {
                  polymat_get_col (subv, Es[i], j);
                  _diag (&row, &col, 2 * k + 1, 2 * j);
                  poly = spolymat_insert_elem (R2primei, row, col);
                  polyvec_dot (poly, subv2, subv);
                }

              /* o(s1)^T*o(Es)^T*Em*m */
              if (Em != NULL && Em[i] != NULL)
                {
                  for (j = 0; j < l; j++)
                    {
                      polymat_get_col (subv, Em[i], j);
                      _diag (&row, &col, 2 * k + 1, 2 * (m1 + Z + j));
                      poly = spolymat_insert_elem (R2primei, row, col);
                      polyvec_dot (poly, subv2, subv);
                    }
                }

              /* o(s1)^T*o(Es)^T*v */
              if (v != NULL && v[i] != NULL)
                {
                  poly = spolyvec_insert_elem (r1primei, 2 * k + 1);
                  polyvec_dot (poly, subv2, v[i]);
                }
            }
        }
      if (Em != NULL && Em[i] != NULL)
        {
          for (k = 0; k < l; k++)
            {
              polymat_get_col (subv, Em[i], k);
              polyvec_auto (subv2, subv);

              if (Es != NULL && Es[i] != NULL)
                {
                  /* o(m)^T*o(Em)^T*Es*s1 */
                  for (j = 0; j < m1; j++)
                    {
                      polymat_get_col (subv, Es[i], j);
                      _diag (&row, &col, 2 * (m1 + Z + k) + 1, 2 * j);
                      poly = spolymat_insert_elem (R2primei, row, col);
                      polyvec_dot (poly, subv2, subv);
                    }
                }

              /* o(m)^T*o(Em)^T*Em*m */
              for (j = 0; j < l; j++)
                {
                  polymat_get_col (subv, Em[i], j);
                  _diag (&row, &col, 2 * (m1 + Z + k) + 1, 2 * (m1 + Z + j));
                  poly = spolymat_insert_elem (R2primei, row, col);
                  polyvec_dot (poly, subv2, subv);
                }

              /* o(m)^T*o(Em)^T*v */
              if (v != NULL && v[i] != NULL)
                {
                  poly = spolyvec_insert_elem (r1primei, 2 * (m1 + Z + k) + 1);
                  polyvec_dot (poly, subv2, v[i]);
                }
            }
        }

      if (v != NULL && v[i] != NULL)
        {
          polyvec_auto (subv2, v[i]);

          if (Es != NULL && Es[i] != NULL)
            {
              /* o(v)^T*Es*s1 */
              polyvec_get_subvec (subv, tmp_polyvec3, 0, m1, 1);
              polyvec_mul2 (subv, subv2, Es[i]);

              _VEC_FOREACH_ELEM (subv, j)
              {
                poly = polyvec_get_elem (subv, j);
                poly2 = spolyvec_insert_elem (r1primei, 0 + 2 * j);
                poly_set (poly2, poly);
              }
            }

          if (Em != NULL && Em[i] != NULL)
            {
              /* o(v)^T*Em*m */
              polyvec_get_subvec (subv, tmp_polyvec3, 0, l, 1);
              polyvec_mul2 (subv, subv2, Em[i]);

              _VEC_FOREACH_ELEM (subv, j)
              {
                poly = polyvec_get_elem (subv, j);
                poly2 = spolyvec_insert_elem (r1primei, 2 * (m1 + Z) + 2 * j);
                poly_set (poly2, poly);
              }
            }
        }

      /* o(bin(Bi^2))*upsiloni */
      poly = spolyvec_insert_elem (r1primei, 2 * (m1 + i));
      int_binexp (NULL, poly, params->l2Bsqr[i]);
      poly_auto_self (poly);

      spolyvec_fromcrt (r1primei);
      spolyvec_sort (r1primei); // compute already sorted
      spolymat_fromcrt (R2primei);
      spolymat_sort (R2primei); // compute already sorted
      ASSERT_ERR (spolymat_is_upperdiag (R2primei));
    }
  polyvec_free (tmp_polyvec3);
}

/*
 * Allocate and init eqs where it can be done w/o knowing the
 * statement.
 */
void
lnp_tbox_eqs_alloc (spolymat_ptr R2[], spolyvec_ptr r1[], poly_ptr r0[],
                    spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                    poly_ptr r0prime[], const lnp_tbox_params_t params)
{
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr quad_eval = quade->quad_eval;
  abdlop_params_srcptr quad = quade->quad_many;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int lambda = quade->lambda;
  /* dimension of quad eqs */
  const unsigned int n = 2 * (tbox->m1 + quad->l);
  /* dimension of eval eqs */
  const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  spolymat_ptr R2i, R2primei;
  spolyvec_ptr r1i, r1primei;
  poly_ptr r0i, r0primei;
  unsigned int i;

  /* alloc quad eqs */

  for (i = 0; i < QUADEQ_BASE_OFF; i++)
    {
      R2i = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2i, Rq, n, n, NELEMS_DIAG (n));
      R2[i] = R2i;

      r1i = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1i, Rq, n, n);
      r1[i] = r1i;

      r0i = _alloc (sizeof (poly_t));
      poly_alloc (r0i, Rq);
      r0[i] = r0i;
    }
  for (i = QUADEQ_BETA3_OFF; i < QUADEQ_BETA3_OFF + QUADEQ_BETA3_N; i++)
    {
      R2i = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2i, Rq, n, n, QUADEQ_BETA3_MAXQ);
      R2[i] = R2i;

      r1[i] = NULL;

      if (r0 != NULL)
        {
          r0i = _alloc (sizeof (poly_t));
          poly_alloc (r0i, Rq);
          r0[i] = r0i;
        }
    }
  for (i = QUADEQ_BETA4_OFF; i < QUADEQ_BETA4_OFF + QUADEQ_BETA4_N; i++)
    {
      R2i = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2i, Rq, n, n, QUADEQ_BETA4_MAXQ);
      R2[i] = R2i;

      r1[i] = NULL;

      if (r0 != NULL)
        {
          r0i = _alloc (sizeof (poly_t));
          poly_alloc (r0i, Rq);
          r0[i] = r0i;
        }
    }

  /* alloc eval eqs*/

  for (i = EVALEQ_BETA3_OFF; i < EVALEQ_BETA3_OFF + EVALEQ_BETA3_N; i++)
    {
      R2prime[i] = NULL;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0prime[i] = NULL;
    }
  for (i = EVALEQ_BETA4_OFF; i < EVALEQ_BETA4_OFF + EVALEQ_BETA4_N; i++)
    {
      R2prime[i] = NULL;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0prime[i] = NULL;
    }
#if 1
  for (i = EVALEQ_Z3_OFF; i < EVALEQ_Z3_OFF + EVALEQ_Z3_N; i++)
    {
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, EVALEQ_Z3_MAXQ);
      R2prime[i] = R2primei;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0primei = _alloc (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime[i] = r0primei;
    }
  for (i = EVALEQ_Z4_OFF; i < EVALEQ_Z4_OFF + EVALEQ_Z4_N; i++)
    {
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, EVALEQ_Z4_MAXQ);
      R2prime[i] = R2primei;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0primei = _alloc (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime[i] = r0primei;
    }
#endif
  for (i = EVALEQ_UPSILON_OFF; i < EVALEQ_UPSILON_OFF + EVALEQ_UPSILON_N; i++)
    {
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, EVALEQ_UPSILON_MAXQ);
      R2prime[i] = R2primei;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0prime[i] = NULL;
    }

  /* init quad eqs */
  _quadeq_beta3 (R2, r0, params);
  _quadeq_beta4 (R2, r0, params);

  /* init eval eqs */
  _evaleq_beta3 (r1prime, params);
  _evaleq_beta4 (r1prime, params);
  _evaleq_upsilon (R2prime, r1prime, params);
}

/*
 * Allocate and initialize eqs where it can be done knowing the
 * binary proof statement.
 */
void
lnp_tbox_eqsbin_init (spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                      poly_ptr r0prime[], polymat_t Ps, polymat_t Pm,
                      polyvec_t f, const lnp_tbox_params_t params)
{
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr quad_eval = quade->quad_eval;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int nbin = params->nbin;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  /* dimension of eval eqs */
  const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
  polyvec_t tmp_polyvec, tmp_polyvec2;
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
  unsigned int i, nelems;

  for (i = EVALEQ_BIN_OFF; i < EVALEQ_BIN_OFF + EVALEQ_BIN_N; i++)
    {
      nelems = 0;

      if (Ps != NULL)
        {
          nelems += m1 * m1;
          if (Pm != NULL)
            {
              nelems += m1 * l;
            }
        }
      if (Pm != NULL)
        {
          nelems += l * l;
          if (Ps != NULL)
            {
              nelems += l * m1;
            }
        }
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, nelems);
      R2prime[i] = R2primei;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0primei = _alloc (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime[i] = r0primei;
    }

  /* scratch space for subroutines */
  polyvec_alloc (tmp_polyvec, Rq, nbin);
  polyvec_alloc (tmp_polyvec2, Rq, nbin);

  _evaleq_bin (R2prime, r1prime, r0prime, Ps, Pm, f, params, tmp_polyvec,
               tmp_polyvec2);

  polyvec_free (tmp_polyvec2);
  polyvec_free (tmp_polyvec);
}

void
lnp_tbox_eqsbin_clear (spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                       poly_ptr r0prime[], const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
  unsigned int i;

  for (i = EVALEQ_BIN_OFF; i < EVALEQ_BIN_OFF + EVALEQ_BIN_N; i++)
    {
      R2primei = R2prime[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0primei = r0prime[i];
      poly_free (r0primei);
      _free (r0primei, sizeof (poly_t));
      r0prime[i] = NULL;
    }
}

/*
 * Allocate and initialize eqs where it can be done knowing the
 * l2-norm proof statement.
 */
void
lnp_tbox_eqsl2_init (spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                     poly_ptr r0prime[], polymat_ptr Es[], polymat_ptr Em[],
                     polyvec_ptr v[], const lnp_tbox_params_t params)
{
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr quad_eval = quade->quad_eval;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int *n = params->n;
  /* dimension of eval eqs */
  const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
  polyvec_t tmp_polyvec;
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
  unsigned int i, nelems;

  for (i = EVALEQ_L2_OFF; i < EVALEQ_L2_OFF + EVALEQ_L2_N; i++)
    {
      nelems = 0;

      if (Es != NULL && Es[i - EVALEQ_L2_OFF] != NULL)
        {
          nelems += m1 * m1;
          if (Em != NULL && Em[i - EVALEQ_L2_OFF] != NULL)
            {
              nelems += m1 * l;
            }
        }
      if (Em != NULL && Em[i - EVALEQ_L2_OFF] != NULL)
        {
          nelems += l * l;
          if (Es != NULL && Es[i - EVALEQ_L2_OFF] != NULL)
            {
              nelems += l * m1;
            }
        }
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, nelems);
      R2prime[i] = R2primei;

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime[i] = r1primei;

      r0primei = _alloc (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime[i] = r0primei;
    }

  nelems = 0;
  for (i = 0; i < Z; i++)
    nelems = MAX (nelems, n[i]);

  /* scratch space for subroutines */
  polyvec_alloc (tmp_polyvec, Rq, nelems);

  _evaleq_l2 (R2prime, r1prime, r0prime, Es, Em, v, params, tmp_polyvec);

  polyvec_free (tmp_polyvec);
}

void
lnp_tbox_eqsl2_clear (spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                      poly_ptr r0prime[], const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  /* dimension of eval eqs */
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
  unsigned int i;

  for (i = EVALEQ_L2_OFF; i < EVALEQ_L2_OFF + EVALEQ_L2_N; i++)
    {
      R2primei = R2prime[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0primei = r0prime[i];
      poly_free (r0primei);
      _free (r0primei, sizeof (poly_t));
      r0prime[i] = NULL;
    }
}

void
lnp_tbox_eqs_free (spolymat_ptr R2[], spolyvec_ptr r1[], poly_ptr r0[],
                   spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                   poly_ptr r0prime[], const lnp_tbox_params_t params)
{
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int lambda = quade->lambda;
  spolymat_ptr R2i, R2primei;
  spolyvec_ptr r1i, r1primei;
  poly_ptr r0i, r0primei;
  unsigned int i;

  /* alloc quad eqs */

  for (i = 0; i < QUADEQ_BASE_OFF; i++)
    {
      R2i = R2[i];
      spolymat_free (R2i);
      _free (R2i, sizeof (spolymat_t));
      R2[i] = NULL;

      r1i = r1[i];
      spolyvec_free (r1i);
      _free (r1i, sizeof (spolyvec_t));
      r1[i] = NULL;

      r0i = r0[i];
      poly_free (r0i);
      _free (r0i, sizeof (poly_t));
      r0[i] = NULL;
    }
  for (i = QUADEQ_BETA3_OFF; i < QUADEQ_BETA3_OFF + QUADEQ_BETA3_N; i++)
    {
      R2i = R2[i];
      spolymat_free (R2i);
      _free (R2i, sizeof (spolymat_t));
      R2[i] = NULL;

      r1[i] = NULL;

      if (r0 != NULL)
        {
          r0i = r0[i];
          poly_free (r0i);
          _free (r0i, sizeof (poly_t));
          r0[i] = NULL;
        }
    }
  for (i = QUADEQ_BETA4_OFF; i < QUADEQ_BETA4_OFF + QUADEQ_BETA4_N; i++)
    {
      R2i = R2[i];
      spolymat_free (R2i);
      _free (R2i, sizeof (spolymat_t));
      R2[i] = NULL;

      r1[i] = NULL;

      if (r0 != NULL)
        {
          r0i = r0[i];
          poly_free (r0i);
          _free (r0i, sizeof (poly_t));
          r0[i] = NULL;
        }
    }

  /* alloc eval eqs*/

  for (i = EVALEQ_BETA3_OFF; i < EVALEQ_BETA3_OFF + EVALEQ_BETA3_N; i++)
    {
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0prime[i] = NULL;
    }
  for (i = EVALEQ_BETA4_OFF; i < EVALEQ_BETA4_OFF + EVALEQ_BETA4_N; i++)
    {
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0prime[i] = NULL;
    }
  for (i = EVALEQ_Z3_OFF; i < EVALEQ_Z3_OFF + EVALEQ_Z3_N; i++)
    {
      R2primei = R2prime[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0primei = r0prime[i];
      poly_free (r0primei);
      _free (r0primei, sizeof (poly_t));
      r0prime[i] = NULL;
    }
  for (i = EVALEQ_Z4_OFF; i < EVALEQ_Z4_OFF + EVALEQ_Z4_N; i++)
    {
      R2primei = R2prime[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0primei = r0prime[i];
      poly_free (r0primei);
      _free (r0primei, sizeof (poly_t));
      r0prime[i] = NULL;
    }
  for (i = EVALEQ_UPSILON_OFF; i < EVALEQ_UPSILON_OFF + EVALEQ_UPSILON_N; i++)
    {
      R2primei = R2prime[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime[i] = NULL;

      r1primei = r1prime[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime[i] = NULL;

      r0prime[i] = NULL;
    }
}

static void
_lnp_tbox_prove (spolymat_ptr R2prime_sz[], spolyvec_ptr r1prime_sz[],
                 poly_ptr r0prime_sz[], uint8_t hash[32], polyvec_t tB,
                 polyvec_t z3, polyvec_t z4, polyvec_t s1, polyvec_t m,
                 polyvec_t s2, polymat_t Bprime, spolymat_ptr R2prime[],
                 spolyvec_ptr r1prime[], poly_ptr r0prime[], unsigned int M,
                 polymat_ptr Es[], polymat_ptr Em[], polyvec_ptr v[],
                 polymat_t Ps, polymat_t Pm, polyvec_t f, polymat_t Ds,
                 polymat_t Dm, polyvec_t u, const uint8_t seed_tbox[32],
                 const lnp_tbox_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  const unsigned int lext = params->tbox->lext;
#endif
  const unsigned int lambda = params->quad_eval->lambda;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int Z = params->Z;
  const unsigned int nprime = params->nprime;
  const unsigned int nex = params->nex;
  const unsigned int nbin = params->nbin;
  const unsigned int *ni = params->n;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int m2 = tbox->m2;
  const unsigned int l = tbox->l;
  const unsigned int kmsis = tbox->kmsis;
  unsigned int i, j;
  INTVEC_T (s3coeffs, nex * d, int_get_nlimbs (q));
  INTVEC_T (s4coeffs, nprime * d, int_get_nlimbs (q));
  INTVEC_T (y3coeffs, 256, int_get_nlimbs (q));
  INTVEC_T (y4coeffs, 256, int_get_nlimbs (q));
  INTVEC_T (z3coeffs, 256, int_get_nlimbs (q));
  INTVEC_T (z4coeffs, 256, int_get_nlimbs (q));
  polyvec_t tmp_polyvec, subv, s1_, m_, s21, y3_, y4_, ty3, ty4, tbeta, beta,
      upsilon, s3, s4, y3, y4, z3_, z4_;
  intvec_ptr coeffs;
#ifdef XXX
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
#endif
  polymat_t By3, By4, Bbeta;
  shake128_state_t hstate;
  coder_state_t cstate;
  rng_state_t rstate_signs;
  rng_state_t rstate_rej;
  uint32_t dom = 0;
  uint8_t rbits;
  unsigned int nrbits, outlen, loff, off;
  uint8_t out[CEIL (256 * 2 * log2q + d * log2q, 8) + 1];
  uint8_t cseed[32]; /* seed for challenge */
  poly_ptr poly;
  int beta3 = 0, beta4 = 0;
  int rej;

  memset (out, 0, CEIL (256 * 2 * log2q + d * log2q, 8) + 1); // XXX

  polyvec_alloc (s3, Rq, nex);
  polyvec_alloc (s4, Rq, nprime);
  polyvec_alloc (y3, Rq, 256 / d);
  polyvec_alloc (y4, Rq, 256 / d);
  polyvec_alloc (z3_, Rq, 256 / d);
  polyvec_alloc (z4_, Rq, 256 / d);
  polyvec_alloc (tmp_polyvec, Rq, MAX (nex, nprime));

  /* s1 = s1_,upsilon, m = m_,y3_,y4_,beta */
  polyvec_get_subvec (s1_, s1, 0, m1, 1);
  if (l > 0)
    polyvec_get_subvec (m_, m, 0, l, 1);
  polyvec_get_subvec (beta, m, l + (256 / d) * 2, 1, 1);
  polyvec_set_zero (beta);

  polyvec_get_subvec (s21, s2, 0, m2 - kmsis, 1);

  /* tB = tB_,ty,tbeta */
  loff = 0;
  if (nex > 0)
    {
      ASSERT_ERR (Ps == NULL || polymat_get_nrows (Ps) == nbin);
      ASSERT_ERR (Ps == NULL || polymat_get_ncols (Ps) == m1);
      ASSERT_ERR (Pm == NULL || polymat_get_nrows (Pm) == nbin);
      ASSERT_ERR (Pm == NULL || polymat_get_ncols (Pm) == l);
      ASSERT_ERR (f == NULL || polyvec_get_nelems (f) == nbin);
      for (i = 0; i < Z; i++)
        {
          ASSERT_ERR (Es == NULL || Es[i] == NULL
                      || polymat_get_nrows (Es[i]) == ni[i]);
          ASSERT_ERR (Es == NULL || Es[i] == NULL
                      || polymat_get_ncols (Es[i]) == m1);
          ASSERT_ERR (Em == NULL || Em[i] == NULL
                      || polymat_get_nrows (Em[i]) == ni[i]);
          ASSERT_ERR (Em == NULL || Em[i] == NULL
                      || polymat_get_ncols (Em[i]) == l);
          ASSERT_ERR (v == NULL || v[i] == NULL
                      || polyvec_get_nelems (v[i]) == ni[i]);
        }

      polyvec_get_subvec (upsilon, s1, m1, Z, 1);
      polyvec_get_subvec (y3_, m, l + loff, 256 / d, 1);
      polyvec_get_subvec (ty3, tB, l + loff, 256 / d, 1);
      polymat_get_submat (By3, Bprime, l + loff, 0, 256 / d, m2 - kmsis, 1, 1);

      polyvec_set_coeffvec2 (s3, s3coeffs);
      polyvec_set_coeffvec2 (y3, y3coeffs);
      polyvec_set_coeffvec2 (z3_, z3coeffs);

      /* s3 */

      off = 0;
      if (nbin > 0)
        {
          polyvec_get_subvec (subv, s3, 0, nbin, 1);
          off += nbin;

          if (f != NULL)
            {
              polyvec_set (subv, f);
              polyvec_fromcrt (subv);
              polyvec_mod (subv, subv);
              polyvec_redc (subv, subv);
            }
          else
            {
              polyvec_set_zero (subv);
            }
          if (Ps != NULL)
            {
              polyvec_addmul (subv, Ps, s1_, 0);
              polyvec_fromcrt (subv);
            }
          if (Pm != NULL)
            {
              polyvec_addmul (subv, Pm, m_, 0);
              polyvec_fromcrt (subv);
            }
        }
      if (Z > 0)
        {
          for (i = 0; i < Z; i++)
            {
              polyvec_get_subvec (subv, s3, off, ni[i], 1);
              off += ni[i];

              if (v != NULL && v[i] != NULL)
                {
                  polyvec_set (subv, v[i]);
                  polyvec_fromcrt (subv);
                  polyvec_mod (subv, subv);
                  polyvec_redc (subv, subv);
                }
              else
                {
                  polyvec_set_zero (subv);
                }
              if (Es != NULL && Es[i] != NULL)
                {
                  polyvec_addmul (subv, Es[i], s1_, 0);
                  polyvec_fromcrt (subv);
                }
              if (Em != NULL && Em[i] != NULL)
                {
                  polyvec_addmul (subv, Em[i], m_, 0);
                }
              polyvec_fromcrt (subv);
            }

          polyvec_get_subvec (subv, s3, off, Z, 1);
          polyvec_set (subv, upsilon);
        }

      loff += 256 / d;
    }
  if (nprime > 0)
    {
      ASSERT_ERR (Ds == NULL || polymat_get_nrows (Ds) == nprime);
      ASSERT_ERR (Ds == NULL || polymat_get_ncols (Ds) == m1);
      ASSERT_ERR (Dm == NULL || polymat_get_nrows (Dm) == nprime);
      ASSERT_ERR (Dm == NULL || polymat_get_ncols (Dm) == l);
      ASSERT_ERR (u == NULL || polyvec_get_nelems (u) == nprime);

      polyvec_get_subvec (y4_, m, l + loff, 256 / d, 1);
      polyvec_get_subvec (ty4, tB, l + loff, 256 / d, 1);
      polymat_get_submat (By4, Bprime, l + loff, 0, 256 / d, m2 - kmsis, 1, 1);

      polyvec_set_coeffvec2 (s4, s4coeffs);
      polyvec_set_coeffvec2 (y4, y4coeffs);
      polyvec_set_coeffvec2 (z4_, z4coeffs);

      /* s4 */

      if (u != NULL)
        {
          polyvec_set (s4, u);
          polyvec_fromcrt (s4);
          polyvec_mod (s4, s4);
          polyvec_redc (s4, s4);
        }
      else
        {
          polyvec_set_zero (s4);
        }
      if (Ds != NULL)
        {
          polyvec_addmul (s4, Ds, s1_, 0);
          polyvec_fromcrt (s4);
        }
      if (Dm != NULL)
        {
          polyvec_addmul (s4, Dm, m_, 0);
          polyvec_fromcrt (s4);
        }

      loff += 256 / d;
    }

  ASSERT_ERR (lext == (loff + 1) + (lambda / 2 + 1));
  polyvec_get_subvec (tbeta, tB, l + loff, 1, 1);
  polymat_get_submat (Bbeta, Bprime, l + loff, 0, 1, m2 - kmsis, 1, 1);

  nrbits = 0;
  rng_init (rstate_rej, seed_tbox, dom++);
  rng_init (rstate_signs, seed_tbox, dom++);

  while (1)
    {
      /* sample signs */
      if (nrbits == 0)
        {
          rng_urandom (rstate_signs, &rbits, 1);
          nrbits = 8;
        }

      if (nex > 0)
        {
          /* y3, append to m  */
          polyvec_grandom (y3, params->log2stdev3, seed_tbox, dom++);
          polyvec_set (y3_, y3);

          /* ty3 */
          polyvec_set (ty3, y3);
          polyvec_addmul (ty3, By3, s21, 0);
          polyvec_mod (ty3, ty3);
          polyvec_redp (ty3, ty3);

          /* beta3  */
          beta3 = (rbits & (1 << (8 - nrbits))) >> (8 - nrbits);
          beta3 = 1 - 2 * beta3; /* {0,1} -> {1,-1} */
          nrbits -= 1;
        }

      if (nprime > 0)
        {
          /* y4, append to m  */
          polyvec_grandom (y4, params->log2stdev4, seed_tbox, dom++);
          polyvec_set (y4_, y4);

          /* ty4 */
          polyvec_set (ty4, y4);
          polyvec_addmul (ty4, By4, s21, 0);
          polyvec_mod (ty4, ty4);
          polyvec_redp (ty4, ty4);

          beta4 = (rbits & (1 << (8 - nrbits + 1))) >> (8 - nrbits + 1);
          beta4 = 1 - 2 * beta4; /* {0,1} -> {1,-1} */
          nrbits -= 1;
        }

      /* tbeta */
      poly = polyvec_get_elem (beta, 0);
      coeffs = poly_get_coeffvec (poly);
      if (nex > 0)
        intvec_set_elem_i64 (coeffs, 0, beta3);
      if (nprime > 0)
        intvec_set_elem_i64 (coeffs, d / 2, beta4);
      polyvec_set (tbeta, beta);
      polyvec_addmul (tbeta, Bbeta, s21, 0);
      polyvec_mod (tbeta, tbeta);
      polyvec_redp (tbeta, tbeta);

      /* encode ty, tbeta, hash of encoding is seed for challenges */

      coder_enc_begin (cstate, out);
      if (nex > 0)
        coder_enc_urandom3 (cstate, ty3, q, log2q);
      if (nprime > 0)
        coder_enc_urandom3 (cstate, ty4, q, log2q);
      coder_enc_urandom3 (cstate, tbeta, q, log2q);
      coder_enc_end (cstate);

      outlen = coder_get_offset (cstate);
      ASSERT_ERR (outlen % 8 == 0);
      ASSERT_ERR (outlen / 8 <= CEIL (256 * 2 * log2q + d * log2q, 8) + 1);
      outlen >>= 3; /* nbits to nbytes */

      shake128_init (hstate);
      shake128_absorb (hstate, hash, 32);
      shake128_absorb (hstate, out, outlen);
      shake128_squeeze (hstate, cseed, 32);

      if (nex > 0)
        {
          INT_T (beta3Rijs3j, int_get_nlimbs (q));
          int8_t Ri[nex * d];
          int_ptr s3coeff, Rs3coeff;

          polyvec_fromcrt (s3);
          polyvec_fromcrt (y3);

          polyvec_set (z3_, y3);
          intvec_set_zero (y3coeffs);

          for (i = 0; i < 256; i++)
            {
              Rs3coeff = intvec_get_elem (y3coeffs, i);

              _expand_R_i (Ri, nex * d, i, cseed);

              for (j = 0; j < nex * d; j++)
                {
                  if (Ri[j] == 0)
                    {
                    }
                  else
                    {
                      ASSERT_ERR (Ri[j] == 1 || Ri[j] == -1);

                      s3coeff = intvec_get_elem (s3coeffs, j);

                      int_set (beta3Rijs3j, s3coeff);
                      int_mul_sgn_self (beta3Rijs3j, Ri[j]);
                      int_add (Rs3coeff, Rs3coeff, beta3Rijs3j);
                    }
                }
            }
          intvec_mul_sgn_self (y3coeffs, beta3);
          intvec_add (z3coeffs, z3coeffs, y3coeffs);
        }

      if (nprime > 0)
        {
          INT_T (beta4Rprimeijs4j, int_get_nlimbs (q));
          int8_t Rprimei[nprime * d];
          int_ptr s4coeff, Rs4coeff;

          polyvec_fromcrt (s4);
          polyvec_fromcrt (y4);

          polyvec_set (z4_, y4);
          intvec_set_zero (y4coeffs);

          for (i = 0; i < 256; i++)
            {
              Rs4coeff = intvec_get_elem (y4coeffs, i);

              _expand_Rprime_i (Rprimei, nprime * d, i, cseed);

              for (j = 0; j < nprime * d; j++)
                {
                  if (Rprimei[j] == 0)
                    {
                    }
                  else
                    {
                      ASSERT_ERR (Rprimei[j] == 1 || Rprimei[j] == -1);

                      s4coeff = intvec_get_elem (s4coeffs, j);

                      int_set (beta4Rprimeijs4j, s4coeff);
                      int_mul_sgn_self (beta4Rprimeijs4j, Rprimei[j]);
                      int_add (Rs4coeff, Rs4coeff, beta4Rprimeijs4j);
                    }
                }
            }
          intvec_mul_sgn_self (y4coeffs, beta4);
          intvec_add (z4coeffs, z4coeffs, y4coeffs);
        }

      /* rejection sampling */

      ASSERT_ERR (params->rej3 == 0 || params->rej3 == 2);
      if (nex > 0 && params->rej3)
        {
          intvec_mul_sgn_self (y3coeffs, beta3); /* revert mul by beta3 */

          rej = rej_bimodal (rstate_rej, z3coeffs, y3coeffs, params->scM3,
                             params->stdev3sqr);
          if (rej)
            continue;
        }
      ASSERT_ERR (params->rej4 == 0 || params->rej4 == 2);
      if (nprime > 0 && params->rej4)
        {
          intvec_mul_sgn_self (y4coeffs, beta4); /* revert mul by beta4 */

          rej = rej_bimodal (rstate_rej, z4coeffs, y4coeffs, params->scM4,
                             params->stdev4sqr);
          if (rej)
            continue;
        }

      break;
    }

  /* additional schwarz-zippel */
  {
    lnp_quad_eval_params_srcptr quade = params->quad_eval;
    abdlop_params_srcptr quad_eval = quade->quad_eval;
    spolymat_t M0, M1;
    spolyvec_t V0, V1;
    INTVEC_T (chal, lambda, q->nlimbs);
    const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
    int_ptr ip;

    spolymat_alloc (M0, Rq, np, np, NELEMS_DIAG (np));
    spolymat_alloc (M1, Rq, np, np, NELEMS_DIAG (np));
    spolyvec_alloc (V0, Rq, np, np);
    spolyvec_alloc (V1, Rq, np, np);

#if 1
    STOPWATCH_START (stopwatch_lnp_tbox_prove_build_jl_eqs,
                     "lnp_tbox_prove_build_jl_eqs");

    _evaleq_z3 (R2prime, r1prime, r0prime, Es, Em, v, Ps, Pm, f, z3coeffs,
                cseed, params, tmp_polyvec, 0);
    _evaleq_z4 (R2prime, r1prime, r0prime, Ds, Dm, u, z4coeffs, cseed, params,
                tmp_polyvec, 0);

    STOPWATCH_STOP (stopwatch_lnp_tbox_prove_build_jl_eqs);
#endif

    STOPWATCH_START (stopwatch_lnp_tbox_prove_schwartz_zippel,
                     "lnp_tbox_prove_schwartz_zippel");

    for (i = 0; i < M; i++)
      {
#if 0
        if (i >= EVALEQ_Z3_OFF && i < EVALEQ_Z3_OFF + EVALEQ_Z3_N)
          {
            R2primei = _alloc (sizeof (spolymat_t));
            spolymat_alloc (R2primei, Rq, np, np, EVALEQ_Z3_MAXQ);
            R2prime[i] = R2primei;

            r1primei = _alloc (sizeof (spolyvec_t));
            spolyvec_alloc (r1primei, Rq, np, np);
            r1prime[i] = r1primei;

            r0primei = _alloc (sizeof (poly_t));
            poly_alloc (r0primei, Rq);
            r0prime[i] = r0primei;

            _evaleq_z3 (R2prime, r1prime, r0prime, Es, Em, v, Ps, Pm, f,
                        z3coeffs, cseed, params, tmp_polyvec,
                        i - EVALEQ_Z3_OFF);
          }
        if (i >= EVALEQ_Z4_OFF && i < EVALEQ_Z4_OFF + EVALEQ_Z4_N)
          {
            R2primei = _alloc (sizeof (spolymat_t));
            spolymat_alloc (R2primei, Rq, np, np, EVALEQ_Z4_MAXQ);
            R2prime[i] = R2primei;

            r1primei = _alloc (sizeof (spolyvec_t));
            spolyvec_alloc (r1primei, Rq, np, np);
            r1prime[i] = r1primei;

            r0primei = _alloc (sizeof (poly_t));
            poly_alloc (r0primei, Rq);
            r0prime[i] = r0primei;

            _evaleq_z4 (R2prime, r1prime, r0prime, Ds, Dm, u, z4coeffs, cseed,
                        params, tmp_polyvec, i - EVALEQ_Z4_OFF);
          }
#endif

        intvec_urandom (chal, q, log2q, cseed, 2 * 256 + i);

        for (j = 0; j < lambda; j++)
          {
            ip = intvec_get_elem (chal, j);

            if (R2prime[i] != NULL)
              {
                spolymat_scale (M0, ip, R2prime[i]);
                spolymat_add (M1, R2prime_sz[j], M0, 0);
                spolymat_set (R2prime_sz[j], M1);
              }

            if (r1prime[i] != NULL)
              {
                spolyvec_scale (V0, ip, r1prime[i]);
                spolyvec_add (V1, r1prime_sz[j], V0, 0);
                spolyvec_set (r1prime_sz[j], V1);
              }

            if (r0prime[i] != NULL)
              poly_addscale (r0prime_sz[j], ip, r0prime[i], 0);
          }

#if 0
        if (i >= EVALEQ_Z3_OFF && i < EVALEQ_Z4_OFF + EVALEQ_Z4_N)
          {
            spolymat_free (R2prime[i]);
            spolyvec_free (r1prime[i]);
            poly_free (r0prime[i]);
          }
#endif
      }

    STOPWATCH_STOP (stopwatch_lnp_tbox_prove_schwartz_zippel);

    spolymat_free (M0);
    spolymat_free (M1);
    spolyvec_free (V0);
    spolyvec_free (V1);
    polyvec_free (tmp_polyvec);
  }

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);

  /* output proof (h,c,z1,z21,hint,z3,z4) */
  polyvec_set (z3, z3_);
  polyvec_set (z4, z4_);

  /* cleanup */
  rng_clear (rstate_signs);
  rng_clear (rstate_rej);
  polyvec_free (s3);
  polyvec_free (s4);
  polyvec_free (y3);
  polyvec_free (y4);
  polyvec_free (z3_);
  polyvec_free (z4_);
}

/*
 * hash hash of tA1, tB
 * s1 = (s1_,upsilon)
 *
 * scratch space:
 * R2,r1,r0 are N+2+lambda/2 arrays
 * R2 is 2*(m1+l)+lambda x 2*(m1+l)+lambda
 */
void
lnp_tbox_prove (uint8_t hash[32], polyvec_t tB, polyvec_t h, poly_t c,
                polyvec_t z1, polyvec_t z21, polyvec_t hint, polyvec_t z3,
                polyvec_t z4, polyvec_t s1, polyvec_t m, polyvec_t s2,
                polyvec_t tA2, polymat_t A1, polymat_t A2prime,
                polymat_t Bprime, spolymat_ptr R2[], spolyvec_ptr r1[],
                unsigned int N, spolymat_ptr R2prime[], spolyvec_ptr r1prime[],
                poly_ptr r0prime[], unsigned int M, polymat_ptr Es[],
                polymat_ptr Em[], polyvec_ptr v[], polymat_t Ps, polymat_t Pm,
                polyvec_t f, polymat_t Ds, polymat_t Dm, polyvec_t u,
                const uint8_t seed[32], const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  abdlop_params_srcptr quad_eval = params->quad_eval->quad_eval;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int Z = params->Z;
  const unsigned int M_ = 2 * (d - 1) + 256 + 256 + Z + 1 + Z;
  const unsigned int lambda = params->quad_eval->lambda;
#if ASSERT == ASSERT_ENABLED
  abdlop_params_srcptr quad = params->quad_eval->quad_many;
  const unsigned int lext = params->tbox->lext;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int m2 = tbox->m2;
  const unsigned int l = tbox->l;
  const unsigned int kmsis = tbox->kmsis;
  const unsigned int nprime = params->nprime;
  const unsigned int nex = params->nex;
  const unsigned int N_ = lambda / 2 + 2;
#endif
  unsigned int i;
  spolymat_ptr R2prime_sz[lambda], R2primei;
  spolyvec_ptr r1prime_sz[lambda], r1primei;
  poly_ptr r0prime_sz[lambda], r0primei;
  shake128_state_t hstate;
  uint8_t expseed[64];
  const uint8_t *seed_tbox = expseed;
  const uint8_t *seed_quad_eval = expseed + 32;
  const unsigned int np = 2 * (tbox->m1 + quad_eval->l);

  STOPWATCH_START (stopwatch_lnp_tbox_prove, "lnp_tbox_prove");

  ASSERT_ERR (nex + nprime
              > 0); /* use quad_eval if no norm proof is needed. */
  ASSERT_ERR (lambda % 2 == 0);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polyvec_get_nelems (h) == lambda / 2);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1 + Z);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_nelems (hint) == kmsis);
  ASSERT_ERR (nex == 0 || polyvec_get_nelems (z3) == 256 / d);
  ASSERT_ERR (nprime == 0 || polyvec_get_nelems (z4) == 256 / d);
  ASSERT_ERR (polyvec_get_nelems (s1) == m1 + Z);
  ASSERT_ERR (polyvec_get_nelems (s2) == m2);
  ASSERT_ERR (polyvec_get_nelems (m) == l + lext);
  ASSERT_ERR (polyvec_get_nelems (tA2) == kmsis);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1 + Z);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == m2 - kmsis);
#if ASSERT == ASSERT_ENABLED
  for (i = N_; i < N_ + N; i++)
    {
      ASSERT_ERR (i >= N || spolymat_is_upperdiag (R2[i]));
      ASSERT_ERR (R2[i] == NULL
                  || spolymat_get_nrows (R2[i]) == 2 * (tbox->m1 + quad->l));
      ASSERT_ERR (R2[i] == NULL
                  || spolymat_get_ncols (R2[i]) == 2 * (tbox->m1 + quad->l));
      ASSERT_ERR (r1[i] == NULL
                  || r1[i]->nelems_max == 2 * (tbox->m1 + quad->l));
    }
  for (i = M_; i < M_ + M; i++)
    {
      ASSERT_ERR (spolymat_is_upperdiag (R2prime[i]));
      ASSERT_ERR (R2prime[i] == NULL
                  || spolymat_get_nrows (R2prime[i])
                         == 2 * (tbox->m1 + quad_eval->l));
      ASSERT_ERR (R2prime[i] == NULL
                  || spolymat_get_ncols (R2prime[i])
                         == 2 * (tbox->m1 + quad_eval->l));
      ASSERT_ERR (r1prime[i] == NULL
                  || r1prime[i]->nelems_max == 2 * (tbox->m1 + quad_eval->l));
    }
#endif

  for (i = 0; i < lambda; i++)
    {
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, NELEMS_DIAG (np));
      R2prime_sz[i] = R2primei;
      spolymat_set_empty (R2prime_sz[i]);

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime_sz[i] = r1primei;
      spolyvec_set_empty (r1prime_sz[i]);

      r0primei = _alloc (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz[i] = r0primei;
      poly_set_zero (r0prime_sz[i]);
    }

  /*
   * Expand input seed into two seeds: one for tbox
   * and one for the sub-protocol quad_eval.
   */
  shake128_init (hstate);
  shake128_absorb (hstate, seed, 32);
  shake128_squeeze (hstate, expseed, sizeof (expseed));
  shake128_clear (hstate);

  _lnp_tbox_prove (R2prime_sz, r1prime_sz, r0prime_sz, hash, tB, z3, z4, s1, m,
                   s2, Bprime, R2prime, r1prime, r0prime, M_ + M, Es, Em, v,
                   Ps, Pm, f, Ds, Dm, u, seed_tbox, params);

  lnp_quad_eval_prove (hash, tB, h, c, z1, z21, hint, s1, m, s2, tA2, A1,
                       A2prime, Bprime, R2, r1, 2 + N, R2prime_sz, r1prime_sz,
                       r0prime_sz, lambda, seed_quad_eval, params->quad_eval);

  for (i = 0; i < lambda; i++)
    {
      R2primei = R2prime_sz[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime_sz[i] = NULL;

      r1primei = r1prime_sz[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime_sz[i] = NULL;

      r0primei = r0prime_sz[i];
      poly_free (r0primei);
      _free (r0primei, sizeof (poly_t));
      r0prime_sz[i] = NULL;
    }
  STOPWATCH_STOP (stopwatch_lnp_tbox_prove);
}

static int
_lnp_tbox_verify (spolymat_ptr R2prime_sz[], spolyvec_ptr r1prime_sz[],
                  poly_ptr r0prime_sz[], uint8_t hash[32], polyvec_t z3,
                  polyvec_t z4, polyvec_t tB, spolymat_ptr R2prime[],
                  spolyvec_ptr r1prime[], poly_ptr r0prime[], unsigned int M,
                  polymat_ptr Es[], polymat_ptr Em[], polyvec_ptr v[],
                  polymat_t Ps, polymat_t Pm, polyvec_t f, polymat_t Ds,
                  polymat_t Dm, polyvec_t u, const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  const unsigned int Z = params->Z;
  const unsigned int lambda = params->quad_eval->lambda;
  UNUSED const unsigned int m1 = tbox->m1 - Z;
#if ASSERT == ASSERT_ENABLED
  const unsigned int lext = params->tbox->lext;
  const unsigned int *ni = params->n;
#endif
  const unsigned int nbin = params->nbin;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int l = tbox->l;
  INT_T (linf, int_get_nlimbs (q));
  INT_T (l2sqr, 2 * int_get_nlimbs (q));
  INTVEC_T (z3coeffs, 256, int_get_nlimbs (q));
  INTVEC_T (z4coeffs, 256, int_get_nlimbs (q));
  unsigned int i, j;
  polyvec_t tmp_polyvec, ty3, ty4, tbeta;
#ifdef XXX
  spolymat_ptr R2primei;
  spolyvec_ptr r1primei;
  poly_ptr r0primei;
#endif
  intvec_t isubv;
  poly_ptr poly;
  intvec_ptr coeffs;
  shake128_state_t hstate;
  coder_state_t cstate;
  unsigned int outlen, loff;
  uint8_t out[CEIL (256 * 2 * log2q + d * log2q, 8) + 1];
  uint8_t cseed[32]; /* seed for challenge */
  int b = 0;

  memset (out, 0, CEIL (256 * 2 * log2q + d * log2q, 8) + 1); // XXX

  loff = 0;
  if (nex > 0)
    {
      ASSERT_ERR (Ps == NULL || polymat_get_nrows (Ps) == nbin);
      ASSERT_ERR (Ps == NULL || polymat_get_ncols (Ps) == m1);
      ASSERT_ERR (Pm == NULL || polymat_get_nrows (Pm) == nbin);
      ASSERT_ERR (Pm == NULL || polymat_get_ncols (Pm) == l);
      ASSERT_ERR (f == NULL || polyvec_get_nelems (f) == nbin);

      polyvec_get_subvec (ty3, tB, l + loff, 256 / d, 1);
      polyvec_mod (ty3, ty3);
      polyvec_redp (ty3, ty3);

      loff += 256 / d;
    }
  if (nprime > 0)
    {
      ASSERT_ERR (Ds == NULL || polymat_get_nrows (Ds) == nprime);
      ASSERT_ERR (Ds == NULL || polymat_get_ncols (Ds) == m1);
      ASSERT_ERR (Dm == NULL || polymat_get_nrows (Dm) == nprime);
      ASSERT_ERR (Dm == NULL || polymat_get_ncols (Dm) == l);
      ASSERT_ERR (u == NULL || polyvec_get_nelems (u) == nprime);
      for (i = 0; i < Z; i++)
        {
          ASSERT_ERR (Es[i] == NULL || polymat_get_nrows (Es[i]) == ni[i]);
          ASSERT_ERR (Es[i] == NULL || polymat_get_ncols (Es[i]) == m1);
          ASSERT_ERR (Em[i] == NULL || polymat_get_nrows (Em[i]) == ni[i]);
          ASSERT_ERR (Em[i] == NULL || polymat_get_ncols (Em[i]) == l);
          ASSERT_ERR (v[i] == NULL || polyvec_get_nelems (v[i]) == ni[i]);
        }

      polyvec_get_subvec (ty4, tB, l + loff, 256 / d, 1);
      polyvec_mod (ty4, ty4);
      polyvec_redp (ty4, ty4);

      loff += 256 / d;
    }
  ASSERT_ERR (lext == (loff + 1) + (lambda / 2 + 1));

  polyvec_alloc (tmp_polyvec, Rq, MAX (nbin, MAX (nex, nprime)));

  polyvec_get_subvec (tbeta, tB, l + loff, 1, 1);
  polyvec_mod (tbeta, tbeta);
  polyvec_redp (tbeta, tbeta);

  /* check  bounds */

  if (nex > 0)
    {
      polyvec_fromcrt (z3);
      polyvec_l2sqr (l2sqr, z3);
      if (int_gt (l2sqr, params->Bz3sqr))
        goto ret;
    }
  if (nprime > 0)
    {
      polyvec_fromcrt (z4);
      polyvec_linf (linf, z4);
      if (int_gt (linf, params->Bz4))
        goto ret;
    }

  /* encode ty, tbeta, hash of encoding is seed for challenges */

  coder_enc_begin (cstate, out);
  if (nex > 0)
    coder_enc_urandom3 (cstate, ty3, q, log2q);
  if (nprime > 0)
    coder_enc_urandom3 (cstate, ty4, q, log2q);
  coder_enc_urandom3 (cstate, tbeta, q, log2q);
  coder_enc_end (cstate);

  outlen = coder_get_offset (cstate);
  ASSERT_ERR (outlen % 8 == 0);
  ASSERT_ERR (outlen / 8 <= CEIL (256 * 2 * log2q + d * log2q, 8) + 1);
  outlen >>= 3; /* nbits to nbytes */

  /* recover challenge */
  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, outlen);
  shake128_squeeze (hstate, cseed, 32);
  shake128_clear (hstate);

  /* get z3 and z4 coefficient vectors */
  for (i = 0; i < 256 / d; i++)
    {
      intvec_get_subvec (isubv, z3coeffs, i * d, d, 1);
      poly = polyvec_get_elem (z3, i);
      coeffs = poly_get_coeffvec (poly);
      intvec_set (isubv, coeffs);

      intvec_get_subvec (isubv, z4coeffs, i * d, d, 1);
      poly = polyvec_get_elem (z4, i);
      coeffs = poly_get_coeffvec (poly);
      intvec_set (isubv, coeffs);
    }

  /* additional schwarz-zippel */
  {
    lnp_quad_eval_params_srcptr quade = params->quad_eval;
    abdlop_params_srcptr quad_eval = quade->quad_eval;
    spolymat_t M0, M1;
    spolyvec_t V0, V1;
    INTVEC_T (chal, lambda, q->nlimbs);
    const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
    int_ptr ip;

    spolymat_alloc (M0, Rq, np, np, NELEMS_DIAG (np));
    spolymat_alloc (M1, Rq, np, np, NELEMS_DIAG (np));
    spolyvec_alloc (V0, Rq, np, np);
    spolyvec_alloc (V1, Rq, np, np);

#if 1
    STOPWATCH_START (stopwatch_lnp_tbox_verify_build_jl_eqs,
                     "lnp_tbox_verify_build_jl_eqs");

    _evaleq_z3 (R2prime, r1prime, r0prime, Es, Em, v, Ps, Pm, f, z3coeffs,
                cseed, params, tmp_polyvec, 0);
    _evaleq_z4 (R2prime, r1prime, r0prime, Ds, Dm, u, z4coeffs, cseed, params,
                tmp_polyvec, 0);

    STOPWATCH_STOP (stopwatch_lnp_tbox_verify_build_jl_eqs);
#endif

    STOPWATCH_START (stopwatch_lnp_tbox_verify_schwartz_zippel,
                     "lnp_tbox_verify_schwartz_zippel");

    for (i = 0; i < M; i++)
      {
#if 0
        if (i >= EVALEQ_Z3_OFF && i < EVALEQ_Z3_OFF + EVALEQ_Z3_N)
          {
            R2primei = _alloc (sizeof (spolymat_t));
            spolymat_alloc (R2primei, Rq, np, np, EVALEQ_Z3_MAXQ);
            R2prime[i] = R2primei;

            r1primei = _alloc (sizeof (spolyvec_t));
            spolyvec_alloc (r1primei, Rq, np, np);
            r1prime[i] = r1primei;

            r0primei = _alloc (sizeof (poly_t));
            poly_alloc (r0primei, Rq);
            r0prime[i] = r0primei;

            _evaleq_z3 (R2prime, r1prime, r0prime, Es, Em, v, Ps, Pm, f,
                        z3coeffs, cseed, params, tmp_polyvec,
                        i - EVALEQ_Z3_OFF);
          }
        if (i >= EVALEQ_Z4_OFF && i < EVALEQ_Z4_OFF + EVALEQ_Z4_N)
          {
            R2primei = _alloc (sizeof (spolymat_t));
            spolymat_alloc (R2primei, Rq, np, np, EVALEQ_Z4_MAXQ);
            R2prime[i] = R2primei;

            r1primei = _alloc (sizeof (spolyvec_t));
            spolyvec_alloc (r1primei, Rq, np, np);
            r1prime[i] = r1primei;

            r0primei = _alloc (sizeof (poly_t));
            poly_alloc (r0primei, Rq);
            r0prime[i] = r0primei;

            _evaleq_z4 (R2prime, r1prime, r0prime, Ds, Dm, u, z4coeffs, cseed,
                        params, tmp_polyvec, i - EVALEQ_Z4_OFF);
          }
#endif

        intvec_urandom (chal, q, log2q, cseed, 2 * 256 + i);

        for (j = 0; j < lambda; j++)
          {
            ip = intvec_get_elem (chal, j);

            if (R2prime[i] != NULL)
              {
                spolymat_scale (M0, ip, R2prime[i]);
                spolymat_add (M1, R2prime_sz[j], M0, 0);
                spolymat_set (R2prime_sz[j], M1);
              }

            if (r1prime[i] != NULL)
              {
                spolyvec_scale (V0, ip, r1prime[i]);
                spolyvec_add (V1, r1prime_sz[j], V0, 0);
                spolyvec_set (r1prime_sz[j], V1);
              }

            if (r0prime[i] != NULL)
              poly_addscale (r0prime_sz[j], ip, r0prime[i], 0);
          }

#if 0
        if (i >= EVALEQ_Z3_OFF && i < EVALEQ_Z4_OFF + EVALEQ_Z4_N)
          {
            spolymat_free (R2prime[i]);
            spolyvec_free (r1prime[i]);
            poly_free (r0prime[i]);
          }
#endif
      }
    STOPWATCH_STOP (stopwatch_lnp_tbox_verify_schwartz_zippel);

    spolymat_free (M0);
    spolymat_free (M1);
    spolyvec_free (V0);
    spolyvec_free (V1);
  }

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);

  b = 1;
ret:
  /* cleanup */
  polyvec_free (tmp_polyvec);
  return b;
}

int
lnp_tbox_verify (uint8_t hash[32], polyvec_t h, poly_t c, polyvec_t z1,
                 polyvec_t z21, polyvec_t hint, polyvec_t z3, polyvec_t z4,
                 polyvec_t tA1, polyvec_t tB, polymat_t A1, polymat_t A2prime,
                 polymat_t Bprime, spolymat_ptr R2[], spolyvec_ptr r1[],
                 poly_ptr r0[], unsigned int N, spolymat_ptr R2prime[],
                 spolyvec_ptr r1prime[], poly_ptr r0prime[], unsigned int M,
                 polymat_ptr Es[], polymat_ptr Em[], polyvec_ptr v[],
                 polymat_t Ps, polymat_t Pm, polyvec_t f, polymat_t Ds,
                 polymat_t Dm, polyvec_t u, const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  abdlop_params_srcptr quad_eval = params->quad_eval->quad_eval;
  const unsigned int Z = params->Z;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int M_ = 2 * (d - 1) + 256 + 256 + Z + 1 + Z;
  const unsigned int np = 2 * (tbox->m1 + quad_eval->l);
  const unsigned int lambda = params->quad_eval->lambda;
#if ASSERT == ASSERT_ENABLED
  const unsigned int lext = params->tbox->lext;
  const unsigned int kmsis = params->tbox->kmsis;
  const unsigned int m2 = params->tbox->m2;
  abdlop_params_srcptr quad = params->quad_eval->quad_many;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int nex = params->nex;
  const unsigned int l = tbox->l;
  const unsigned int nprime = params->nprime;
  const unsigned int N_ = lambda / 2 + 2;
#endif
  spolymat_ptr R2prime_sz[lambda], R2primei;
  spolyvec_ptr r1prime_sz[lambda], r1primei;
  poly_ptr r0prime_sz[lambda], r0primei;
  unsigned int i;
  int b;

  STOPWATCH_START (stopwatch_lnp_tbox_verify, "lnp_tbox_verify");

  ASSERT_ERR (lambda % 2 == 0);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polyvec_get_nelems (h) == lambda / 2);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1 + Z);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_nelems (hint) == kmsis);
  ASSERT_ERR (nex == 0 || polyvec_get_nelems (z3) == 256 / d);
  ASSERT_ERR (nprime == 0 || polyvec_get_nelems (z4) == 256 / d);
  ASSERT_ERR (polyvec_get_nelems (tA1) == kmsis);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1 + Z);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == m2 - kmsis);
#if ASSERT == ASSERT_ENABLED
  for (i = N_; i < N_ + N; i++)
    {
      ASSERT_ERR (i >= N || spolymat_is_upperdiag (R2[i]));
      ASSERT_ERR (R2[i] == NULL
                  || spolymat_get_nrows (R2[i]) == 2 * (tbox->m1 + quad->l));
      ASSERT_ERR (R2[i] == NULL
                  || spolymat_get_ncols (R2[i]) == 2 * (tbox->m1 + quad->l));
      ASSERT_ERR (r1[i] == NULL
                  || r1[i]->nelems_max == 2 * (tbox->m1 + quad->l));
    }
  for (i = M_; i < M_ + M; i++)
    {
      ASSERT_ERR (spolymat_is_upperdiag (R2prime[i]));
      ASSERT_ERR (R2prime[i] == NULL
                  || spolymat_get_nrows (R2prime[i])
                         == 2 * (tbox->m1 + quad_eval->l));
      ASSERT_ERR (R2prime[i] == NULL
                  || spolymat_get_ncols (R2prime[i])
                         == 2 * (tbox->m1 + quad_eval->l));
      ASSERT_ERR (r1prime[i] == NULL
                  || r1prime[i]->nelems_max == 2 * (tbox->m1 + quad_eval->l));
    }
#endif

  for (i = 0; i < lambda; i++)
    {
      R2primei = _alloc (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np, np, NELEMS_DIAG (np));
      R2prime_sz[i] = R2primei;
      spolymat_set_empty (R2prime_sz[i]);

      r1primei = _alloc (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np, np);
      r1prime_sz[i] = r1primei;
      spolyvec_set_empty (r1prime_sz[i]);

      r0primei = _alloc (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz[i] = r0primei;
      poly_set_zero (r0prime_sz[i]);
    }

  b = _lnp_tbox_verify (R2prime_sz, r1prime_sz, r0prime_sz, hash, z3, z4, tB,
                        R2prime, r1prime, r0prime, M_ + M, Es, Em, v, Ps, Pm,
                        f, Ds, Dm, u, params);
  if (b != 1)
    goto ret;

  /* verify proof (h,c,z1,z21,hint,z3,z4) */
  b = lnp_quad_eval_verify (hash, h, c, z1, z21, hint, tA1, tB, A1, A2prime,
                            Bprime, R2, r1, r0, 2 + N, R2prime_sz, r1prime_sz,
                            r0prime_sz, lambda, params->quad_eval);
  if (b != 1)
    goto ret;

  b = 1;
ret:
  for (i = 0; i < lambda; i++)
    {
      R2primei = R2prime_sz[i];
      spolymat_free (R2primei);
      _free (R2primei, sizeof (spolymat_t));
      R2prime_sz[i] = NULL;

      r1primei = r1prime_sz[i];
      spolyvec_free (r1primei);
      _free (r1primei, sizeof (spolyvec_t));
      r1prime_sz[i] = NULL;

      r0primei = r0prime_sz[i];
      poly_free (r0primei);
      _free (r0primei, sizeof (poly_t));
      r0prime_sz[i] = NULL;
    }
  STOPWATCH_STOP (stopwatch_lnp_tbox_verify);
  return b;
}
