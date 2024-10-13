#include "lazer.h"
#include "lnp-tbox.h"
#include "stopwatch.h"

/*
 * XXX make these cases work (less eqs etc)
 * nex == 0 => no l2 and binary proof
 * nprime == 0 => no arp
 */

#ifdef XXX
/* R2, R2_ may not overlap */
static void
_scatter_mat (polymat_ptr R2, polymat_ptr R2_, unsigned int m1, unsigned int Z,
              unsigned int l)
{
  polymat_t subm, subm2;

  polymat_get_submat (subm2, R2_, 0, 0, 2 * m1, 2 * m1, 1, 1);
  polymat_get_submat (subm, R2, 0, 0, 2 * m1, 2 * m1, 1, 1);
  polymat_set (subm, subm2);

  polymat_get_submat (subm2, R2_, 0, 2 * m1, 2 * m1, 2 * l, 1, 1);
  polymat_get_submat (subm, R2, 0, 2 * (m1 + Z), 2 * m1, 2 * l, 1, 1);
  polymat_set (subm, subm2);

  polymat_get_submat (subm2, R2_, 2 * m1, 2 * m1, 2 * l, 2 * l, 1, 1);
  polymat_get_submat (subm, R2, 2 * (m1 + Z), 2 * (m1 + Z), 2 * l, 2 * l, 1,
                      1);
  polymat_set (subm, subm2);
}
#endif

/* R2 != R2_ */
static void
_scatter_smat (spolymat_ptr R2, spolymat_ptr R2_, unsigned int m1,
               unsigned int Z, unsigned int l)
{
  const unsigned int nelems = R2_->nelems;
  unsigned int i, row, col;
  poly_ptr poly, poly2;

  ASSERT_ERR (R2->nelems_max >= R2_->nelems);
  ASSERT_ERR (spolymat_is_upperdiag (R2_));

  (void)l; /* XXX unused */

  spolymat_set_empty (R2);

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

#ifdef XXX
/* r1, r1_ may not overlap */
static void
_scatter_vec (polyvec_ptr r1, polyvec_ptr r1_, unsigned int m1, unsigned int Z,
              unsigned int l)
{
  polyvec_t subv, subv2;

  polyvec_get_subvec (subv2, r1_, 0, 2 * m1, 1);
  polyvec_get_subvec (subv, r1, 0, 2 * m1, 1);
  polyvec_set (subv, subv2);

  polyvec_get_subvec (subv2, r1_, 2 * m1, 2 * l, 1);
  polyvec_get_subvec (subv, r1, 2 * (m1 + Z), 2 * l, 1);
  polyvec_set (subv, subv2);
}
#endif
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
_lnp_hash_statement_quadeqs (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int n = 2 * (m1 + l);
  const unsigned int lambda = params->quad_eval->lambda;
  const size_t buflen
      = CEIL ((((n * n - n) / 2 + n) + (n) + (1)) * d * log2q, 8) + 1;
  unsigned int i;
  shake128_state_t hstate;
  coder_state_t cstate;
  size_t len;
  uint8_t *buf;

  buf = _alloc (buflen);

  shake128_init (hstate);

  for (i = QUADEQ_INPUT_OFF; i < QUADEQ_INPUT_OFF + state->N; i++)
    {
      ASSERT_ERR (spolymat_get_nrows (state->R2[i]) == n);
      ASSERT_ERR (spolymat_get_ncols (state->R2[i]) == n);
      ASSERT_ERR (state->r1[i]->nelems_max == n);
      ASSERT_ERR (spolymat_is_upperdiag (state->R2[i]));

      coder_enc_begin (cstate, buf);
      if (state->R2[i] != NULL)
        {
          spolymat_redp (state->R2[i]);
          coder_enc_urandom5 (cstate, state->R2[i], q, log2q);
        }
      if (state->r1[i] != NULL)
        {
          spolyvec_redp (state->r1[i]);
          coder_enc_urandom6 (cstate, state->r1[i], q, log2q);
        }
      if (state->r0[i] != NULL)
        {
          poly_redp (state->r0[i], state->r0[i]);
          coder_enc_urandom2 (cstate, state->r0[i], q, log2q);
        }
      coder_enc_end (cstate);
      len = coder_get_offset (cstate);
      ASSERT_ERR (len % 8 == 0);
      ASSERT_ERR (len / 8 <= buflen);
      len >>= 3; /* nbits to nbytes */

      shake128_absorb (hstate, buf, len);
    }

  shake128_squeeze (hstate, state->hash_quadeqs, 32);
  shake128_clear (hstate);

  _free (buf, buflen);
}

static void
_lnp_hash_statement_evaleqs (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int n = 2 * (m1 + l);
  const size_t buflen
      = CEIL ((((n * n - n) / 2 + n) + (n) + (1)) * d * log2q, 8) + 1;
  unsigned int i;
  shake128_state_t hstate;
  coder_state_t cstate;
  size_t len;
  uint8_t *buf;

  buf = _alloc (buflen);

  shake128_init (hstate);

  for (i = EVALEQ_INPUT_OFF; i < EVALEQ_INPUT_OFF + state->M; i++)
    {
      ASSERT_ERR (spolymat_get_nrows (state->R2prime[i]) == n);
      ASSERT_ERR (spolymat_get_ncols (state->R2prime[i]) == n);
      ASSERT_ERR (state->r1prime[i]->nelems_max == n);
      ASSERT_ERR (spolymat_is_upperdiag (state->R2prime[i]));

      coder_enc_begin (cstate, buf);
      if (state->R2prime[i] != NULL)
        {
          spolymat_redp (state->R2prime[i]);
          coder_enc_urandom5 (cstate, state->R2prime[i], q, log2q);
        }
      if (state->r1prime[i] != NULL)
        {
          spolyvec_redp (state->r1prime[i]);
          coder_enc_urandom6 (cstate, state->r1prime[i], q, log2q);
        }
      if (state->r0prime[i] != NULL)
        {
          poly_redp (state->r0prime[i], state->r0prime[i]);
          coder_enc_urandom2 (cstate, state->r0prime[i], q, log2q);
        }
      coder_enc_end (cstate);
      len = coder_get_offset (cstate);
      ASSERT_ERR (len % 8 == 0);
      ASSERT_ERR (len / 8 <= buflen);
      len >>= 3; /* nbits to nbytes */

      shake128_absorb (hstate, buf, len);
    }

  shake128_squeeze (hstate, state->hash_evaleqs, 32);
  shake128_clear (hstate);

  _free (buf, buflen);
}

static void
_lnp_hash_statement_l2 (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int *n = params->n;
  unsigned int nelems;
  size_t buflen;
  uint8_t *buf;
  unsigned int i;
  shake128_state_t hstate;
  coder_state_t cstate;
  size_t len;

  nelems = 0;
  for (i = 0; i < Z; i++)
    nelems = MAX (nelems, n[i]);

  buflen = CEIL ((nelems * m1 + nelems * l + nelems) * d * log2q, 8) + 1;
  buf = _alloc (buflen);

  shake128_init (hstate);

  for (i = 0; i < Z; i++)
    {
      coder_enc_begin (cstate, buf);
      if (state->Es != NULL && !(state->Es[i] == NULL))
        {
          polymat_redp (state->Es[i], state->Es[i]);
          coder_enc_urandom4 (cstate, state->Es[i], q, log2q);
        }
      if (state->Em != NULL && !(state->Em[i] == NULL))
        {
          polymat_redp (state->Em[i], state->Em[i]);
          coder_enc_urandom4 (cstate, state->Em[i], q, log2q);
        }
      if (state->v != NULL && !(state->v[i] == NULL))
        {
          polyvec_redp (state->v[i], state->v[i]);
          coder_enc_urandom3 (cstate, state->v[i], q, log2q);
        }
      coder_enc_end (cstate);
      len = coder_get_offset (cstate);
      ASSERT_ERR (len % 8 == 0);
      ASSERT_ERR (len / 8 <= buflen);
      len >>= 3; /* nbits to nbytes */

      shake128_absorb (hstate, buf, len);
    }

  shake128_squeeze (hstate, state->hash_l2, 32);
  shake128_clear (hstate);

  _free (buf, buflen);
}

static void
_lnp_hash_statement_bin (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nbin = params->nbin;
  const size_t buflen
      = CEIL ((m1 * nbin + l * nbin + nbin) * d * log2q, 8) + 1;
  shake128_state_t hstate;
  coder_state_t cstate;
  size_t len;
  uint8_t *buf;

  buf = _alloc (buflen);

  coder_enc_begin (cstate, buf);
  if (state->Ps != NULL)
    {
      polymat_redp (state->Ps, state->Ps);
      coder_enc_urandom4 (cstate, state->Ps, q, log2q);
    }
  if (state->Pm != NULL)
    {
      polymat_redp (state->Pm, state->Pm);
      coder_enc_urandom4 (cstate, state->Pm, q, log2q);
    }
  if (state->f != NULL)
    {
      polyvec_redp (state->f, state->f);
      coder_enc_urandom3 (cstate, state->f, q, log2q);
    }
  coder_enc_end (cstate);
  len = coder_get_offset (cstate);
  ASSERT_ERR (len % 8 == 0);
  ASSERT_ERR (len / 8 <= buflen);
  len >>= 3; /* nbits to nbytes */

  shake128_init (hstate);
  shake128_absorb (hstate, buf, len);
  shake128_squeeze (hstate, state->hash_bin, 32);
  shake128_clear (hstate);

  _free (buf, buflen);
}

static void
_lnp_hash_statement_arp (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nprime = params->nprime;
  const size_t buflen
      = CEIL ((m1 * nprime + l * nprime + nprime) * d * log2q, 8) + 1;
  shake128_state_t hstate;
  coder_state_t cstate;
  size_t len;
  uint8_t *buf;

  buf = _alloc (buflen);

  coder_enc_begin (cstate, buf);
  if (state->Ds != NULL)
    {
      polymat_redp (state->Ds, state->Ds);
      coder_enc_urandom4 (cstate, state->Ds, q, log2q);
    }
  if (state->Dm != NULL)
    {
      polymat_redp (state->Dm, state->Dm);
      coder_enc_urandom4 (cstate, state->Dm, q, log2q);
    }
  if (state->u != NULL)
    {
      polyvec_redp (state->u, state->u);
      coder_enc_urandom3 (cstate, state->u, q, log2q);
    }
  coder_enc_end (cstate);
  len = coder_get_offset (cstate);
  ASSERT_ERR (len % 8 == 0);
  ASSERT_ERR (len / 8 <= buflen);
  len >>= 3; /* nbits to nbytes */

  shake128_init (hstate);
  shake128_absorb (hstate, buf, len);
  shake128_squeeze (hstate, state->hash_arp, 32);
  shake128_clear (hstate);

  _free (buf, buflen);
}

static void
_lnp_init (_lnp_state_t state, const uint8_t ppseed[32],
           const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int Z = params->Z;
  const unsigned int d = Rq->d;
  const unsigned int kmsis = tbox->kmsis;
  const unsigned int lambda = params->quad_eval->lambda;

  memset (state, 0, sizeof (_lnp_state_t));

  state->params = params;
  memcpy (state->ppseed, ppseed, 32);

  polymat_alloc (state->A1, Rq, kmsis, tbox->m1);
  polymat_alloc (state->A2prime, Rq, kmsis, tbox->m2 - kmsis);
  polymat_alloc (state->Bprime, Rq, tbox->l + tbox->lext, tbox->m2 - kmsis);
  polyvec_alloc (state->tA1, Rq, kmsis);
  polyvec_alloc (state->tA2, Rq, kmsis);
  polyvec_alloc (state->tB, Rq, tbox->l + tbox->lext);
  polyvec_alloc (state->h, Rq, quade->lambda / 2);
  polyvec_alloc (state->hint, Rq, kmsis);
  polyvec_alloc (state->z1, Rq, tbox->m1);
  polyvec_alloc (state->z21, Rq, tbox->m2 - kmsis);
  polyvec_alloc (state->z3, Rq, 256 / d);
  polyvec_alloc (state->z4, Rq, 256 / d);
  poly_alloc (state->c, Rq);

  state->R2 = _alloc (sizeof (spolymat_ptr) * QUADEQ_INPUT_OFF);
  state->r1 = _alloc (sizeof (spolyvec_ptr) * QUADEQ_INPUT_OFF);
  state->r0 = _alloc (sizeof (poly_ptr) * QUADEQ_INPUT_OFF);
  memset (state->R2, 0, sizeof (spolymat_ptr) * QUADEQ_INPUT_OFF);
  memset (state->r1, 0, sizeof (spolyvec_ptr) * QUADEQ_INPUT_OFF);
  memset (state->r0, 0, sizeof (poly_ptr) * QUADEQ_INPUT_OFF);

  state->R2prime = _alloc (sizeof (spolymat_ptr) * EVALEQ_INPUT_OFF);
  state->r1prime = _alloc (sizeof (spolyvec_ptr) * EVALEQ_INPUT_OFF);
  state->r0prime = _alloc (sizeof (poly_ptr) * EVALEQ_INPUT_OFF);
  memset (state->R2prime, 0, sizeof (spolymat_ptr) * EVALEQ_INPUT_OFF);
  memset (state->r1prime, 0, sizeof (spolyvec_ptr) * EVALEQ_INPUT_OFF);
  memset (state->r0prime, 0, sizeof (poly_ptr) * EVALEQ_INPUT_OFF);

  if (Z > 0)
    {
      state->Es = _alloc (sizeof (polymat_ptr) * Z);
      state->Em = _alloc (sizeof (polymat_ptr) * Z);
      state->v = _alloc (sizeof (polyvec_ptr) * Z);
      memset (state->Es, 0, sizeof (polymat_ptr) * Z);
      memset (state->Em, 0, sizeof (polymat_ptr) * Z);
      memset (state->v, 0, sizeof (polyvec_ptr) * Z);
    }

  /* expand public randomness */
  abdlop_keygen (state->A1, state->A2prime, state->Bprime, ppseed, tbox);
}

void
lnp_prover_init (lnp_prover_state_t state, const uint8_t ppseed[32],
                 const lnp_tbox_params_t params)
{
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  _lnp_init (state->state, ppseed, params);

  polyvec_alloc (state->s1, Rq, tbox->m1);
  polyvec_alloc (state->s2, Rq, tbox->m2);
  polyvec_alloc (state->m, Rq, tbox->l + tbox->lext);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_init (lnp_verifier_state_t state, const uint8_t ppseed[32],
                   const lnp_tbox_params_t params)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_init (state->state, ppseed, params);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

/*
 * compute upsilon and append to s1
 * (upsilon is part of the witness, but can only be set
 *  after the l2-norm statement has been set, since it depends on v).
 */
static void
_lnp_prover_set_witness_upsilon (lnp_prover_state_t state)
{
  lnp_tbox_params_srcptr params = state->state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int *n = params->n;
  INT_T (l2sqr, 2 * Rq->q->nlimbs);
  INT_T (l2sqr2, 2 * Rq->q->nlimbs);
  polyvec_t upsilon, subv, s1_, m_, tmp;
  polyvec_ptr v_;
  polymat_ptr Es_, Em_;
  unsigned int i, nelems;

  ASSERT_ERR (state->witness_set == 1);
  ASSERT_ERR (state->state->statement_l2_set == 1);
  ASSERT_ERR (Z > 0);

  nelems = 0;
  for (i = 0; i < Z; i++)
    nelems = MAX (nelems, n[i]);

  polyvec_alloc (tmp, Rq, nelems);

  polyvec_get_subvec (m_, state->m, 0, l, 1);
  polyvec_get_subvec (s1_, state->s1, 0, m1, 1);
  polyvec_get_subvec (upsilon, state->s1, m1, Z, 1);
  for (i = 0; i < Z; i++)
    {
      Es_ = state->state->Es != NULL ? state->state->Es[i] : NULL;
      Em_ = state->state->Em != NULL ? state->state->Em[i] : NULL;
      v_ = state->state->v != NULL ? state->state->v[i] : NULL;

      polyvec_get_subvec (subv, tmp, 0, params->n[i], 1);
      int_set (l2sqr2, params->l2Bsqr[i]);

      if (v_ != NULL)
        polyvec_set (subv, v_);
      else
        polyvec_set_zero (subv);

      if (Em_ != NULL)
        polyvec_addmul (subv, Em_, m_, 0);

      if (Es_ != NULL)
        polyvec_addmul (subv, Es_, s1_, 0);

      polyvec_l2sqr (l2sqr, subv);
      int_sub (l2sqr2, l2sqr2, l2sqr);

      int_binexp (polyvec_get_elem (upsilon, i), NULL, l2sqr2);
    }

  polyvec_free (tmp);
}

void
lnp_prover_set_witness (lnp_prover_state_t state, polyvec_t s1, polyvec_t m)
{
  lnp_tbox_params_srcptr params = state->state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  polyvec_t subv;
  polyring_srcptr Rprime;
  const unsigned int l = tbox->l;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  unsigned int k;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  polyvec_get_subvec (subv, state->s1, 0, m1, 1);
  Rprime = polyvec_get_ring (s1);

  k = Rprime->d / d;
  ASSERT_ERR (k * d == Rprime->d);
  ASSERT_ERR (polyvec_get_nelems (s1) == m1 / k);

  if (k == 1)
    polyvec_set (subv, s1);
  else
    polyvec_toisoring (subv, s1);

  if (l > 0)
    {
      polyvec_get_subvec (subv, state->m, 0, l, 1);
      Rprime = polyvec_get_ring (m);
      k = Rprime->d / d;
      ASSERT_ERR (k * d == Rprime->d);
      ASSERT_ERR (polyvec_get_nelems (m) == l / k);
      if (k == 1)
        polyvec_set (subv, m);
      else
        polyvec_toisoring (subv, m);
    }

  state->witness_set = 1;

  if (Z > 0 && state->state->statement_l2_set)
    _lnp_prover_set_witness_upsilon (state);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

static void
_lnp_quadeqs_clear (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  const unsigned int lambda = params->quad_eval->lambda;
  const unsigned int N = state->N;
  unsigned int i;

  for (i = QUADEQ_INPUT_OFF; i < QUADEQ_INPUT_OFF + N; i++)
    {
      if (state->R2[i] != NULL)
        {
          spolymat_free (state->R2[i]);
          state->R2[i] = NULL;
        }
      if (state->r1[i] != NULL)
        {
          spolyvec_free (state->r1[i]);
          state->r1[i] = NULL;
        }
      if (state->r0[i] != NULL)
        {
          poly_free (state->r0[i]);
          state->r0[i] = NULL;
        }
    }
}

static void
_lnp_evaleqs_clear (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int M = state->M;
  unsigned int i;

  for (i = EVALEQ_INPUT_OFF; i < EVALEQ_INPUT_OFF + M; i++)
    {
      if (state->R2prime[i] != NULL)
        {
          spolymat_free (state->R2prime[i]);
          state->R2prime[i] = NULL;
        }
      if (state->r1prime[i] != NULL)
        {
          spolyvec_free (state->r1prime[i]);
          state->r1prime[i] = NULL;
        }
      if (state->r0prime[i] != NULL)
        {
          poly_free (state->r0prime[i]);
          state->r0prime[i] = NULL;
        }
    }
}

static void
_lnp_set_statement_quadeqs (_lnp_state_t state, spolymat_ptr R2[],
                            spolyvec_ptr r1[], poly_ptr r0[], unsigned int N)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr quad = quade->quad_many;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int lambda = params->quad_eval->lambda;
  const unsigned int nN = QUADEQ_INPUT_OFF + N;
  const unsigned int oN = QUADEQ_INPUT_OFF + state->N;
  unsigned int i;

  ASSERT_ERR (state->N > 0);

  _lnp_quadeqs_clear (state);

  if (N != state->N)
    {
      state->R2 = _realloc (state->R2, sizeof (spolymat_ptr) * oN,
                            sizeof (spolymat_ptr) * nN);
      state->r1 = _realloc (state->r1, sizeof (spolyvec_ptr) * oN,
                            sizeof (spolyvec_ptr) * nN);
      state->r0 = _realloc (state->r0, sizeof (poly_ptr) * oN,
                            sizeof (poly_ptr) * nN);
    }

  memset (state->R2 + QUADEQ_INPUT_OFF, 0, sizeof (spolymat_ptr) * N);
  memset (state->r1 + QUADEQ_INPUT_OFF, 0, sizeof (spolyvec_ptr) * N);
  memset (state->r0 + QUADEQ_INPUT_OFF, 0, sizeof (poly_ptr) * N);

  for (i = QUADEQ_INPUT_OFF; i < QUADEQ_INPUT_OFF + N; i++)
    {
      ASSERT_ERR (R2 == NULL || R2[i] == NULL
                  || spolymat_get_ring (R2[i]) == Rq);
      ASSERT_ERR (R2 == NULL || R2[i] == NULL
                  || spolymat_get_nrows (R2[i]) == 2 * (m1 + l));
      ASSERT_ERR (R2 == NULL || R2[i] == NULL
                  || spolymat_get_ncols (R2[i]) == 2 * (m1 + l));
      ASSERT_ERR (r1 == NULL || r1[i] == NULL
                  || spolyvec_get_ring (r1[i]) == Rq);
      ASSERT_ERR (r1 == NULL || r1[i] == NULL
                  || r1[i]->nelems_max == 2 * (m1 + l));
      ASSERT_ERR (r0 == NULL || r0[i] == NULL || poly_get_ring (r0[i]) == Rq);

      if (R2[i] != NULL)
        {
          state->R2[i] = _alloc (sizeof (spolymat_t));
          spolymat_alloc (state->R2[i], Rq, 2 * (tbox->m1 + quad->l),
                          2 * (tbox->m1 + quad->l), R2[i]->nelems);
          _scatter_smat (state->R2[i], R2[i], m1, Z, l);
        }
      if (r1[i] != NULL)
        {
          state->r1[i] = _alloc (sizeof (spolyvec_t));
          spolyvec_alloc (state->r1[i], Rq, 2 * (tbox->m1 + quad->l),
                          2 * (tbox->m1 + quad->l));
          _scatter_vec (state->r1[i], r1[i], m1, Z);
        }
      if (r0[i] != NULL)
        {
          state->r0[i] = _alloc (sizeof (poly_t));
          poly_alloc (state->r0[i], Rq);
          poly_set (state->r0[i], r0[i]);
        }
    }

  state->N = N;
  _lnp_hash_statement_quadeqs (state);
}

void
lnp_prover_set_statement_quadeqs (lnp_prover_state_t state, spolymat_ptr R2[],
                                  spolyvec_ptr r1[], poly_ptr r0[],
                                  unsigned int N)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_quadeqs (state->state, R2, r1, r0, N);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_set_statement_quadeqs (lnp_verifier_state_t state,
                                    spolymat_ptr R2[], spolyvec_ptr r1[],
                                    poly_ptr r0[], unsigned int N)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_quadeqs (state->state, R2, r1, r0, N);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

static void
_lnp_set_statement_evaleqs (_lnp_state_t state, spolymat_ptr R2prime[],
                            spolyvec_ptr r1prime[], poly_ptr r0prime[],
                            unsigned int M)
{
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  lnp_quad_eval_params_srcptr quade = params->quad_eval;
  abdlop_params_srcptr quad = quade->quad_many;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nM = EVALEQ_INPUT_OFF + M;
  const unsigned int oM = EVALEQ_INPUT_OFF + state->M;
  unsigned int i;

  ASSERT_ERR (state->M > 0);

  _lnp_evaleqs_clear (state);

  if (M != state->M)
    {
      state->R2prime = _realloc (state->R2prime, sizeof (spolymat_ptr) * oM,
                                 sizeof (spolymat_ptr) * nM);
      state->r1prime = _realloc (state->r1prime, sizeof (spolyvec_ptr) * oM,
                                 sizeof (spolyvec_ptr) * nM);
      state->r0prime = _realloc (state->r0prime, sizeof (poly_ptr) * oM,
                                 sizeof (poly_ptr) * nM);
    }

  memset (state->R2prime + EVALEQ_INPUT_OFF, 0, sizeof (spolymat_ptr) * M);
  memset (state->r1prime + EVALEQ_INPUT_OFF, 0, sizeof (spolyvec_ptr) * M);
  memset (state->r0prime + EVALEQ_INPUT_OFF, 0, sizeof (poly_ptr) * M);

  for (i = EVALEQ_INPUT_OFF; i < EVALEQ_INPUT_OFF + M; i++)
    {
      ASSERT_ERR (R2prime == NULL || R2prime[i] == NULL
                  || spolymat_get_ring (R2prime[i]) == Rq);
      ASSERT_ERR (R2prime == NULL || R2prime[i] == NULL
                  || spolymat_get_nrows (R2prime[i]) == 2 * (m1 + l));
      ASSERT_ERR (R2prime == NULL || R2prime[i] == NULL
                  || spolymat_get_ncols (R2prime[i]) == 2 * (m1 + l));
      ASSERT_ERR (r1prime == NULL || r1prime[i] == NULL
                  || spolyvec_get_ring (r1prime[i]) == Rq);
      ASSERT_ERR (r1prime == NULL || r1prime[i] == NULL
                  || r1prime[i]->nelems_max == 2 * (m1 + l));
      ASSERT_ERR (r0prime == NULL || r0prime[i] == NULL
                  || poly_get_ring (r0prime[i]) == Rq);

      if (R2prime[i] != NULL)
        {
          state->R2prime[i] = _alloc (sizeof (spolymat_t));
          spolymat_alloc (state->R2prime[i], Rq, 2 * (tbox->m1 + quad->l),
                          2 * (tbox->m1 + quad->l), R2prime[i]->nelems);
          _scatter_smat (state->R2prime[i], R2prime[i], m1, Z, l);
        }
      if (r1prime[i] != NULL)
        {
          state->r1prime[i] = _alloc (sizeof (spolyvec_t));
          spolyvec_alloc (state->r1prime[i], Rq, 2 * (tbox->m1 + quad->l),
                          2 * (tbox->m1 + quad->l));
          _scatter_vec (state->r1prime[i], r1prime[i], m1, Z);
        }
      if (r0prime[i] != NULL)
        {
          state->r0prime[i] = _alloc (sizeof (poly_t));
          poly_alloc (state->r0prime[i], Rq);
          poly_set (state->r0prime[i], r0prime[i]);
        }
    }

  state->M = M;
  _lnp_hash_statement_evaleqs (state);
}

void
lnp_prover_set_statement_evaleqs (lnp_prover_state_t state,
                                  spolymat_ptr R2prime[],
                                  spolyvec_ptr r1prime[], poly_ptr r0prime[],
                                  unsigned int M)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_evaleqs (state->state, R2prime, r1prime, r0prime, M);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_set_statement_evaleqs (lnp_verifier_state_t state,
                                    spolymat_ptr R2prime[],
                                    spolyvec_ptr r1prime[], poly_ptr r0prime[],
                                    unsigned int M)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_evaleqs (state->state, R2prime, r1prime, r0prime, M);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

static void
_lnp_set_statement_l2 (_lnp_state_t state, polymat_ptr Es[], polymat_ptr Em[],
                       polyvec_ptr v[])
{
  lnp_tbox_params_srcptr params = state->params;
  const unsigned int Z = params->Z;
#if ASSERT == ASSERT_ENABLED
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int *n = params->n;
#endif
  unsigned int i;

  ASSERT_ERR (Z > 0);

  if (state->statement_l2_set == 1)
    state->statement_l2_set = 0;

  for (i = 0; i < Z; i++)
    {
      // XXX assert at least 1 of Esi, Emi != NULL
      ASSERT_ERR (Es == NULL || Es[i] == NULL
                  || polymat_get_ring (Es[i]) == Rq);
      ASSERT_ERR (Es == NULL || Es[i] == NULL
                  || polymat_get_nrows (Es[i]) == n[i]);
      ASSERT_ERR (Es == NULL || Es[i] == NULL
                  || polymat_get_ncols (Es[i]) == m1);
      ASSERT_ERR (Em == NULL || Em[i] == NULL
                  || polymat_get_ring (Em[i]) == Rq);
      ASSERT_ERR (Em == NULL || Em[i] == NULL
                  || polymat_get_nrows (Em[i]) == n[i]);
      ASSERT_ERR (Em == NULL || Em[i] == NULL
                  || polymat_get_ncols (Em[i]) == l);
      ASSERT_ERR (v == NULL || v[i] == NULL || polyvec_get_ring (v[i]) == Rq);
      ASSERT_ERR (v == NULL || v[i] == NULL
                  || polyvec_get_nelems (v[i]) == n[i]);

      state->Es[i] = (Es == NULL ? NULL : Es[i]);
      state->Em[i] = (Em == NULL ? NULL : Em[i]);
      state->v[i] = (v == NULL ? NULL : v[i]);
    }

  _lnp_hash_statement_l2 (state);
  state->statement_l2_set = 1;
}

void
lnp_prover_set_statement_l2 (lnp_prover_state_t state, polymat_ptr Es[],
                             polymat_ptr Em[], polyvec_ptr v[])
{
  lnp_tbox_params_srcptr params = state->state->params;
  const unsigned int Z = params->Z;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_l2 (state->state, Es, Em, v);

  if (Z > 0 && state->witness_set)
    _lnp_prover_set_witness_upsilon (state);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_set_statement_l2 (lnp_verifier_state_t state, polymat_ptr Es[],
                               polymat_ptr Em[], polyvec_ptr v[])
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_l2 (state->state, Es, Em, v);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

static void
_lnp_set_statement_bin (_lnp_state_t state, polymat_t Ps, polymat_t Pm,
                        polyvec_t f)
{
#if ASSERT == ASSERT_ENABLED
  lnp_tbox_params_srcptr params = state->params;
  const unsigned int Z = params->Z;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nbin = params->nbin;
#endif

  ASSERT_ERR (nbin > 0);
  ASSERT_ERR (Pm != NULL || Ps != NULL);
  ASSERT_ERR (Ps == NULL || polymat_get_ring (Ps) == Rq);
  ASSERT_ERR (Ps == NULL || polymat_get_nrows (Ps) == nbin);
  ASSERT_ERR (Ps == NULL || polymat_get_ncols (Ps) == m1);
  ASSERT_ERR (Pm == NULL || polymat_get_ring (Pm) == Rq);
  ASSERT_ERR (Pm == NULL || polymat_get_nrows (Pm) == nbin);
  ASSERT_ERR (Pm == NULL || polymat_get_ncols (Pm) == l);
  ASSERT_ERR (f == NULL || polyvec_get_ring (f) == Rq);
  ASSERT_ERR (f == NULL || polyvec_get_nelems (f) == nbin);

  if (state->statement_bin_set == 1)
    state->statement_bin_set = 0;

  state->Ps = Ps;
  state->Pm = Pm;
  state->f = f;

  _lnp_hash_statement_bin (state);
  state->statement_bin_set = 1;
}

void
lnp_prover_set_statement_bin (lnp_prover_state_t state, polymat_t Ps,
                              polymat_t Pm, polyvec_t f)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_bin (state->state, Ps, Pm, f);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_set_statement_bin (lnp_verifier_state_t state, polymat_t Ps,
                                polymat_t Pm, polyvec_t f)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_bin (state->state, Ps, Pm, f);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

static void
_lnp_set_statement_arp (_lnp_state_t state, polymat_t Ds, polymat_t Dm,
                        polyvec_t u)
{
#if ASSERT == ASSERT_ENABLED
  lnp_tbox_params_srcptr params = state->params;
  const unsigned int Z = params->Z;
  abdlop_params_srcptr tbox = params->tbox;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nprime = params->nprime;
#endif

  ASSERT_ERR (nprime > 0);
  ASSERT_ERR (Ds != NULL || Dm != NULL);
  ASSERT_ERR (Ds == NULL || polymat_get_ring (Ds) == Rq);
  ASSERT_ERR (Ds == NULL || polymat_get_nrows (Ds) == nprime);
  ASSERT_ERR (Ds == NULL || polymat_get_ncols (Ds) == m1);
  ASSERT_ERR (Dm == NULL || polymat_get_ring (Dm) == Rq);
  ASSERT_ERR (Dm == NULL || polymat_get_nrows (Dm) == nprime);
  ASSERT_ERR (Dm == NULL || polymat_get_ncols (Dm) == l);
  ASSERT_ERR (u == NULL || polyvec_get_ring (u) == Rq);
  ASSERT_ERR (u == NULL || polyvec_get_nelems (u) == nprime);

  state->Ds = Ds;
  state->Dm = Dm;
  state->u = u;

  _lnp_hash_statement_arp (state);
  state->statement_arp_set = 1;
}

void
lnp_prover_set_statement_arp (lnp_prover_state_t state, polymat_t Ds,
                              polymat_t Dm, polyvec_t u)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_arp (state->state, Ds, Dm, u);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_set_statement_arp (lnp_verifier_state_t state, polymat_t Ds,
                                polymat_t Dm, polyvec_t u)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_set_statement_arp (state->state, Ds, Dm, u);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

static void
_lnp_clear (_lnp_state_t state)
{
  lnp_tbox_params_srcptr params = state->params;
  polyring_srcptr Rq = params->tbox->ring;
  const unsigned int d = Rq->d;
  const unsigned int lambda = params->quad_eval->lambda;
  const unsigned int Z = params->Z;
  const unsigned int N = state->N;
  const unsigned int M = state->M;

  polymat_free (state->A1);
  polymat_free (state->A2prime);
  polymat_free (state->Bprime);
  polyvec_free (state->tA1);
  polyvec_free (state->tA2);
  polyvec_free (state->tB);
  polyvec_free (state->h);
  polyvec_free (state->hint);
  polyvec_free (state->z1);
  polyvec_free (state->z21);
  polyvec_free (state->z3);
  polyvec_free (state->z4);
  poly_free (state->c);

  _lnp_quadeqs_clear (state);
  _lnp_evaleqs_clear (state);

  _free (state->R2, sizeof (spolymat_ptr) * (QUADEQ_INPUT_OFF + N));
  _free (state->r1, sizeof (spolyvec_ptr) * (QUADEQ_INPUT_OFF + N));
  _free (state->r0, sizeof (poly_ptr) * (QUADEQ_INPUT_OFF + N));

  _free (state->R2prime, sizeof (spolymat_ptr) * (EVALEQ_INPUT_OFF + M));
  _free (state->r1prime, sizeof (spolyvec_ptr) * (EVALEQ_INPUT_OFF + M));
  _free (state->r0prime, sizeof (poly_ptr) * (EVALEQ_INPUT_OFF + M));

  _free (state->Es, sizeof (polymat_ptr) * Z);
  _free (state->Em, sizeof (spolymat_ptr) * Z);
  _free (state->v, sizeof (polyvec_ptr) * Z);
}

static void
_lnp_tbox_encproof (uint8_t *out, size_t *len, polyvec_t tA1, polyvec_t tB,
                    polyvec_t h, poly_t c, polyvec_t z1, polyvec_t z21,
                    polyvec_t hint, polyvec_t z3, polyvec_t z4,
                    const lnp_tbox_params_t params)
{
  polyring_srcptr Rq = params->tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int d = Rq->d;
  INTVEC_T (coeffs, d, q->nlimbs);
  const unsigned int log2q = Rq->log2q;
  const unsigned int log2omega = params->quad_eval->quad_many->log2omega;
  const unsigned int D = params->quad_eval->quad_many->dcompress->D;
  const int64_t omega = params->quad_eval->quad_many->omega;
  coder_state_t cstate;
  size_t prooflen;
  INT_T (mod, q->nlimbs);

  coder_enc_begin (cstate, out);

  /* full-sized elements (log2q bits) */
  polyvec_fromcrt (tB);
  polyvec_mod (tB, tB);
  polyvec_redp (tB, tB);
  coder_enc_urandom3 (cstate, tB, q, log2q);

  polyvec_fromcrt (h);
  polyvec_mod (h, h);
  polyvec_redp (h, h);
  coder_enc_urandom3 (cstate, h, q, log2q);

  /* compressed elements (log2q - D bits) */
  int_set_one (mod);
  int_lshift (mod, mod, log2q - D);
  polyvec_fromcrt (tA1);
  polyvec_mod (tA1, tA1);
  polyvec_redp (tA1, tA1);
  coder_enc_urandom3 (cstate, tA1, mod, log2q - D);

  /* challenge (log2omega bits) */
  int_set_i64 (mod, 2 * omega + 1);
  poly_fromcrt (c);
  intvec_set (coeffs, poly_get_coeffvec (c));
  intvec_redp (coeffs, coeffs, mod);
  coder_enc_urandom (cstate, coeffs, mod, log2omega);

  /* hints */
  coder_enc_ghint3 (cstate, hint);

  /* gaussian elements */
  polyvec_fromcrt (z1);
  polyvec_fromcrt (z21);
  polyvec_fromcrt (z3);
  polyvec_fromcrt (z4);
  polyvec_mod (z1, z1);
  polyvec_mod (z21, z21);
  polyvec_mod (z3, z3);
  polyvec_mod (z4, z4);
  polyvec_redc (z1, z1);
  polyvec_redc (z21, z21);
  polyvec_redc (z3, z3);
  polyvec_redc (z4, z4);
  coder_enc_grandom3 (cstate, z1, params->quad_eval->quad_many->log2stdev1);
  coder_enc_grandom3 (cstate, z21, params->quad_eval->quad_many->log2stdev2);
  coder_enc_grandom3 (cstate, z3, params->log2stdev3);
  coder_enc_grandom3 (cstate, z4, params->log2stdev4);

  coder_enc_end (cstate);

  prooflen = coder_get_offset (cstate);
  ASSERT_ERR (prooflen % 8 == 0);
  (prooflen) >>= 3; /* nbits to nbytes */

  if (len != NULL)
    *len = prooflen;
}

static int
_lnp_tbox_decproof (size_t *len, const uint8_t *in, polyvec_t tA1,
                    polyvec_t tB, polyvec_t h, poly_t c, polyvec_t z1,
                    polyvec_t z21, polyvec_t hint, polyvec_t z3, polyvec_t z4,
                    const lnp_tbox_params_t params)
{
  polyring_srcptr Rq = params->tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int log2omega = params->quad_eval->quad_many->log2omega;
  const unsigned int D = params->quad_eval->quad_many->dcompress->D;
  const int64_t omega = params->quad_eval->quad_many->omega;
  coder_state_t cstate;
  intvec_ptr coeffs;
  size_t prooflen = 0;
  int succ = 0, rc;
  INT_T (mod, q->nlimbs);

  coder_dec_begin (cstate, in);

  /* full-sized elements (log2q bits) */
  rc = coder_dec_urandom3 (cstate, tB, q, log2q);
  if (rc != 0)
    goto ret;

  rc = coder_dec_urandom3 (cstate, h, q, log2q);
  if (rc != 0)
    goto ret;

  /* compressed elements (log2q - D bits) */
  int_set_one (mod);
  int_lshift (mod, mod, log2q - D);
  rc = coder_dec_urandom3 (cstate, tA1, mod, log2q - D);
  if (rc != 0)
    goto ret;

  /* challenge (log2omega bits) */
  int_set_i64 (mod, 2 * omega + 1);
  rc = coder_dec_urandom2 (cstate, c, mod, log2omega);
  if (rc != 0)
    goto ret;

  coeffs = poly_get_coeffvec (c);
  intvec_redc (coeffs, coeffs, mod);

  /* hints */
  coder_dec_ghint3 (cstate, hint);

  /* gaussian elements */
  coder_dec_grandom3 (cstate, z1, params->quad_eval->quad_many->log2stdev1);
  coder_dec_grandom3 (cstate, z21, params->quad_eval->quad_many->log2stdev2);
  coder_dec_grandom3 (cstate, z3, params->log2stdev3);
  coder_dec_grandom3 (cstate, z4, params->log2stdev4);

  rc = coder_dec_end (cstate);
  if (rc != 1)
    goto ret;

  prooflen = coder_get_offset (cstate);
  ASSERT_ERR (prooflen % 8 == 0);
  (prooflen) >>= 3; /* nbits to nbytes */

  succ = 1;
ret:
  if (len != NULL)
    *len = prooflen;

  return succ;
}

/*
 * Hash the seed of the public parameters and the substatements to
 * obtain a hash of the public parameters and the statement.
 */
static void
_lnp_hash_pp_and_statement (_lnp_state_t state, uint8_t hash[32])
{
  lnp_tbox_params_srcptr params = state->params;
  const unsigned int N = state->N;
  const unsigned int M = state->M;
  const unsigned int Z = params->Z;
  const unsigned int nbin = params->nbin;
  const unsigned int nprime = params->nprime;
  shake128_state_t hstate;

  ASSERT_ERR (Z == 0 || state->statement_l2_set == 1);
  ASSERT_ERR (nbin == 0 || state->statement_bin_set == 1);
  ASSERT_ERR (nprime == 0 || state->statement_arp_set == 1);

  shake128_init (hstate);
  shake128_absorb (hstate, state->ppseed, 32);
  if (N > 0)
    shake128_absorb (hstate, state->hash_quadeqs, 32);
  if (M > 0)
    shake128_absorb (hstate, state->hash_evaleqs, 32);
  if (Z > 0)
    shake128_absorb (hstate, state->hash_l2, 32);
  if (nbin > 0)
    shake128_absorb (hstate, state->hash_bin, 32);
  if (nprime > 0)
    shake128_absorb (hstate, state->hash_arp, 32);
  shake128_squeeze (hstate, hash, 32);
  shake128_clear (hstate);
}

#if ASSERT == ASSERT_ENABLED
static void
_verify_statement (lnp_prover_state_t state_)
{
  _lnp_state_srcptr state = state_->state;
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  const unsigned int lambda = params->quad_eval->lambda;
  const unsigned int N = state->N;
  const unsigned int M = state->M;
  const unsigned int Z = params->Z;
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int l = tbox->l;
  const unsigned int nex = params->nex;
  const unsigned int nprime = params->nprime;
  const unsigned int nbin = params->nbin;
  const unsigned int *n = params->n;
  polyring_srcptr Rq = params->tbox->ring;
  polyvec_t s, s1_, m_, sub, tmp, t;
  int_srcptr q = Rq->q;
  poly_t tmp2, zero;
  INT_T (l2sqr, 2 * q->nlimbs);
  INT_T (l2sqr_, q->nlimbs);
  const unsigned int d = Rq->d;
  unsigned int i, j;
  int64_t coeff;
  poly_ptr poly;

  poly_alloc (tmp2, Rq);
  poly_alloc (zero, Rq);
  polyvec_alloc (s, Rq, 2 * (m1 + l));
  polyvec_alloc (tmp, Rq, MAX (nex, MAX (nprime, 2 * (m1 + l))));
  poly_set_zero (zero);

  /* s = (<s1>,<m>) */

  polyvec_get_subvec (s1_, state_->s1, 0, m1, 1);
  polyvec_get_subvec (sub, s, 0, m1, 2);
  polyvec_set (sub, s1_);
  polyvec_get_subvec (sub, s, 1, m1, 2);
  polyvec_auto (sub, s1_);
  if (l > 0)
    {
      polyvec_get_subvec (m_, state_->m, 0, l, 1);
      polyvec_get_subvec (sub, s, m1 * 2, l, 2);
      polyvec_set (sub, m_);
      polyvec_get_subvec (sub, s, m1 * 2 + 1, l, 2);
      polyvec_auto (sub, m_);
    }

  for (i = QUADEQ_INPUT_OFF; i < QUADEQ_INPUT_OFF + N; i++)
    {
      polyvec_get_subvec (t, tmp, 0, 2 * (m1 + l), 1);

      polyvec_dot2 (tmp2, state->r1[i], s);
      polyvec_mulsparse (t, state->R2[i], s);
      polyvec_fromcrt (t);
      poly_adddot (tmp2, s, t, 0);
      ASSERT_ERR (poly_eq (tmp2, zero) == 1);
    }
  for (i = EVALEQ_INPUT_OFF; i < EVALEQ_INPUT_OFF + M; i++)
    {
      polyvec_get_subvec (t, tmp, 0, 2 * (m1 + l), 1);

      polyvec_dot2 (tmp2, state->r1[i], s);
      polyvec_mulsparse (t, state->R2[i], s);
      polyvec_fromcrt (t);
      poly_adddot (tmp2, s, t, 0);
      for (j = 0; j < d; j++)
        ASSERT_ERR (int_eqzero (poly_get_coeff (tmp2, j)) == 1);
    }
  for (i = 0; i < Z; i++)
    {
      polyvec_get_subvec (t, tmp, 0, n[i], 1);

      if (state->v != NULL && state->v[i] != NULL)
        polyvec_set (t, state->v[i]);
      else
        polyvec_set_zero (t);
      if (state->Es != NULL && state->Es[i] != NULL)
        polyvec_addmul (t, state->Es[i], s1_, 0);
      if (state->Em != NULL && state->Em[i] != NULL)
        polyvec_addmul (t, state->Em[i], m_, 0);
      polyvec_fromcrt (t);
      polyvec_mod (t, t);
      polyvec_redc (t, t);

      polyvec_l2sqr (l2sqr, t);
      int_mod (l2sqr_, l2sqr, q);
      int_redc (l2sqr_, l2sqr_, q);
      ASSERT_ERR (int_le (l2sqr_, params->l2Bsqr[i]) == 1);
    }

  polyvec_get_subvec (t, tmp, 0, nprime, 1);

  if (state->u != NULL)
    polyvec_set (t, state->u);
  else
    polyvec_set_zero (t);
  if (state->Ds != NULL)
    polyvec_addmul (t, state->Ds, s1_, 0);
  if (state->Dm != NULL)
    polyvec_addmul (t, state->Dm, m_, 0);
  polyvec_fromcrt (t);
  polyvec_mod (t, t);
  polyvec_redc (t, t);

  polyvec_l2sqr (l2sqr, t);
  int_mod (l2sqr_, l2sqr, q);
  int_redc (l2sqr_, l2sqr_, q);
  // XXX ASSERT_ERR (l2sqr_, params->Bprime);

  polyvec_get_subvec (t, tmp, 0, nbin, 1);

  if (state->f != NULL)
    polyvec_set (t, state->f);
  else
    polyvec_set_zero (t);
  if (state->Ps != NULL)
    polyvec_addmul (t, state->Ps, s1_, 0);
  if (state->Pm != NULL)
    polyvec_addmul (t, state->Pm, m_, 0);
  polyvec_fromcrt (t);
  polyvec_mod (t, t);
  polyvec_redc (t, t);
  for (i = 0; i < polyvec_get_nelems (t); i++)
    {
      poly = polyvec_get_elem (t, i);
      for (j = 0; j < d; j++)
        {
          coeff = int_get_i64 (poly_get_coeff (poly, j));
          ASSERT_ERR (coeff == 0 || coeff == 1);
        }
    }

  polyvec_free (s);
  polyvec_free (tmp);
  poly_free (zero);
  poly_free (tmp2);
}
#endif

void
lnp_prover_prove (lnp_prover_state_t state_, uint8_t *proof, size_t *len,
                  const uint8_t seed[32])
{
  _lnp_state_ptr state = state_->state;
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  uint8_t hash[32], expseed[64];
  const uint8_t *seedproto = expseed;
  const uint8_t *seedsubproto = expseed + 32;
  shake128_state_t hstate;
  size_t prooflen;
  INT_T (lo, 1);
  INT_T (hi, 1);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  STOPWATCH_START (stopwatch_lnp_prover_prove, "lnp_prover_prove");

#if ASSERT == ASSERT_ENABLED
  _verify_statement (state_);
#endif

  /*
   * expand prover randomness to obtain seeds
   * for the protocol and the subprotocol.
   */
  shake128_init (hstate);
  shake128_absorb (hstate, seed, 32);
  shake128_squeeze (hstate, expseed, 64);
  shake128_clear (hstate);

  /* expand prover randomness */
  int_set_i64 (lo, -tbox->nu);
  int_set_i64 (hi, tbox->nu);
  polyvec_urandom_bnd (state_->s2, lo, hi, seedproto, 0);

  /* hash public parameters and statement */
  _lnp_hash_pp_and_statement (state, hash);

  /* commit */
  abdlop_commit (state->tA1, state->tA2, state->tB, state_->s1, state_->m,
                 state_->s2, state->A1, state->A2prime, state->Bprime, tbox);

  /* hash in commitment */
  abdlop_hashcomm (hash, state->tA1, state->tB, tbox);

  /* generate proof */
  lnp_tbox_prove (
      hash, state->tB, state->h, state->c, state->z1, state->z21, state->hint,
      state->z3, state->z4, state_->s1, state_->m, state_->s2, state->tA2,
      state->A1, state->A2prime, state->Bprime, state->R2, state->r1, state->N,
      state->R2prime, state->r1prime, state->r0prime, state->M, state->Es,
      state->Em, state->v, state->Ps, state->Pm, state->f, state->Ds,
      state->Dm, state->u, seedsubproto, params);

  /* encode commitment and proof */
  _lnp_tbox_encproof (proof, &prooflen, state->tA1, state->tB, state->h,
                      state->c, state->z1, state->z21, state->hint, state->z3,
                      state->z4, params);

  if (len != NULL)
    *len = prooflen;

  STOPWATCH_STOP (stopwatch_lnp_prover_prove);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

int
lnp_verifier_verify (lnp_verifier_state_t state_, const uint8_t *proof,
                     size_t *len)
{
  _lnp_state_ptr state = state_->state;
  lnp_tbox_params_srcptr params = state->params;
  abdlop_params_srcptr tbox = params->tbox;
  uint8_t hash[32];
  size_t prooflen;
  int b;

  STOPWATCH_START (stopwatch_lnp_verifier_verify, "lnp_verifier_verify");
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  /* decode commitment and proof */
  b = _lnp_tbox_decproof (&prooflen, proof, state->tA1, state->tB, state->h,
                          state->c, state->z1, state->z21, state->hint,
                          state->z3, state->z4, params);
  if (len != NULL)
    *len = prooflen;
  if (b != 1)
    goto ret;

  /* reduce inputs XXX required ? */
  // XXX hint
  polyvec_mod (state->h, state->h);
  polyvec_redc (state->h, state->h);
  poly_mod (state->c, state->c);
  poly_redc (state->c, state->c);
  polyvec_mod (state->z1, state->z1);
  polyvec_redc (state->z1, state->z1);
  polyvec_mod (state->z21, state->z21);
  polyvec_redc (state->z21, state->z21);
  polyvec_mod (state->z3, state->z3);
  polyvec_redc (state->z3, state->z3);
  polyvec_mod (state->z4, state->z4);
  polyvec_redc (state->z4, state->z4);

  /* hash public parameters and statement */
  _lnp_hash_pp_and_statement (state, hash);

  /* hash in commitment */
  abdlop_hashcomm (hash, state->tA1, state->tB, tbox);

  /* gverify proof */
  b = lnp_tbox_verify (
      hash, state->h, state->c, state->z1, state->z21, state->hint, state->z3,
      state->z4, state->tA1, state->tB, state->A1, state->A2prime,
      state->Bprime, state->R2, state->r1, state->r0, state->N, state->R2prime,
      state->r1prime, state->r0prime, state->M, state->Es, state->Em, state->v,
      state->Ps, state->Pm, state->f, state->Ds, state->Dm, state->u, params);
  if (b != 1)
    goto ret;

ret:
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
  STOPWATCH_STOP (stopwatch_lnp_verifier_verify);
  return b;
}

void
lnp_prover_clear (lnp_prover_state_t state)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  polyvec_free (state->s1);
  polyvec_free (state->s2);
  polyvec_free (state->m);
  _lnp_clear (state->state);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

void
lnp_verifier_clear (lnp_verifier_state_t state)
{
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  _lnp_clear (state->state);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}
