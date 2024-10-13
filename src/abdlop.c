#include "lazer.h"

/*
 * A1 uniform in Rq^(kmsis x m1)
 *
 * A2 = (A2prime, 1)
 * A2prime uniform in Rq^(kmsis x (m2-kmsis))
 * 1 in Rq^(kmsis x kmsis)
 *
 * B = (Bprime, 0)
 * Bprime uniform in Rq^(l x (m2-kmsis))
 * 0 in Rq^(l x kmsis)
 */

/*
 * Compute public key (A1,A2prime,Bprime) from seed.
 *
 * Expand uniformly random A1,A2prime,Bprime
 * from seed||0, seed||1, seed||2 respectively:
 * A1 uniform in Rq^(kmsis x m1),
 * A2prime uniform in Rq^(kmsis x (m2-kmsis))
 * Bprime = (B,Bext) uniform in Rq^((l+lext) x (m2-kmsis))
 *
 * Caller must allocate A1,A2prime,Bprime with dimensions and ring Rq
 * given by params.
 */
void
abdlop_keygen (polymat_t A1, polymat_t A2prime, polymat_t Bprime,
               const uint8_t seed[32], const abdlop_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  const unsigned int kmsis = params->kmsis;
  const unsigned int m2 = params->m2;
#endif
  polyring_srcptr Rq = params->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int m1 = params->m1;
  const unsigned int l = params->l;
  const unsigned int lext = params->lext;
  const unsigned int l_ = l + lext;

  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (l_ == 0 || polymat_get_ring (Bprime) == Rq);
  ASSERT_ERR (l_ == 0 || polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (l_ == 0 || polymat_get_ncols (Bprime) == m2 - kmsis);

  if (m1 > 0)
    {
      polymat_urandom (A1, q, log2q, seed, 0);
      polymat_urandom (A2prime, q, log2q, seed, 1);
    }

  if (l_ > 0)
    polymat_urandom (Bprime, q, log2q, seed, 2);
}

/*
 * Compute commitment (tA1,tB) to "small" message s1 and message m_
 * from randomness s2 and public key (A1,A2prime,Bprime).
 * Also computes tA2 which is needed to compute an opening proof.
 *
 * s1 in Rq^m1, l2(s1) <= alpha
 * m_ in Rq^l
 * s2 = (s21,s22) in uniform in [-nu,nu]^(m2-kmsis) x [-nu,nu]^kmsis
 * A1 uniform in Rq^(kmsis x m1),
 * A2prime uniform in Rq^(kmsis x (m2-kmsis))
 * Bprime = (Bprime_,Bext) uniform in Rq^((l+lext) x (m2-kmsis))
 *
 * tA = A1*s1 + A2prime * s21 + s22 in Rq^kmsis
 * tA1 = Power2Round(tA,D) in Rq^kmsis
 * tA2 = tA - 2^D * tA1 in Rq^kmsis
 *
 * m = (m_,mext)
 * tB = (tB_,tBext)
 * tB_ = B * s21 + m_ in Rq^l
 *
 * Caller must allocate tA1,tA2,tB,s1,m,s2,A1,A2prime,Bprime
 * with dimensions and ring Rq given by params.
 */
void
abdlop_commit (polyvec_t tA1, polyvec_t tA2, polyvec_t tB, polyvec_t s1,
               polyvec_t m, polyvec_t s2, polymat_t A1, polymat_t A2prime,
               polymat_t Bprime, const abdlop_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  polyring_srcptr Rq = params->ring;
  const unsigned int lext = params->lext;
#endif
  const unsigned int D = params->dcompress->D;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int m2 = params->m2;
  const unsigned int l = params->l;
  dcompress_params_srcptr dcomp_param = params->dcompress;
  polyvec_t s21, s22, tB_, m_;
  polymat_t Bprime_;

  ASSERT_ERR (polyvec_get_ring (tA1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA1) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tA2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA2) == kmsis);
  ASSERT_ERR (l + lext == 0 || polyvec_get_ring (tB) == Rq);
  ASSERT_ERR (l + lext == 0 || polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polyvec_get_ring (s1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (s1) == m1);
  ASSERT_ERR (l + lext == 0 || polyvec_get_ring (m) == Rq);
  ASSERT_ERR (l + lext == 0 || polyvec_get_nelems (m) == l + lext);
  ASSERT_ERR (polyvec_get_ring (s2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (s2) == m2);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (l + lext == 0 || polymat_get_ring (Bprime) == Rq);
  ASSERT_ERR (l + lext == 0 || polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (l + lext == 0 || polymat_get_ncols (Bprime) == m2 - kmsis);

  polyvec_get_subvec (s21, s2, 0, m2 - kmsis, 1);
  polyvec_get_subvec (s22, s2, m2 - kmsis, kmsis, 1);

  if (m1 > 0)
    {
      polyvec_set (tA2, s22);
      polyvec_addmul (tA2, A1, s1, 0);
      polyvec_addmul (tA2, A2prime, s21, 0);

      polyvec_dcompress_power2round (tA1, tA2, dcomp_param);
      polyvec_sublshift (tA2, tA1, D);

      polyvec_mod (tA1, tA1);
      polyvec_mod (tA2, tA2);
    }

  if (l > 0)
    {
      polyvec_get_subvec (m_, m, 0, l, 1);
      polyvec_get_subvec (tB_, tB, 0, l, 1);
      polymat_get_submat (Bprime_, Bprime, 0, 0, l, m2 - kmsis, 1, 1);

      polyvec_set (tB_, m_);
      polyvec_addmul (tB_, Bprime_, s21, 0);
      polyvec_mod (tB_, tB_);
    }
}

void
abdlop_enccomm (uint8_t *buf, size_t *buflen, polyvec_t tA1, polyvec_t tB,
                const abdlop_params_t params)
{
  polyring_srcptr Rq = params->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int D = params->dcompress->D;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int l = params->l;
#if ASSERT == ASSERT_ENABLED
  size_t outlen;
#endif
  coder_state_t cstate;
  polyvec_t tB_;
  const unsigned int len
      = CEIL (kmsis * d * (log2q - D) + l * d * log2q, 8) + 1;

  if (buflen != NULL)
    *buflen = len;

  if (buf == NULL)
    return;

  coder_enc_begin (cstate, buf);
  if (m1 > 0)
    {
      INT_T (mod, q->nlimbs);

      int_set_one (mod);
      int_lshift (mod, mod, log2q - D);

      polyvec_redp (tA1, tA1);
      coder_enc_urandom3 (cstate, tA1, mod, log2q - D);
    }
  if (l > 0)
    {
      polyvec_get_subvec (tB_, tB, 0, l, 1);
      polyvec_redp (tB_, tB_);
      coder_enc_urandom3 (cstate, tB_, q, log2q);
    }
  coder_enc_end (cstate);

#if ASSERT == ASSERT_ENABLED
  outlen =
#endif
      coder_get_offset (cstate);
  ASSERT_ERR (outlen % 8 == 0);
  ASSERT_ERR (outlen / 8 <= len);
}

/*
 * hash = H(hash||tA1||tB).
 * Input hash is usually the seed for the public parameters
 * used by abdlop_keygen.
 */
void
abdlop_hashcomm (uint8_t hash[32], polyvec_t tA1, polyvec_t tB,
                 const abdlop_params_t params)
{
  polyring_srcptr Rq = params->ring;
  const unsigned int log2q = Rq->log2q;
  const unsigned int d = Rq->d;
  const unsigned int D = params->dcompress->D;
  const unsigned int kmsis = params->kmsis;
  const unsigned int l = params->l;
  shake128_state_t hstate;
  const size_t outlen = CEIL (kmsis * d * (log2q - D) + l * d * log2q, 8) + 1;
  uint8_t out[outlen];

  abdlop_enccomm (out, NULL, tA1, tB, params);

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, outlen);
  shake128_squeeze (hstate, hash, 32);
  shake128_clear (hstate);
}

/*
 * Compute compressed opening proof (c,z1,z21,h) from hash of transcript,
 * "short" message s1 and randomness s2, tA2, the public key
 * parts (A1,A2prime) and a seed.
 * Also update hash of transcript by hashing c into it.
 *
 * tA2 = Rq^kmsis
 * s1 in Rq^m1, l2(s1) <= alpha
 * s2 = (s21,s22) in uniform in [-nu,nu]^(m2-kmsis) x [-nu,nu]^kmsis
 * A1 uniform in Rq^(kmsis x m1),
 * A2prime uniform in Rq^(kmsis x (m2-kmsis))

 * c in Rq from challenge space defined by params,
 * z1 in Rq^m1,
 * z21 in Rq^(m2-kmsis),
 * h hint in Rq^(kmsis),
 */
void
abdlop_prove (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
              polyvec_t h, polyvec_t tA2, polyvec_t s1, polyvec_t s2,
              polymat_t A1, polymat_t A2prime, const uint8_t seed[32],
              const abdlop_params_t params)
{
  polyring_srcptr Rq = params->ring;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int m2 = params->m2;
  const unsigned int d = polyring_get_deg (Rq);
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int nlimbs = int_get_nlimbs (q);
  dcompress_params_srcptr dcomp_param = params->dcompress;
  int_srcptr gamma = dcompress_get_gamma (dcomp_param);
  int_srcptr m = dcompress_get_m (dcomp_param);
  const unsigned int log2m = dcompress_get_log2m (dcomp_param);
  polyvec_t y21, y22, s21, s22, y1, y2, cs1, cs2, w, w1, w0;
  INTVEC_T (z1coeffs, d * m1, nlimbs);
  INTVEC_T (z2coeffs, d * m2, nlimbs);
  INTVEC_T (cs1coeffs, d * m1, nlimbs);
  INTVEC_T (cs2coeffs, d * m2, nlimbs);
  INT_T (norm, nlimbs * 2);
  shake128_state_t hstate;
  coder_state_t cstate;
  rng_state_t rngstate;
  uint32_t dom;
  /* buff for encoding of w1 */
  uint8_t out[CEIL (log2m * d * kmsis, 8) + 1];
  unsigned int outlen;
  uint8_t cseed[32]; /* challenge seed */
  uint8_t yseed[32]; /* mask seed */
  int rej;

  ASSERT_ERR (poly_get_ring (c) == Rq);
  ASSERT_ERR (polyvec_get_ring (z1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1);
  ASSERT_ERR (polyvec_get_ring (z21) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_ring (h) == Rq);
  ASSERT_ERR (polyvec_get_nelems (h) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tA2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA2) == kmsis);
  ASSERT_ERR (polyvec_get_ring (s1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (s1) == m1);
  ASSERT_ERR (polyvec_get_ring (s2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (s2) == m2);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);

  polyvec_alloc (y1, Rq, m1);
  polyvec_alloc (y2, Rq, m2);
  polyvec_alloc (cs1, Rq, m1);
  polyvec_alloc (cs2, Rq, m2);
  polyvec_alloc (w, Rq, kmsis);
  polyvec_alloc (w1, Rq, kmsis);
  polyvec_alloc (w0, Rq, kmsis);

  /* reuse y1, y2=(y21,y22) as z1, z2=(z21,z22) later */
  polyvec_set_coeffvec2 (y1, z1coeffs);
  polyvec_set_coeffvec2 (y2, z2coeffs);
  polyvec_set_coeffvec2 (cs1, cs1coeffs);
  polyvec_set_coeffvec2 (cs2, cs2coeffs);

  polyvec_get_subvec (s21, s2, 0, m2 - kmsis, 1);
  polyvec_get_subvec (s22, s2, m2 - kmsis, kmsis, 1);

  polyvec_get_subvec (y21, y2, 0, m2 - kmsis, 1);
  polyvec_get_subvec (y22, y2, m2 - kmsis, kmsis, 1);

  rng_init (rngstate, seed, 0);
  rng_urandom (rngstate, yseed, 32);

  dom = 0;
  while (1)
    {
      polyvec_grandom (y1, params->log2stdev1, yseed, dom);
      dom++;
      polyvec_grandom (y2, params->log2stdev2, yseed, dom);
      dom++;

      polyvec_set (w, y22);
      polyvec_addmul (w, A1, y1, 0);
      polyvec_addmul (w, A2prime, y21, 0); // XXX w correct, 7 rejs

      polyvec_dcompress_decompose (
          w1, w0, w, dcomp_param); // w1*gamma+w0 == w XXX correct

      coder_enc_begin (cstate, out);
      coder_enc_urandom3 (cstate, w1, m, log2m);
      coder_enc_end (cstate);

      outlen = coder_get_offset (cstate);
      ASSERT_ERR (outlen % 8 == 0);
      ASSERT_ERR (outlen / 8 <= CEIL (log2m * d * kmsis, 8) + 1);
      outlen >>= 3; /* nbits to nbytes */

      shake128_init (hstate);
      shake128_absorb (hstate, hash, 32);
      shake128_absorb (hstate, out, outlen);
      shake128_squeeze (hstate, cseed, 32);

      poly_urandom_autostable (c, params->omega, params->log2omega, cseed, 0);

      polyvec_scale2 (cs1, c, s1);
      polyvec_scale2 (cs2, c, s2); // XXX correct
      polyvec_add (y1, y1, cs1, 0);
      polyvec_add (y2, y2, cs2,
                   0); // XXX correct, XXX sample sign with bimodal ?

      if (params->rej1)
        {
          polyvec_fromcrt (y1);
          polyvec_fromcrt (cs1);
          rej = rej_standard (rngstate, z1coeffs, cs1coeffs, params->scM1,
                              params->stdev1sqr);
          if (rej)
            continue;
        }
      if (params->rej2)
        {
          polyvec_fromcrt (y2);
          polyvec_fromcrt (cs2);
          rej = rej_standard (rngstate, z2coeffs, cs2coeffs, params->scM2,
                              params->stdev2sqr);
          if (rej)
            continue;
        }

      polyvec_subscale2 (y22, c, tA2, 0);
      polyvec_sub (y22, y22, w0, 0); // XXX correct

      polyvec_l2sqr (norm, y2);

      rej = int_gt (norm, params->Bsqr);
      if (rej)
        continue;

      break;
    }

  polyvec_scale (w1, gamma, w1);
  polyvec_sub (w1, w1, y22, 0); // correct XXX

  /* output proof (z1,z21,h,c) */
  polyvec_set (z1, y1);
  polyvec_set (z21, y21);

  polyvec_mod (w1, w1);
  polyvec_redc (w1, w1); // XXX
  polyvec_dcompress_make_ghint (h, y22, w1, dcomp_param);

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);

  /* cleanup */
  shake128_clear (hstate);
  rng_clear (rngstate);
  polyvec_free (y1);
  polyvec_free (y2);
  polyvec_free (cs1);
  polyvec_free (cs2);
  polyvec_free (w);
  polyvec_free (w1);
  polyvec_free (w0);
}

/*
 * Verify opening proof (c,z1,z21,h) from hash of transcript, commitment
 * part tA1 and public key parts (A1,A2prime).
 * Also update hash of transcript by hashing c into it.
 * Outputs 1 on successful verification.
 *
 * c in Rq from challenge space defined by params,
 * z1 in Rq^m1,
 * z21 in Rq^(m2-kmsis),
 * h hint in Rq^(kmsis),
 * tA1 in Rq^kmsis
 * A1 uniform in Rq^(kmsis x m1),
 * A2prime uniform in Rq^(kmsis x (m2-kmsis))
 */
int
abdlop_verify (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
               polyvec_t h, polyvec_t tA1, polymat_t A1, polymat_t A2prime,
               const abdlop_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  const unsigned int m2 = params->m2;
#endif
  polyring_srcptr Rq = params->ring;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int d = polyring_get_deg (Rq);
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int nlimbs = int_get_nlimbs (q);
  dcompress_params_srcptr dcomp_param = params->dcompress;
  int_srcptr gamma = dcompress_get_gamma (dcomp_param);
  int_srcptr m = dcompress_get_m (dcomp_param);
  const unsigned int D = dcompress_get_d (dcomp_param);
  const unsigned int log2m = dcompress_get_log2m (dcomp_param);
  polyvec_t w1, tmp1;
  poly_t tmp2, c2;
  INT_T (l2sqr, nlimbs * 2);
  INT_T (bnd, nlimbs * 2);
  INT_T (linf, nlimbs);
  INT_T (tmp, 1);
  unsigned int outlen;
  shake128_state_t hstate;
  coder_state_t cstate;
  uint8_t out[CEIL (log2m * d * kmsis, 8) + 1];
  uint8_t cseed[32];
  int b, accept = 0;

  ASSERT_ERR (poly_get_ring (c) == Rq);
  ASSERT_ERR (polyvec_get_ring (z1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1);
  ASSERT_ERR (polyvec_get_ring (z21) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_ring (h) == Rq);
  ASSERT_ERR (polyvec_get_nelems (h) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tA1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA1) == kmsis);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);

  polyvec_alloc (w1, Rq, kmsis);
  polyvec_alloc (tmp1, Rq, kmsis);
  poly_alloc (tmp2, Rq);
  poly_alloc (c2, Rq);

  /* recover w1 */

  polyvec_mul (tmp1, A1, z1);
  polyvec_addmul (tmp1, A2prime, z21, 0);
  poly_lshift (tmp2, c, D);
  polyvec_subscale2 (tmp1, tmp2, tA1, 0);

  polyvec_mod (tmp1, tmp1); // XXX
  polyvec_dcompress_use_ghint (w1, h, tmp1, dcomp_param);

  /* recover challenge from w1 */

  coder_enc_begin (cstate, out);
  coder_enc_urandom3 (cstate, w1, m, log2m);
  coder_enc_end (cstate);

  outlen = coder_get_offset (cstate);
  ASSERT_ERR (outlen % 8 == 0);
  ASSERT_ERR (outlen / 8 <= CEIL (log2m * d * kmsis, 8) + 1);
  outlen >>= 3; /* nbits to nbytes */

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, outlen);
  shake128_squeeze (hstate, cseed, 32);

  poly_urandom_autostable (c2, params->omega, params->log2omega, cseed, 0);
  b = poly_eq (c, c2);
  if (!b)
    goto ret;

  /* check bounds */

  /*
   * l2(z1)^2 <= sigma1^2 * 2 * m1 *d = 1.55^2 * 2^(log2sigma1*2) * 2 * m1 * d
   * 1.55^2=2.4025, scale by 400: 1.55^2*400=962
   */
  int_set_i64 (tmp, 400);
  int_set_i64 (bnd, 2 * m1 * d * 962);
  int_lshift (bnd, bnd, 2 * params->log2stdev1);
  int_div (bnd, tmp, bnd, tmp); /* div by 400, overlap ok ? XXX */

  polyvec_l2sqr (l2sqr, z1);
  b = int_le (l2sqr, bnd);
  if (!b)
    goto ret;

  polyvec_subscale (tmp1, gamma, w1, 0);
  polyvec_l2sqr (l2sqr, tmp1);
  b = int_le (l2sqr, params->Bsqr);
  if (!b)
    goto ret;

  /* 2*linf(h) <= m */
  polyvec_linf (linf, h);
  b = int_le (linf, m);
  if (!b)
    goto ret;

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);
  accept = 1;
ret:
  /* cleanup */
  shake128_clear (hstate);
  polyvec_free (w1);
  polyvec_free (tmp1);
  poly_free (tmp2);
  poly_free (c2);
  return accept;
}
