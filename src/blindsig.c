#include "blindsig.h"
#include "stopwatch.h"

void
falcon_redc (int16_t c[DEG_])
{
  unsigned int i;

  for (i = 0; i < DEG_; i++)
    {
      if (c[i] > (Q_ - 1) / 2)
        c[i] -= Q_;
      else if (c[i] < -(Q_ - 1) / 2)
        c[i] += Q_;
    }
}

void
falcon_add (int16_t c[DEG_], const int16_t a[DEG_], const int16_t b[DEG_])
{
  unsigned int i;

  for (i = 0; i < DEG_; i++)
    c[i] = (a[i] + b[i]) % Q_;
}

void
falcon_mul (int16_t c[DEG_], const int16_t a[DEG_], const int16_t b[DEG_])
{
  int32_t t[2 * DEG_] = { 0 };
  unsigned int i, j;

  for (i = 0; i < DEG_; i++)
    for (j = 0; j < DEG_; j++)
      t[i + j] += (int32_t)a[i] * b[j] % Q_;

  for (i = 0; i < DEG_; i++)
    c[i] = (t[i] - t[DEG_ + i]) % Q_;
}

/* generate compressed secret and public key. */
void
falcon_keygen (uint8_t sk[PRIVKEYLEN], uint8_t pk[PUBKEYLEN])
{
  shake256_context rng;
  uint8_t tmpkg[TMPKGLEN];
  int r;

  shake256_init_prng_from_system (&rng);

  r = falcon_keygen_make (&rng, SIGNER_LOGN, (void *)sk, PRIVKEYLEN,
                          (void *)pk, PUBKEYLEN, tmpkg, TMPKGLEN);
  if (r != 0)
    {
      fprintf (stderr, "falcon keygen failed: %d.\n", r);
      exit (EXIT_FAILURE);
    }
}

/* compressed public key to coefficient representation */
void
falcon_decode_pubkey (int16_t h[DEG_], const uint8_t pk[PUBKEYLEN])
{
  if (Zf (modq_decode) ((uint16_t *)h, SIGNER_LOGN, pk + 1, PUBKEYLEN - 1)
      != PUBKEYLEN - 1)
    {
      fprintf (stderr, "falcon decoding of pubkey failed.\n");
      exit (EXIT_FAILURE);
    }
}

/* find (s1,s2) s.t.: (1,h) * (s1,s2)^T = t*/
void
falcon_preimage_sample (int16_t s1[DEG_], int16_t s2[DEG_],
                        const int16_t t[DEG_], const uint8_t sk[PRIVKEYLEN])
{
  __attribute__ ((aligned (8))) uint8_t tmp[72 * DEG_];
  int8_t f[DEG_], g[DEG_], F[DEG_], G[DEG_];
  uint16_t h[DEG_], tu[DEG_];
  shake256_context rng;
  unsigned oldcw;
  int u, v;

  shake256_init_prng_from_system (&rng);

  /* decode private key elements */
  u = 1;
  v = Zf (trim_i8_decode) (f, SIGNER_LOGN, Zf (max_fg_bits)[SIGNER_LOGN],
                           sk + u, SKBYTES - u);
  if (v == 0)
    goto err;

  u += v;
  v = Zf (trim_i8_decode) (g, SIGNER_LOGN, Zf (max_fg_bits)[SIGNER_LOGN],
                           sk + u, SKBYTES - u);
  if (v == 0)
    goto err;

  u += v;
  v = Zf (trim_i8_decode) (F, SIGNER_LOGN, Zf (max_FG_bits)[SIGNER_LOGN],
                           sk + u, SKBYTES - u);
  if (v == 0)
    goto err;

  u += v;
  if (u != SKBYTES)
    goto err;

  /* complete private key */
  if (!Zf (complete_private) (G, f, g, F, SIGNER_LOGN, tmp))
    goto err;

  for (u = 0; u < DEG_; u++)
    tu[u] = t[u] + ((t[u] >> 15) & Q_);

  oldcw = set_fpu_cw (2);
  Zf (sign_dyn) (s2, (inner_shake256_context *)&rng, f, g, F, G, tu,
                 SIGNER_LOGN, tmp);
  set_fpu_cw (oldcw);

  Zf (compute_public) (h, f, g, SIGNER_LOGN, tmp);
  Zf (to_ntt_monty) (h, SIGNER_LOGN);
  if (!Zf (reconstruct_s1) (s1, tu, s2, h, SIGNER_LOGN, tmp))
    goto err;

  return;
err:
  fprintf (stderr, "falcon preimage sampling failed.\n");
  exit (EXIT_FAILURE);
}

void
signer_clear (signer_state_t state)
{
  lin_verifier_clear (state->p1);
}

void
signer_keygen (uint8_t sk[PRIVKEYLEN], uint8_t pk[PUBKEYLEN])
{
  falcon_keygen (sk, pk);
}

void
signer_init (signer_state_t state, const uint8_t pubkey[PUBKEYLEN],
             const uint8_t privkey[PRIVKEYLEN])
{
  polyring_srcptr Rq1 = p1_params->tbox_params->tbox->ring;
  int_srcptr q1 = polyring_get_mod (Rq1);
  POLYRING_T (Rprime1, q1, 512);
  shake128_state_t hstate;
  uint8_t expseed[64];
  const uint8_t *signerppseed = expseed;
  const uint8_t *p1ppseed = expseed + 32;
  INTVEC_T (coeffvec, 512, 1);
  polymat_t A;
  poly_ptr poly;
  unsigned int i;
  INT_T (lo, 1);
  INT_T (hi, 1);

  polymat_alloc (A, Rprime1, 1, 3);

  memset (state, 0, sizeof (signer_state_t));

  memcpy (state->privkey, privkey, PRIVKEYLEN);
  falcon_decode_pubkey (state->pubkey, pubkey);

  int_set_i64 (lo, -(FALCON_P - 1) / 2);
  int_set_i64 (hi, (FALCON_P - 1) / 2);

  /* expand public parameter seed to obtain seeds for protocol and
   * sub-protocols */
  shake128_init (hstate);
  shake128_absorb (hstate, public_randomness, 32);
  shake128_squeeze (hstate, expseed, sizeof (expseed));
  shake128_clear (hstate);

  /* Ar1 */
  intvec_urandom_bnd (coeffvec, lo, hi, signerppseed, 0);
  for (i = 0; i < 512; i++)
    state->Ar1[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Ar2 */
  intvec_urandom_bnd (coeffvec, lo, hi, signerppseed, 1);
  for (i = 0; i < 512; i++)
    state->Ar2[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Am */
  intvec_urandom_bnd (coeffvec, lo, hi, signerppseed, 2);
  for (i = 0; i < 512; i++)
    state->Am[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Atau */
  intvec_urandom_bnd (coeffvec, lo, hi, signerppseed, 3);
  for (i = 0; i < 512; i++)
    state->Atau[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* A = [Ar1,Ar2,Am] in Rprime^(1x3) */
  poly = polymat_get_elem (A, 0, 0);
  poly_set_coeffvec_i16 (poly, state->Ar1);
  poly = polymat_get_elem (A, 0, 1);
  poly_set_coeffvec_i16 (poly, state->Ar2);
  poly = polymat_get_elem (A, 0, 2);
  poly_set_coeffvec_i16 (poly, state->Am);

  /*
   * P1:
   * As+t=0 over Rprime
   * A=(Ar1,Ar2,Am) in Rprime^(1x3)
   * s=(r1,r2,m) in Rprime^3
   * t=-t in Rprime
   */
  lin_verifier_init (state->p1, p1ppseed, p1_params);
  lin_verifier_set_statement_A (state->p1, A);

  polymat_free (A);
}

static void
_signer_sign (signer_state_t state, uint8_t *blindsig, size_t *blindsiglen,
              const int16_t t[512])
{
  INTVEC_T (s1_, 512, 1);
  INTVEC_T (s2_, 512, 1);
  int16_t tau[512];
  int16_t s1[512];
  int16_t s2[512];
  uint8_t tau_[512 / 8];
  int16_t t_[512];
  coder_state_t cstate;
  unsigned int i, len;
  INT_T (mod, 1);

  /* l2(s1,s2) <= sqrt(BOUND), so linf(s1,s2) <= floor(sqrt(BOUND)) = 5833 */
  int_set_i64 (mod, 2 * 5833 + 1); /* [-5833,5833] -> [0,2*5833+1) */
  // unsigned int log2mod = 14;                    /* ceil(log(2*5833+1)) */

  bytes_urandom (tau_, sizeof (tau_));

  /* expand tau_ */
  for (i = 0; i < 512; i++)
    {
      const unsigned int q = i >> 3;
      const unsigned int r = i - (q << 3);

      tau[i] = (tau_[q] & (1 << r)) >> r;
    }

  /* t_ = t+Atau*tau */
  falcon_mul (t_, state->Atau, tau);
  falcon_add (t_, t_, t);
  /* sample (s1,s2) s.t.: (1,pk)*(s1,s2)^T = t_ */
  falcon_preimage_sample (s1, s2, t_, state->privkey);

  /* encode blindsig*/

  intvec_set_i16 (s1_, s1);
  intvec_set_i16 (s2_, s2);

  coder_enc_begin (cstate, blindsig);
  coder_enc_bytes (cstate, tau_, sizeof (tau_));
  coder_enc_grandom (cstate, s1_, 7); // 7 = ceil(log(165/1.55,2))
  coder_enc_grandom (cstate, s2_, 7);
  coder_enc_end (cstate);
  len = coder_get_offset (cstate);
  ASSERT_ERR (len % 8 == 0);
  len >>= 3; /* nbits to nbytes */
  *blindsiglen = len;
}

int
signer_sign (signer_state_t state, uint8_t *blindsig, size_t *blindsiglen,
             const uint8_t *masked_msg, UNUSED size_t masked_msglen)
{
  polyring_srcptr Rq1 = p1_params->tbox_params->tbox->ring;
  int_srcptr q1 = polyring_get_mod (Rq1);
  POLYRING_T (Rprime1, q1, 512);
  const unsigned int d1 = polyring_get_deg (Rq1);
  const unsigned int k1 = 512 / d1;
  INT_T (p, 1);
  INT_T (pinv1, 1);
  poly_ptr poly;
  size_t len;
  int b = 0, rc;
  int16_t t[512];
  coder_state_t cstate;
  unsigned int log2p, off;
  INTVEC_T (t_, 512, 1);
  polyvec_t tvec, u;

  STOPWATCH_START (stopwatch_blindsig_signer_sign, "signer_sign");

  polyvec_alloc (tvec, Rprime1, 1);
  polyvec_alloc (u, Rq1, 1 * k1);

  int_set_i64 (p, FALCON_P);
  log2p = 14; /* ceil(log(12289)) */

  /* decode masked message */

  coder_dec_begin (cstate, masked_msg);
  rc = coder_dec_urandom (cstate, t_, p, log2p);
  if (rc != 0)
    goto ret; /* invalid encoding */
  rc = coder_dec_end (cstate);
  if (rc != 1)
    goto ret; /* invalid padding */
  off = coder_get_offset (cstate);
  ASSERT_ERR (off % 8 == 0);
  off >>= 3; /* nbits to nbytes */

  intvec_redc (t_, t_, p);
  intvec_get_i16 (t, t_);

  /* verify P1 */

  poly = polyvec_get_elem (tvec, 0);
  poly_set_coeffvec (poly, t_);
  poly_neg_self (poly);

  lin_verifier_set_statement_t (state->p1, tvec);

  /* verify proof appended to the encoding of t */
  b = lin_verifier_verify (state->p1, masked_msg + off, &len);
  ASSERT_ERR (off + len == masked_msglen);
  if (b != 1)
    goto ret;

  /* create blind signature */

  _signer_sign (state, blindsig, blindsiglen, t);

ret:
  polyvec_free (tvec);
  polyvec_free (u);

  STOPWATCH_STOP (stopwatch_blindsig_signer_sign);
  return b;
}

#ifdef XXX
void
signer_test (void)
{
  int16_t lhs[512], rhs[512];
  uint8_t seed[32];
  uint64_t dom;
  pubkey_t pubkey;
  blindsig_t blindsig;
  masked_msg_t masked_msg;
  signer_state_t state;
  rng_state_t rng;
  int32_t sqrnorm;
  unsigned int i;

  /* create some random t and ppseed for testing */

  bytes_urandom (seed, sizeof (seed));

  dom = 0;
  rng_init (rng, seed, dom);

  rng_urandom (rng, (uint8_t *)masked_msg->t, sizeof (masked_msg->t));
  for (i = 0; i < DEG_; i++)
    masked_msg->t[i] %= Q_;

  signer_init (state, pubkey);
  _signer_sign (state, blindsig, masked_msg);

  falcon_mul (lhs, pubkey->pk, blindsig->s2);
  falcon_add (lhs, lhs, blindsig->s1);

  falcon_mul (rhs, state->Atau, blindsig->tau);
  falcon_add (rhs, rhs, masked_msg->t);

  /* check norm (s1,s2) */
  sqrnorm = 0;
  for (i = 0; i < DEG_; i++)
    {
      sqrnorm += (int32_t)blindsig->s1[i] * blindsig->s1[i];
      sqrnorm += (int32_t)blindsig->s2[i] * blindsig->s2[i];
    }
  if (sqrnorm > BOUND)
    {
      fprintf (stderr, "signer test failed: (s1,s2) not small.\n");
      exit (EXIT_FAILURE);
    }

  /* check lhs == rhs mod Q */
  for (i = 0; i < DEG_; i++)
    {
      if ((lhs[i] - rhs[i]) % Q_ != 0)
        {
          fprintf (stderr,
                   "signer test failed: (1,pk)*(s1,s2)^T != t+Atau*tau.\n");
          exit (EXIT_FAILURE);
        }
    }

  lnp_verifier_clear (state->p1);
}
#endif

// Falcon-512
#define Q 12289

void
verifier_clear (verifier_state_t state)
{
  lin_verifier_clear (state->p2);
}

void
verifier_init (verifier_state_t state, const uint8_t pubkey[PUBKEYLEN])
{
  polyring_srcptr Rq2 = p2_params->tbox_params->tbox->ring;
  int_srcptr q2 = polyring_get_mod (Rq2);
  POLYRING_T (Rprime2, q2, 512);
  INTVEC_T (coeffvec, 512, 1);
  shake128_state_t hstate;
  unsigned int i;
  uint8_t expseed[32 * 3];
  const uint8_t *verifierppseed = expseed;
  const uint8_t *p2ppseed = expseed + 64;
  INT_T (lo, 1);
  INT_T (hi, 1);
  polymat_t A2;
  polyvec_t tvec2;
  poly_ptr poly;
  int_ptr coeff;

  polymat_alloc (A2, Rprime2, 1, 5);
  polyvec_alloc (tvec2, Rprime2, 1);

  memset (state, 0, sizeof (verifier_state_t));

  falcon_decode_pubkey (state->pubkey, pubkey);

  int_set_i64 (lo, -(FALCON_P - 1) / 2);
  int_set_i64 (hi, (FALCON_P - 1) / 2);

  /* expand public parameter seed to obtain seeds for protocol and
   * sub-protocols */
  shake128_init (hstate);
  shake128_absorb (hstate, public_randomness, 32);
  shake128_squeeze (hstate, expseed, sizeof (expseed));
  shake128_clear (hstate);

  /* Ar1 */
  intvec_urandom_bnd (coeffvec, lo, hi, verifierppseed, 0);
  for (i = 0; i < 512; i++)
    state->Ar1[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Ar2 */
  intvec_urandom_bnd (coeffvec, lo, hi, verifierppseed, 1);
  for (i = 0; i < 512; i++)
    state->Ar2[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Am */
  intvec_urandom_bnd (coeffvec, lo, hi, verifierppseed, 2);
  for (i = 0; i < 512; i++)
    state->Am[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Atau */
  intvec_urandom_bnd (coeffvec, lo, hi, verifierppseed, 3);
  for (i = 0; i < 512; i++)
    state->Atau[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* A = [Ar1,Ar2,Atau,-1,-B2] in Rprime^(1x5) */
  poly = polymat_get_elem (A2, 0, 0);
  poly_set_coeffvec_i16 (poly, state->Ar1);
  poly = polymat_get_elem (A2, 0, 1);
  poly_set_coeffvec_i16 (poly, state->Ar2);
  poly = polymat_get_elem (A2, 0, 2);
  poly_set_coeffvec_i16 (poly, state->Atau);
  poly = polymat_get_elem (A2, 0, 3);
  poly_set_zero (poly);
  int_set_i64 (poly_get_coeff (poly, 0), -1);
  poly = polymat_get_elem (A2, 0, 4);
  for (i = 0; i < 512; i++)
    {
      coeff = poly_get_coeff (poly, i);
      int_set_i64 (coeff, -state->pubkey[i]);
    }

  lin_verifier_init (state->p2, p2ppseed, p2_params);
  lin_verifier_set_statement_A (state->p2, A2);

  /*
   * P2:
   * As*s+Am*m+t=0 over Rprime
   * As=(Ar1,Ar2,Tau) in Rprime^(1x3)
   * s=(r1,r2,tau) in Rprime^3
   * Am=(-B1,-B2) in Rprime^2
   * m=(s1,s2) in Rprime
   * t=Am*m in Rprime
   */

  polyvec_free (tvec2);

  polymat_free (A2);
}

int
verifier_vrfy (verifier_state_t state, const uint8_t msg[512 / 8],
               const uint8_t *sig, UNUSED size_t siglen)
{
  polyring_srcptr Rq2 = p2_params->tbox_params->tbox->ring;
  int_srcptr q2 = polyring_get_mod (Rq2);
  POLYRING_T (Rprime2, q2, 512);
  polyvec_t u_, tvec;
  size_t len;
  int b = 0;
  int16_t m[512];
  unsigned int i;
  poly_ptr poly;

  STOPWATCH_START (stopwatch_blindsig_verifier_vrfy, "verifier_vrfy");

  polyvec_alloc (tvec, Rprime2, 1);

  /* expand message */
  for (i = 0; i < 512; i++)
    {
      const unsigned int quot = i >> 3;
      const unsigned int rem = i - (quot << 3);

      m[i] = (msg[quot] & (1 << rem)) >> rem;
    }

  /* m <- Am*m */
  falcon_mul (m, state->Am, m);
  falcon_redc (m);

  polyvec_get_subvec (u_, tvec, 0, 1, 1);
  poly = polyvec_get_elem (u_, 0);
  poly_set_coeffvec_i16 (poly, m);

  lin_verifier_set_statement_t (state->p2, u_);

  b = lin_verifier_verify (state->p2, sig, &len);
  ASSERT_ERR (len == siglen);
  if (b != 1)
    goto ret;

ret:
  polyvec_free (tvec);
  STOPWATCH_STOP (stopwatch_blindsig_verifier_vrfy);
  return b;
}

void
user_clear (user_state_t state)
{
  lin_prover_clear (state->p1);
  lin_prover_clear (state->p2);
}

void
user_init (user_state_t state, const uint8_t pubkey[PUBKEYLEN])
{
  polyring_srcptr Rq1 = p1_params->tbox_params->tbox->ring;
  polyring_srcptr Rq2 = p2_params->tbox_params->tbox->ring;
  int_srcptr q1 = polyring_get_mod (Rq1);
  int_srcptr q2 = polyring_get_mod (Rq2);
  POLYRING_T (Rprime1, q1, 512);
  POLYRING_T (Rprime2, q2, 512);
  shake128_state_t hstate;
  uint8_t expseed[32 * 3];
  const uint8_t *userppseed = expseed;
  const uint8_t *p1ppseed = expseed + 32;
  const uint8_t *p2ppseed = expseed + 64;
  INTVEC_T (coeffvec, 512, 1);
  INT_T (lo, 1);
  INT_T (hi, 1);
  polymat_t A1, A2;
  poly_ptr poly;
  unsigned int i;
  int_ptr coeff;

  polymat_alloc (A1, Rprime1, 1, 3);
  polymat_alloc (A2, Rprime2, 1, 5);

  memset (state, 0, sizeof (user_state_t));

  falcon_decode_pubkey (state->pubkey, pubkey);

  int_set_i64 (lo, -(FALCON_P - 1) / 2);
  int_set_i64 (hi, (FALCON_P - 1) / 2);

  /* expand public parameter seed to obtain seeds for protocol and
   * sub-protocols */
  shake128_init (hstate);
  shake128_absorb (hstate, public_randomness, 32);
  shake128_squeeze (hstate, expseed, sizeof (expseed));
  shake128_clear (hstate);

  /* Ar1 */
  intvec_urandom_bnd (coeffvec, lo, hi, userppseed, 0);
  for (i = 0; i < 512; i++)
    state->Ar1[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Ar2 */
  intvec_urandom_bnd (coeffvec, lo, hi, userppseed, 1);
  for (i = 0; i < 512; i++)
    state->Ar2[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Am */
  intvec_urandom_bnd (coeffvec, lo, hi, userppseed, 2);
  for (i = 0; i < 512; i++)
    state->Am[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* Atau */
  intvec_urandom_bnd (coeffvec, lo, hi, userppseed, 3);
  for (i = 0; i < 512; i++)
    state->Atau[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);

  /* A = [Ar1,Ar2,Am] in Rprime^(1x3) */
  poly = polymat_get_elem (A1, 0, 0);
  poly_set_coeffvec_i16 (poly, state->Ar1);
  poly = polymat_get_elem (A1, 0, 1);
  poly_set_coeffvec_i16 (poly, state->Ar2);
  poly = polymat_get_elem (A1, 0, 2);
  poly_set_coeffvec_i16 (poly, state->Am);

  /*
   * P1:
   * As+t=0 over Rprime
   * A=(Ar1,Ar2,Am) in Rprime^(1x3)
   * s=(r1,r2,m) in Rprime^3
   * t=-t in Rprime
   */
  lin_prover_init (state->p1, p1ppseed, p1_params);
  lin_prover_set_statement_A (state->p1, A1);

  /*
   * P2:
   * As*s+Am*m+t=0 over Rprime
   * As=(Ar1,Ar2,Tau) in Rprime^(1x3)
   * s=(r1,r2,tau) in Rprime^3
   * Am=(-B1,-B2) in Rprime^2
   * m=(s1,s2) in Rprime
   * t=Am*m in Rprime
   */

  /* A = [Ar1,Ar2,Atau,-B1,-B2] in Rprime^(1x5) */
  poly = polymat_get_elem (A2, 0, 0);
  poly_set_coeffvec_i16 (poly, state->Ar1);
  poly = polymat_get_elem (A2, 0, 1);
  poly_set_coeffvec_i16 (poly, state->Ar2);
  poly = polymat_get_elem (A2, 0, 2);
  poly_set_coeffvec_i16 (poly, state->Atau);
  poly = polymat_get_elem (A2, 0, 3);
  poly_set_zero (poly);
  coeff = poly_get_coeff (poly, 0);
  int_set_i64 (coeff, -1);
  poly = polymat_get_elem (A2, 0, 4);
  for (i = 0; i < 512; i++)
    {
      coeff = poly_get_coeff (poly, i);
      int_set_i64 (coeff, -state->pubkey[i]);
    }

  /* init prover state */

  lin_prover_init (state->p2, p2ppseed, p2_params);
  lin_prover_set_statement_A (state->p2, A2);

  polymat_free (A1);
  polymat_free (A2);
}

void
user_maskmsg (user_state_t state, uint8_t *masked_msg, size_t *masked_msglen,
              const uint8_t msg[512 / 8])
{
  polyring_srcptr Rq1 = p1_params->tbox_params->tbox->ring;
  int_srcptr q1 = polyring_get_mod (Rq1);
  POLYRING_T (Rprime1, q1, 512);
  INT_T (p, 1);
  polyvec_t u_, s_;
  INTVEC_T (coeffvec, 512, 1);
  int16_t tmp[512];
  unsigned int i;
  int16_t m[512];
  uint8_t dom, seed[64];
  const uint8_t *protoseed = seed;
  const uint8_t *subprotoseed = seed + 32;
  int64_t l2sqr;
  int32_t coeff = 0;
  INTVEC_T (t_, 512, 1);
  int16_t t[512];
  coder_state_t cstate;
  unsigned int log2p, len;
  polyvec_t tvec;
  poly_ptr poly;

  STOPWATCH_START (stopwatch_blindsig_user_maskmsg, "user_maskmsg");

  polyvec_alloc (tvec, Rprime1, 3);

  memcpy (state->m, msg, 512 / 8);

  int_set_i64 (p, FALCON_P);
  log2p = 14; /* ceil(log(12289)) */

  /* randomness for (sub)protocols */
  bytes_urandom (seed, sizeof (seed));

  /* expand message */
  for (i = 0; i < 512; i++)
    {
      const unsigned int quot = i >> 3;
      const unsigned int rem = i - (quot << 3);

      m[i] = (state->m[quot] & (1 << rem)) >> rem;
    }

  /* sample gaussian (r1,r2) with sigma_r = 1.55*2 = 3.1 */
  dom = 0;
  do
    {
      l2sqr = 0;

      intvec_grandom (coeffvec, 1, protoseed, dom++);
      for (i = 0; i < 512; i++)
        {
          coeff = (int32_t)intvec_get_elem_i64 (coeffvec, i);
          state->r1[i] = coeff;
          l2sqr += (int64_t)coeff * coeff;
        }

      intvec_grandom (coeffvec, 1, protoseed, dom++);
      for (i = 0; i < 512; i++)
        {
          state->r2[i] = (int32_t)intvec_get_elem_i64 (coeffvec, i);
          l2sqr += (int64_t)coeff * coeff;
        }
    }
  while (l2sqr > 109 * 109);

  /* t = Ar1*r1 + Ar2*r2 + Am*m */
  falcon_mul (t, state->Ar1, state->r1);
  falcon_mul (tmp, state->Ar2, state->r2);
  falcon_add (t, t, tmp);
  falcon_mul (tmp, state->Am, m);
  falcon_add (t, t, tmp);

  /* encode t */

  intvec_set_i16 (t_, t);
  intvec_redp (t_, t_, p);

  coder_enc_begin (cstate, masked_msg);
  coder_enc_urandom (cstate, t_, p, log2p);
  coder_enc_end (cstate);
  len = coder_get_offset (cstate);
  ASSERT_ERR (len % 8 == 0);
  len >>= 3; /* nbits to nbytes */

  /*
   * P1:
   * As+t=0 over Rprime
   * A=(Ar1,Ar2,Am) in Rprime^(1x3)
   * s=(r1,r2,m) in Rprime^3
   * t=-t in Rprime
   */

  intvec_redc (t_, t_, p);
  polyvec_get_subvec (u_, tvec, 0, 1, 1);
  poly = polyvec_get_elem (u_, 0);
  poly_set_coeffvec (poly, t_);
  poly_neg_self (poly);

  lin_prover_set_statement_t (state->p1, u_);

  polyvec_get_subvec (s_, tvec, 0, 3, 1);
  poly = polyvec_get_elem (s_, 0);
  poly_set_coeffvec_i16 (poly, state->r1);
  poly = polyvec_get_elem (s_, 1);
  poly_set_coeffvec_i16 (poly, state->r2);
  poly = polyvec_get_elem (s_, 2);
  poly_set_coeffvec_i16 (poly, m);

  lin_prover_set_witness (state->p1, s_);

  /* append P1 to encoding of t */
  lin_prover_prove (state->p1, masked_msg + len, masked_msglen, subprotoseed);
  *masked_msglen += len;

  polyvec_free (tvec);
  STOPWATCH_STOP (stopwatch_blindsig_user_maskmsg);
}

int
user_sign (user_state_t state, uint8_t *sig, size_t *siglen,
           const uint8_t *blindsig, UNUSED size_t blindsiglen)
{
  polyring_srcptr Rq2 = p2_params->tbox_params->tbox->ring;
  int_srcptr q2 = polyring_get_mod (Rq2);
  POLYRING_T (Rprime2, q2, 512);
  poly_ptr poly;
  uint8_t seed[32];
  polyvec_t u_, tvec;
  coder_state_t cstate;
  int16_t msg[512];
  unsigned int i;
  UNUSED unsigned int len;
  INTVEC_T (s1_, 512, 1);
  INTVEC_T (s2_, 512, 1);
  int16_t tau[512];
  int16_t s1[512];
  int16_t s2[512];
  uint8_t tau_[512 / 8];
  INT_T (mod, 1);
  int rc = 0;
  int rv = -1;

  STOPWATCH_START (stopwatch_blindsig_user_sign, "user_sign");

  polyvec_alloc (tvec, Rprime2, 5);

  /* decode blindsig */

  /* l2(s1,s2) <= sqrt(BOUND), so linf(s1,s2) <= floor(sqrt(BOUND)) = 5833 */
  int_set_i64 (mod, 2 * 5833 + 1); /* [-5833,5833] -> [0,2*5833+1) */
  // unsigned int log2mod = 14;                    /* ceil(log(2*5833+1)) */

  coder_dec_begin (cstate, blindsig);
  rc = coder_dec_bytes (cstate, tau_, sizeof (tau_));
  coder_dec_grandom (cstate, s1_, 7); // 7 = ceil(log(165/1.55,2))
  coder_dec_grandom (cstate, s2_, 7);
  if (rc != 0)
    goto ret; /* decoding failed */
  rc = coder_dec_end (cstate);
  if (rc != 1)
    goto ret; /* invalid padding */
  len = coder_get_offset (cstate);
  ASSERT_ERR (len % 8 == 0);
  len >>= 3; /* nbits to nbytes */
  ASSERT_ERR (len == blindsiglen);

  /* expand tau_ */
  for (i = 0; i < 512; i++)
    {
      const unsigned int q = i >> 3;
      const unsigned int r = i - (q << 3);

      tau[i] = (tau_[q] & (1 << r)) >> r;
    }

  intvec_get_i16 (s1, s1_);
  intvec_get_i16 (s2, s2_);

  /* generate sig */

  /* randomness for (sub)protocols */
  bytes_urandom (seed, sizeof (seed));

  /* expand message */
  for (i = 0; i < 512; i++)
    {
      const unsigned int quot = i >> 3;
      const unsigned int rem = i - (quot << 3);

      msg[i] = (state->m[quot] & (1 << rem)) >> rem;
    }

  /* m <- Am*m */
  falcon_mul (msg, state->Am, msg);
  falcon_redc (msg);

  polyvec_get_subvec (u_, tvec, 0, 1, 1);
  poly = polyvec_get_elem (u_, 0);
  poly_set_coeffvec_i16 (poly, msg);

  lin_prover_set_statement_t (state->p2, u_);

  poly = polyvec_get_elem (tvec, 0);
  poly_set_coeffvec_i16 (poly, state->r1);
  poly = polyvec_get_elem (tvec, 1);
  poly_set_coeffvec_i16 (poly, state->r2);
  poly = polyvec_get_elem (tvec, 2);
  poly_set_coeffvec_i16 (poly, tau);

  poly = polyvec_get_elem (tvec, 3);
  poly_set_coeffvec_i16 (poly, s1);
  poly = polyvec_get_elem (tvec, 4);
  poly_set_coeffvec_i16 (poly, s2);

  lin_prover_set_witness (state->p2, tvec);
  lin_prover_prove (state->p2, sig, siglen, seed);

  rv = 1;
ret:
  polyvec_free (tvec);

  STOPWATCH_STOP (stopwatch_blindsig_user_sign);
  return rv;
}
