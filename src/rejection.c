#include "lazer.h"

#include <mpfr.h>

/* bits of M * 2^128. assume M < 2^6*/
#define NBITS_M (128 + 6)

/*
 * rejection sampling algorithms
 *
 * The sigma2 parameter is the precopmuted square of the standard deviation.
 * All functions returns 1 for reject, 0 for accept.
 */

/* int z is destroyed, XXX what if z=0 ? */
static inline void
_int2mpfr (mpfr_t f, int_t z, unsigned int prec)
{
  mpfr_exp_t exp = z->nlimbs * NBITS_LIMB;
  const int sgn = int_sgn (z);

  ASSERT_ERR (z->nlimbs * NBITS_LIMB >= prec);

  if (int_eqzero(z) != 0) {
    mpfr_custom_init_set (f, sgn * MPFR_ZERO_KIND, exp, prec, z->limbs);
    return;
  }

  while (z->limbs[z->nlimbs - 1] == 0 && exp > 0)
    {
      int_lshift (z, z, NBITS_LIMB);
      exp -= NBITS_LIMB;
    }
  while ((z->limbs[z->nlimbs - 1] & ((limb_t)1 << (NBITS_LIMB - 1))) == 0
         && exp > 0)
    {
      int_lshift (z, z, 1);
      exp--;
    }

  if (z->nlimbs * NBITS_LIMB > prec)
    int_rshift (z, z, z->nlimbs * NBITS_LIMB - prec);

  mpfr_custom_init_set (f, sgn * MPFR_REGULAR_KIND, exp, prec, z->limbs);
}

/*
 * standard rejection sampling
 *
 * u <- [0,1)
 * if u > 1/M * exp( (-2<z,v> + <v,v>) / (2*sigma^2) )
 *   reject
 * else
 *   accept
 *
 * we compute instead
 *
 * u <- {0, ..., 2^128 - 1}
 * scM <- round(M * 2^128)
 * if M*u > 2^256 * exp( (-2*<z,v> + <v,v>) / (2*sigma^2) )
 *   reject
 * else
 *   accept
 *
 * For |v| <= T, kappa security param,
 * sigma = gamma * T
 * M = exp(sqrt(2(kappa+1)/log(e)) * 1/gamma + 1/(2*gamma^2)):
 * Distribution of (v,z) conditioned on accept is within statistical
 * distance of 2^(-kappa) of distribution of v x discrete gaussian
 * of standard deviation sigma.
 */
int
rej_standard (rng_state_t state, const intvec_t z, const intvec_t v,
              const int_t scM, const int_t sigma2)
{
  INT_T (sigma2dbl, sigma2->nlimbs); /* XXX potential overflow */
  INT_T (t1, z->nlimbs << 1);
  INT_T (t2, v->nlimbs << 1);
  INT_T (u, CEIL (NBITS_M, NBITS_LIMB));
  INT_T (mu, CEIL (NBITS_M, NBITS_LIMB) << 1);
  mpfr_t nom, denom, t3, t4;
  int reject;

  ASSERT_ERR (z->nelems == v->nelems);
  ASSERT_ERR (z->nlimbs == v->nlimbs);
  ASSERT_ERR (sigma2dbl->nlimbs * NBITS_LIMB >= 128);
  ASSERT_ERR (t2->nlimbs * NBITS_LIMB >= 128);
  ASSERT_ERR (mu->nlimbs * NBITS_LIMB >= 128);
  ASSERT_ERR (scM->nlimbs == CEIL (NBITS_M, NBITS_LIMB));

  mpfr_init2 (t3, 128);

  /* u <- {0, ..., 2^128 - 1} */
  int_set_i64 (u, 0);
  rng_urandom (state, (uint8_t *)u->limbs, 128 / 8);

  int_mul (mu, scM, u); /* M * u */

  intvec_dot (t1, z, v);  /* <z,v> */
  int_lshift (t1, t1, 1); /* 2*<z,v> */
  intvec_dot (t2, v, v);  /* <v,v>*/
  int_sub (t2, t2, t1);   /* <v,v> - 2*<z,v> */

  int_lshift (sigma2dbl, sigma2, 1); /* 2*sigma^2 */

  _int2mpfr (nom, t2, 128);
  _int2mpfr (denom, sigma2dbl, 128);
  _int2mpfr (t4, mu, 128);

  mpfr_div (t3, nom, denom, MPFR_RNDN); /* (-2<z,v> + <v,v>) / (2*sigma^2) */
  mpfr_exp (t3, t3, MPFR_RNDN); /* exp((-2<z,v> + <v,v>) / (2*sigma^2)) */

  mpfr_mul_2ui (t3, t3, 256, MPFR_RNDN); /* 2^128 * exp(...) */
  if (mpfr_cmp (t4, t3) > 0)
    reject = 1;
  else
    reject = 0;

  mpfr_clear (t3);
  return reject;
}

/*
 * bimodal rejection sampling
 *
 * u <- [0,1)
 * if u > 1 / ( M * exp( -<v,v> / (2*sigma^2) ) * cosh ( <z,v> / (sigma^2) ) )
 *   reject
 * else
 *   accept
 *
 * we compute instead
 *
 * u <- {0, ..., 2^128 - 1}
 * scM <- round(M * 2^128)
 * if ( M * exp(-<v,v>/(2*sigma^2)) * cosh(<z,v>/(sigma^2))) * u > 2^256
 *   reject
 * else
 *   accept
 *
 * For |v| <= T,
 * sigma = gamma*T,
 * M = exp(1/(2*gamma^2)):
 * Distribution of (v,z) conditioned on accept is within statistical
 * distance of 2^(-kappa) of distribution of v x discrete gaussian
 * of standard deviation sigma.
 */
int
rej_bimodal (rng_state_t state, const intvec_t z, const intvec_t v,
             const int_t scM, const int_t sigma2)
{
  INT_T (sigma2dbl, sigma2->nlimbs); /* XXX potential overflow */
  INT_T (_sigma2, sigma2->nlimbs);   /* XXX potential overflow */
  INT_T (t1, z->nlimbs << 1);
  INT_T (t2, v->nlimbs << 1);
  INT_T (u, CEIL (NBITS_M, NBITS_LIMB));
  INT_T (Mu, CEIL (NBITS_M, NBITS_LIMB) << 1);
  mpfr_t nom, denom, t3, t4, t5, t6;
  int reject;

  ASSERT_ERR (z->nelems == v->nelems);
  ASSERT_ERR (z->nlimbs == v->nlimbs);
  ASSERT_ERR (scM->nlimbs == CEIL (NBITS_M, NBITS_LIMB));

  mpfr_init2 (t3, 128);
  mpfr_init2 (t4, 128);
  mpfr_init2 (t6, 128);

  /* u <- {0, ..., 2^128 - 1} */
  int_set_i64 (u, 0);
  rng_urandom (state, (uint8_t *)u->limbs,
               128 / 8); /* XXX not endian neutral */

  int_mul (Mu, scM, u); /* M * u */

  intvec_dot (t1, z, v); /* <z,v> */
  intvec_dot (t2, v, v); /* <v,v>*/
  int_neg (t2, t2);      /* -<v,v>*/

  int_lshift (sigma2dbl, sigma2, 1); /* 2*sigma^2 */
  int_set (_sigma2, sigma2);

  _int2mpfr (nom, t2, 128);
  _int2mpfr (denom, sigma2dbl, sigma2dbl->nlimbs * NBITS_LIMB);
  mpfr_div (t3, nom, denom, MPFR_RNDN); /* (-<v,v>) / (2*sigma^2) */
  mpfr_exp (t3, t3, MPFR_RNDN);         /* exp((-<v,v>) / (2*sigma^2)) */

  _int2mpfr (nom, t1, 128);
  _int2mpfr (denom, _sigma2, sigma2->nlimbs * NBITS_LIMB);
  mpfr_div (t4, nom, denom, MPFR_RNDN); /* (<z,v>) / (sigma^2) */
  mpfr_cosh (t4, t4, MPFR_RNDN);        /* cosh((<z,v>) / (sigma^2)) */

  _int2mpfr (t5, Mu, 128);
  mpfr_mul (t5, t5, t3, MPFR_RNDN);
  mpfr_mul (t5, t5, t4, MPFR_RNDN);

  mpfr_set_ui_2exp (t6, 1, 256, MPFR_RNDN);

  if (mpfr_cmp (t5, t6) > 0)
    reject = 1;
  else
    reject = 0;

  mpfr_clear (t3);
  mpfr_clear (t4);
  mpfr_clear (t6);
  return reject;
}
