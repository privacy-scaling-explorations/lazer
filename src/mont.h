#ifndef MONT_H
#define MONT_H
#include "lazer.h"
#include <stdint.h>

/*
 * reduce product
 * c in [-2^31*p,2^31*p]
 * returns t = c*2^(-32) mod p, t in (-p,p)
 */
static inline crtcoeff_t
_mont_rprod (crtcoeff_dbl_t c, crtcoeff_t p, crtcoeff_t pinv)
{
  crtcoeff_t t;

  ASSERT_ERR (c >= -((crtcoeff_dbl_t)1 << (CRTCOEFF_NBITS - 1)) * p);
  ASSERT_ERR (c <= ((crtcoeff_dbl_t)1 << (CRTCOEFF_NBITS - 1)) * p);

  t = (crtcoeff_dbl_t)(crtcoeff_t)c * pinv;
  t = (c - (crtcoeff_dbl_t)t * p) >> CRTCOEFF_NBITS;

  ASSERT_ERR (t > -p);
  ASSERT_ERR (t < p);
  return t;
}

/*
 * reduce sum
 * c in (-2p,2p)
 * returns t = c mod p, t in (-p,p)
 */
static inline crtcoeff_t
_mont_rsum (crtcoeff_t c, crtcoeff_t p)
{
  crtcoeff_t abs_c, t, maskc, maskt;

  ASSERT_ERR (c > -((crtcoeff_dbl_t)2 * p));
  ASSERT_ERR (c < ((crtcoeff_dbl_t)2 * p));

  maskc = c >> (CRTCOEFF_NBITS - 1);
  abs_c = (maskc & -c) | ((~maskc) & c);

  t = abs_c - p;
  maskt = t >> (CRTCOEFF_NBITS - 1);

  t = ((maskt & abs_c) | ((~maskt) & t));
  t = (maskc & -t) | ((~maskc) & t);

  ASSERT_ERR (t == c % p);
  return t;
}

#endif
