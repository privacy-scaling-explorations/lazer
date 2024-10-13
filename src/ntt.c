#include "lazer.h"
#include "mont.h"
#include <stdint.h>

static void
_ntt (crtcoeff_t c[], const modulus_t _p, const unsigned int d)
{
  const crtcoeff_t *roots = _p->roots;
  const crtcoeff_t pinv = _p->mont_pinv;
  const crtcoeff_t rr = _p->mont_redr;
  const crtcoeff_t p = _p->p;
  unsigned int len, start, j, k;
  crtcoeff_t root, t;

  /* to montgomery domain */
  for (j = 0; j < d; ++j)
    c[j] = _mont_rprod ((crtcoeff_dbl_t)rr * c[j], p, pinv);

  k = 0;
  for (len = (d >> 1); len > 0; len >>= 1)
    {
      for (start = 0; start < d; start = j + len)
        {
          root = roots[++k];
          for (j = start; j < start + len; ++j)
            {
              t = _mont_rprod ((crtcoeff_dbl_t)root * c[j + len], p, pinv);
              c[j + len] = c[j] - t;
              c[j] = c[j] + t;
            }
        }
    }

  /* XXX reduce to (-p,p) */
  for (j = 0; j < d; ++j)
    c[j] = _mont_rprod ((crtcoeff_dbl_t)roots[0] * c[j], p, pinv);
}

/*
 * Input and output coeffs in (-p,p).
 */
static void
_intt (crtcoeff_t c[], const modulus_t _p, const unsigned int d)
{
  const crtcoeff_t *roots = _p->roots;
  const crtcoeff_t pinv = _p->mont_pinv;
  const crtcoeff_t p = _p->p;
  const crtcoeff_t intt_const = _p->intt_const;
  unsigned int start, len, j, k;
  crtcoeff_t t, root;

#if ASSERT == ASSERT_ENABLED
  {
    unsigned int i;

    for (i = 0; i < d; i++)
      {
        ASSERT_ERR (c[i] > -p);
        ASSERT_ERR (c[i] < p);
      }
  }
#endif

  k = d;
  for (len = 1; len < d; len <<= 1)
    {
      for (start = 0; start < d; start = j + len)
        {
          root = -roots[--k];
          for (j = start; j < start + len; ++j)
            {
              t = c[j];
              c[j] = t + c[j + len];
              c[j + len] = t - c[j + len];
              c[j + len] = _mont_rprod ((int64_t)root * c[j + len], p, pinv);
            }
        }
    }

  /* multiply by 1/deg and return from montgomery domain */
  for (j = 0; j < d; ++j)
    c[j] = _mont_rprod ((crtcoeff_dbl_t)intt_const * c[j], p, pinv);

#if ASSERT == ASSERT_ENABLED
  {
    unsigned int i;

    for (i = 0; i < d; i++)
      {
        ASSERT_ERR (c[i] > -p);
        ASSERT_ERR (c[i] < p);
      }
  }
#endif
}
