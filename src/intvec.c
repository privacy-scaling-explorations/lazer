#include "intvec.h"
#include "brandom.h"
#include "grandom.h"
#include "lazer.h"
#include "memory.h"
#include "urandom.h"

#include <string.h>

static union
{
  unsigned int uint;
  unsigned char little;
} _endian = { .uint = 1 };

/* little-endian increment */
static inline void
_inc32 (uint32_t *ctr)
{
  if (_endian.little)
    {
      (*ctr)++;
    }
  else
    {
      unsigned int i;
      uint8_t *ptr = (uint8_t *)ctr;
      uint16_t c = 1;

      for (i = 0; i < 4; i++)
        {
          c += ptr[i];
          ptr[i] = (uint8_t)c;
          c >>= 8;
        }
    }
}

void
intvec_alloc (intvec_ptr r, unsigned int nelems, unsigned int nlimbs)
{
  void *mem;

  mem = _alloc (_sizeof_intvec_data (nelems, nlimbs));

  _intvec_init (r, nelems, nlimbs, mem);
}

void
intvec_free (intvec_ptr r)
{
  if (r == NULL)
    return;

  _free (r->bytes, _sizeof_intvec_data (r->nelems, r->nlimbs));
}

void
intvec_get_subvec (intvec_t subvec, const intvec_t vec, unsigned int elem,
                   unsigned int nelems, unsigned int stride)
{
  ASSERT_ERR (1 + FLOOR (vec->nelems - elem, stride) >= nelems);

  subvec->bytes = vec->bytes;
  subvec->nbytes = vec->nbytes;

  subvec->elems = intvec_get_elem (vec, elem);

  subvec->nelems = nelems;
  subvec->stride_elems = stride * vec->stride_elems;
  subvec->nlimbs = vec->nlimbs;
}

int
intvec_eq (const intvec_t a, const intvec_t b)
{
  int_srcptr aptr, bptr;
  unsigned int neq = 0;
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);

  _VEC_FOREACH_ELEM (a, i)
  {
    aptr = intvec_get_elem_src (a, i);
    bptr = intvec_get_elem_src (b, i);

    neq |= (1 ^ int_eq (aptr, bptr));
  }
  return 1 ^ neq;
}

int
intvec_lt (const intvec_t a, const int_t m)
{
  int_srcptr aptr;
  unsigned int ge = 0;
  unsigned int i;

  _VEC_FOREACH_ELEM (a, i)
  {
    aptr = intvec_get_elem_src (a, i);

    ge |= (1 ^ int_lt (aptr, m));
  }
  return 1 ^ ge;
}

int
intvec_gt (const intvec_t a, const int_t m)
{
  int_srcptr aptr;
  unsigned int le = 0;
  unsigned int i;

  _VEC_FOREACH_ELEM (a, i)
  {
    aptr = intvec_get_elem_src (a, i);

    le |= (1 ^ int_gt (aptr, m));
  }
  return 1 ^ le;
}

int
intvec_le (const intvec_t a, const int_t m)
{
  int_srcptr aptr;
  unsigned int gt = 0;
  unsigned int i;

  _VEC_FOREACH_ELEM (a, i)
  {
    aptr = intvec_get_elem_src (a, i);

    gt |= (1 ^ int_le (aptr, m));
  }
  return 1 ^ gt;
}

int
intvec_ge (const intvec_t a, const int_t m)
{
  int_srcptr aptr;
  unsigned int lt = 0;
  unsigned int i;

  _VEC_FOREACH_ELEM (a, i)
  {
    aptr = intvec_get_elem_src (a, i);

    lt |= (1 ^ int_ge (aptr, m));
  }
  return 1 ^ lt;
}

void
intvec_rshift (intvec_t r, const intvec_t a, unsigned int n)
{
  unsigned int i;

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_rshift (intvec_get_elem (r, i), intvec_get_elem_src (a, i), n);
  }
}

void
intvec_lshift (intvec_t r, const intvec_t a, unsigned int n)
{
  unsigned int i;

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_lshift (intvec_get_elem (r, i), intvec_get_elem_src (a, i), n);
  }
}

void
intvec_rrot (intvec_t r, const intvec_t a, unsigned int n)
{
  const unsigned int nelems = r->nelems;
  unsigned int i;
  int_ptr t;
  INTVEC_T (tmp, r->nelems, r->nlimbs);

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);
  ASSERT_ERR (n < nelems);

  for (i = 1; i <= n; i++)
    {
      t = intvec_get_elem (tmp, nelems - i);
      int_set (t, intvec_get_elem_src (a, n - i));
      int_neg (t, t);
    }
  for (i = 0; i < nelems - n; i++)
    {
      intvec_set_elem (tmp, i, intvec_get_elem_src (a, i + n));
    }

  intvec_set (r, tmp);
}

void
intvec_lrot (intvec_t r, const intvec_t a, unsigned int n)
{
  const unsigned int nelems = r->nelems;
  unsigned int i;
  int_ptr t;
  INTVEC_T (tmp, r->nelems, r->nlimbs);

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);
  ASSERT_ERR (n < nelems);

  for (i = 1; i <= n; i++)
    {
      t = intvec_get_elem (tmp, n - i);
      int_set (t, intvec_get_elem_src (a, nelems - i));
      int_neg_self (t);
    }
  for (i = n; i < nelems; i++)
    {
      intvec_set_elem (tmp, i, intvec_get_elem_src (a, i - n));
    }

  intvec_set (r, tmp);
}

void
intvec_add (intvec_t r, const intvec_t a, const intvec_t b)
{
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);
  ASSERT_ERR (r->nelems == a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_add (intvec_get_elem (r, i), intvec_get_elem_src (a, i),
             intvec_get_elem_src (b, i));
  }
}

void
intvec_sub (intvec_t r, const intvec_t a, const intvec_t b)
{
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);
  ASSERT_ERR (r->nelems == a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_sub (intvec_get_elem (r, i), intvec_get_elem_src (a, i),
             intvec_get_elem_src (b, i));
  }
}

void
intvec_mul (intvec_t r, const intvec_t a, const intvec_t b)
{
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);
  ASSERT_ERR (r->nelems == 2 * a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_mul (intvec_get_elem (r, i), intvec_get_elem_src (a, i),
             intvec_get_elem_src (b, i));
  }
}

void
intvec_mul_matvec (intvec_t r, const intmat_t mat, const intvec_t vec)
{
  unsigned int i, j;
  int_srcptr mp, vp;
  int_ptr rp;
  INT_T (prod, r->nlimbs);

  ASSERT_ERR (vec->nlimbs == mat->nlimbs);
  ASSERT_ERR (r->nlimbs == 2 * mat->nlimbs);
  ASSERT_ERR (r->nelems == mat->nrows);
  ASSERT_ERR (vec->nelems == mat->ncols);

  _MAT_FOREACH_ROW (mat, i)
  {
    rp = intvec_get_elem (r, i);
    int_set_i64 (rp, 0);

    _MAT_FOREACH_COL (mat, j)
    {
      mp = intmat_get_elem_src (mat, i, j);
      vp = intvec_get_elem_src (vec, j);

      int_mul (prod, mp, vp);
      int_add (rp, rp, prod);
    }
  }
}

void
intvec_div (intvec_t rq, intvec_t rr, const intvec_t a, const intvec_t b)
{
  unsigned int i;

  ASSERT_ERR (rq->nelems == rr->nelems);
  ASSERT_ERR (rr->nelems == a->nelems);
  ASSERT_ERR (a->nelems == b->nelems);

  _VEC_FOREACH_ELEM (rq, i)
  {
    int_div (intvec_get_elem (rq, i), intvec_get_elem (rr, i),
             intvec_get_elem_src (a, i), intvec_get_elem (b, i));
  }
}

void
intvec_scale (intvec_t r, const int_t a, const intvec_t b)
{
  unsigned int i;

  ASSERT_ERR (a->nlimbs == b->nlimbs);
  ASSERT_ERR (r->nlimbs == 2 * a->nlimbs);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_mul (intvec_get_elem (r, i), a, intvec_get_elem_src (b, i));
  }
}

void
intvec_mod (intvec_t r, const intvec_t a, const int_t m)
{
  unsigned int i;

  ASSERT_ERR (r->nlimbs == m->nlimbs);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_mod (intvec_get_elem (r, i), intvec_get_elem_src (a, i), m);
  }
}

void
intvec_redc (intvec_t r, const intvec_t a, const int_t m)
{
  unsigned int i;
  int_ptr ri;
  int_srcptr ai;

  ASSERT_ERR (a->nlimbs == m->nlimbs);
  ASSERT_ERR (r->nlimbs == a->nlimbs);

  _VEC_FOREACH_ELEM (r, i)
  {
    ai = intvec_get_elem_src (a, i);
    ri = intvec_get_elem (r, i);
    int_redc (ri, ai, m);
  }
}

void
intvec_redp (intvec_t r, const intvec_t a, const int_t m)
{
  unsigned int i;

  ASSERT_ERR (a->nlimbs == m->nlimbs);
  ASSERT_ERR (r->nlimbs == a->nlimbs);

  _VEC_FOREACH_ELEM (r, i)
  {
    int_redp (intvec_get_elem (r, i), intvec_get_elem_src (a, i), m);
  }
}

void
intvec_auto_self (intvec_t r)
{
  INT_T (tmp, r->nlimbs);
  int_ptr ptr1, ptr2;
  unsigned int i;

  ASSERT_ERR (r->nelems % 2 == 0);

  ptr1 = intvec_get_elem (r, r->nelems / 2);
  int_neg_self (ptr1);

  for (i = 1; i < r->nelems / 2; i++)
    {
      ptr1 = intvec_get_elem (r, i);
      ptr2 = intvec_get_elem (r, r->nelems - i);

      int_set (tmp, ptr1);
      int_neg (ptr1, ptr2);
      int_neg (ptr2, tmp);
    }
}

void
intvec_mul_sgn_self (intvec_t r, int sgn)
{
  int_ptr ptr;
  unsigned int i;

  ASSERT_ERR (sgn == 1 || sgn == -1);

  _VEC_FOREACH_ELEM (r, i)
  {
    ptr = intvec_get_elem (r, i);
    int_mul_sgn_self (ptr, sgn);
  }
}

void
intvec_auto (intvec_t r, const intvec_t a)
{
  int_srcptr src1, src2;
  int_ptr dst1, dst2;
  unsigned int i;

  ASSERT_ERR (r->elems != a->elems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);
  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nelems % 2 == 0);

  src1 = intvec_get_elem_src (a, 0);
  intvec_set_elem (r, 0, src1);

  src1 = intvec_get_elem_src (a, r->nelems / 2);
  dst1 = intvec_get_elem (r, r->nelems / 2);
  int_neg (dst1, src1);

  for (i = 1; i < r->nelems / 2; i++)
    {
      src1 = intvec_get_elem_src (a, i);
      dst1 = intvec_get_elem (r, r->nelems - i);

      src2 = intvec_get_elem_src (a, r->nelems - i);
      dst2 = intvec_get_elem (r, i);

      int_neg (dst1, src1);
      int_neg (dst2, src2);
    }
}

void
intvec_dot (int_t r, const intvec_t a, const intvec_t b)
{
  int_srcptr ai, bi;
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);
  ASSERT_ERR (r->nlimbs == 2 * a->nlimbs);

  int_set_i64 (r, 0);

  _VEC_FOREACH_ELEM (a, i)
  {
    ai = intvec_get_elem_src (a, i);
    bi = intvec_get_elem_src (b, i);
    int_addmul (r, ai, bi);
  }
}

void
intvec_l2sqr (int_t r, const intvec_t a)
{
  int_srcptr ai;
  unsigned int i;

  ASSERT_ERR (r->nlimbs == 2 * a->nlimbs);

  int_set_i64 (r, 0);

  _VEC_FOREACH_ELEM (a, i)
  {
    ai = intvec_get_elem_src (a, i);
    int_addsqr (r, ai);
  }
}

/*
 * not ctime.
 */
void
intvec_linf (int_t r, const intvec_t a)
{
  int_srcptr ai, max;
  unsigned int i;

  ASSERT_ERR (r->nlimbs == a->nlimbs);

  max = intvec_get_elem_src (a, 0);

  for (i = 1; i < a->nelems; i++)
    {
      ai = intvec_get_elem_src (a, i);
      if (int_absgt (ai, max))
        max = ai;
    }

  int_set (r, max);
  r->neg = 0;
}

void
intvec_urandom_autostable (intvec_t r, int64_t bnd, unsigned int log2,
                           const uint8_t seed[32], uint32_t dom)
{
  _intvec_urandom_autostable (r, bnd, log2, seed, dom);
}

void
intvec_grandom (intvec_t r, unsigned int log2o, const uint8_t seed[32],
                uint32_t dom)
{
  _intvec_grandom (r, log2o, seed, dom);
}

void
intvec_brandom (intvec_t r, unsigned int k, const uint8_t seed[32],
                uint32_t dom)
{
  _intvec_brandom (r, k, seed, dom);
}

void
intvec_urandom (intvec_t r, const int_t mod, unsigned int log2mod,
                const uint8_t seed[32], uint32_t dom)
{
  _intvec_urandom (r, mod, log2mod, seed, dom);
}

void
intvec_urandom_bnd (intvec_t r, const int_t lo, const int_t hi,
                    const uint8_t seed[32], uint32_t dom)
{
  _urandom_bnd (r, lo, hi, seed, dom);
}

size_t
intvec_out_str (FILE *stream, int base, const intvec_t a)
{
  int_srcptr ptr;
  unsigned int i;
  size_t nbytes = 0;

  fprintf (stream, "(");
  nbytes += 1;

  _VEC_FOREACH_ELEM (a, i)
  {
    ptr = intvec_get_elem_src (a, i);
    nbytes += int_out_str (stream, base, ptr);

    if (i + 1 < a->nelems)
      {
        fprintf (stream, ",");
        nbytes += 1;
      }
  }

  fprintf (stream, ")");
  nbytes += 1;

  return nbytes;
}

void
intvec_dump (intvec_t vec)
{
  intvec_out_str (stdout, 10, vec);
  fprintf (stdout, "\n");
  fflush (stdout);
}

void
intvec_clear (intvec_t r)
{
  explicit_bzero (r->bytes, r->nbytes);
  memset (r, 0, sizeof (intvec_t));
}

static void
_intvec_urandom_autostable (intvec_t r, int64_t bnd, unsigned int log2,
                            const uint8_t seed[32], uint64_t dom)
{
  const unsigned int nchal_coeffs = r->nelems / 2;
  unsigned int i;
  int64_t chal_coeffs[nchal_coeffs];
  int64_t mod;

  ASSERT_ERR (r->nelems % 2 == 0);

  mod = 2 * bnd + 1;

  /* coeffs in [0, 2*bnd] */
  _urandom_i64 (chal_coeffs, nchal_coeffs, mod, log2, seed, dom);

  for (i = 0; i < nchal_coeffs; i++)
    {
      /* coeffs in [-bnd, bnd] */
      chal_coeffs[i] -= bnd;
      intvec_set_elem_i64 (r, i, chal_coeffs[i]);
    }
  intvec_set_elem_i64 (r, nchal_coeffs, 0);
  for (i = nchal_coeffs + 1; i < r->nelems; i++)
    intvec_set_elem_i64 (r, i, -chal_coeffs[2 * nchal_coeffs - i]);
}

static void
_intvec_urandom (intvec_t r, const int_t mod, unsigned int log2mod,
                 const uint8_t seed[32], uint64_t dom)
{
  _urandom (r, mod, log2mod, seed, dom);
}

static void
_intvec_brandom (intvec_t r, unsigned int k, const uint8_t seed[32],
                 uint64_t dom)
{
  int8_t vec[r->nelems];

  _brandom (vec, r->nelems, k, seed, dom);
  intvec_set_i8 (r, vec);
}

static void
_intvec_grandom (intvec_t r, unsigned int log2o, const uint8_t seed[32],
                 uint64_t dom)
{
  int32_t samples[r->nelems];

  if (LIKELY (log2o < 24))
    {
      _grandom_sample_i32 (samples, r->nelems, log2o, seed, dom);
      intvec_set_i32 (r, samples);
    }
  else
    {
      unsigned int i;

      _VEC_FOREACH_ELEM (r, i)
      {
        _grandom_sample (intvec_get_elem (r, i), log2o, seed, dom);
      }
    }
}
