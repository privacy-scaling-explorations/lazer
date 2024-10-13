#include "brandom.h"
#include "dom.h"
#include "grandom.h"
#include "lazer.h"
#include "memory.h"
#include "poly.h"
#include "urandom.h"

#include <string.h>

void
polymat_alloc (polymat_ptr r, const polyring_t ring, unsigned int nrows,
               unsigned int ncols)
{
  void *mem;

  mem = _alloc (_sizeof_polymat_data (ring, nrows, ncols));

  _polymat_init (r, ring, nrows, ncols, mem);
  r->mem = mem;
}

void
polymat_free (polymat_ptr r)
{
  unsigned int i, j;
  poly_ptr poly;

  if (r == NULL)
    return;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    poly = polymat_get_elem (r, i, j);
    _free (poly->crtrep, _sizeof_crtrep_data (poly->ring));
  }

  _free (r->mem, _sizeof_polymat_data (r->ring, r->nrows, r->ncols));
}

int
polymat_is_upperdiag (polymat_t a)
{
  if (a == NULL)
    return 1;

  unsigned int i, j;
  poly_ptr elem;
  polyring_srcptr Rq = polymat_get_ring (a);
  const unsigned int nrows = polymat_get_nrows (a);
  const unsigned int ncols = polymat_get_ncols (a);
  poly_t zero;

  poly_alloc (zero, Rq);

  poly_set_zero (zero);
  polymat_fromcrt (a);

  for (i = 1; i < nrows; i++)
    {
      for (j = 0; j < MIN (i, ncols); j++)
        {
          elem = polymat_get_elem (a, i, j);
          if (!(elem->flags & FZERO) && !poly_eq (elem, zero))
            return 0;
        }
    }

  poly_free (zero);
  return 1;
}

void
polymat_subdiags_set_zero (polymat_t r)
{
  const unsigned int nrows = polymat_get_nrows (r);
  polyvec_t subv;
  unsigned int diag;

  polymat_fromcrt (r);

  for (diag = 1; diag < nrows; diag++)
    {
      polymat_get_diag (subv, r, -diag);
      polyvec_set_zero (subv);
    }
  ASSERT_ERR (polymat_is_upperdiag (r));
}

void
polymat_auto (polymat_t r, polymat_t a)
{
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    poly_auto (polymat_get_elem (r, i, j), polymat_get_elem (a, i, j));
  }
}

void
polymat_tocrt (polymat_t r)
{
  unsigned int i, j;
  poly_ptr ri;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    poly_tocrt (ri);
  }
}

void
polymat_tocrtdiag (polymat_t r)
{
  unsigned int i, j;
  poly_ptr ri;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    poly_tocrt (ri);
  }
}

void
polymat_fromcrt (polymat_t r)
{
  unsigned int i, j;
  poly_ptr ri;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    poly_fromcrt (ri);
  }
}

void
polymat_fromcrtdiag (polymat_t r)
{
  unsigned int i, j;
  poly_ptr ri;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    poly_fromcrt (ri);
  }
}

void
polymat_add (polymat_t r, polymat_t a, polymat_t b, int crt)
{
  poly_ptr pr, pa, pb;
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    pr = polymat_get_elem (r, i, j);
    pa = polymat_get_elem (a, i, j);
    pb = polymat_get_elem (b, i, j);
    poly_add (pr, pa, pb, crt);
  }
}

void
polymat_adddiag (polymat_t r, polymat_t a, polymat_t b, int crt)
{
  poly_ptr pr, pa, pb;
  unsigned int i, j;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    pr = polymat_get_elem (r, i, j);
    pa = polymat_get_elem (a, i, j);
    pb = polymat_get_elem (b, i, j);
    poly_add (pr, pa, pb, crt);
  }
}

void
polymat_sub (polymat_t r, polymat_t a, polymat_t b, int crt)
{
  poly_ptr pr, pa, pb;
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    pr = polymat_get_elem (r, i, j);
    pa = polymat_get_elem (a, i, j);
    pb = polymat_get_elem (b, i, j);
    poly_sub (pr, pa, pb, crt);
  }
}

void
polymat_subdiag (polymat_t r, polymat_t a, polymat_t b, int crt)
{
  poly_ptr pr, pa, pb;
  unsigned int i, j;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    pr = polymat_get_elem (r, i, j);
    pa = polymat_get_elem (a, i, j);
    pb = polymat_get_elem (b, i, j);
    poly_add (pr, pa, pb, crt);
  }
}

void
polymat_scale (polymat_t r, const int_t a, polymat_t b)
{
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    poly_scale (polymat_get_elem (r, i, j), a, polymat_get_elem (b, i, j));
  }
}

void
polymat_scalediag (polymat_t r, const int_t a, polymat_t b)
{
  unsigned int i, j;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    poly_scale (polymat_get_elem (r, i, j), a, polymat_get_elem (b, i, j));
  }
}

void
polymat_addscale (polymat_t r, const int_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scale (tmp, a, b);
  polymat_add (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_addscalediag (polymat_t r, const int_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scalediag (tmp, a, b);
  polymat_adddiag (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_subscale (polymat_t r, const int_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scale (tmp, a, b);
  polymat_sub (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_subscalediag (polymat_t r, const int_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scalediag (tmp, a, b);
  polymat_subdiag (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_scale2 (polymat_t r, poly_t a, polymat_t b)
{
  poly_ptr ri, bi;
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    bi = polymat_get_elem (b, i, j);
    poly_mul (ri, a, bi);
  }
}

void
polymat_scalediag2 (polymat_t r, poly_t a, polymat_t b)
{
  poly_ptr ri, bi;
  unsigned int i, j;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    bi = polymat_get_elem (b, i, j);
    poly_mul (ri, a, bi);
  }
}

void
polymat_addscale2 (polymat_t r, poly_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scale2 (tmp, a, b);
  polymat_add (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_addscalediag2 (polymat_t r, poly_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scalediag2 (tmp, a, b);
  polymat_adddiag (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_subscale2 (polymat_t r, poly_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scale2 (tmp, a, b);
  polymat_sub (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_subscalediag2 (polymat_t r, poly_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polymat_get_ring (r);
  const unsigned int nrows = polymat_get_nrows (r);
  const unsigned int ncols = polymat_get_ncols (r);
  polymat_t tmp;

  polymat_alloc (tmp, ring, nrows, ncols);

  polymat_scalediag2 (tmp, a, b);
  polymat_subdiag (r, r, tmp, crt);

  polymat_free (tmp);
}

void
polymat_rrot (polymat_t r, polymat_t a, unsigned int n)
{
  poly_ptr ap, rp;
  unsigned int i, j;

  ASSERT_ERR (r->nrows == a->nrows);
  ASSERT_ERR (r->ncols == a->ncols);
  ASSERT_ERR (n < r->ring->d);

  _MAT_FOREACH_ELEM (r, i, j)
  {
    rp = polymat_get_elem (r, i, j);
    ap = polymat_get_elem (a, i, j);
    poly_rrot (rp, ap, n);
  }
}

void
polymat_rrotdiag (polymat_t r, polymat_t a, unsigned int n)
{
  poly_ptr ap, rp;
  unsigned int i, j;

  ASSERT_ERR (r->nrows == a->nrows);
  ASSERT_ERR (r->ncols == a->ncols);
  ASSERT_ERR (n < r->ring->d);

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    rp = polymat_get_elem (r, i, j);
    ap = polymat_get_elem (a, i, j);
    poly_rrot (rp, ap, n);
  }
}

void
polymat_lrot (polymat_t r, polymat_t a, unsigned int n)
{
  poly_ptr ap, rp;
  unsigned int i, j;

  ASSERT_ERR (r->nrows == a->nrows);
  ASSERT_ERR (r->ncols == a->ncols);
  ASSERT_ERR (n < r->ring->d);

  _MAT_FOREACH_ELEM (r, i, j)
  {
    rp = polymat_get_elem (r, i, j);
    ap = polymat_get_elem (a, i, j);
    poly_lrot (rp, ap, n);
  }
}

void
polymat_lrotdiag (polymat_t r, polymat_t a, unsigned int n)
{
  poly_ptr ap, rp;
  unsigned int i, j;

  ASSERT_ERR (r->nrows == a->nrows);
  ASSERT_ERR (r->ncols == a->ncols);
  ASSERT_ERR (n < r->ring->d);

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    rp = polymat_get_elem (r, i, j);
    ap = polymat_get_elem (a, i, j);
    poly_lrot (rp, ap, n);
  }
}

void
polymat_urandom (polymat_t r, const int_t mod, unsigned int log2mod,
                 const uint8_t seed[32], uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  poly_ptr ptr;
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ptr = polymat_get_elem (r, i, j);
    _dom.d32[0] = i * r->ncols + j;
    _poly_urandom (ptr, mod, log2mod, seed, _dom.d64);
  }
}

void
polymat_brandom (polymat_t r, unsigned int k, const uint8_t seed[32],
                 uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  poly_ptr ptr;
  unsigned int i, j;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ptr = polymat_get_elem (r, i, j);
    _dom.d32[0] = i * r->ncols + j;
    _poly_brandom (ptr, k, seed, _dom.d64);
  }
}

void
polymat_mod (polymat_t r, polymat_t a)
{
  unsigned int i, j;
  poly_ptr ri, ai;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    ai = polymat_get_elem (a, i, j);
    poly_mod (ri, ai);
  }
}

void
polymat_moddiag (polymat_t r, polymat_t a)
{
  unsigned int i, j;
  poly_ptr ri, ai;

  _MAT_FOREACH_ELEM_UPPER (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    ai = polymat_get_elem (a, i, j);
    poly_mod (ri, ai);
  }
}

void
polymat_redp (polymat_t r, polymat_t a)
{
  unsigned int i, j;
  poly_ptr ri, ai;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    ai = polymat_get_elem (a, i, j);
    poly_redp (ri, ai);
  }
}

void
polymat_redc (polymat_t r, polymat_t a)
{
  unsigned int i, j;
  poly_ptr ri, ai;

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ri = polymat_get_elem (r, i, j);
    ai = polymat_get_elem (a, i, j);
    poly_redc (ri, ai);
  }
}

size_t
polymat_out_str (FILE *stream, int base, const polymat_t a)
{
  poly_ptr ptr;
  unsigned int i, j;
  size_t nbytes = 0;

  fprintf (stream, "[");
  nbytes += 1;

  _MAT_FOREACH_ROW (a, i)
  {
    fprintf (stream, "[");
    nbytes += 1;

    _MAT_FOREACH_COL (a, j)
    {
      ptr = polymat_get_elem (a, i, j);
      nbytes += poly_out_str (stream, base, ptr);

      if (j + 1 < a->ncols)
        {
          fprintf (stream, ",");
          nbytes += 1;
        }
    }

    fprintf (stream, "]");
    nbytes += 1;

    if (i + 1 < a->nrows)
      {
        fprintf (stream, ",");
        nbytes += 1;
      }
  }

  fprintf (stream, "]");
  nbytes += 1;

  return nbytes;
}

void
polymat_dump (polymat_t mat)
{
  polymat_out_str (stdout, 10, mat);
  fprintf (stdout, "\n");
  fflush (stdout);
}
