#include "brandom.h"
#include "dom.h"
#include "grandom.h"
#include "lazer.h"
#include "memory.h"
#include "poly.h"
#include "urandom.h"

#include <string.h>

void
polyvec_alloc (polyvec_ptr r, const polyring_t ring, unsigned int nelems)
{
  void *mem;

  mem = _alloc (_sizeof_polyvec_data (ring, nelems));

  _polyvec_init (r, ring, nelems, mem);
  r->mem = mem;
}

void
polyvec_free (polyvec_ptr r)
{
  poly_ptr poly;
  unsigned int i;

  if (r == NULL)
    return;

  _VEC_FOREACH_ELEM (r, i)
  {
    poly = polyvec_get_elem (r, i);
    _free (poly->crtrep, _sizeof_crtrep_data (poly->ring));
  }

  _free (r->mem, _sizeof_polyvec_data (r->ring, r->nelems));
}

void
polyvec_get_subvec (polyvec_t subvec, const polyvec_t vec, unsigned int elem,
                    unsigned int nelems, unsigned int stride)
{
  ASSERT_ERR (1 + FLOOR (vec->nelems - elem, stride) >= nelems);

  subvec->elems = polyvec_get_elem (vec, elem);

  subvec->nelems = nelems;
  subvec->stride_elems = stride * vec->stride_elems;
  subvec->ring = vec->ring;
  subvec->mem = NULL;
}

int
polyvec_eq (polyvec_t a, polyvec_t b)
{
  poly_ptr aptr, bptr;
  unsigned int neq = 0;
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);

  _VEC_FOREACH_ELEM (a, i)
  {
    aptr = polyvec_get_elem (a, i);
    bptr = polyvec_get_elem (b, i);

    neq |= (1 ^ poly_eq (aptr, bptr));
  }
  return 1 ^ neq;
}

void
polyvec_rshift (polyvec_t r, polyvec_t a, unsigned int n)
{
  unsigned int i;

  ASSERT_ERR (r->nelems == a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_rshift (polyvec_get_elem (r, i), polyvec_get_elem (a, i), n);
  }
}

void
polyvec_lshift (polyvec_t r, polyvec_t a, unsigned int n)
{
  unsigned int i;

  ASSERT_ERR (r->nelems == a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_lshift (polyvec_get_elem (r, i), polyvec_get_elem (a, i), n);
  }
}

void
polyvec_rrot (polyvec_t r, polyvec_t a, unsigned int n)
{
  poly_ptr ap, rp;
  unsigned int i;

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (n < r->ring->d);

  _VEC_FOREACH_ELEM (r, i)
  {
    rp = polyvec_get_elem (r, i);
    ap = polyvec_get_elem (a, i);
    poly_rrot (rp, ap, n);
  }
}

void
polyvec_lrot (polyvec_t r, polyvec_t a, unsigned int n)
{
  poly_ptr ap, rp;
  unsigned int i;

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (n < r->ring->d);

  _VEC_FOREACH_ELEM (r, i)
  {
    rp = polyvec_get_elem (r, i);
    ap = polyvec_get_elem (a, i);
    poly_lrot (rp, ap, n);
  }
}

void
polyvec_add (polyvec_t r, polyvec_t a, polyvec_t b, int crt)
{
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);
  ASSERT_ERR (r->nelems == a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_add (polyvec_get_elem (r, i), polyvec_get_elem (a, i),
              polyvec_get_elem (b, i), crt);
  }
}

void
polyvec_sub (polyvec_t r, polyvec_t a, polyvec_t b, int crt)
{
  unsigned int i;

  ASSERT_ERR (a->nelems == b->nelems);
  ASSERT_ERR (r->nelems == a->nelems);

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_sub (polyvec_get_elem (r, i), polyvec_get_elem (a, i),
              polyvec_get_elem (b, i), crt);
  }
}

void
polyvec_mul (polyvec_t r, polymat_t a, polyvec_t b)
{
  polyvec_t rowi;
  poly_ptr ri;
  unsigned int i;

  ASSERT_ERR (polyvec_get_ring (r) == polyvec_get_ring (b));
  ASSERT_ERR (polyvec_get_ring (r) == polymat_get_ring (a));
  ASSERT_ERR (polyvec_get_nelems (r) == polymat_get_nrows (a));
  ASSERT_ERR (polyvec_get_nelems (b) == polymat_get_ncols (a));

  _MAT_FOREACH_ROW (a, i)
  {
    ri = polyvec_get_elem (r, i);

    polymat_get_row (rowi, a, i);
    polyvec_dot (ri, rowi, b);
  }
}

void
polyvec_mulsparse (polyvec_t r, spolymat_t a, polyvec_t b)
{
  poly_ptr rp, ap, bp;
  unsigned int i, row, col;

  ASSERT_ERR (polyvec_get_ring (r) == polyvec_get_ring (b));
  ASSERT_ERR (polyvec_get_ring (r) == spolymat_get_ring (a));

  _VEC_FOREACH_ELEM (r, i)
  {
    rp = polyvec_get_elem (r, i);

    if (rp->crtrep == NULL)
      {
        ASSERT_ERR (rp->crt == 0);
        rp->crtrep = _alloc (_sizeof_crtrep_data (rp->ring));
      }
    rp->crt = 1;
    memset (rp->crtrep, 0, _sizeof_crtrep_data (rp->ring));
  }

  _SMAT_FOREACH_ELEM (a, i)
  {
    ap = spolymat_get_elem (a, i);
    row = spolymat_get_row (a, i);
    col = spolymat_get_col (a, i);

    bp = polyvec_get_elem (b, col);
    rp = polyvec_get_elem (r, row);

    poly_addmul (rp, ap, bp, 1);
  }
}

void
polyvec_muldiag (polyvec_t r, polymat_t diag, polyvec_t b)
{
  poly_ptr ri, diagi, bi;
  unsigned int i, j;
  const unsigned int ncols = polymat_get_ncols (diag);

  ASSERT_ERR (polyvec_get_ring (r) == polyvec_get_ring (b));
  ASSERT_ERR (polyvec_get_ring (r) == polymat_get_ring (diag));
  ASSERT_ERR (polyvec_get_nelems (r) == polymat_get_nrows (diag));
  ASSERT_ERR (polyvec_get_nelems (b) == polymat_get_ncols (diag));

  _VEC_FOREACH_ELEM (r, i)
  {
    ri = polyvec_get_elem (r, i);

    if (ri->crtrep == NULL)
      {
        ASSERT_ERR (ri->crt == 0);
        ri->crtrep = _alloc (_sizeof_crtrep_data (ri->ring));
      }
    ri->crt = 1;
    memset (ri->crtrep, 0, _sizeof_crtrep_data (r->ring));

    for (j = i; j < ncols; j++)
      {
        diagi = polymat_get_elem (diag, i, j);
        bi = polyvec_get_elem (b, j);

        poly_addmul (ri, diagi, bi, 1);
      }
  }
}

void
polyvec_mul2 (polyvec_t r, polyvec_t a, polymat_t b)
{
  polyvec_t coli;
  poly_ptr ri;
  unsigned int i;

  ASSERT_ERR (polyvec_get_ring (r) == polyvec_get_ring (a));
  ASSERT_ERR (polyvec_get_ring (r) == polymat_get_ring (b));
  ASSERT_ERR (polyvec_get_nelems (r) == polymat_get_ncols (b));
  ASSERT_ERR (polyvec_get_nelems (a) == polymat_get_nrows (b));

  _MAT_FOREACH_COL (b, i)
  {
    ri = polyvec_get_elem (r, i);

    polymat_get_col (coli, b, i);
    polyvec_dot (ri, a, coli);
  }
}

void
polyvec_muldiag2 (polyvec_t r, polyvec_t a, polymat_t diag)
{
  poly_ptr ri, diagi, ai;
  unsigned int i, j;
  const unsigned int nrows = polymat_get_nrows (diag);

  ASSERT_ERR (polyvec_get_ring (r) == polyvec_get_ring (a));
  ASSERT_ERR (polyvec_get_ring (r) == polymat_get_ring (diag));
  ASSERT_ERR (polyvec_get_nelems (r) == polymat_get_ncols (diag));
  ASSERT_ERR (polyvec_get_nelems (a) == polymat_get_nrows (diag));

  _VEC_FOREACH_ELEM (r, i)
  {
    ri = polyvec_get_elem (r, i);

    if (ri->crtrep == NULL)
      {
        ASSERT_ERR (ri->crt == 0);
        ri->crtrep = _alloc (_sizeof_crtrep_data (ri->ring));
      }
    ri->crt = 1;
    memset (ri->crtrep, 0, _sizeof_crtrep_data (r->ring));

    for (j = 0; j < MIN (i + 1, nrows); j++)
      {
        ai = polyvec_get_elem (a, j);
        diagi = polymat_get_elem (diag, j, i);

        poly_addmul (ri, ai, diagi, 1);
      }
  }
}

void
polyvec_addmul (polyvec_t r, polymat_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_mul (tmp, a, b);
  polyvec_add (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_addmul2 (polyvec_t r, polyvec_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_mul2 (tmp, a, b);
  polyvec_add (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_submul (polyvec_t r, polymat_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_mul (tmp, a, b);
  polyvec_sub (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_submul2 (polyvec_t r, polyvec_t a, polymat_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_mul2 (tmp, a, b);
  polyvec_sub (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_addrshift (polyvec_t r, polyvec_t a, unsigned int n)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_rshift (tmp, a, n);
  polyvec_add (r, r, tmp, 0);

  polyvec_free (tmp);
}

void
polyvec_subrshift (polyvec_t r, polyvec_t a, unsigned int n)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_rshift (tmp, a, n);
  polyvec_sub (r, r, tmp, 0);

  polyvec_free (tmp);
}

void
polyvec_addlshift (polyvec_t r, polyvec_t a, unsigned int n)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_lshift (tmp, a, n);
  polyvec_add (r, r, tmp, 0);

  polyvec_free (tmp);
}

void
polyvec_sublshift (polyvec_t r, polyvec_t a, unsigned int n)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_lshift (tmp, a, n);
  polyvec_sub (r, r, tmp, 0);

  polyvec_free (tmp);
}

void
polyvec_scale (polyvec_t r, const int_t a, polyvec_t b)
{
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_scale (polyvec_get_elem (r, i), a, polyvec_get_elem (b, i));
  }
}

void
polyvec_addscale (polyvec_t r, const int_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_scale (tmp, a, b);
  polyvec_add (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_subscale (polyvec_t r, const int_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_scale (tmp, a, b);
  polyvec_sub (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_scale2 (polyvec_t r, poly_t a, polyvec_t b)
{
  poly_ptr ri, bi;
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    ri = polyvec_get_elem (r, i);
    bi = polyvec_get_elem (b, i);
    poly_mul (ri, a, bi);
  }
}

void
polyvec_addscale2 (polyvec_t r, poly_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_scale2 (tmp, a, b);
  polyvec_add (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_subscale2 (polyvec_t r, poly_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = polyvec_get_ring (r);
  const unsigned int nelems = polyvec_get_nelems (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, nelems);

  polyvec_scale2 (tmp, a, b);
  polyvec_sub (r, r, tmp, crt);

  polyvec_free (tmp);
}

void
polyvec_mod (polyvec_t r, polyvec_t a)
{
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_mod (polyvec_get_elem (r, i), polyvec_get_elem (a, i));
  }
}

void
polyvec_redc (polyvec_t r, polyvec_t a)
{
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_redc (polyvec_get_elem (r, i), polyvec_get_elem (a, i));
  }
}

void
polyvec_redp (polyvec_t r, polyvec_t a)
{
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_redp (polyvec_get_elem (r, i), polyvec_get_elem (a, i));
  }
}

void
polyvec_tocrt (polyvec_t r)
{
  unsigned int i;
  poly_ptr ri;

  _VEC_FOREACH_ELEM (r, i)
  {
    ri = polyvec_get_elem (r, i);
    poly_tocrt (ri);
  }
}

void
polyvec_fromcrt (polyvec_t r)
{
  unsigned int i;
  poly_ptr ri;

  _VEC_FOREACH_ELEM (r, i)
  {
    ri = polyvec_get_elem (r, i);
#if DEBUGINFO == DEBUGINFO_ENABLED
    if (ri->crt == 1)
      DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt polvec_fromcrt");
#endif
    poly_fromcrt (ri);
  }
}

void
polyvec_toisoring (polyvec_t vec, polyvec_t a)
{
  polyring_srcptr Rprime = polyvec_get_ring (a);
  polyring_srcptr R = polyvec_get_ring (vec);
  const unsigned int dprime = polyring_get_deg (Rprime);
  const unsigned int d = polyring_get_deg (R);
  const unsigned int k = dprime / d;
  unsigned int i;
  polyvec_t subv;
  poly_ptr poly;

  ASSERT_ERR (d * k == dprime);
  ASSERT_ERR (polyvec_get_nelems (vec) == k * polyvec_get_nelems (a));
  // XXXASSERT_ERR (int_eq (polyring_get_mod (Rprime), polyring_get_mod (R)) ==
  // 1);

  _VEC_FOREACH_ELEM (a, i)
  {
    poly = polyvec_get_elem (a, i);
    polyvec_get_subvec (subv, vec, i * k, k, 1);
    poly_toisoring (subv, poly);
  }
}

void
polyvec_fromisoring (polyvec_t a, polyvec_t vec)
{
  polyring_srcptr Rprime = polyvec_get_ring (a);
  polyring_srcptr R = polyvec_get_ring (vec);
  const unsigned int dprime = polyring_get_deg (Rprime);
  const unsigned int d = polyring_get_deg (R);
  const unsigned int k = dprime / d;
  unsigned int i;
  polyvec_t subv;
  poly_ptr poly;

  ASSERT_ERR (d * k == dprime);
  ASSERT_ERR (polyvec_get_nelems (vec) == k * polyvec_get_nelems (a));
  // XXXASSERT_ERR (int_eq (polyring_get_mod (Rprime), polyring_get_mod (R)) ==
  // 1);

  _VEC_FOREACH_ELEM (a, i)
  {
    poly = polyvec_get_elem (a, i);
    polyvec_get_subvec (subv, vec, i * k, k, 1);
    poly_fromisoring (poly, subv);
  }
}

void
polyvec_auto_self (polyvec_t r)
{
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i) { poly_auto_self (polyvec_get_elem (r, i)); }
}

void
polyvec_auto (polyvec_t r, polyvec_t a)
{
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    poly_auto (polyvec_get_elem (r, i), polyvec_get_elem (a, i));
  }
}

void
polyvec_urandom_autostable (polyvec_t r, int64_t bnd, unsigned int log2,
                            const uint8_t seed[32], uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    _dom.d32[0]++;
    _poly_urandom_autostable (polyvec_get_elem (r, i), bnd, log2, seed,
                              _dom.d64);
  }
}

void
polyvec_elem_mul (polyvec_t r, polyvec_t a, polyvec_t b)
{
  poly_ptr rp, ap, bp;
  unsigned int i;

  ASSERT_ERR (r->ring == b->ring);
  ASSERT_ERR (a->ring == b->ring);

  _VEC_FOREACH_ELEM (a, i)
  {
    rp = polyvec_get_elem (r, i);
    ap = polyvec_get_elem (a, i);
    bp = polyvec_get_elem (b, i);

    if (rp->crtrep == NULL)
      {
        ASSERT_ERR (rp->crt == 0);
        rp->crtrep = _alloc (_sizeof_crtrep_data (r->ring));
      }
    rp->crt = 1;
    memset (rp->crtrep, 0, _sizeof_crtrep_data (r->ring));

    poly_mul (rp, ap, bp);
  }
}

void
polyvec_dot (poly_t r, polyvec_t a, polyvec_t b)
{
  poly_ptr ap, bp;
  unsigned int i;
  poly_t prod;

  poly_alloc (prod, a->ring);

  ASSERT_ERR (a->ring == b->ring);

  if (r->crtrep == NULL)
    {
      ASSERT_ERR (r->crt == 0);
      r->crtrep = _alloc (_sizeof_crtrep_data (r->ring));
    }
  r->crt = 1;
  memset (r->crtrep, 0, _sizeof_crtrep_data (r->ring));

  _VEC_FOREACH_ELEM (a, i)
  {
    ap = polyvec_get_elem (a, i);
    bp = polyvec_get_elem (b, i);
    poly_mul (prod, ap, bp);
    poly_add (r, r, prod, 1);
  }

  poly_free (prod);
}

void
polyvec_dot2 (poly_t r, spolyvec_t a, polyvec_t b)
{
  poly_ptr ap, bp;
  unsigned int i, elem;
  poly_t prod;

  poly_alloc (prod, a->ring);

  ASSERT_ERR (a->ring == b->ring);

  if (r->crtrep == NULL)
    {
      ASSERT_ERR (r->crt == 0);
      r->crtrep = _alloc (_sizeof_crtrep_data (r->ring));
    }
  r->crt = 1;
  memset (r->crtrep, 0, _sizeof_crtrep_data (r->ring));

  _SVEC_FOREACH_ELEM (a, i)
  {
    ap = spolyvec_get_elem (a, i);
    elem = spolyvec_get_elem_ (a, i);
    bp = polyvec_get_elem (b, elem);
    poly_mul (prod, ap, bp);
    poly_add (r, r, prod, 1);
  }

  poly_free (prod);
}

void
polyvec_dcompress_power2round (polyvec_t r, polyvec_t a,
                               const dcompress_params_t params)
{
  unsigned int i;
  poly_ptr rp, ap;

  ASSERT_ERR (r->ring == a->ring);

  _VEC_FOREACH_ELEM (a, i)
  {
    rp = polyvec_get_elem (r, i);
    ap = polyvec_get_elem (a, i);
    poly_dcompress_power2round (rp, ap, params);
  }
}

void
polyvec_dcompress_decompose (polyvec_t r1, polyvec_t r0, polyvec_t r,
                             const dcompress_params_t params)
{
  unsigned int i;
  poly_ptr r1p, r0p, rp;

  ASSERT_ERR (r->ring == r1->ring);
  ASSERT_ERR (r->ring == r0->ring);

  _VEC_FOREACH_ELEM (r1, i)
  {
    r1p = polyvec_get_elem (r1, i);
    r0p = polyvec_get_elem (r0, i);
    rp = polyvec_get_elem (r, i);
    poly_dcompress_decompose (r1p, r0p, rp, params);
  }
}

void
polyvec_dcompress_use_ghint (polyvec_t ret, polyvec_t y, polyvec_t r,
                             const dcompress_params_t params)
{
  unsigned int i;
  poly_ptr retp, yp, rp;

  ASSERT_ERR (ret->ring == y->ring);
  ASSERT_ERR (ret->ring == r->ring);

  _VEC_FOREACH_ELEM (ret, i)
  {
    retp = polyvec_get_elem (ret, i);
    yp = polyvec_get_elem (y, i);
    rp = polyvec_get_elem (r, i);
    poly_dcompress_use_ghint (retp, yp, rp, params);
  }
}

void
polyvec_dcompress_make_ghint (polyvec_t ret, polyvec_t z, polyvec_t r,
                              const dcompress_params_t params)
{
  unsigned int i;
  poly_ptr retp, zp, rp;

  ASSERT_ERR (ret->ring == z->ring);
  ASSERT_ERR (ret->ring == r->ring);

  _VEC_FOREACH_ELEM (ret, i)
  {
    retp = polyvec_get_elem (ret, i);
    rp = polyvec_get_elem (r, i);
    zp = polyvec_get_elem (z, i);
    poly_dcompress_make_ghint (retp, zp, rp, params);
  }
}

/*
 * not ctime.
 */
void
polyvec_linf (int_t r, polyvec_t a)
{
  INT_T (tmp, r->nlimbs);
  poly_ptr ap;
  unsigned int i;

  ASSERT_ERR (r->nlimbs == a->ring->q->nlimbs);

  int_set_i64 (r, 0);

  _VEC_FOREACH_ELEM (a, i)
  {
    ap = polyvec_get_elem (a, i);
    poly_linf (tmp, ap);
    if (int_absgt (tmp, r))
      int_set (r, tmp);
  }
}

void
polyvec_l2sqr (int_t r, polyvec_t a)
{
  INT_T (tmp, r->nlimbs);
  poly_ptr ap;
  unsigned int i;

  ASSERT_ERR (r->nlimbs == 2 * a->ring->q->nlimbs);

  int_set_i64 (r, 0);

  _VEC_FOREACH_ELEM (a, i)
  {
    ap = polyvec_get_elem (a, i);
    poly_l2sqr (tmp, ap);
    int_add (r, r, tmp);
  }
}

void
polyvec_grandom (polyvec_t r, unsigned int log2o, const uint8_t seed[32],
                 uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  unsigned int i;
  poly_ptr ri;

  _VEC_FOREACH_ELEM (r, i)
  {
    _dom.d32[0]++;
    ri = polyvec_get_elem (r, i);
    _poly_grandom (ri, log2o, seed, _dom.d64);
  }
}

void
polyvec_brandom (polyvec_t r, unsigned int k, const uint8_t seed[32],
                 uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    _dom.d32[0]++;
    _poly_brandom (polyvec_get_elem (r, i), k, seed, _dom.d64);
  }
}

void
polyvec_urandom (polyvec_t r, const int_t mod, unsigned int log2mod,
                 const uint8_t seed[32], uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    _dom.d32[0]++;
    _poly_urandom (polyvec_get_elem (r, i), mod, log2mod, seed, _dom.d64);
  }
}

void
polyvec_urandom_bnd (polyvec_t r, const int_t lo, const int_t hi,
                     const uint8_t seed[32], uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  unsigned int i;

  _VEC_FOREACH_ELEM (r, i)
  {
    _dom.d32[0]++;
    _poly_urandom_bnd (polyvec_get_elem (r, i), lo, hi, seed, _dom.d64);
  }
}

size_t
polyvec_out_str (FILE *stream, int base, polyvec_t a)
{
  poly_ptr ptr;
  unsigned int i;
  size_t nbytes = 0;

  fprintf (stream, "(");
  nbytes += 1;

  _VEC_FOREACH_ELEM (a, i)
  {
    ptr = polyvec_get_elem (a, i);
    nbytes += poly_out_str (stream, base, ptr);

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
polyvec_dump (polyvec_t vec)
{
  polyvec_out_str (stdout, 10, vec);
  fprintf (stdout, "\n");
  fflush (stdout);
}
