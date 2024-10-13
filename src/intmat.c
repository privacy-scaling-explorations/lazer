#include "brandom.h"
#include "dom.h"
#include "lazer.h"
#include "memory.h"
#include "rng.h"
#include "urandom.h"

#include <string.h>

void
intmat_alloc (intmat_ptr r, unsigned int nrows, unsigned int ncols,
              unsigned int nlimbs)
{
  void *mem;

  mem = _alloc (_sizeof_intmat_data (nrows, ncols, nlimbs));

  _intmat_init (r, nrows, ncols, nlimbs, mem);
}

void
intmat_free (intmat_ptr r)
{
  if (r == NULL)
    return;

  _free (r->bytes, _sizeof_intmat_data (r->nrows, r->ncols, r->nlimbs));
}

int
intmat_eq (const intmat_t a, const intmat_t b)
{
  int_srcptr aptr, bptr;
  unsigned int eq = 0;
  unsigned int i, j;

  ASSERT_ERR (a->nrows == b->nrows);
  ASSERT_ERR (a->ncols == b->ncols);

  _MAT_FOREACH_ELEM (a, i, j)
  {
    aptr = intmat_get_elem_src (a, i, j);
    bptr = intmat_get_elem_src (b, i, j);

    eq |= (1 ^ int_eq (aptr, bptr));
  }
  return 1 ^ eq;
}

void
intmat_mul_sgn_self (intmat_t r, int sgn)
{
  int_ptr ptr;
  unsigned int i, j;

  ASSERT_ERR (sgn == 1 || sgn == -1);

  _MAT_FOREACH_ELEM (r, i, j)
  {
    ptr = intmat_get_elem (r, i, j);
    int_mul_sgn_self (ptr, sgn);
  }
}

void
intmat_brandom (intmat_t m, unsigned int k, const uint8_t seed[32],
                uint32_t dom)
{
  const unsigned int nelems = m->nrows * m->ncols;
  int8_t vec[nelems];

  _brandom (vec, nelems, k, seed, dom);
  intmat_set_i8 (m, vec);
}

void
intmat_urandom (intmat_t r, const int_t mod, unsigned int log2mod,
                const uint8_t seed[32], uint32_t dom)
{
  union dom _dom = { { 0, dom } };
  intvec_t rowi;
  unsigned int i;

  _MAT_FOREACH_ROW (r, i)
  {
    intmat_get_row (rowi, r, i);
    _dom.d32[0] = i * r->ncols;

    _urandom (rowi, mod, log2mod, seed, _dom.d64);
  }
}

size_t
intmat_out_str (FILE *stream, int base, const intmat_t a)
{
  int_srcptr ptr;
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
      ptr = intmat_get_elem_src (a, i, j);
      nbytes += int_out_str (stream, base, ptr);

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
intmat_dump (intmat_t mat)
{
  intmat_out_str (stdout, 10, mat);
  fprintf (stdout, "\n");
  fflush (stdout);
}

void
intmat_clear (intmat_t r)
{
  explicit_bzero (r->bytes, r->nbytes);
  memset (r, 0, sizeof (intmat_t));
}
