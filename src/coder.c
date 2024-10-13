#include "lazer.h"

static inline void _inc_idx (unsigned int *byte, unsigned int *bit);
static inline void _inc_idx_zero (uint8_t *buf, unsigned int *byte,
                                  unsigned int *bit);
static unsigned int _uencode (uint8_t **byte, unsigned int *bit,
                              const intvec_t v, UNUSED const int_t m,
                              unsigned int mbits);
static unsigned int _udecode (intvec_t v, const uint8_t **byte,
                              unsigned int *bit, const int_t m,
                              unsigned int mbits);

void
coder_enc_begin (coder_state_t state, uint8_t *out)
{
  state->out = out;
  state->in = NULL;
  state->byte_off = 0;
  state->bit_off = 0;
}

void
coder_dec_begin (coder_state_t state, const uint8_t *in)
{
  state->out = NULL;
  state->in = in;
  state->byte_off = 0;
  state->bit_off = 0;
}

unsigned int
coder_get_offset (coder_state_t state)
{
  return (state->byte_off << 3) + state->bit_off;
}

void
coder_enc_bytes (coder_state_t state, const uint8_t *bytes,
                 unsigned int nbytes)
{
  /* encoding bytes first avoids bit shifting */
  ASSERT_ERR (state->bit_off == 0);

  memcpy (state->out, bytes, nbytes);
  state->byte_off += nbytes;
  state->out += nbytes;
}

void
coder_enc_ghint (coder_state_t state, const intvec_t ghint)
{
  int64_t elem;
  unsigned int i, j, _bit;
  unsigned int nbytes = 0;
  uint8_t *_byte;

  ASSERT_ERR (state->out != NULL);
  ASSERT_ERR (state->in == NULL);

  _byte = state->out;
  _bit = state->bit_off;

  /* zero unset bits in first byte */
  _byte[0] &= ~((uint8_t)(~0) << _bit);

  _VEC_FOREACH_ELEM (ghint, i)
  {
    elem = int_get_i64 (intvec_get_elem_src (ghint, i));
    if (elem == 0)
      {
        _inc_idx_zero (_byte, &nbytes, &_bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
    else if (elem == 1)
      {
        _inc_idx_zero (_byte, &nbytes, &_bit);
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
    else if (elem == -1)
      {
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
    else if (elem >= 2)
      {
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
        for (j = 0; j < 2 * elem - 4; j++)
          {
            _inc_idx_zero (_byte, &nbytes, &_bit);
          }
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
    else
      {
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
        for (j = 0; j < 2 * (-elem) - 3; j++)
          {
            _inc_idx_zero (_byte, &nbytes, &_bit);
          }
        _byte[nbytes] |= (1 << _bit);
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
  }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->out += nbytes;
}

void
coder_dec_ghint (coder_state_t state, intvec_t ghint)
{
  int_ptr elem;
  unsigned int i, _bit;
  unsigned int nbytes = 0;
  const uint8_t *_byte;
  int bit0, bit1, j;

  _byte = state->in;
  _bit = state->bit_off;

  _VEC_FOREACH_ELEM (ghint, i)
  {
    elem = intvec_get_elem (ghint, i);

    bit0 = !!(_byte[nbytes] & (1 << _bit));
    _inc_idx (&nbytes, &_bit);
    bit1 = !!(_byte[nbytes] & (1 << _bit));
    _inc_idx (&nbytes, &_bit);

    if (bit0 == 0 && bit1 == 0)
      {
        int_set_i64 (elem, 0);
      }
    else if (bit0 == 0 && bit1 == 1)
      {
        int_set_i64 (elem, 1);
      }
    else if (bit0 == 1 && bit1 == 0)
      {
        int_set_i64 (elem, -1);
      }
    else
      {
        j = 0;
        while (!(_byte[nbytes] & (1 << _bit)))
          {
            j++;
            _inc_idx (&nbytes, &_bit);
          }
        _inc_idx (&nbytes, &_bit);

        if (j % 2 == 0)
          int_set_i64 (elem, (j + 4) / 2);
        else
          int_set_i64 (elem, -(j + 3) / 2);
      }
  }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->in += nbytes;
}

void
coder_enc_urandom (coder_state_t state, const intvec_t v, const int_t m,
                   unsigned int mbits)
{
  unsigned int nbits;

  ASSERT_ERR (state->out != NULL);
  ASSERT_ERR (state->in == NULL);

  nbits = _uencode (&(state->out), &(state->bit_off), v, m, mbits);
  state->byte_off += (nbits >> 3);
}

/* clang-format off */
/*
 * for split: z = a*sigma*z1 + z0
 * z0 := z mod a*sigma in [-a*sigma/2, a*sigma/2)
 * z1 := (z-z0)/(a*sigma)
 *
 * a   | weighted avg huff bits | bin bits ceil(log(ao)) | weighted avg bits total | code unitary
 * --------------------------------------------------------------------------------------------
 * 8   |                   1.00 |       3 + ceil(log(o)) | 4.00 + ceil(log(o))     | yes
 * 4   |                   1.07 |       2 + ceil(log(o)) | 3.07 + ceil(log(o))     | yes
 * 2   |                   1.48 |       1 + ceil(log(o)) | 2.48 + ceil(log(o))     | yes
 * 1   |                   2.22 |       0 + ceil(log(o)) | 2.22 + ceil(log(o))     | yes
 * 1/2 |                   3.12 |      -1 + ceil(log(o)) | 2.12 + ceil(log(o))     | outside of [-4,4]
 * 1/4 |                   4.09 |      -2 + ceil(log(o)) | 2.09 + ceil(log(o))     | outside of [-17,17]
 * 1/8 |                   5.09 |      -3 + ceil(log(o)) | 2.09 + ceil(log(o))     | outside of [-70,70]
 *
 * Out sigmas are of the form 1.55*2^log2o, so we chose a = 2/1.55 ~ 1.29.
 * For this 1 < a < 2, the huffman encoding of z1 is unitary, so we dont need
 * to store a table.
 * Moreover, the division by a*sigma = 2/1.55 * 1.55*2^log2o = 2^(log2o+1)
 * is just a right-shift by log2o+1 bits.
 *
 * Assume v is centered and its elements are single-precision.
 */
/* clang-format on */
void
coder_enc_grandom (coder_state_t state, const intvec_t v, unsigned int log2o)
{
  const unsigned int log2ao = log2o + 1;
  const long ao = (limb_t)1 << log2ao;
  unsigned int i, j, nones, _bit;
  unsigned int nbytes = 0;
  long z, z0, z1;
  int_srcptr elem;
  uint8_t *_byte;

  ASSERT_ERR (state->out != NULL);
  ASSERT_ERR (state->in == NULL);
  ASSERT_ERR (v->nlimbs == 1);

  _byte = state->out;
  _bit = state->bit_off;

  /* zero unset bits in first byte */
  _byte[0] &= ~((uint8_t)(~0) << _bit);

  _VEC_FOREACH_ELEM (v, i)
  {
    elem = intvec_get_elem_src (v, i);
    z = int_get_i64 (elem);

    z0 = z - ((z & (~(limb_t)0 << log2ao)));
    // z0 = z - ((z >> log2ao) << log2ao);
    /* from (-a*sigma, a*sigma) to [-a*sigma/2, a*sigma/2) */
    if (z0 >= ao / 2)
      z0 -= ao;
    else if (z0 < -ao / 2)
      z0 += ao;

    z1 = (z - z0) >> log2ao;

    /* z1: unitary encoding */

    nones = (z1 <= 0 ? -2 * z1 : 2 * z1 - 1);
    for (j = 0; j < nones; j++)
      {
        _byte[nbytes] |= (1 << _bit); /* ones */
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
    _inc_idx_zero (_byte, &nbytes, &_bit); /* final zero */

    /* z0: binary encoding (two's complement) */

    for (j = 0; j < log2ao; j++)
      {
        _byte[nbytes] |= ((z0 & 1) << _bit);

        z0 >>= 1;
        _inc_idx_zero (_byte, &nbytes, &_bit);
      }
  }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->out += nbytes;
}

void
coder_dec_grandom (coder_state_t state, intvec_t v, unsigned int log2o)
{
  const unsigned int log2ao = log2o + 1;
  const long ao = (limb_t)1 << log2ao;
  long z, z0, z1;
  const uint8_t *_byte;
  unsigned int nbytes = 0;
  unsigned int i, j, _bit, nones;
  int_ptr elem;

  ASSERT_ERR (state->out == NULL);
  ASSERT_ERR (state->in != NULL);

  _byte = state->in;
  _bit = state->bit_off;

  _VEC_FOREACH_ELEM (v, i)
  {
    /* z1: unitary encoding */

    for (nones = 0; (_byte[nbytes] & (1 << _bit)) > 0; nones++)
      _inc_idx (&nbytes, &_bit);
    _inc_idx (&nbytes, &_bit); /* final zero */

    z1 = (nones % 2 == 0 ? -(long)nones / 2 : ((long)nones + 1) / 2);

    /* z0: binary encoding (two's complement) */
    z0 = 0;
    for (j = 0; j < log2ao; j++)
      {
        z0 |= ((((long)_byte[nbytes] & ((long)1 << _bit)) >> _bit) << j);

        _inc_idx (&nbytes, &_bit);
      }

    if ((z0 & (1 << (log2ao - 1))) > 0)
      z0 |= (((unsigned long)~0) << log2ao); /* sign-extend */

    z = ao * z1 + z0;

    elem = intvec_get_elem (v, i);
    int_set_i64 (elem, z);
  }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->in += nbytes;
}

/* Pad with a one. Fill zeros to next byte boundary. */
void
coder_enc_end (coder_state_t state)
{
  unsigned int _bit;
  uint8_t *_byte;

  ASSERT_ERR (state->out != NULL);
  ASSERT_ERR (state->in == NULL);

  _byte = state->out;
  _bit = state->bit_off;

  /* zero unset bits in first byte */
  _byte[0] &= ~((uint8_t)(~0) << _bit);
  _byte[0] |= (1 << _bit);

  /* we are done, so no need to wrap mod 8 and carry to bytes_off */
  state->out = NULL;
  state->bit_off = 0;
  state->byte_off++;
}

int
coder_dec_bytes (coder_state_t state, uint8_t *bytes, unsigned int nbytes)
{
  if (UNLIKELY (state->bit_off != 0))
    return 1;

  memcpy (bytes, state->in, nbytes);
  state->byte_off += nbytes;
  state->in += nbytes;
  return 0;
}

int
coder_dec_urandom (coder_state_t state, intvec_t v, const int_t m,
                   unsigned int mbits)
{
  unsigned int nbits;

  ASSERT_ERR (state->out == NULL);
  ASSERT_ERR (state->in != NULL);

  nbits = _udecode (v, &(state->in), &(state->bit_off), m, mbits);
  state->byte_off += (nbits >> 3);
  return nbits == 0 ? 1 : 0;
}

/* check padding. */
int
coder_dec_end (coder_state_t state)
{

  unsigned int _bit;
  const uint8_t *_byte;
  uint8_t succ = 0;

  ASSERT_ERR (state->out == NULL);
  ASSERT_ERR (state->in != NULL);

  _byte = state->in;
  _bit = state->bit_off;

  /* check final one and possibly trailing zeros */
  if ((_byte[0] & ((uint8_t)(~0) << _bit)) == ((uint8_t)1 << _bit))
    {
      succ = 1;
      state->bit_off = 0;
      state->byte_off++;
    }

  state->in = NULL;
  return succ;
}

static inline void
_inc_idx (unsigned int *byte, unsigned int *bit)
{
  (*bit)++;
  if (*bit > 7)
    {
      *bit = 0;
      (*byte)++;
    }
}

static inline void
_inc_idx_zero (uint8_t *buf, unsigned int *byte, unsigned int *bit)
{
  (*bit)++;
  if (*bit > 7)
    {
      *bit = 0;
      (*byte)++;
      buf[*byte] = 0;
    }
}

/*
 * Encode vector v uniform in {0,...,m-1}.
 */
static unsigned int
_uencode (uint8_t **byte, unsigned int *bit, const intvec_t v,
          UNUSED const int_t m, unsigned int mbits)
{
  unsigned int i, j, k, _bit, nbits = 0, _mbits;
  unsigned int nbytes = 0;
  int_srcptr elem;
  uint8_t *_byte;
  limb_t limb;

#if ASSERT == ASSERT_ENABLED
  {
    INT_T (zero, v->nlimbs);

    int_set_i64 (zero, 0);

    ASSERT_ERR (*bit <= 7);
    ASSERT_ERR (intvec_ge (v, zero) == 1);
    ASSERT_ERR (intvec_lt (v, m) == 1);
  }
#endif

  _byte = *byte;
  _bit = *bit;

  /* zero unset bits in first byte */
  _byte[0] &= ~((uint8_t)(~0) << _bit);

  _VEC_FOREACH_ELEM (v, i)
  {
    elem = intvec_get_elem_src (v, i);
    _mbits = mbits;

    for (j = 0; j < v->nlimbs; j++)
      {
        limb = elem->limbs[j];

        for (k = 0; k < MIN (NBITS_LIMB, _mbits); k++)
          {
            _byte[nbytes] |= ((limb & 1) << _bit);

            limb >>= 1;
            _inc_idx_zero (_byte, &nbytes, &_bit);
            nbits++;
          }

        _mbits -= k;
        if (_mbits == 0)
          break;
      }
  }

  *byte = _byte + nbytes;
  *bit = _bit;
  return nbits;
}

static unsigned int
_udecode (intvec_t v, const uint8_t **byte, unsigned int *bit, const int_t m,
          unsigned int mbits)
{
  INT_T (zero, 1);
  unsigned int i, j, k, _bit, nbits = 0, _mbits;
  unsigned int nbytes = 0;
  const uint8_t *_byte;
  int_ptr elem;
  limb_t limb;

  ASSERT_ERR (*bit <= 7);

  _byte = *byte;
  _bit = *bit;
  int_set_i64 (zero, 0);

  _VEC_FOREACH_ELEM (v, i)
  {
    elem = intvec_get_elem (v, i);
    _mbits = mbits;
    elem->neg = 0;

    for (j = 0; j < v->nlimbs; j++)
      {
        limb = 0;

        for (k = 0; k < MIN (NBITS_LIMB, _mbits); k++)
          {
            limb |= ((((limb_t)_byte[nbytes] & ((limb_t)1 << _bit)) >> _bit)
                     << k);

            _inc_idx (&nbytes, &_bit);
            nbits++;
          }

        _mbits -= k;
        elem->limbs[j] = limb;
      }
  }

  if (!intvec_lt (v, m))
    return 0; /* decoding failed */

  *byte = _byte + nbytes;
  *bit = _bit;
  return nbits;
}

void
coder_enc_urandom2 (coder_state_t state, poly_t v, const int_t m,
                    unsigned int mbits)
{
  intvec_ptr coeffvec;

  coeffvec = poly_get_coeffvec (v);
  coder_enc_urandom (state, coeffvec, m, mbits);
}

int
coder_dec_urandom2 (coder_state_t state, poly_t v, const int_t m,
                    unsigned int mbits)
{
  intvec_ptr coeffvec;

  coeffvec = poly_get_coeffvec (v);
  return coder_dec_urandom (state, coeffvec, m, mbits);
}

void
coder_enc_grandom2 (coder_state_t state, poly_t v, unsigned int log2o)
{
  intvec_ptr coeffvec;

  coeffvec = poly_get_coeffvec (v);
  coder_enc_grandom (state, coeffvec, log2o);
}

void
coder_enc_ghint2 (coder_state_t state, poly_t ghint)
{
  intvec_ptr coeffvec;

  coeffvec = poly_get_coeffvec (ghint);
  coder_enc_ghint (state, coeffvec);
}

void
coder_dec_grandom2 (coder_state_t state, poly_t v, unsigned int log2o)
{
  intvec_ptr coeffvec;

  coeffvec = poly_get_coeffvec (v);
  coder_dec_grandom (state, coeffvec, log2o);
}

void
coder_dec_ghint2 (coder_state_t state, poly_t ghint)
{
  intvec_ptr coeffvec;

  coeffvec = poly_get_coeffvec (ghint);
  coder_dec_ghint (state, coeffvec);
}

void
coder_enc_urandom3 (coder_state_t state, polyvec_t v, const int_t m,
                    unsigned int mbits)
{
  poly_ptr poly;
  unsigned int i;

  _VEC_FOREACH_ELEM (v, i)
  {
    poly = polyvec_get_elem (v, i);
    coder_enc_urandom2 (state, poly, m, mbits);
  }
}

void
coder_enc_urandom4 (coder_state_t state, polymat_t v, const int_t m,
                    unsigned int mbits)
{
  poly_ptr poly;
  unsigned int i, j;

  _MAT_FOREACH_ELEM (v, i, j)
  {
    poly = polymat_get_elem (v, i, j);
    coder_enc_urandom2 (state, poly, m, mbits);
  }
}

void
coder_enc_urandom4diag (coder_state_t state, polymat_t v, const int_t m,
                        unsigned int mbits)
{
  poly_ptr poly;
  unsigned int i, j;

  ASSERT_ERR (polymat_is_upperdiag (v));

  _MAT_FOREACH_ELEM_UPPER (v, i, j)
  {
    poly = polymat_get_elem (v, i, j);
    coder_enc_urandom2 (state, poly, m, mbits);
  }
}

/*
 * encoding of (sorted) sparse matrix:
 * 16 bit number of elements
 * || (per element: 16 bit row || 16 bit col)
 * || (per element: encoding of poly)
 * 16 bit integers are little-endian.
 */
void
coder_enc_urandom5 (coder_state_t state, spolymat_t v, const int_t m,
                    unsigned int mbits)
{
  const unsigned int len = 1 + v->nelems * 2;
  uint16_t buf[len], row, col;
  uint8_t *ptr;
  poly_ptr poly;
  unsigned int i, _bit;
  unsigned int nbytes = 0;
  uint8_t *_byte;

  ASSERT_ERR (state->out != NULL);
  ASSERT_ERR (state->in == NULL);
  ASSERT_ERR (v->nelems > 0);
  ASSERT_ERR (v->sorted);
  ASSERT_ERR (v->nelems <= UINT16_MAX);

  _byte = state->out;
  _bit = state->bit_off;

  buf[0] = htole16 (v->nelems);

  _SMAT_FOREACH_ELEM (v, i)
  {
    row = spolymat_get_row (v, i);
    col = spolymat_get_col (v, i);

    buf[1 + 2 * i] = htole16 (row);
    buf[1 + 2 * i + 1] = htole16 (col);
  }

  ptr = (uint8_t *)buf;
  for (i = 0; i < len * 2 * 8; i++)
    {
      if (ptr[i / 8] & (1 << (i % 8)))
        _byte[nbytes] |= (1 << _bit);

      _inc_idx_zero (_byte, &nbytes, &_bit);
    }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->out += nbytes;

  _SMAT_FOREACH_ELEM (v, i)
  {
    poly = spolymat_get_elem (v, i);

    coder_enc_urandom2 (state, poly, m, mbits);
  }
}

/*
 * encoding of (sorted) sparse vector:
 * 16 bit number of elements
 * || (per element: 16 bit idx)
 * || (per element: encoding of poly)
 * 16 bit integers are little-endian.
 */
void
coder_enc_urandom6 (coder_state_t state, spolyvec_t v, const int_t m,
                    unsigned int mbits)
{
  const unsigned int len = 1 + v->nelems;
  uint16_t buf[len], idx;
  uint8_t *ptr;
  poly_ptr poly;
  unsigned int i, _bit;
  unsigned int nbytes = 0;
  uint8_t *_byte;

  ASSERT_ERR (state->out != NULL);
  ASSERT_ERR (state->in == NULL);
  ASSERT_ERR (v->nelems > 0);
  ASSERT_ERR (v->sorted);
  ASSERT_ERR (v->nelems <= UINT16_MAX);

  _byte = state->out;
  _bit = state->bit_off;

  buf[0] = htole16 (v->nelems);

  _SMAT_FOREACH_ELEM (v, i)
  {
    idx = spolyvec_get_elem_ (v, i);

    buf[1 + i] = htole16 (idx);
  }

  ptr = (uint8_t *)buf;
  for (i = 0; i < len * 2 * 8; i++)
    {
      if (ptr[i / 8] & (1 << (i % 8)))
        _byte[nbytes] |= (1 << _bit);

      _inc_idx_zero (_byte, &nbytes, &_bit);
    }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->out += nbytes;

  _SMAT_FOREACH_ELEM (v, i)
  {
    poly = spolyvec_get_elem (v, i);

    coder_enc_urandom2 (state, poly, m, mbits);
  }
}

void
coder_dec_urandom5 (coder_state_t state, spolymat_t v, const int_t m,
                    unsigned int mbits)
{
  uint16_t header, row, col;
  uint8_t *ptr;
  unsigned int nelems, len;
  poly_ptr poly;
  unsigned int i, _bit;
  unsigned int nbytes = 0;
  const uint8_t *_byte;

  ASSERT_ERR (state->out == NULL);
  ASSERT_ERR (state->in != NULL);

  spolymat_set_empty (v);

  _byte = state->in;
  _bit = state->bit_off;

  header = 0;
  ptr = (uint8_t *)&header;
  for (i = 0; i < 16; i++)
    {
      if (_byte[nbytes] & (1 << _bit))
        ptr[i / 8] |= (1 << (i % 8));

      _inc_idx (&nbytes, &_bit);
    }
  nelems = le16toh (header);

  ASSERT_ERR (nelems <= v->nelems_max);

  len = nelems * 2;
  uint16_t buf[len];

  ptr = (uint8_t *)buf;
  for (i = 0; i < len * 2 * 8; i++)
    {
      if (_byte[nbytes] & (1 << _bit))
        ptr[i / 8] |= (1 << (i % 8));

      _inc_idx (&nbytes, &_bit);
    }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->in += nbytes;

  _SMAT_FOREACH_ELEM (v, i)
  {
    row = le16toh (buf[2 * i]);
    col = le16toh (buf[2 * i + 1]);
    poly = spolymat_insert_elem (v, row, col);

    coder_dec_urandom2 (state, poly, m, mbits);
  }
}

void
coder_dec_urandom6 (coder_state_t state, spolyvec_t v, const int_t m,
                    unsigned int mbits)
{
  uint16_t header, idx;
  uint8_t *ptr;
  unsigned int nelems, len;
  poly_ptr poly;
  unsigned int i, _bit;
  unsigned int nbytes = 0;
  const uint8_t *_byte;

  ASSERT_ERR (state->out == NULL);
  ASSERT_ERR (state->in != NULL);

  spolyvec_set_empty (v);

  _byte = state->in;
  _bit = state->bit_off;

  header = 0;
  ptr = (uint8_t *)&header;
  for (i = 0; i < 16; i++)
    {
      if (_byte[nbytes] & (1 << _bit))
        ptr[i / 8] |= (1 << (i % 8));

      _inc_idx (&nbytes, &_bit);
    }
  nelems = le16toh (header);

  ASSERT_ERR (nelems <= v->nelems_max);

  len = nelems;
  uint16_t buf[len];

  ptr = (uint8_t *)buf;
  for (i = 0; i < len * 2 * 8; i++)
    {
      if (_byte[nbytes] & (1 << _bit))
        ptr[i / 8] |= (1 << (i % 8));

      _inc_idx (&nbytes, &_bit);
    }

  state->byte_off += nbytes;
  state->bit_off = _bit;
  state->in += nbytes;

  _SMAT_FOREACH_ELEM (v, i)
  {
    idx = le16toh (buf[2 * i]);
    poly = spolyvec_insert_elem (v, idx);

    coder_dec_urandom2 (state, poly, m, mbits);
  }
}

int
coder_dec_urandom3 (coder_state_t state, polyvec_t v, const int_t m,
                    unsigned int mbits)
{
  poly_ptr poly;
  unsigned int i;
  int rc;

  _VEC_FOREACH_ELEM (v, i)
  {
    poly = polyvec_get_elem (v, i);
    rc = coder_dec_urandom2 (state, poly, m, mbits);
    if (rc != 0)
      return 1;
  }
  return 0;
}

int
coder_dec_urandom4 (coder_state_t state, polymat_t v, const int_t m,
                    unsigned int mbits)
{
  poly_ptr poly;
  unsigned int i, j;
  int rc;

  _MAT_FOREACH_ELEM_UPPER (v, i, j)
  {
    poly = polymat_get_elem (v, i, j);
    rc = coder_dec_urandom2 (state, poly, m, mbits);
    if (rc != 0)
      return 1;
  }
  return 0;
}

int
coder_dec_urandom4diag (coder_state_t state, polymat_t v, const int_t m,
                        unsigned int mbits)
{
  poly_ptr poly;
  unsigned int i, j;
  int rc;

  ASSERT_ERR (polymat_is_upperdiag (v));

  _MAT_FOREACH_ELEM_UPPER (v, i, j)
  {
    poly = polymat_get_elem (v, i, j);
    rc = coder_dec_urandom2 (state, poly, m, mbits);
    if (rc != 0)
      return 1;
  }
  return 0;
}

void
coder_enc_grandom3 (coder_state_t state, polyvec_t v, unsigned int log2o)
{
  poly_ptr poly;
  unsigned int i;

  _VEC_FOREACH_ELEM (v, i)
  {
    poly = polyvec_get_elem (v, i);
    coder_enc_grandom2 (state, poly, log2o);
  }
}

void
coder_enc_ghint3 (coder_state_t state, polyvec_t ghint)
{
  poly_ptr poly;
  unsigned int i;

  _VEC_FOREACH_ELEM (ghint, i)
  {
    poly = polyvec_get_elem (ghint, i);
    coder_enc_ghint2 (state, poly);
  }
}

void
coder_dec_grandom3 (coder_state_t state, polyvec_t v, unsigned int log2o)
{
  poly_ptr poly;
  unsigned int i;

  _VEC_FOREACH_ELEM (v, i)
  {
    poly = polyvec_get_elem (v, i);
    coder_dec_grandom2 (state, poly, log2o);
  }
}

void
coder_dec_ghint3 (coder_state_t state, polyvec_t ghint)
{
  poly_ptr poly;
  unsigned int i;

  _VEC_FOREACH_ELEM (ghint, i)
  {
    poly = polyvec_get_elem (ghint, i);
    coder_dec_ghint2 (state, poly);
  }
}
