#include "rng.h"
#include "aes256ctr.h"
#include "lazer.h"
#include "shake128.h"

#include <string.h>

void
rng_init (rng_state_t state, const uint8_t seed[32], uint64_t dom)
{
  _rng_init (state, seed, dom);
}

void
rng_urandom (rng_state_t state, uint8_t *out, size_t outlen)
{
  _rng_urandom (state, out, outlen);
}

void
rng_clear (rng_state_t state)
{
  _rng_clear (state);
}

static void
_rng_init (rng_state_t state, const uint8_t seed[32], uint64_t dom)
{
  dom = htole64 (dom);

#if RNG == RNG_SHAKE128
  _shake128_init (state->state);
  _shake128_absorb (state->state, seed, 32);
  _shake128_absorb (state->state, (uint8_t *)&dom, sizeof (dom));
#elif RNG == RNG_AES256CTR
  {
    uint64_t nonce[2] = { dom, 0 };

    _aes256ctr_init (state->state, seed, (uint8_t *)&nonce);
  }
#endif
}

static void
_rng_urandom (rng_state_t state, uint8_t *out, size_t outlen)
{
  ASSERT_ERR (outlen > 0);

#if RNG == RNG_SHAKE128
  _shake128_squeeze (state->state, out, outlen);
#elif RNG == RNG_AES256CTR
  _aes256ctr_stream (state->state, out, outlen);
#endif
}

static void
_rng_clear (rng_state_t state)
{
#if RNG == RNG_SHAKE128
  _shake128_clear (state->state);
#elif RNG == RNG_AES256CTR
  _aes256ctr_clear (state->state);
#endif
}
