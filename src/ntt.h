#include "lazer.h"
#include <stdint.h>

static void _ntt (crtcoeff_t c[], const modulus_t _p, const unsigned int d);
static void _intt (crtcoeff_t c[], const modulus_t _p, const unsigned int d);
