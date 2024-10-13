#ifndef BLINDSIG_H
#define BLINDSIG_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define FALCON_P 12289

#include "blindsig-p1-params.h"
#include "blindsig-p2-params.h"
#include "falcon.h"
#include "inner.h"

// Falcon-512
#define Q_ 12289
#define DEG_ 512
#define SKBYTES 1281
#define BOUND 34034726
#define SIGNER_LOGN 9
#define PUBKEYLEN FALCON_PUBKEY_SIZE (SIGNER_LOGN)
#define PRIVKEYLEN FALCON_PRIVKEY_SIZE (SIGNER_LOGN)
#define TMPKGLEN FALCON_TMPSIZE_KEYGEN (SIGNER_LOGN)

// XXX moved to lazer.h
#ifdef XXX
void falcon_keygen (uint8_t sk[PRIVKEYLEN], uint8_t pk[PUBKEYLEN]);
void falcon_preimage_sample (int16_t s1[DEG_], int16_t s2[DEG_],
                             const int16_t t[DEG_],
                             const uint8_t sk[PRIVKEYLEN]);
void falcon_decode_pubkey (int16_t h[DEG_], const uint8_t pk[PUBKEYLEN]);
void falcon_redc (int16_t c[DEG_]);
void falcon_add (int16_t c[DEG_], const int16_t a[DEG_], const int16_t b[DEG_]);
void falcon_mul (int16_t c[DEG_], const int16_t a[DEG_], const int16_t b[DEG_]);
#endif

/* SHA3-256(0x00) */
static const uint8_t public_randomness[32]
    = { 0x5d, 0x53, 0x46, 0x9f, 0x20, 0xfe, 0xf4, 0xf8, 0xea, 0xb5, 0x2b,
        0x88, 0x04, 0x4e, 0xde, 0x69, 0xc7, 0x7a, 0x6a, 0x68, 0xa6, 0x07,
        0x28, 0x60, 0x9f, 0xc4, 0xa6, 0x5f, 0xf5, 0x31, 0xe7, 0xd0 };

#endif
