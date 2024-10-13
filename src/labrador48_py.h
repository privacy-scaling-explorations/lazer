#ifndef LABRADOR48_PY_H
#define LABRADOR48_PY_H

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

#define LABRADOR48_N 64
#define LABRADOR48_K 8
#define LABRADOR48_L 4

typedef union {
  __m512i v[LABRADOR48_N/32];
  int16_t c[LABRADOR48_N];
} labrador48_vecn;

typedef struct {
  labrador48_vecn vec[1];
} labrador48_poly;

typedef struct {
  labrador48_poly vec[LABRADOR48_K];
} labrador48_polx;

typedef struct {
  labrador48_vecn limbs[LABRADOR48_L];
} labrador48_polz;

typedef struct {
  size_t len;
  size_t *rows;
  size_t *cols;
  labrador48_polx *coeffs;
} labrador48_sparsemat;

typedef struct {
  size_t deg;
  labrador48_sparsemat a[1];
  size_t nz;
  size_t *idx;
  size_t *off;
  size_t *len;
  size_t *mult;
  labrador48_polx **phi;
  labrador48_polx *b;
} labrador48_sparsecnst;

typedef struct {
  size_t r;
  size_t *n;
  uint64_t *normsq;
  labrador48_poly **s;
} labrador48_witness;

typedef struct {
  size_t deg;
  size_t nz;
  size_t *idx;
  labrador48_polx **phi;
  labrador48_polx *b;
} labrador48_lincnst;

typedef struct {
  size_t r;
  size_t *n;
  size_t fu;
  size_t bu;
  size_t kappa;
  size_t kappa1;
  labrador48_polz *u;
  labrador48_polx *alpha;
} labrador48_commitment;

typedef struct {
  size_t f;       // amortized opening decomposition parts
  size_t fu;      // uniform decomposition parts
  size_t fg;      // quadratic garbage decomposition parts
  size_t b;       // amortized opening decomposition basis
  size_t bu;      // uniform decomposition basis
  size_t bg;      // quadratic garbage decomposition basis
  size_t kappa;   // inner commitment rank
  size_t kappa1;  // outer commitment rank
  size_t u1len;
  size_t u2len;
} labrador48_comparams;

typedef struct {
  size_t deg;  // extension degree
  labrador48_sparsemat a[1];
  labrador48_polx *phi;
  labrador48_polx *b;
} labrador48_constraint;

typedef struct {
  size_t r;
  size_t *n;
  size_t k;
  labrador48_sparsecnst *cnst;
  uint64_t *betasq;
  uint8_t h[16];
  labrador48_polx *u0;
  labrador48_polx *alpha;
} labrador48_smplstmnt;

typedef struct {
  size_t r;          // input witness multiplicity
  size_t *n;         // input witness ranks (r)
  size_t *nu;        // input witness decomposition parts (r)
  int tail;
  labrador48_comparams cpp[1];  // commitment parameters
  labrador48_polz *u1;          // outer commitment 1 (kappa1)
  size_t jlnonce;    // JL matrix nonce
  int32_t p[256];    // JL projection
  labrador48_polz *bb;          // int to pol extensions (4)
  labrador48_polz *u2;          // outer commitment 2 (kappa1)
  uint64_t normsq;   // output witness norm
} labrador48_proof;

typedef struct {
  size_t l;
  double size;
  labrador48_proof *pi[16];
  labrador48_witness owt;
} labrador48_composite;

void labrador48_init_witness_raw(labrador48_witness *wt, size_t r, size_t n[]);
int labrador48_set_witness_vector_raw(labrador48_witness *wt, size_t i, size_t n, size_t deg, const int64_t s[]);
void labrador48_init_smplstmnt_raw(labrador48_smplstmnt *st, size_t r, size_t n[], uint64_t betasq[], size_t k);
int labrador48_set_smplstmnt_lincnst_raw(labrador48_smplstmnt *st, size_t i, size_t nz, size_t idx[], size_t n[],size_t deg, int64_t *phi, int64_t *b);
int labrador48_simple_verify(const labrador48_smplstmnt *st, const labrador48_witness *wt);
void labrador48_init_comkey(size_t n);
void labrador48_free_comkey();
void labrador48_free_commitment(labrador48_commitment *com);
void labrador48_free_witness(labrador48_witness *wt);
void labrador48_free_composite(labrador48_composite *proof);
void labrador48_free_smplstmnt(labrador48_smplstmnt *st);

int labrador48_composite_prove_simple(labrador48_composite *proof, labrador48_commitment *com, const labrador48_smplstmnt *st, const labrador48_witness *wt);
int labrador48_composite_verify_simple(const labrador48_composite *proof, const labrador48_commitment *com, const labrador48_smplstmnt *st);

#endif
