#ifndef LABRADOR24_PY_H
#define LABRADOR24_PY_H

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

#define LABRADOR24_N 64
#define LABRADOR24_K 5
#define LABRADOR24_L 2

typedef union {
  __m512i v[LABRADOR24_N/32];
  int16_t c[LABRADOR24_N];
} labrador24_vecn;

typedef struct {
  labrador24_vecn vec[1];
} labrador24_poly;

typedef struct {
  labrador24_poly vec[LABRADOR24_K];
} labrador24_polx;

typedef struct {
  labrador24_vecn limbs[LABRADOR24_L];
} labrador24_polz;

typedef struct {
  size_t len;
  size_t *rows;
  size_t *cols;
  labrador24_polx *coeffs;
} labrador24_sparsemat;

typedef struct {
  size_t deg;
  labrador24_sparsemat a[1];
  size_t nz;
  size_t *idx;
  size_t *off;
  size_t *len;
  size_t *mult;
  labrador24_polx **phi;
  labrador24_polx *b;
} labrador24_sparsecnst;

typedef struct {
  size_t r;
  size_t *n;
  uint64_t *normsq;
  labrador24_poly **s;
} labrador24_witness;

typedef struct {
  size_t deg;
  size_t nz;
  size_t *idx;
  labrador24_polx **phi;
  labrador24_polx *b;
} labrador24_lincnst;

typedef struct {
  size_t r;
  size_t *n;
  size_t fu;
  size_t bu;
  size_t kappa;
  size_t kappa1;
  labrador24_polz *u;
  labrador24_polx *alpha;
} labrador24_commitment;

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
} labrador24_comparams;

typedef struct {
  size_t deg;  // extension degree
  labrador24_sparsemat a[1];
  labrador24_polx *phi;
  labrador24_polx *b;
} labrador24_constraint;

typedef struct {
  size_t r;
  size_t *n;
  size_t k;
  labrador24_sparsecnst *cnst;
  uint64_t *betasq;
  uint8_t h[16];
  labrador24_polx *u0;
  labrador24_polx *alpha;
} labrador24_smplstmnt;

typedef struct {
  size_t r;          // input witness multiplicity
  size_t *n;         // input witness ranks (r)
  size_t *nu;        // input witness decomposition parts (r)
  int tail;
  labrador24_comparams cpp[1];  // commitment parameters
  labrador24_polz *u1;          // outer commitment 1 (kappa1)
  size_t jlnonce;    // JL matrix nonce
  int32_t p[256];    // JL projection
  labrador24_polz *bb;          // int to pol extensions (4)
  labrador24_polz *u2;          // outer commitment 2 (kappa1)
  uint64_t normsq;   // output witness norm
} labrador24_proof;

typedef struct {
  size_t l;
  double size;
  labrador24_proof *pi[16];
  labrador24_witness owt;
} labrador24_composite;

void labrador24_init_witness_raw(labrador24_witness *wt, size_t r, size_t n[]);
int labrador24_set_witness_vector_raw(labrador24_witness *wt, size_t i, size_t n, size_t deg, const int64_t s[]);
void labrador24_init_smplstmnt_raw(labrador24_smplstmnt *st, size_t r, size_t n[], uint64_t betasq[], size_t k);
int labrador24_set_smplstmnt_lincnst_raw(labrador24_smplstmnt *st, size_t i, size_t nz, size_t idx[], size_t n[],size_t deg, int64_t *phi, int64_t *b);
int labrador24_simple_verify(const labrador24_smplstmnt *st, const labrador24_witness *wt);
void labrador24_init_comkey(size_t n);
void labrador24_free_comkey();
void labrador24_free_commitment(labrador24_commitment *com);
void labrador24_free_witness(labrador24_witness *wt);
void labrador24_free_composite(labrador24_composite *proof);
void labrador24_free_smplstmnt(labrador24_smplstmnt *st);

int labrador24_composite_prove_simple(labrador24_composite *proof, labrador24_commitment *com, const labrador24_smplstmnt *st, const labrador24_witness *wt);
int labrador24_composite_verify_simple(const labrador24_composite *proof, const labrador24_commitment *com, const labrador24_smplstmnt *st);

#endif
