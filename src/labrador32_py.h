#ifndef LABRADOR32_PY_H
#define LABRADOR32_PY_H

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

#define LABRADOR32_N 64
#define LABRADOR32_K 6
#define LABRADOR32_L 3

typedef union {
  __m512i v[LABRADOR32_N/32];
  int16_t c[LABRADOR32_N];
} labrador32_vecn;

typedef struct {
  labrador32_vecn vec[1];
} labrador32_poly;

typedef struct {
  labrador32_poly vec[LABRADOR32_K];
} labrador32_polx;

typedef struct {
  labrador32_vecn limbs[LABRADOR32_L];
} labrador32_polz;

typedef struct {
  size_t len;
  size_t *rows;
  size_t *cols;
  labrador32_polx *coeffs;
} labrador32_sparsemat;

typedef struct {
  size_t deg;
  labrador32_sparsemat a[1];
  size_t nz;
  size_t *idx;
  size_t *off;
  size_t *len;
  size_t *mult;
  labrador32_polx **phi;
  labrador32_polx *b;
} labrador32_sparsecnst;

typedef struct {
  size_t r;
  size_t *n;
  uint64_t *normsq;
  labrador32_poly **s;
} labrador32_witness;

typedef struct {
  size_t deg;
  size_t nz;
  size_t *idx;
  labrador32_polx **phi;
  labrador32_polx *b;
} labrador32_lincnst;

typedef struct {
  size_t r;
  size_t *n;
  size_t fu;
  size_t bu;
  size_t kappa;
  size_t kappa1;
  labrador32_polz *u;
  labrador32_polx *alpha;
} labrador32_commitment;

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
} labrador32_comparams;

typedef struct {
  size_t deg;  // extension degree
  labrador32_sparsemat a[1];
  labrador32_polx *phi;
  labrador32_polx *b;
} labrador32_constraint;

typedef struct {
  size_t r;
  size_t *n;
  size_t k;
  labrador32_sparsecnst *cnst;
  uint64_t *betasq;
  uint8_t h[16];
  labrador32_polx *u0;
  labrador32_polx *alpha;
} labrador32_smplstmnt;

typedef struct {
  size_t r;          // input witness multiplicity
  size_t *n;         // input witness ranks (r)
  size_t *nu;        // input witness decomposition parts (r)
  int tail;
  labrador32_comparams cpp[1];  // commitment parameters
  labrador32_polz *u1;          // outer commitment 1 (kappa1)
  size_t jlnonce;    // JL matrix nonce
  int32_t p[256];    // JL projection
  labrador32_polz *bb;          // int to pol extensions (4)
  labrador32_polz *u2;          // outer commitment 2 (kappa1)
  uint64_t normsq;   // output witness norm
} labrador32_proof;

typedef struct {
  size_t l;
  double size;
  labrador32_proof *pi[16];
  labrador32_witness owt;
} labrador32_composite;

void labrador32_init_witness_raw(labrador32_witness *wt, size_t r, size_t n[]);
int labrador32_set_witness_vector_raw(labrador32_witness *wt, size_t i, size_t n, size_t deg, const int64_t s[]);
void labrador32_init_smplstmnt_raw(labrador32_smplstmnt *st, size_t r, size_t n[], uint64_t betasq[], size_t k);
int labrador32_set_smplstmnt_lincnst_raw(labrador32_smplstmnt *st, size_t i, size_t nz, size_t idx[], size_t n[],size_t deg, int64_t *phi, int64_t *b);
int labrador32_simple_verify(const labrador32_smplstmnt *st, const labrador32_witness *wt);
void labrador32_init_comkey(size_t n);
void labrador32_free_comkey();
void labrador32_free_commitment(labrador32_commitment *com);
void labrador32_free_witness(labrador32_witness *wt);
void labrador32_free_composite(labrador32_composite *proof);
void labrador32_free_smplstmnt(labrador32_smplstmnt *st);

int labrador32_composite_prove_simple(labrador32_composite *proof, labrador32_commitment *com, const labrador32_smplstmnt *st, const labrador32_witness *wt);
int labrador32_composite_verify_simple(const labrador32_composite *proof, const labrador32_commitment *com, const labrador32_smplstmnt *st);

#endif
