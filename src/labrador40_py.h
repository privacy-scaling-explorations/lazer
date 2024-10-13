#ifndef LABRADOR40_PY_H
#define LABRADOR40_PY_H

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

#define LABRADOR40_N 64
#define LABRADOR40_K 7
#define LABRADOR40_L 3

typedef union {
  __m512i v[LABRADOR40_N/32];
  int16_t c[LABRADOR40_N];
} labrador40_vecn;

typedef struct {
  labrador40_vecn vec[1];
} labrador40_poly;

typedef struct {
  labrador40_poly vec[LABRADOR40_K];
} labrador40_polx;

typedef struct {
  labrador40_vecn limbs[LABRADOR40_L];
} labrador40_polz;

typedef struct {
  size_t len;
  size_t *rows;
  size_t *cols;
  labrador40_polx *coeffs;
} labrador40_sparsemat;

typedef struct {
  size_t deg;
  labrador40_sparsemat a[1];
  size_t nz;
  size_t *idx;
  size_t *off;
  size_t *len;
  size_t *mult;
  labrador40_polx **phi;
  labrador40_polx *b;
} labrador40_sparsecnst;

typedef struct {
  size_t r;
  size_t *n;
  uint64_t *normsq;
  labrador40_poly **s;
} labrador40_witness;

typedef struct {
  size_t deg;
  size_t nz;
  size_t *idx;
  labrador40_polx **phi;
  labrador40_polx *b;
} labrador40_lincnst;

typedef struct {
  size_t r;
  size_t *n;
  size_t fu;
  size_t bu;
  size_t kappa;
  size_t kappa1;
  labrador40_polz *u;
  labrador40_polx *alpha;
} labrador40_commitment;

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
} labrador40_comparams;

typedef struct {
  size_t deg;  // extension degree
  labrador40_sparsemat a[1];
  labrador40_polx *phi;
  labrador40_polx *b;
} labrador40_constraint;

typedef struct {
  size_t r;
  size_t *n;
  size_t k;
  labrador40_sparsecnst *cnst;
  uint64_t *betasq;
  uint8_t h[16];
  labrador40_polx *u0;
  labrador40_polx *alpha;
} labrador40_smplstmnt;

typedef struct {
  size_t r;          // input witness multiplicity
  size_t *n;         // input witness ranks (r)
  size_t *nu;        // input witness decomposition parts (r)
  int tail;
  labrador40_comparams cpp[1];  // commitment parameters
  labrador40_polz *u1;          // outer commitment 1 (kappa1)
  size_t jlnonce;    // JL matrix nonce
  int32_t p[256];    // JL projection
  labrador40_polz *bb;          // int to pol extensions (4)
  labrador40_polz *u2;          // outer commitment 2 (kappa1)
  uint64_t normsq;   // output witness norm
} labrador40_proof;

typedef struct {
  size_t l;
  double size;
  labrador40_proof *pi[16];
  labrador40_witness owt;
} labrador40_composite;

void labrador40_init_witness_raw(labrador40_witness *wt, size_t r, size_t n[]);
int labrador40_set_witness_vector_raw(labrador40_witness *wt, size_t i, size_t n, size_t deg, const int64_t s[]);
void labrador40_init_smplstmnt_raw(labrador40_smplstmnt *st, size_t r, size_t n[], uint64_t betasq[], size_t k);
int labrador40_set_smplstmnt_lincnst_raw(labrador40_smplstmnt *st, size_t i, size_t nz, size_t idx[], size_t n[],size_t deg, int64_t *phi, int64_t *b);
int labrador40_simple_verify(const labrador40_smplstmnt *st, const labrador40_witness *wt);
void labrador40_init_comkey(size_t n);
void labrador40_free_comkey();
void labrador40_free_commitment(labrador40_commitment *com);
void labrador40_free_witness(labrador40_witness *wt);
void labrador40_free_composite(labrador40_composite *proof);
void labrador40_free_smplstmnt(labrador40_smplstmnt *st);

int labrador40_composite_prove_simple(labrador40_composite *proof, labrador40_commitment *com, const labrador40_smplstmnt *st, const labrador40_witness *wt);
int labrador40_composite_verify_simple(const labrador40_composite *proof, const labrador40_commitment *com, const labrador40_smplstmnt *st);

#endif
