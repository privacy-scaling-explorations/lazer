#include "lazer.h"
#include "test.h"
#include <math.h>

#define NSAMPLES 2000000
#define LOG2O 10

static void standard_rej_test (void);
static void bimodal_rej_test (void);

int
main (void)
{
  lazer_init();

  standard_rej_test ();
  bimodal_rej_test ();
  TEST_PASS ();
}

static void
standard_rej_test (void)
{
  rng_state_t state;
  const long double sigma = (long double)1.55 * (1 << LOG2O);
  uint8_t seed[32] = { 0 };
  uint64_t dom = 0;
  long sum = 0, sum_sqr = 0;
  long double empiric_mean, empiric_var;
  unsigned int M;
  unsigned int nsamples;
  INT_T (sample, 1);
  INT_T (sample_sqr, 2);
  INT_T (T, 1);
  INT_T (scM, 3);
  INT_T (o2, 2);
  INTVEC_T (v, 1, 1);
  INTVEC_T (y, 1, 1);
  INTVEC_T (z, 1, 1);
  unsigned int i;
  int64_t tmp;
  unsigned int nrejs = 0;

  /*
   * o=gama*T,
   * o=1.55*2^10, M=4, gamma=9.6836868201367992, T=163.904515860580
   */
  int_set_i64 (T, 163);
  M = 4;
  /* M * 2^128 */
  int_set_i64 (scM, 0);
  scM->limbs[2] = M;
  int_set_i64 (o2, 2519203); /* o^2 = 2519203.84 */

  bytes_urandom (seed, sizeof (seed));
  rng_init (state, seed, dom);
  bytes_urandom (seed, sizeof (seed));

  printf ("standard rejection test\n\n");

  for (i = 0; i < NSAMPLES; i++)
    {
      intvec_urandom (v, T, 8, seed, dom);
      dom++;
      intvec_grandom (y, LOG2O, seed, dom);
      dom++;

      intvec_add (z, y, v);

      if (rej_standard (state, z, v, scM, o2) == 0)
        {

          tmp = intvec_get_elem_i64 (z, 0);
          sum += tmp;
          sum_sqr += tmp * tmp;
        }
      else
        {
          nrejs++;
        }
    }

  rng_clear (state);

  nsamples = NSAMPLES - nrejs;

  empiric_mean = (long double)sum / nsamples;
  empiric_var = (long double)sum_sqr / nsamples;

  printf ("expected rejection rate: %0.2f\n", (float)1 - (float)1 / M);
  printf ("empiric rejection rate:  %0.2f\n", (float)nrejs / NSAMPLES);

  printf ("expected mean:           0\n");
  printf ("empiric mean:            %Lf\n", empiric_mean);

  printf ("expected variance:       %Lf\n", sigma * sigma);
  printf ("empiric variance:        %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (empiric_mean) < 10);
  TEST_EXPECT (fabsl (empiric_var - sigma * sigma) < 20000);
}

static void
bimodal_rej_test (void)
{
  rng_state_t state;
  const long double sigma = (long double)1.55 * (1 << LOG2O);
  uint8_t seed[32] = { 0 };
  uint64_t dom = 0;
  long sum = 0, sum_sqr = 0;
  long double empiric_mean, empiric_var;
  unsigned int M;
  unsigned int nsamples;
  INT_T (sample, 1);
  INT_T (sample_sqr, 2);
  INT_T (T, 1);
  INT_T (scM, 3);
  INT_T (o2, 2);
  INTVEC_T (v, 1, 1);
  INTVEC_T (y, 1, 1);
  INTVEC_T (z, 1, 1);
  unsigned int i;
  int64_t tmp;
  unsigned int nrejs = 0;
  int8_t sign;

  /*
   * o=gama*T,
   * o=1.55*2^10, M=4, gamma=, T=1868.78518773656
   */
  int_set_i64 (T, 1868);
  M = 2;
  /* M * 2^128 */
  int_set_i64 (scM, 0);
  scM->limbs[2] = M;
  int_set_i64 (o2, 2519203); /* o^2 = 2519203.84 */

  bytes_urandom (seed, sizeof (seed));
  rng_init (state, seed, dom);
  bytes_urandom (seed, sizeof (seed));

  printf ("bimodal rejection test\n\n");

  for (i = 0; i < NSAMPLES; i++)
    {
      intvec_urandom (v, T, 11, seed, dom);
      dom++;
      intvec_grandom (y, LOG2O, seed, dom);
      dom++;

      rng_urandom (state, (uint8_t *)&sign, sizeof (sign));
      if ((sign & 1) == 0)
        sign = -1;
      else
        sign = 1;
      intvec_mul_sgn_self (v, sign);

      intvec_add (z, y, v);

      if (rej_bimodal (state, z, v, scM, o2) == 0)
        {

          tmp = intvec_get_elem_i64 (z, 0);
          sum += tmp;
          sum_sqr += tmp * tmp;
        }
      else
        {
          nrejs++;
        }
    }

  rng_clear (state);

  nsamples = NSAMPLES - nrejs;

  empiric_mean = (long double)sum / nsamples;
  empiric_var = (long double)sum_sqr / nsamples;

  printf ("expected rejection rate: %0.2f\n", (float)1 - (float)1 / M);
  printf ("empiric rejection rate:  %0.2f\n", (float)nrejs / NSAMPLES);

  printf ("expected mean:           0\n");
  printf ("empiric mean:            %Lf\n", empiric_mean);

  printf ("expected variance:       %Lf\n", sigma * sigma);
  printf ("empiric variance:        %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (empiric_mean) < 10);
  TEST_EXPECT (fabsl (empiric_var - sigma * sigma) < 10000);
}
