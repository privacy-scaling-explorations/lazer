#include "test.h"
#include <math.h>

#define SP_NSAMPLES 5000000
#define SP_LOG2O 10

#define MP_LOG2NSAMPLES 22
#define MP_LOG2O 128

#define SAMPLEVEC_DIM 4

static void single_prec (void);
static void multiple_prec (void);

int
main (void)
{
  lazer_init();

  single_prec ();
  multiple_prec ();
  TEST_PASS ();
}

static void
single_prec (void)
{
  const long double sigma = (long double)1.55 * (1 << SP_LOG2O);
  uint8_t seed[32] = { 0 };
  uint64_t dom = 0;
  uint64_t *seed_ptr = (uint64_t *)&(seed[0]);
  long sum = 0, sum_sqr = 0;
  long double empiric_mean, empiric_var;
  INT_T (sample, 1);
  INT_T (sample_sqr, 2);
  INTVEC_T (samplevec, SAMPLEVEC_DIM, 1);
  unsigned int i;
  int64_t tmp;

  bytes_urandom (seed, sizeof (seed));

  printf ("single precision test\n\n");

  for (*seed_ptr = 0; *seed_ptr < SP_NSAMPLES; (*seed_ptr)++)
    {
      dom = 0;
      int_grandom (sample, SP_LOG2O, seed, dom);
      dom = 1;
      intvec_grandom (samplevec, SP_LOG2O, seed, dom);

      tmp = int_get_i64 (sample);
      sum += tmp;
      sum_sqr += tmp * tmp;

      _VEC_FOREACH_ELEM (samplevec, i)
      {
        tmp = intvec_get_elem_i64 (samplevec, i);
        sum += tmp;
        sum_sqr += tmp * tmp;
      }
    }

  empiric_mean = (long double)sum / (SP_NSAMPLES * (1 + SAMPLEVEC_DIM));
  empiric_var = (long double)sum_sqr / (SP_NSAMPLES * (1 + SAMPLEVEC_DIM));

  printf ("expected mean:     0\n");
  printf ("empiric mean:      %Lf\n", empiric_mean);

  printf ("expected variance: %Lf\n", sigma * sigma);
  printf ("empiric variance:  %Lf\n", empiric_var);
  printf ("\n");

  TEST_EXPECT (fabsl (empiric_mean) < 1);
  TEST_EXPECT (fabsl (empiric_var - sigma * sigma) < 10000);
}

static void
multiple_prec (void)
{
  uint8_t seed[32] = { 0 };
  uint32_t dom = 0;
  uint64_t *seed_ptr = (uint64_t *)&(seed[0]);
  INT_T (sample, 4);
  INT_T (sum, 4);
  INT_T (sample_sqr, 8);
  INT_T (sum_sqr, 8);

  bytes_urandom (seed, sizeof (seed));

  printf ("multiple precision test\n\n");

  int_set_i64 (sum, 0);
  int_set_i64 (sum_sqr, 0);

  for (*seed_ptr = 0; *seed_ptr < (1 << MP_LOG2NSAMPLES); (*seed_ptr)++)
    {
      int_grandom (sample, MP_LOG2O, seed, dom);

      int_mul (sample_sqr, sample, sample);
      int_add (sum, sum, sample);
      int_add (sum_sqr, sum_sqr, sample_sqr);
    }

  /* divide by number of samples */
  int_rshift (sum, sum, MP_LOG2NSAMPLES);
  int_rshift (sum_sqr, sum_sqr, MP_LOG2NSAMPLES);

  printf ("expected mean:     0\n");
  printf ("empiric mean:      ");
  int_out_str (stdout, 10, sum);
  printf ("\n");

  printf ("expected variance: "
          "2781904943926521595051292914833726986174811381592014551047968455790"
          "11293959946.24\n");
  printf ("empiric variance:  ");
  int_out_str (stdout, 10, sum_sqr);
  printf ("\n\n");
}
