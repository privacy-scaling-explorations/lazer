#include "test.h"
#include <math.h>

#define NSAMPLES 1000000
#define NROWS 3
#define NCOLS 3

int
main (void)
{
  uint8_t seed[32] = { 0 };
  unsigned int i, j, o, k, l, n;
  int64_t tmp;
  uint32_t dom = 0;
  long double mean, empiric_mean;
  long double var, empiric_var;

  lazer_init();

  bytes_urandom (seed, sizeof(seed));

  for (k = 1; k <= 3; k++)
    {
      for (i = 1; i <= NROWS; i++)
        {
          for (j = 1; j <= NCOLS; j++)
            {
              INT_T (z, 2);
              INTVEC_T (v, j, 2);
              INTMAT_T (m, i, j, 2);
              mean = 0;
              var = 0;

              for (o = 0; o < NSAMPLES; o++)
                {
                  int_brandom (z, k, seed, dom);
                  dom++;
                  tmp = int_get_i64 (z);
                  mean += tmp;
                  var += (tmp * tmp);

                  intvec_brandom (v, k, seed, dom);
                  dom++;
                  _VEC_FOREACH_ELEM (v, l)
                  {
                    tmp = intvec_get_elem_i64 (v, l);
                    mean += tmp;
                    var += (tmp * tmp);
                  }

                  intmat_brandom (m, k, seed, dom);
                  dom++;
                  _MAT_FOREACH_ELEM (m, l, n)
                  {
                    tmp = intmat_get_elem_i64 (m, l, n);
                    mean += tmp;
                    var += (tmp * tmp);
                  }
                }

              empiric_mean = mean / (NSAMPLES * (i * j + j + 1));
              empiric_var = var / (NSAMPLES * (i * j + j + 1));

              printf ("bin_%d\n", k);
              printf ("expected mean:     0\n");
              printf ("empiric mean:      %Lf\n", empiric_mean);
              printf ("expected variance: %.2f\n", (float)k / 2);
              printf ("empiric variance:  %Lf\n", empiric_var);
              printf ("\n");

              TEST_EXPECT (fabsl (empiric_mean) < 0.01);
              TEST_EXPECT (fabsl (empiric_var - (float)k / 2) < 0.01);
            }
        }
    }

  TEST_PASS ();
}
