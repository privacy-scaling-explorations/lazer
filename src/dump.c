#include <stdio.h>
#include <mpfr.h>

void _dump_mpfr (mpfr_srcptr op);

void _dump_mpfr (mpfr_srcptr op)
{
  mpfr_out_str (stdout, 10, 0, op, MPFR_RNDN);
  fprintf (stdout, "\n");
  fflush (stdout);
}
