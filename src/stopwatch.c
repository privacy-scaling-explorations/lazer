#include "stopwatch.h"
#include "lazer.h"

#define INCINDENT 2

STOPWATCH_T (stopwatch_blindsig_user_maskmsg);
STOPWATCH_T (stopwatch_blindsig_signer_sign);
STOPWATCH_T (stopwatch_blindsig_user_sign);
STOPWATCH_T (stopwatch_blindsig_verifier_vrfy);

STOPWATCH_T (stopwatch_lnp_prover_prove);
STOPWATCH_T (stopwatch_lnp_verifier_verify);

STOPWATCH_T (stopwatch_lnp_tbox_prove);
STOPWATCH_T (stopwatch_lnp_tbox_prove_tg);
STOPWATCH_T (stopwatch_lnp_tbox_prove_z34);
STOPWATCH_T (stopwatch_lnp_tbox_prove_auto);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_beta3);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_beta4);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_upsilon);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_bin);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_l2);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_z4);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_z3);
STOPWATCH_T (stopwatch_lnp_tbox_prove_sz_auto);
STOPWATCH_T (stopwatch_lnp_tbox_prove_hi);

STOPWATCH_T (stopwatch_lnp_tbox_verify);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_beta3);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_beta4);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_upsilon);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_bin);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_l2);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_z4);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_z3);
STOPWATCH_T (stopwatch_lnp_tbox_verify_sz_auto);

STOPWATCH_T (stopwatch_lnp_tbox_prove);
STOPWATCH_T (stopwatch_lnp_tbox_verify);

STOPWATCH_T (stopwatch_lnp_quad_many_prove);
STOPWATCH_T (stopwatch_lnp_quad_many_verify);

STOPWATCH_T (stopwatch_lnp_quad_prove);
STOPWATCH_T (stopwatch_lnp_quad_verify);

void
print_stopwatch_user_maskmsg (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_blindsig_user_maskmsg, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_prover_prove (indent + INCINDENT);
}

void
print_stopwatch_signer_sign (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_blindsig_signer_sign, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_verifier_verify (indent + INCINDENT);
}

void
print_stopwatch_user_sign (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_blindsig_user_sign, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_prover_prove (indent + INCINDENT);
}

void
print_stopwatch_verifier_vrfy (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_blindsig_verifier_vrfy, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_verifier_verify (indent + INCINDENT);
}

void
print_stopwatch_lnp_prover_prove (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_prover_prove, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_tbox_prove (indent + INCINDENT);
}

void
print_stopwatch_lnp_verifier_verify (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_verifier_verify, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_tbox_verify (indent + INCINDENT);
}

void
print_stopwatch_lnp_tbox_prove (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove, STOPWATCH_MSEC, indent);

  print_stopwatch_lnp_tbox_prove_tg (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_z34 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_auto (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_beta3 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_beta4 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_upsilon (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_bin (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_l2 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_z4 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_z3 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_sz_auto (indent + INCINDENT);
  print_stopwatch_lnp_tbox_prove_hi (indent + INCINDENT);

  print_stopwatch_lnp_quad_many_prove (indent + INCINDENT);
}

void
print_stopwatch_lnp_tbox_verify (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify, STOPWATCH_MSEC, indent);

  print_stopwatch_lnp_tbox_verify_sz_beta3 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_beta4 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_upsilon (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_bin (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_l2 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_z4 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_z3 (indent + INCINDENT);
  print_stopwatch_lnp_tbox_verify_sz_auto (indent + INCINDENT);

  print_stopwatch_lnp_quad_many_verify (indent + INCINDENT);
}

void
print_stopwatch_lnp_quad_many_prove (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_quad_many_prove, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_quad_prove (indent + INCINDENT);
}

void
print_stopwatch_lnp_quad_many_verify (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_quad_many_verify, STOPWATCH_MSEC, indent);
  print_stopwatch_lnp_quad_verify (indent + INCINDENT);
}

void
print_stopwatch_lnp_quad_prove (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_quad_prove, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_quad_verify (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_quad_verify, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_tg (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_tg, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_z34 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_z34, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_auto (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_auto, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_beta3 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_beta3, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_beta4 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_beta4, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_upsilon (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_upsilon, STOPWATCH_MSEC,
                   indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_bin (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_bin, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_l2 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_l2, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_z4 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_z4, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_z3 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_z3, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_sz_auto (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_sz_auto, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_prove_hi (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_prove_hi, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_beta3 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_beta3, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_beta4 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_beta4, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_upsilon (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_upsilon, STOPWATCH_MSEC,
                   indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_bin (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_bin, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_l2 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_l2, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_z4 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_z4, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_z3 (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_z3, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_lnp_tbox_verify_sz_auto (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_tbox_verify_sz_auto, STOPWATCH_MSEC, indent);
}

/* quad_eval not used anymore by tbox */

STOPWATCH_T (stopwatch_lnp_quad_eval_prove);
STOPWATCH_T (stopwatch_lnp_quad_eval_prove_compute_h);
STOPWATCH_T (stopwatch_lnp_quad_eval_verify);
STOPWATCH_T (stopwatch_lnp_quad_eval_schwartz_zippel_quad);
STOPWATCH_T (stopwatch_lnp_quad_eval_schwartz_zippel_lin);
STOPWATCH_T (stopwatch_lnp_quad_eval_schwartz_zippel_const);

void
print_stopwatch_lnp_quad_eval_prove (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_prove, STOPWATCH_MSEC, indent);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_prove_compute_h, STOPWATCH_MSEC,
                   indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_schwartz_zippel_quad,
                   STOPWATCH_MSEC, indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_schwartz_zippel_lin, STOPWATCH_MSEC,
                   indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_schwartz_zippel_const,
                   STOPWATCH_MSEC, indent + INCINDENT);
  print_stopwatch_lnp_quad_many_prove (indent + INCINDENT);
}

void
print_stopwatch_lnp_quad_eval_verify (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_verify, STOPWATCH_MSEC, indent);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_schwartz_zippel_quad,
                   STOPWATCH_MSEC, indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_schwartz_zippel_lin, STOPWATCH_MSEC,
                   indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_lnp_quad_eval_schwartz_zippel_const,
                   STOPWATCH_MSEC, indent + INCINDENT);
  print_stopwatch_lnp_quad_many_verify (indent + INCINDENT);
}
