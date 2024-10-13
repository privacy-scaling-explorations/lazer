#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>

#include "blindsig-demo.h"
#include "lazer.h"

/* use the keypair example from blindsig-demo.h */
#define USE_STATIC_KEYPAIR 1

/* run this many iterations */
#define ITERATIONS 1

int
main (void)
{
  uint8_t privkey[1281];
  uint8_t pubkey[897];

  uint8_t blindsig[3000];
  size_t blindsiglen;

  uint8_t sig[32000];
  size_t siglen;

  uint8_t masked_msg[22000];
  size_t masked_msglen;

  uint8_t message[512 / 8];

  verifier_state_t verifier;
  signer_state_t signer;
  user_state_t user;

  int i, rc;

  lazer_init ();

  // signer_test (); // XXX

  printf ("lazer blind-signature demo\n");
  printf ("--------------------------\n\n");

#if USE_STATIC_KEYPAIR == 1
  memcpy (pubkey, static_pubkey, 897);
  memcpy (privkey, static_privkey, 1281);
#else
  printf ("Generate a random a public/private-keypair ... ");
  fflush (stdout);
  signer_keygen (privkey, pubkey);
  printf ("[OK]\n");
  printf ("keypair (pubkey,privkey): %d bytes\n\n", 1281 + 897);
#endif

  printf ("Initialize user with public key ... ");
  fflush (stdout);
  user_init (user, pubkey);
  printf ("[OK]\n\n");

  printf ("Initialize signer with public and private key ... ");
  fflush (stdout);
  signer_init (signer, pubkey, privkey);
  printf ("[OK]\n\n");

  printf ("Initialize verifier with public key ... ");
  fflush (stdout);
  verifier_init (verifier, pubkey);
  printf ("[OK]\n\n");

  for (i = 0; i < ITERATIONS; i++)
    {
      printf ("User outputs masked message (including a proof of its "
              "well-formedness) ... ");
      fflush (stdout);
      bytes_urandom (message,
                     512 / 8); /* generate a random message for testing */
      user_maskmsg (user, masked_msg, &masked_msglen, message);
      printf ("[OK]\n");
      print_stopwatch_user_maskmsg (0);

      printf ("masked message (t,P1): %lu bytes\n\n", masked_msglen);

      printf ("Signer checks the proof and if it verifies outputs a blind "
              "signature ... ");
      fflush (stdout);
      rc = signer_sign (signer, blindsig, &blindsiglen, masked_msg,
                        masked_msglen);
      if (rc != 1)
        {
          printf ("masked message is invalid.\n");
          exit (EXIT_FAILURE);
        }
      printf ("[OK]\n");

      print_stopwatch_signer_sign (0);

      printf ("blind signature (tau,s1,s2): %lu bytes\n\n", blindsiglen);

      printf ("User outputs a signature on the message ... ");
      fflush (stdout);
      rc = user_sign (user, sig, &siglen, blindsig, blindsiglen);
      if (rc != 1)
        {
          printf ("decoding blindsig failed.\n");
          exit (EXIT_FAILURE);
        }
      printf ("[OK]\n");
      print_stopwatch_user_sign (0);

      printf ("signature (P2): %lu bytes\n\n", siglen);

      printf ("Verfifier verifies the signature on the message ... ");
      fflush (stdout);
      rc = verifier_vrfy (verifier, message, sig, siglen);
      if (rc != 1)
        {
          printf ("verification failed.\n");
          exit (EXIT_FAILURE);
        }
      printf ("[OK]\n\n");
      print_stopwatch_verifier_vrfy (0);
    }

  user_clear (user);
  signer_clear (signer);
  verifier_clear (verifier);
  mpfr_free_cache ();
  return 0;
}
