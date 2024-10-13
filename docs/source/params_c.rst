.. _genparams_c:

Generating parameters
=====================

Most proof system parameters are constants that are fixed at compile time.

This section explains how to ...
 - specify your proof system parameters in python, using an easy syntax.
 - use sage scripts to generate C code from these specifications.

The generated C code is basically a single structure collecting all constants corresponding to the specification.


Linear relations and norms
--------------------------

:ref:`linrel_ex` example can be generalized to specify proof system parameters for proving linear relations with norms.

The lin-codegen.sage script can then be used to generate a C header file demo-params.h from the specification in demo-params.py:

.. code-block:: console

    cd scripts
    sage lin-codegen.sage demo-params.py > demo-params.h

The header file defines a variable of type lin_params_t with the name chosen in the specification ("param" in the example).

Here is an example C program, that initializes a prover and a verifier from the example's parameter set:

.. code-block:: C

    #include "lazer.h"
    #include "demo-params.h" // variable "params" defined here

    int main(void)
    {
      uint8_t pp[32] = { /* public randomness */ };
      lin_prover_state_t prover;
      lin_verifier_state_t verifier;

      lin_prover_init (prover, pp, params); 
      lin_verifier_init (verifier, pp, params); 

      /* do stuff ... */

      lin_prover_clear (prover); 
      lin_verifier_clear (verifier); 
      return 0;
    }



