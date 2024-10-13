Generating parameters
=====================

..
    The previous :ref:`genparams_c` section explained how to generate C code for proof system parameters from a specification.

First we create a python parameter file and run a sage script to create a C header file that will hold 
the optimized parameters used for the linear ZK proof.

If, for instance, we have the parameter file for :ref:`linrel_ex` example, then we can create a C header file as

.. code-block:: console

    cd scripts
    sage lin-codegen.sage demo-params.py > demo-params.h

The header file defines a variable of type lin_params_t with the name chosen in the specification (“param” in the example).

The params_cffi_build.py script can then be used to compile the generated C code into a python module:

.. code-block:: console

    cd python
    python3 params_cffi_build.py demo_params.h

This will create a python module _demo_params_cffi exporting a lib object that contains the parameter set.
The lib object can contain more than one parameter set (e.g., if there are multiple parameter sets defined in demo_params.h).
Therefore it has a get_params method that can be used to query parameter sets by name.

Here is a python script, that initializes a prover and a verifier from the example's parameter set:


.. code-block:: python

    from lazer import *                 # import lazer python module
    from _demo_params_cffi import lib   # import parameter sets

    seed     = # public randomness
    prover   = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))   

