Linear relations and norms
==========================

In this example, found in lazer/python/demo, we prove the simple norm-bounded linear relation :math:`As+t=0` with
:math:`\|s\|\leq \sqrt{2048}` over some polynomial ring.

.. _linrel_ex:

The demo_params.py file
-----------------------

We first give the parameter set a name, which we will import in the actual proof system file 

.. code-block:: python

    vname = "param"               

We then set the polynomial ring to :math:`R_p=\mathbb{Z}_p[X]/(X^{256}+1)`, for :math:`p=2^{32}-4607`,
and set the dimensions of the matrix :math:`A` to be :math:`4\times 8` over the ring :math:`R_p`

.. code-block:: python

    deg   = 256                  
    mod   = 2**32 - 4607          
    dim   = (4,8)                 

One then needs to explain what one would like to prove about the witness :math:`s` that is supposed to 
satisfy :math:`As+t=0`. In this example, we will prove something about the entire witness. In particular,
we would like to show that :math:`\|s\|\leq \sqrt{2048}`. For this we set the partition list of lists -- wpart, to be all 
8 polynomials of :math:`s`, and set the wl2 list to :math:`\sqrt{2048}`.  Since we do not wish to prove that the
coefficients of :math:`s` are binary, we set the wbin element to 0.

.. code-block:: python

    wpart = [ list(range(0,8)) ]  
    wl2   = [ sqrt(2048)       ]  
    wbin  = [ 0                ]

Suppose that we wanted to do something more complicated. Say that :math:`s` is viewed concatenation of 3 vectors,
:math:`s_1 s_2 s_3`, each of 4,3, and 1 polynomials, respectively, and we want to 
prove that :math:`\|s_1\|\leq \sqrt{1024}`, that :math:`s_2` is binary, and that :math:`\|s_3\|\leq \sqrt{64}`.
We would then set the lists wpart,wl2, and wbin as follows:  

.. code-block:: python

    wpart = [ list(range(0,4)), list(range(5,7)),[7] ]  
    wl2   = [ sqrt(1024), 0, sqrt(64)       ]  
    wbin  = [ 0, 1, 0                ]

.. Going back to the demo_params.py file, the next line states whether rejection sampling is to be used for 
.. proving the corresponding secret. One does not need to use rejection sampling in the case that the secret
.. is some vector which is chosen fresh every time in the protocol - for example the randomness of a commitment.
.. Rejection sampling should, on the other hand, be used if the vector is some static secret. The coefficients
.. of the list wrej are set to 1 if rejection sampling is to be performed and 0 otherwise. 

.. .. code-block:: python

..     wrej  = [ 1                ]

Finally, we have the optional variable wlinf. If one has a bound on the :math:`\ell_\infty`-norm of the witness,
then the parameters of the ZK proof can be set slightly tighter and the proof will be a bit shorter. In the case of
the demo.py example, we will be setting the coefficients of the witness to be -1,0, or 1, and so we set  

.. code-block:: python

    wlinf  = 1

We would like to point out that this *does not* imply that we will be additionally proving that the :math:`\ell_\infty`-norm
is at most 1. It's just information that allows the parameters of the proof to be a little smaller.

The demo.py file
----------------

We need to import the "lazer" module 

.. code-block:: python

    from lazer import * 

and if it's not in the same directory, also set the appropriate path

.. code-block:: python
    
    import sys
    sys.path.append('..') 

The parameters from the demo_params.h file are then imported -- the modulus, ring degree, and the matrix dimensions

.. code-block:: python

    from demo_params import mod, deg, dim
    d, p, m, n = deg, mod, dim[0], dim[1]

The variable "seed" is a 32-byte public seed that is shared between the prover and verifier. It is used
to create things like the random matrix used in the underlying commitment scheme. One can 
choose it randomly from the system randomness, if one wishes.  We then import the automatically-generated
file _demo_params_cffi as lib, and then initialize the prover and the verifier using the "params"
variables that were specified in the demo_params.py file. This is the only time that the user 
needs to touch the C interface.

.. code-block:: python

    seed = b'\0' * 32
    from _demo_params_cffi import lib
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))

The following line declares Rp to be the polynomial ring :math:`R_p=\mathbb{Z}_p[X]/(X^d+1)`

.. code-block:: python

    Rp = polyring_t(d, p)

We then create a random matrix :math:`A\in R_p^{m\times n}`

.. code-block:: python

    A = polymat_t(Rp, m, n)
    A.urandom(p, seed, 0)

One similarly creates a random vector :math:`s\in R_p^n` with coefficients in :math:`\{-1,0,1\}`,
having Bernoulli distribution, and then sets :math:`t=-A*s` 

.. code-block:: python
    
    s = polyvec_t(Rp, n)
    s.brandom(1, seed, 0)
    t = -A*s

The prover then sets the statement that he would like to prove: :math:`As+t=0`, the witness :math:`s`,
and creates the ZK proof of knowledge specified in the demo_params.py file

.. code-block:: python

    prover.set_statement(A, t)
    prover.set_witness(s)
    proof = prover.prove()

The verifier sets the same statement and checks whether the proof output by the prover is valid

.. code-block:: python

    verifier.set_statement(A, t)
    verifier.verify(proof)
