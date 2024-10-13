Kyber1024
=========

We show how to prove knowledge of a secret Kyber1024 key.

Define the ring :math:`R_p=\mathbb{Z}_p[X]/(X^{256}+1)` for :math:`p=3329`.
The Kyber-1024 secret key consists of a vector :math:`\vec s\in R_p^8` which is distributed as

.. math::
    \vec s=\vec s_1+\vec s_2 - \vec s_3 -\vec s_4
    :label: skgen

where all the coefficients of each :math:`\vec s_i`
are chosen uniformly and independently from :math:`\{0,1\}`.

The Kyber-1024 public key consists of a matrix :math:`A=[~A_1~ | ~I~]`, where :math:`A_1` is chosen uniformly at
random from :math:`R_p^{4\times 4}` (or more specifically, it is expanded using a cryptographic hash
function from some public seed :math:`\rho`), :math:`I\in R_p^{4\times 4}` is the identity matrix, and

.. math::
    \vec t=A\vec s.
    :label: pk

Suppose that we would like to prove that the pair :math:`(A,\vec t)` is a valid Kyber-1024 public key. More specifically,
we would like to show that there exists a vector :math:`\vec s` distributed as in :eq:`skgen` satisfying :eq:`pk`.
The reason that one may want to perform such a proof is to avoid attacks that can be mounted in protocols where
some party uses a secret key (or ciphertext randomness, which requires proving a very similar relation) that is too large

As a simple example where such a proof is useful is in a non-interactive key exchange where two users,
having secret keys :math:`\vec s_1,\vec s_2` can use (something very similar to) their Kyber
public keys (and secret keys) to exchange a secret random string (c.f. :cite:`swoosh`). The
probability that this procedure fails is directly related to the
size of the coefficients of :math:`\langle \vec s_1,\vec s_2\rangle` -- i.e. the larger the
coefficeints of the inner product of their secret keys, the larger the probability that the procedure fails is.
Furthermore, the security of each user depends on the procedure never failing. An adversary who tries to attack
user 1's key, may generate an invalid public key by making his secret :math:`\vec s_2` longer and thus make the
probability that

.. math::
    \langle \vec s_1,\vec s_2\rangle
    :label: errors

has large coefficients. A very similar attack can
occur in other scenarios where the adversary may want to create an invalid ciphertext whose randomness is too
large; and avoiding such an attack requires one to prove that the randomness used in creating the ciphertext
is not large.

The question is now how should "large" be defined. Naively, one might just want to prove that :math:`\vec s` is a
"valid" Kyber secret key -- that is, it has coefficients between -2 and 2. But proving that
:math:`\|\vec s\|_\infty \leq 2` is, for most purposes however, sub-optimal or even wrong. The reason is that, assuming
that one of the vectors in :eq:`errors` is honestly distributed as in :eq:`skgen`, then the size of the coefficients in
:eq:`errors` depends on the :math:`\bf{\ell_2}`-norm of :math:`\vec s_1,\vec s_2`. Furthermore, the :math:`\ell_2`-norm
of :math:`\vec s` chosen as in :eq:`skgen` is very tightly concentrated around :math:`\sqrt{n}`, where :math:`n` is The
number of integer coefficients making up :math:`\vec s` (so in the case of Kyber-1024, it is tightly concentrated around
:math:`\sqrt{2048}`).

So if the goal of the adversary is to choose his :math:`\vec s_i` to make the coefficients of :eq:`errors` too large,
he should choose it with a large :math:`\ell_2`-norm. Observe that if we only require that the :math:`\ell_\infty`
norm of :math:`\vec s_i` is at most 2, the adversary can set all the coefficients to 2, making the :math:`\ell_2`-norm
double what one would expect if the secret key were honestly generated (the probability that an honestly-generated
vector has such a large norm is :math:`2^{-3\cdot 2048}`). In short, the :math:`\ell_\infty`-norm
doesn't matter at all -- in fact, it will even be detrimentary to adversary trying to maximize the coefficients
in :eq:`errors` to choose a vector that has a large :math:`\ell_\infty`-norm (when the :math:`\ell_2`-norm is fixed).

So in our proof of a correctly-generated Kyber-1024 key, we will be proving that the size of the coefficients of :eq:`errors`
has a small :math:`\ell_2`-norm. The structure of the Python code is very similar to that of the "demo_params.py" and "demo.py" files,
which are also given as examples.



The kyber1024_params.py file
----------------------------

We first give the parameter set a name, which we will import in the actual proof system file

.. code-block:: python

    vname = "param"

We then set the polynomial ring to :math:`R_p=\mathbb{Z}_p[X]/(X^{256}+1)`, for :math:`p=2^{32}-4607`,
and set the dimensions of the matrix :math:`A` to be :math:`4\times 8` over the ring :math:`R_p`

.. code-block:: python

    deg   = 256
    mod   = 3329
    dim   = (4,8)

One then needs to explain what one would like to prove about the witness :math:`\vec s` that is supposed to
satisfy :math:`A\vec s-\vec t=0`. In case of Kyber-1024, we will prove something about the entire witness. As
discussed above, the :math:`\ell_2`-norm of :math:`\vec s` is tightly concentrated around
:math:`\sqrt{2048}`
and with probability negligibly-close to 1, it will be less than :math:`1.2\sqrt{2048}` (depending on the application,
one can set this value to even be exactly :math:`\sqrt{2048}`, and then simply require that the party generating
:math:`\vec s` tries a few times until it is small-enough -- this does not meaningfully change the security of the scheme).

We set the partition list of lists -- wpart, to be all
8 polynomials of :math:`\vec s`, and set the wl2 list to :math:`1.2*\sqrt{2048}`.  Since we do not wish to prove that the
coefficients of :math:`\vec s` are binary, we set the wbin element to 0.

.. code-block:: python

    wpart = [ list(range(0,8)) ]
    wl2   = [ 1.2*sqrt(2048)       ]
    wbin  = [ 0                ]

Finally, we have the optional variable wlinf. If one has a bound on the :math:`\ell_\infty`-norm of the witness,
then the parameters of the ZK proof can be set slightly tighter and the proof will be a bit shorter. In the case of
Kyber-1024, all the coefficients are at most 2 in magnitude, and so we set

.. code-block:: python

    wlinf  = 2

We would like to point out that this *does not* imply that we will be additionally proving that the :math:`\ell_\infty`-norm
is at most 2. It's just information that allows the parameters of the proof to be a little smaller.

The kyber1024.py file
---------------------

We start by declaring public randomness that is used to derive the seeds
for generating the Kyber key and the public randomness used in the
commitment scheme.  We point out that we do not derive the kyber public
key in the same way as in the official kyber standard (bu from the same
distribution).

.. code-block:: python

    shake128 = hashlib.shake_128(bytes.fromhex("00"))
    KYBERPP = shake128.digest(32)   # kyber public randmoness
    shake128 = hashlib.shake_128(bytes.fromhex("01"))
    P1PP = shake128.digest(32)      # proof system public radomness

We then import the parameters from the parameter file.

.. code-block:: python

    from kyber1024_params import mod, deg, m, n     # import kyber parameters
    from _kyber1024_params_cffi import lib          # import proof system parameters
    prover = lin_prover_state_t(P1PP, lib.get_params("param"))
    verifier = lin_verifier_state_t(P1PP, lib.get_params("param"))

We then create the matrix :math:`A=[~A_1~ | ~I~]`, the secret key :math:`sk`
according to the binomial distribution as in :eq:`skgen`, and the second part
of the public key as in :eq:`pk`.

.. code-block:: python

    R = polyring_t(deg, mod)
    A1 = polymat_t.urandom_static(R, m, m, mod, KYBERPP, 0)
    A2 = polymat_t.identity(R, m)
    A = polymat_t(R, m, n, [A1, A2])
    sk = polyvec_t.brandom_static(R, n, 2, secrets.token_bytes(32), 0)
    pk = A*sk

The prover is then initialized to prove the linear relation :math:`A*sk+pk=0`
corresponding to :eq:`pk`. And then generate the ZK proof byte array.

.. code-block:: python

    prover.set_statement(A, -pk)
    prover.set_witness(sk)
    proof = prover.prove()

Then the verifier is initialized with the same public parameters and tries
to verify the ZK proof.

.. code-block:: python

    try:
        verifier.verify(proof)
    except VerificationError:
        print("reject")
    else:
        print("accept")
