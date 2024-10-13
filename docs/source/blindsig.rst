Blind signature
===============

We show how to implement a blind signature from :cite:`BLNS23`.

The scheme works over the polynomial ring :math:`R_p=\mathbb{Z}_p[X]/(X^{512}+1)` for :math:`p=12289`. 
The public parameters consist of random polynomials :math:`a_m,a_\tau\in R_p` and a vector :math:`\vec a_r\in R_p^2`.

The blind signature uses a pre-image sampling procedure which is identical to the one
used in the FALCON digital signature scheme :cite:`falcon`. In order to perform this
pre-image sampling,  The signer additionally creates a FALCON  public key :math:`b\in R_p` and a corresponding secret basis 
with which he is able to, for any :math:`t\in R_p`, create short preimages :math:`s_1,s_2\in R_p`
satisfying 

.. math::
    bs_1+s_2 &= t\\
    \|(s_1,s_2)\| &\leq\sqrt{34034726}

As in the FALCON-512 digital signature scheme, the polynomials :math:`s_i` are sampled from a 
discrete Gaussian distribution with standard deviation approximately 165 (see Table 3.3 of :cite:`falcon`). 

In the first step of the protocol, the user who wants to get a message :math:`m\in\{0,1\}^{512}` 
(viewed as a binary polynomial in :math:`R_p`) signed blindly, picks a short vector :math:`\vec r\in R_p^2` from
some distribution (it will be a discrete Gaussian with standard deviation 1.55*2) and outputs 

.. math::
    t=\vec a_r\cdot \vec r + a_m m
    :label: zk1

together with a Zero-Knowledge proof of knowledge of a binary :math:`m` and a vector :math:`\vec r` such that 
:math:`\|\vec r\|\leq 109`.

Upon receiving the :math:`t` and the zero-knowledge proof, the signer verifies the proof of :eq:`zk1` and if
it is satisfied, generates a binary polynomial :math:`\tau\in R_p` and creates a FALCON pre-image :math:`s_1,s_2` 
satisfying 

.. math::
    bs_1+s_2=t+a_\tau \tau. 
    :label: sigout

The values :math:`\tau,s_1,s_2` are then sent to the user.
The user's digital signature of a message :math:`m` consists of the message and a zero-knowledge proof of knowledge of
:math:`\vec r,s_1,s_2,\tau` satisfying

.. math::
    bs_1+s_2 &= \vec a_r\cdot \vec r + a_m m + a_\tau \tau\\
    \|(s_1,s_2)\| &\leq\sqrt{34034726}\\
    \|\vec r\| &\leq 109\\
    \tau &\in\{0,1\}^{512}
    :label: zk2

---------------

We now move on to see how this scheme is implemented using the LaZer library.
Because there are two zero-knowledge proofs, we will need two parameter files. 
The file setting up the proof in :eq:`zk1` is in the blindsig_p1_params.py file and 
the file setting up the proof in :eq:`zk2` is in the blindsig_p2_params.py file.


The blindsig_p1_params.py file
-------------------------------

We first give the parameter set a name:

.. code-block:: python

    vname = "p1_param"        

We then set the the polynomial ring to :math:`R_p=\mathbb{Z}_p[X]/(X^{512}+1)`, for :math:`p=12289`. The linear
relation involving the witness in :eq:`zk1` can be written as a product of a :math:`1\times 3` matrix and a :math:`3\times 1`
vector. We thus set the dimensions of the public matrix to be :math:`1\times 3` over the ring :math:`R_p`

.. code-block:: python

    deg   = 512                  
    mod   = 12289          
    dim   = (1,3)

The witness satisfying qquation :eq:`zk1` is logically broken down into 2 parts - the vector :math:`\vec r` and the message :math:`m`.
The vector :math:`\vec r` consists of two polynomials satisfying :math:`\|\vec r\|\leq 109`, and :math:`m` is a binary polynomial.
These properties of the witness are expressed as    

.. code-block:: python

    wpart = [ [0,1],   [2] ]  
    wl2   = [   109,     0 ]  
    wbin  = [     0,     1 ]  

.. Notice that the vector :math:`vec r` is chosen fresh every time in the protocol, and so it is not a problem if some information
.. about it is leaked. On the other hand, the message :math:`m` might be some static value that we may want to sign many times
.. and therefore we do not want a small part of it leaking out with every protocol run. We thus do not need to do rejection
.. sampling on :math:`\vec r` and do want to do it on :math:`m`. 

.. .. code-block:: python

..     wrej  = [     0,     1 ]

.. If one is, however using this blind signature in an application where the message is also chosen fresh every time 
.. (e.g. in a digital cash scheme), then perhaps one can indicate that rejection sampling on the message is also unnecessary,
.. and one could make the proof slightly shorter.

The final optional argument, wlinf, indicates the maximum norm in the entire witness. With very high probability, we 
will not have :math:`\|\vec r\|_\infty > 19`, and so we write 

.. code-block:: python

    wlinf=19

The blindsig_p2_params.py file
---------------------------------

The parameters for the second proof are named p2_param, and as before, we are working over the ring :math:`R_p=\mathbb{Z}_p[X]/(X^{512}+1)`, for :math:`p=12289`.

.. code-block:: python
    
    vname = "p2_param"                      

    deg   = 512                            
    mod   = 12289                       

The linear part of equation in :eq:`zk2` can be written as a product of a :math:`1\times 5` 
matrix (which consists of :math:`\vec a, a_\tau, -b, -1`) with a :math:`5-dimensional` witness over :math:`R_p`.
The :math:`\ell_2`-norm of the witness being multiplied with :math:`\vec r` is :math:`109`, the 
witness being multiplied by :math:`a_\tau` is a binary polynomial and the :math:`\ell_2`-norm of 
:math:`(s_1,s_2)` is :math:`\sqrt{34034726}`. 

.. code-block:: python

    wpart = [ [0,1], [2],          [3,4] ]  
    wl2   = [   109,   0, sqrt(34034726) ]  
    wbin  = [     0,   1,              0 ] 

The optional wlinf parameter, which is the bound on the :math:`\ell_\infty`-norm of the witness 
reduces the proof by a few bytes.

.. code-block:: python

    wlinf = 5833 



The blindsig.py file
----------------------

We first import the C functions that will allow us to obtain the autogenerated proof system parameters.
We also import some values from the parameters file 

.. code-block:: python

    from _blindsig_params_cffi import lib
    from blindsig_p1_params import mod, deg, wl2

We then set up the public randomness which will be used to generate the public parameters of the 
blind signature and of the two ZK proof systems.

.. code-block:: python

    shake128 = hashlib.shake_128(bytes.fromhex("00"))
    BLINDSIGPP = shake128.digest(32)
    shake128 = hashlib.shake_128(bytes.fromhex("01"))
    P1PP = shake128.digest(32)
    shake128 = hashlib.shake_128(bytes.fromhex("02"))
    P2PP = shake128.digest(32)

The ring is set to :math:`R_p=R_p=\mathbb{Z}_p[X]/(X^{512}+1)` and then uniformly-random polynomials are created
in the ring.  The polynomial B1 is set as the identity.

.. code-block:: python

    RING = polyring_t(deg, mod) # falcon ring
    BND = int((mod-1)/2)    
    AR1, AR2, AM, ATAU = poly_t(RING), poly_t(RING), poly_t(RING), poly_t(RING)
    AR1.urandom_bnd(-BND, BND, BLINDSIGPP, 0)
    AR2.urandom_bnd(-BND, BND, BLINDSIGPP, 1)
    AM.urandom_bnd(-BND, BND, BLINDSIGPP, 2)
    ATAU.urandom_bnd(-BND, BND, BLINDSIGPP, 3)
    B1 = poly_t(RING, {0: 1})

We then define the "user" class. The user is initialized with the public key of the signer, which he 
converts to the polynomial B2. The first and the second proof systems are then initialized with the public 
randomness and their (auto-generated) parameters.

.. code-block:: python

    def __init__(self, pk: falcon_pkenc):
        #self.B2 = poly_t(RING, falcon_decode_pk(pk))
        self.B2 = falcon_decode_pk(pk)
        self.p1_prover = lin_prover_state_t(P1PP, lib.get_params("p1_param"))
        self.p2_prover = lin_prover_state_t(P2PP, lib.get_params("p2_param"))

The maskmsg function performs the user's first move leading up to :eq:`zk1`.  He takes a 64-byte message, 
converts it into a poly_t, and also creates two polynomials :math:`r_1,r_2` from a discrete Gaussian distribution
with standard deviation :math:`1.55*2^{logsigma}` such that :math:`\|(r_1,r_2)\|\leq wl2[0]=109`. He then computes the 
polynomial :math:`t` as in :eq:`zk1`.

.. code-block:: python

    def maskmsg(self, msg: bytes):
        if len(msg) != 64:
            raise ValueError("msg must be 32 bytes.")

        # sample (r1,r1)
        r1, r2 = poly_t(RING), poly_t(RING)
        m = poly_t(RING, msg)
        seed = secrets.token_bytes(32)  # internal coins
        logsigma = 1                    # sigma = 1.55*2^logsigma
        l2sqr_bnd = wl2[0] * wl2[0]     # l2(r1,r2)^2
        ctr = 0
        while (True):
            r1.grandom(logsigma, seed, ctr)
            ctr += 1
            r2.grandom(logsigma, seed, ctr)
            ctr += 1

            l2sqr = r1.l2sq() + r2.l2sq()
            if (l2sqr <= l2sqr_bnd):
                break

        t = AR1*r1 + AR2*r2 + AM*m

The user then creates the linear equation and the witness with which to produce a zero-kowledge proof of
knowledge of the witness.

.. code-block:: python

    A = polymat_t(RING, 1, 3, [AR1, AR2, AM])
    u = polyvec_t(RING, 1, [-t])
    w = polyvec_t(RING, 3, [r1, r2, m])

Note that :math:`A*w+u=0`, and thus the proof statement, witness, and the proof are generated as follows:

.. code-block:: python

    self.p1_prover.set_statement(A, u)
    self.p1_prover.set_witness(w)
    proof = self.p1_prover.prove()

The user then prepares an encoding of the variable t which he needs to send to the signer. He first instantiates
a coder_t object,

.. code-block:: python

    coder = coder_t()

then initializes it with the maximum number of bytes that are needed (this is something that the implementer of the
scheme needs to compute or approximate), then encodes the variable t in it with the hint to the encoder that t is 
uniformly random modulo mod, and then closes the encoder creating a byte stric on the encoding of t.

.. code-block:: python

    coder = coder_t()
    coder.enc_begin(22000)
    coder.enc_urandom(mod, t)
    tenc = coder.enc_end()

Because the variables :math:`r_1,r_2,m` are going to be needed in the next move of the user (i.e. when he produces
the ZK proof for :eq:`zk2`), we save these variables for later in the user object.

.. code-block:: python

    self.r1 = r1
    self.r2 = r2
    self.m = m

The function then returns the concatenation of the encoding of t with the ZK proof

.. code-block:: python

    return tenc + proof

The sign function of the user corresponds to his second move in which he ends up producing the proof in :eq:`zk2`. 
In this move, he receives from the signer a polynomial :math:`\tau` and polynomials :math:`\tau,s_1,s_2` which satisfy :eq:`sigout`. Because
these polynomials are distributed according to a discrete Gaussian distribution with standard deviation approximately
165, they are encoded using the enc_grandom function, and thus decoded using the dec_grandom function. The polynomial 
:math:`\tau` is binary and so its representation is most naturally encoded as a byte string. The user thus
decodes the byte string :math:`\tau` and the polynomials :math:`s_1,s_2`, and then converts the byte string to a 
polynomial variable :math:`tau\_`.

.. code-block:: python

    try:
        coder = coder_t()
        coder.dec_begin(blindsig)
        coder.dec_bytes(tau)
        coder.dec_grandom(165, s1)
        coder.dec_grandom(165, s2)
        coder.dec_end()
    except DecodingError:
        raise InvalidMaskedMsg

    tau_ = poly_t(RING, tau)

Now the user needs to set up the linear equation for which he would like to generate a ZK proof.  He creates 
the matrix :math:`A` and vectors :math:`\vec u,\vec w` such that :math:`A\vec w+\vec u=0` corresponding to 
:eq:`zk2`. The user then sets the statement and witness, and creates the ZK proof.

.. code-block:: python

    A = polymat_t(RING, 1, 5, [AR1, AR2, ATAU, -B1, -self.B2])
    u = polyvec_t(RING, 1, [AM * self.m])
    w = polyvec_t(RING, 5, [self.r1, self.r2, tau_, s1, s2])
    self.p2_prover.set_statement(A, u)
    self.p2_prover.set_witness(w)
    proof = self.p2_prover.prove()

We now move on to the signer class. The signer is initialized with a FALCON secret key, which he stores. Since the
signer will need to verify the first ZK proof from :eq:`zk1`, he also sets up a verification state using the
public parameter P1PP and the parameters from the blindsig_p1_params.py file.

.. code-block:: python

    def __init__(self, sk: falcon_skenc):
        self.sk = sk
        self.p1_verifier = lin_verifier_state_t(P1PP, lib.get_params("p1_param"))

In the sign function, the signer needs to verify the ZK from :eq:`zk1` and then output the :math:`\tau,s_1,s_2` 
satisfying :eq:`sigout`.  The input to the sign function is the encoding of :math:`t` and the proof produced by the user in the
maskmsg function. The signer performs the reverse of the encoding procedure to get :math:`t`
and also creates an integer tlen, which gives the length of the encoding. 

.. code-block:: python

    t = poly_t(RING)
    try:
        coder = coder_t()
        coder.dec_begin(masked_msg)
        coder.dec_urandom(mod, t)
        tlen = coder.dec_end()
    except DecodingError:
        raise InvalidMaskedMsg

Next the signer creates :math:`A,\vec u`, which are the public part of the linear equation :math:`A\vec w+\vec u=0`.
He also extracts the proof from the masked_msg variable (since masked_msg was a concatenation of the encoding of 
:math:`t` and the proof, and we know that tlen is the number of bytes in :math:`t`, we can recover the proof). 
The signer then sets the statement and tries to verify the proof.

.. code-block:: python

    A = polymat_t(RING, 1, 3, [AR1, AR2, AM])
    u = polyvec_t(RING, 1, [-t])
    proof = masked_msg[tlen:]

    self.p1_verifier.set_statement(A, u)
    try:
        self.p1_verifier.verify(proof)
    except VerificationError:
        raise InvalidMaskedMsg("Masked message invalid.")

Assuming that verification passed, the signer creates a binary polynomial :math:`\tau`
and then sampling a short pre-image using the FALCON secret key sk :math:`s_1,s_2` satisfying
:eq:`sigout`.

.. code-block:: python

    tau_ = secrets.token_bytes(64)
    tau = poly_t(RING, tau_)

    s1, s2 = falcon_preimage_sample(self.sk, ATAU * tau + t)

The byte-string :math:`tau\_`, which is used to create the polynomial :math:`\tau` is then
encoded along with :math:`s_1,s_2` and outputted. Since :math:`s_1,s_2` are discrete Gaussians
with standard deviation around 165, we specify that we would like to use an encoding for
Gaussians (one could use the uniform encoding, but it will be somewhat longer).

.. code-block:: python

    coder = coder_t()
    coder.enc_begin(2000)
    coder.enc_bytes(tau_)
    coder.enc_grandom(165, s1)
    coder.enc_grandom(165, s2)
    blindsig = coder.enc_end()

We now get to the verifier who is going to be verifying the final signature from :eq:`zk2`.
The verifier is initialized with a FALCON public key and he also creates the public 
parameters for verifying :eq:`zk2`.   

.. code-block:: python

    def __init__(self, pk: falcon_pkenc):
        self.B2 = falcon_decode_pk(pk)        
        self.p2_verifier = lin_verifier_state_t(P2PP, lib.get_params("p2_param"))

The verification function takes as input the message byte string and a byte string representing the 
ZK proof. It converts the byte string to a poly_t type and then sets up the matrix 
:math:`A` and vector :math:`\vec u` such that :math:`A\vec w+\vec u=0` for the witness in 
:eq:`zk2`. It then runs the verification to check that the proof is correct. 

.. code-block:: python 

    m = poly_t(RING, msg)

    # verify proof for PoK(w): Aw + u = 0
    A = polymat_t(RING, 1, 5, [AR1, AR2, ATAU, -B1, -self.B2])
    u = polyvec_t(RING, 1, [AM * m])

    self.p2_verifier.set_statement(A, u)
    try:
        self.p2_verifier.verify(sig)
    except VerificationError:
        raise InvalidSignature("Signature invalid.")

We now briefly look at the main() function which controls the flow of the blind signature
protocol. First, a random message is chosen and we sample the secret and public FALCON keys.
In practice, these keys would be sampled by the signer.

.. code-block:: python

    sk, pk, _ = falcon_keygen()

We then create the user, signer, and verifier. 

.. code-block:: python

    user = user_t(pk)
    signer = signer_t(sk)
    verifier = verifier_t(pk)

The user's first move in the blind signature is then executed and the masked_msg output produced.

.. code-block:: python

    masked_msg = user.maskmsg(msg)

The signer then performs his move and creates the output in :eq:`sigout` and encodes it into
the variable blindsig.

.. code-block:: python

    blindsig = signer.sign(masked_msg)

The user then takes this output and performs his second move by creating the ZK proof 
in :eq:`zk2`, which the verifier verifies.

.. code-block:: python

    sig = user.sign(blindsig)
    verifier.verify(msg, sig)

