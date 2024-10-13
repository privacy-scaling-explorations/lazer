Aggregate signature
===================

We implement an aggregate signature scheme, which combines the lattice-based FALCON :cite:`falcon`
hash-and-sign signature scheme together with a succinct proof of knowledge based on LaBRADOR :cite:`labrador`
implemented in the LaZer library. The construction is simple: user i has a public key polynomial 
:math:`a^{(i)}\in R_p=\mathbb{Z}_{12289}[X]/(X^{512}+1)` and a secret trapdoor basis that allows him, for any :math:`t^{(i)}\in R_p`,
to sample two polynomials :math:`s_1^{(i)},s_2^{(i)}\in R_p` with small coefficients satisfying 

.. math::
    a^{(i)}*s_2^{(i)}+s_1^{(i)} = t^{(i)} \bmod{p}
    :label: falcon

The :math:`t^{(i)}` are the (hashes of the) messages and the polynomials :math:`s_1^{(i)},s_2^{(i)}` are the
signature (actually, only :math:`s_2^{(i)}\in R_p` is the FALCON signature, but :math:`s_1^{(i)}` can be trvially
computed from the other known values in :eq:`falcon`). The aggregator, who sees n equations of the form 
:eq:`falcon` creates a proof that for all :math:`1\leq i\leq n` there exist short :math:`s_1^{(i)},s_2^{(i)}`
satisfying :eq:`falcon`. A verifier who has access to all the public keys and messages is then able to verify
the proof.

Because the succinct proof part of the LaZer library is only optimized to work modulo certain primes (24,32,40, and 48-bit),
it is up to the protocol designer to change the equations that he would like to prove into equivalent
ones over the special primes.  The standard way to convert equations :eq:`falcon` over :math:`R_p` to equivalent ones
over :math:`R_q=\mathbb{Z}_q[X]/(X^{512}+1)` for :math:`q\gg p` is to first rewrite the equation in :eq:`falcon` 
over the ring  :math:`R=\mathbb{Z}[X]/(X^{512}+1)` as   

.. math::
    a*s_2+s_1 = t \bmod{p} \Leftrightarrow \exists v\in R\,\, s.t. a*s_2+s_1 + v*p = t
    :label: falconZ

where the polynomial :math:`v` has coefficients bounded by :math:`\approx \sqrt{d}*\|s_2\|` (this is because
:math:`\|a\| \leq p*\sqrt{d}`, and so by Cauchy-Schwarz, :math:`\|a*s_2\|_\infty \leq p*\sqrt{d}*\|s_2\|`,
with :math:`s_1` playing a negligible role, and so we can ignore it). If we additionally prove that the
norms of :math:`s_1,s_2,v` are small-enough, then 

.. math::
     a*s_2+s_1 + v*p = t \bmod{q} \Rightarrow  a*s_2+s_1 + v*p = t 
    :label: falconBig

because for a large-enough :math:`q`, no coefficient wrap-around occurs.  Thus instead of proving :eq:`falcon`,
the aggregator will be proving 

.. math::
    a^{(i)}*s_2^{(i)}+s_1^{(i)} + p*v^{(i)} = t^{(i)} \bmod{q}
    :label: falconReal

It turns out that a :math:`40`-bit prime is sufficient for ensuring that no wraparound occurs.
We will now go through the aggregate signature example to give intuition for how one would use the succinct
proof part of the LaZer library.

The agg_sig.py file
--------------------

We start by importing the lazer and labrador files (and any other functions that are needed)

.. code-block:: python
    
    from lazer import *
    from labrador import *


We then specify the number of signatures that will be aggregated as well as the :math:`\ell_2^2`-norms of
the polynomials :math:`s_2^{(i)},s_1^{(i)},v^{(i)}` in :eq:`falconReal`. Note that in the actual FALCON signature,
a norm-bound on :math:`(s_1,s_2)` is given, whereas we chose to bound the individual polynomials by dividing that 
bound by :math:`\sqrt{2}`.  The security of FALCON remains exactly the same, but we slightly increase the probability
that the norms are too big, and so we have to try and re-sign.  The norm on :math:`v` is obtained by computing 
its expected value and adjusting slightly slightly so that we don't have to resign more than :math:`1\%` of 
the time. 

.. code-block:: python
    
    sig_num=1024
    norms=[17017363,17017363,round(1248245003*.75)]

We then create a polyring_t class for the FALCON ring :math:`R_p` and also for the modulus that we will use in the 
succinct proof. The modulus needs to correspond to the modulus in one of
the (currently) supported rings, which are specified in the labrador.py file and they are:

LAB_RING_24=polyring_t(64,2**24-3),
LAB_RING_32=polyring_t(64,2**32-99),
LAB_RING_40=polyring_t(64,2**40-195),
LAB_RING_48=polyring_t(64,2**48-59),

which are rings of the form :math:`\mathbb{Z}_{q}[X]/(X^{64}+1)`, for :math:`q=2^{24}-3,2^{32}-99,2^{40}-195,2^{48}-59`.
Note that the ring we will be working with has dimension :math:`512`, but as long as the modulus matches one
of the supported rings, the proof system will use standard algebra to map operations over :math:`\mathbb{Z}_q[X]/(X^{512}+1)`
to those over :math:`\mathbb{Z}_q[X]/(X^{64}+1)`. 

For our aggregate signature example, the :math:`40`-bit modulus is sufficient, and so we select it as the modulus
of the large-modulus ring we will be working with. We will also need the string "40" to be passed to some
functions, and so we store it as PRIMESIZE. Note that a smaller-than-needed modulus is chosen, the proof system
will print out some error messages.

.. code-block:: python
    
    mod=12289
    deg=512
    FALCON_RING=polyring_t(deg,mod)
    BIGMOD_RING=polyring_t(deg,LAB_RING_40.mod)
    PRIMESIZE=str(math.ceil(math.log2(BIGMOD_RING.mod)))

We did not optimize the FACLON key generation procedure, and this is by far the most expensive computation (taking
about :math:`10`X the time of the proof itself), and so for benchmarking purposes, one can just choose to use
the same key for all the signers by setting the SAME_KEY flag to True. We next set up the randomness and create
the :math:`1` polynomial and compute the inverse of :math:`12289` modulo the 40-bit prime :math:`q`.

.. code-block:: python

    shake128 = hashlib.shake_128(bytes.fromhex("00"))
    TARGPP = shake128.digest(32)
    ID=int_to_poly(1,BIGMOD_RING)
    inv_fal_mod=_invmod(mod,BIGMOD_RING.mod)


We now move to setting up the parameters for the proof system. To initialize the proof system, we need
to specify the structure of the entire witness vector :math:`\vec w` (which is the concatenation of all the 
:math:`s_2^{(i)},s_1^{(i)},v^{(i)}`). Every element of :math:`\vec w` needs
to be either a polynomial or a polynomial vector. We need to specify the degree of the ring in which each of
these polynomials (or polynomial vectors) is in. In the case of the aggregate signature scheme, each equation
:eq:`falconReal` consists of three polynomials of degree :math:`512`. Thus the entire witness consists of
3*sig_num witnesses of degree :math:`512`.

.. code-block:: python

    deg_list=[deg]*(3*sig_num)

Next we specify the number of polynomials that each witness has. In our case, each witness consists of one
polynomial, and so we define the num_pols_list variable as a list of 3*sig_num ones.

.. code-block:: python

    num_pols_list=[1]*(3*sig_num)

We then specify the norm bounds (:math:`\ell_2`-squared) that we would like to prove for each of the witnesses. Because we plan to
store the witnesses as

.. math::
    [s_2^{(1)},s_1^{(1)},v^{(1)},s_2^{(2)},s_1^{(2)},v^{(2)},\ldots,s_2^{(sig_num)},s_1^{(sig_num)},v^{(sig_num)}]

the norm bounds are going to be 

.. math::
    [17017363,17017363,round(1248245003*.75)]
    :label: norms

repeated sig_num times, and so we define this list as 

.. code-block:: python

    norm_list=norms*sig_num

We finally define the number of constraints (which is the number of signatures) and initialize the proof statement

.. code-block:: python

    num_constraints=sig_num
    PS=proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)

We now create the public/secret key pairs for the signature schemes. If we decided to use the same key
for all the signatures (for testing purposes), then we just generate one key pair over the ring :math:`R_p`
and then lift it to the ring :math:`R_q`

.. code-block:: python

    skenc,pkenc,pkpol=falcon_keygen()
    l_pk=pkpol.lift(BIGMOD_RING) 

If we are not using the same key, then we simply create a list of sig_num public and private keys, where
the pulic keys are lifted in the same manner as above.

Next we start to create signatures and add them to the proof statement.

First, we create the target (i.e. hash of the message) polynomial :math:`t^{(i)}` in :eq:`falcon` and 
lift it to the ring :math:`R_q` as in :eq:`falconReal`. Normally, :math:`t` would be a hash of a message
(and a nonce), but how the signatures are created is not really relevant to the functionality of the 
succinct proof system and so we ignore this part.

.. code-block:: python

    f_t=poly_t.urandom_static(FALCON_RING,FALCON_RING.mod,TARGPP,0)
    l_t=f_t.lift(BIGMOD_RING)

We then pre-image sample the :math:`s_2^{(i)},s_1^{(i)}` as in :eq:`falcon` and then similarly lift them
to :math:`R_q`.

.. code-block:: python

    l_s1, l_s2 = falcon_preimage_sample(skenc, f_t)
    l_s1=l_s1.lift(BIGMOD_RING)
    l_s2=l_s2.lift(BIGMOD_RING)

Finally, we compute the polynomial :math:`v^{(i)}` in :eq:`falconReal` and centralize it.

.. code-block:: python

    l_v=poly_t(BIGMOD_RING)
    v=(l_t-l_s1-l_pk*l_s2)*inv_fal_mod
    v.redc()

We then check that the norms are less than what they are required to be in :eq:`norms`, and resample if they are not.
We also perform a few "samity check" assertions to make sure that the witnesses satisfy the succinct
proof requirement of having at most :math:`14`-bit coefficients and that the modulus of the ring in 
:eq:`falconReal` is large enough so that the implication in :eq:`falconBig` holds.

We next set up an equation that we will add to the statement. 

.. code-block:: python

    stat_left=[l_pk,ID,ID*mod]
    wit=[l_s2,l_s1,v]
    PS.fresh_statement(stat_left,wit,l_t)

The above means that l_pk*l_s2+1*l_s1+p*v = t in the ring BIGMOD_RING, which is exactly the equation from 
:eq:`falconReal`.  We then loop through adding sig_num equations to the statement.  At the end, we output
the statement stmnt=PS.output_statement(). (In reality, the verifier should compute this independently,
but it does not matter for benchmarking).  

After the loop, we perform a sanity check to make sure that the equation we entered actually satisfies what 
we want to prove by running PS.smpl_verify(). We then produce the zero-knowledge proof and verify it

.. code-block:: python

    proof = PS.pack_prove()
    pack_verify(proof,stmnt,PRIMESIZE)







