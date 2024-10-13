Anonymous credentials
=====================

We show how to implement an anonymous credentials system from :cite:`BLNS23`, which is similar to the blind
signature from the same paper, and so we urge to look at the blind signature implementation (and explanation) first.
One of the main differences in the implementation is the fact that instead of all operations of the scheme
being perfomred over the ring :math:`R_p=\mathbb{Z}_p[X]/(X^{512}+1)`, they are now 
done over the ring :math:`R_p=\mathbb{Z}_p[X]/(X^{64}+1)`. The reason for this change is that in an anonymous
credential scheme, the number of credentials may be large, so the length of the "message" might be larger than
in the blind signature, which would increase the size of the proof. But each credential does not require a 
512-dimensional vector to store -- it (or its hash) can be represented as binary vectors that can fit into one (or a few) 
64-degree polynomials. Thus it makes sense to work over smaller-degree rings. In the instantiation, this "simply" requires
us to use polynomial vectors of degree :math:`512/64=8` in order to achieve the same security of the underlying
lattice problem.

We will now go over the anonymous credential scheme in :cite:`BLNS23` and then explain the relevant changes that 
need to be made in the previously-described blind signature scheme.

The scheme works over the polynomial ring :math:`R_p=\mathbb{Z}_p[X]/(X^{64}+1)` for :math:`p=12289`. 
The public parameters consist of random matrices :math:`A_m,A_\tau\in R_p^{8\times 8}` and a 
matrix :math:`A_r\in R_p^{8\times 16}`.

The anonymous credential scheme uses a pre-image sampling procedure which is identical to the one
used in the FALCON digital signature scheme :cite:`falcon`. In order to perform this
pre-image sampling,  The signer additionally creates a FALCON public key :math:`b\in \mathbb{Z}_p[X]/(X^{512}+1)` and a corresponding secret basis 
with which he is able to, for any :math:`t\in \mathbb{Z}_p[X]/(X^{512}+1)`, create short preimages :math:`s_1,s_2\in \mathbb{Z}_p[X]/(X^{512}+1)`
satisfying 

.. math::
    bs_1+s_2 &= t\\
    \|(s_1,s_2)\| &\leq\sqrt{34034726}

Note that now we have two rings -- the ring :math:`R_p`, and the FALCON ring :math:`R_F=\mathbb{Z}_p[X]/(X^{512}+1)`. 
Because :math:`R_p` is a subring of the FALCON ring, one can define a (length-preserving and invertible) ring homomorphism mapping
polynomials from :math:`R_F` to dimension-8 vectors over :math:`R_p`. While the exact details 
are not needed to use the LaZer library, the standard such homomorphism is desrcibed in Section 2.8 of :cite:`L21`.
The important thing to know is that addition in the ring over :math:`R_p^8` is the ususal coordinate-wise addition,
while multiplication of two elements :math:`f,g\in R_p^8` can be written as the usual polynomial 
matrix-vector product :math:`M_f*g` over :math:`R_p`, where :math:`M_f\in R_p^{8\times 8}` and is derived from :math:`f`
in some particular way (again see :cite:`L21`).  

The LaZer library can output the FALCON public keys either as polynomials over :math:`R_F` or 
as matrices over :math:`R_p^{k\times k}`. It can have the pre-image sampling function similarly take inputs and produce
outputs that are vectors over :math:`R_p^k`.  In particular, suppose that :math:`sigma` is the above-mentioned ring
homomorphism :math:`\sigma: R_F \rightarrow R_p^k`.  We can then use the falcon_decode_pk function to take 
the falcon public key and produce a matrix :math:`B\in R+p^k\times k`; and we can then tell the pre-image sampling function that the target
:math:`t\in R_p` is in the ring :math:`R_p`, and to produce pre-images :math:`\vec s_1,\vec s_2\in R_p^k` such that 
:math:`B\vec s_1+\vec s_2=\vec t`. If the FALCON public key over :math:`R_f` is :math:`b`, then 

.. math::
    b\sigma^{-1}(\vec s_1)+\sigma^{-1}(\vec s_2) = \sigma^{-1}(\vec t)

over the ring :math:`R_F`. 

In the first step of the protocol, the user who wants to get 8 credentials :math:`m\in R_p^{8}` 
(each of which is a binary polynomial in :math:`R_p`) signed blindly, picks a short vector :math:`\vec r\in R_p^2` from
some distribution (it will be a discrete Gaussian with standard deviation 1.55*2) and outputs 

.. math::
    \vec t=A_r\cdot \vec r + A_m \vec m
    :label: zk1_anon

together with a Zero-Knowledge proof of knowledge of a binary :math:`m` and a vector :math:`\vec r` such that 
:math:`\|\vec r\|\leq 109`.

Upon receiving the :math:`\vec t` and the zero-knowledge proof, the signer verifies the proof of :eq:`zk1` and if
it is satisfied, generates a binary polynomial vector :math:`\vec \tau\in R_p^8` and creates a FALCON pre-image 
:math:`\vec s_1,\vec s_2` (using the above-described homomorphism) satisfying 

.. math::
    B\vec s_1+\vec s_2=\vec t+A_\tau \vec \tau. 
    :label: sigout_anon

The values :math:`\vec \tau,\vec s_1,\vec s_2` are then sent to the user.
If the user would like to open some :math:`\alpha` credentials that are part of :math:`m\in R_p^8`, then he can reveal 
the corresponding part (call it :math:`m_1\in R_p^\alpha`), and 
give a zero-knowledge proof of knowledge of
:math:`\vec r,\vec s_1,\vec s_2,\vec \tau,\vec m_2` satisfying

.. math::
    B\vec s_1+\vec s_2 &= A_r\cdot \vec r + A_{m_1} \vec m_1 + A_{m_2} \vec m_2+ A_\tau \vec \tau\\
    \|(\vec s_1,\vec s_2)\| &\leq\sqrt{34034726}\\
    \|\vec r\| &\leq 109\\
    \vec \tau &\in\{0,1\}^{512}
    :label: zk2_anon

where :math:`\vec m_2` is the part of :math:`\vec m` that remains hidden and :math:`A_{m_i}` are the corresponding 
parts of :math:`A_m`.

We will now go through the part of the anon_cred.py python program that are different to those in blindsig.py
due to the relationship between the FALCON ring and the ring :math:`R_p` over which the rest of the protocol uses.

In the signer class, we would like to pre-image sample vectors :math:`\vec s_1,\vec s_2` for an image 
:math:`A_\tau\vec \tau +\vec t`, which is done by calling 
the following function:

.. code-block:: python

    s1, s2 = falcon_preimage_sample(self.sk, ATAU * tau + t,RING)

Notice that ATAU * tau + t is an element of :math:`R_p^8`, and so the function needs to first 
convert ATAU * tau + t to an element in the FALCON ring as :math:`t'=\sigma^{-1}(A_\tau\vec \tau+\vec t)`, then 
do pre-image sampling to obtain :math:`s_1',s_2'\in R_F`, and then finally output 
:math:`\vec s_1 = \sigma(s_1'),\vec s_2=\sigma(s_2')`. The fact that we would like the last conversion done is 
specified by the optional argument RING, which is set to be :math:`R_p`.

The only other major difference between the blind signature scheme and the anonymous credential one is the 
fact that the "message" :math:`\vec m` is not revealed in its entirety in the latter. As in :eq:`zk2_anon`, the message
is split into the revealed public part :math:`\vec m_1` and the hidden part :math:`\vec m_2`.  The hidden part :math:`\vec m_2`
becomes part of the witness; but because we need to declare the witness size ahead of time in the parameters file,
we cannot have :math:`\vec m_2` be a vector of some variable number of dimensions. This is why we instead zero-out the
columns of the matrix :math:`A_m` that correspond to the positions of the public vector 

.. code-block:: python

    AM_priv=AM.zero_out_cols(pub_mvec)

and then simply set :math:`\vec m_2=\vec m` (or set  :math:`\vec m_2` to be :math:`\vec m` with the public positions zeroed-out).



