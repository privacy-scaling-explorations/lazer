# Set proof system parameters for proving knowledge
# of an MLWE secret (s,e) over Rp:
#
#   t = A*s + e
#   l2(s,e) <= B
#
# A in Rp^(n x m), t,e in Rp^n, s in Rp^m,
# Rp = Zp[X]/(X^d + 1)

name = "falcon_params"       # variable name

dprime = 512                     # Kyber ring degree
p = 12289                    # Kyber modulus

dimn = 1                       # (2,3,4) for Kyber(512,768,1024)
dimm = 1                       # (2,3,4) for Kyber(512,768,1024)

S = 1                       # linf(s) <= S, (3,2,2) for Kyber(512,768,1024)
E = 1                     # linf(e) <= E, (3,2,2) for Kyber(512,768,1024)

# not affected by ring iso: dimensions n and m are multiplied by
# k but degree dprime is divided by k
#BOUND = sqrt(S ** 2 * dimm * dprime + E ** 2 * dimn * dprime) # XXX tigther bound 1.1*expected, binomial secrets
BOUND = 32

# Set LNP proof system parameters to prove knowledge of (s,v)
# over Rq:
#
#   l2(s, t-As-pv) <= B i.e., linf([[1,0],[-A,-p1]]*(s,v) + (0,t)) <= B
#   linf(v) <= psi*Bv
#
# where q is set large enough so that computations in R
# do not wrap mod q, such that above statement holds over
# the integers, which again implies the original statement
# over Rp.
#
# Notation:
#   R = Z[X]/(X^d+1), Rq = Zq[X]/(X^d+1) polynomial rings,
#    d is either 64 or 128, q is either q1 or q1*q2,
#    q1, q2 odd primes congruent to 5 mod 8, q1 > q2.
#
#   o in Aut(R), automorphism defined by o(X) = X^-1.
#
#   <x> = poly||o(x), x in R, works elementwise on vectors and
#    matrices of polynomials.

# Prove knowledge of polyvec (s1,m) in Rq^(m1+l), l2(s1) <= alpha,
# and prove statements (1-5).

d = 64                  # set ring degree (either 64 or 128)
k = dprime/d

PSI = 2238
P = (p-1)/2
# not affected by k since dimension m is multiplied by k
# but degree dprime is divided by k
# set ring modulus q's bit-size
log2q = ceil(log(2 * (1+PSI)*(P + dimm*P*dprime*S + E) + 1, 2))   # XXX
# log2q1 = 0            # if q = q1*q2, also set q1's bit-size
m1 = k*(dimm+dimn+dimn)              # set length of vector s1
alpha = BOUND        # set l2-norm bound on vector s1
l = 0                  # set length of vector m

# 1. Quadratic relations over Rq with automorphisms / inner products over Rq:
#  (s1,m) and its o-automorphisms satisfiy N quadratic equations over Rq.
#
#  For i in [1,N]:
#
#  o<s1,m>^T * R2[i] * o<s1,m> + r1[i] * o<s1,m> + r0[i] = 0
#
#  polymat R2[i] in Rq^(2(m1+l) x 2(m1+l))
#  polyvec r1[i] in Rq^(2(m1+l))
#  poly    r0[i] in Rq

# 2. Quadratic relations over Zq with automorphisms / inner products over Zq:
#  (s1,m) and its o-automorphisms satisfy M quadratic equations over Zq.
#
#  For i in [1,M]:
#
#  constant coefficient of
#  <s1,m>^T * R2prime[i] * <s1,m> + r1prime[i] * <s1,m> +r0prime[i] = 0
#
#  polymat R2prime[i] in Rq^(2(m1+l) x 2(m1+l))
#  polyvec r1prime[i] in Rq^(2(m1+l))
#  poly    r0prime[i] in Rq

# 3. Exact infinity norm:
#  The evaluation of a linear function Rq^(m1+l) -> Rq^nbin at (s1,m)
#  has binary coefficients.
#
#  coefficients of
#  Ps * s1 + Pm * m + f in {0,1}^(nbin*d)
#
#  polymat Ps in Rq^(nbin x m1)
#  polymat Pm in Rq^(nbin x l)
#  polyvec f in Rq^(nbin)

# nbin =                # set length of vector with binary coefficients

# 4. Exact euclidean norm:
#  The evaluations of a linear functions Rq^(m1+l) -> Rq^n[i] at (s1,m)
#  has l2 norm bounded by B[i], B[i] <= sqrt(q), i=1,...,Z.
#
#  For i in [1,Z]:
#
#  l2(Es[i] * s1 + Em[i] * m + v[i]) <= B[i]
#
#  polymat Es[i] on Rq^(n[i] x m1)
#  polymat Em[i] on Rq^(n[i] x l)
#  polyvec v[i] on Rq^(n[i])

n = [k*(dimm+dimn)]                # set lengths of vectors bounded in l2 norm
B = [BOUND]                  # set l2 norm bounds

# 5. Approximate infinity norm:
#  The evaluation of a linear function Rq^(m1+l) -> Rq^nprime at (s1,m)
#  whose linf norm is bounded by Bprime has linf norm bounded by psi * Bprime,
#  psi >= 1.
#
#  Given
#   linf(Ds * s1 + Dm * m + u) <= Bprime,
#  prove that
#   linf(Ds * s1 + Dm * m + u) <= psi * Bprime
#
#  polymat Ds in Rq^(nprime x m1)
#  polymat Dm in Rq^(nprime x l)
#  polyvec u in Rq^(nprime)

nprime = k*dimn             # set length of vector bounded in linf norm
Bprime = (P+dimm*P*dprime*S+E)/p             # set linf norm bound
