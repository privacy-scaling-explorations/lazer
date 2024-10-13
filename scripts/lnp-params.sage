# Set LNP proof system parameters.
#
# Notation:
#   R = Z[X]\(X^d+1), Rq = Zq[X]\(X^d+1) polynomial rings,
#    d is either 64 or 128, q is either q1 or q1*q2,
#    q1, q2 odd primes congruent to 5 mod 8, q1 > q2.
#
#   o in Aut(R), automorphism defined by o(X) = X^-1.
#
#   <x> = poly||o(x), x in R, works elementwise on vectors and
#    matrices of polynomials.

name = "params"         # parameter set name

# Prove knowledge of polyvec (s1,m) in Rq^(m1+l), l2(s1) <= alpha,
# and prove statements (1-5).

d = 64                  # set ring degree (either 64 or 128)
log2q = 32              # set ring modulus q's bit-size
# log2q1 = 0            # if q = q1*q2, also set q1's bit-size
m1 = 16                 # set length of vector s1
alpha = sqrt(1024)      # set l2-norm bound on vector s1
l = 0                   # set length of vector m

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

n = [32]                # set lengths of vectors bounded in l2 norm
B = [sqrt(2048)]        # set l2 norm bounds

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

# nprime =              # set length of vector bounded in linf norm
# Bprime =              # set linf norm bound
