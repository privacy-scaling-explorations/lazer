# Set proof system parameters for proving knowledge
# of an MLWE secret (s,e) over Rp:
#
#   t = A*s + e
#   l2(s,e) <= B
#
# A in Rp^(n x m), t,e in Rp^n, s in Rp^m,
# Rp = Zp[X]/(X^d + 1)

name = "kyber_params"       # variable name

dprime = 256                     # Kyber ring degree
p = 3329                    # Kyber modulus

dimn = 2                       # (2,3,4) for Kyber(512,768,1024)
dimm = 2                       # (2,3,4) for Kyber(512,768,1024)

S = 3                       # linf(s) <= S, (3,2,2) for Kyber(512,768,1024)
E = 3                       # linf(e) <= E, (3,2,2) for Kyber(512,768,1024)

# mu and var per bernoulli trial B
mu_ber = 0.5 * 0 + 0.5 * 1
var_ber = 0.5 * (0 - 0.5)^2 + 0.5 * (1 - 0.5)^2
# Bin_S = sum((-1)^i * B) for i = 0,...,2S-1
var = 2*S*var_ber

# not affected by ring iso: dimensions n and m are multiplied by
# k but degree dprime is divided by k
BOUND = 1.1 * sqrt((dimn + dimm) * dprime * sqrt(var)^2) # 1.1*expected norm

d = 64                  # set ring degree (either 64 or 128)
k = dprime/d

PSI = 3771
P = (p-1)/2
# not affected by k since dimension m is multiplied by k
# but degree dprime is divided by k
# set ring modulus q's bit-size
log2q = ceil(log(2 * (1+PSI)*(P + dimm*P*dprime*S + E) + 1, 2))
# log2q1 = 0            # if q = q1*q2, also set q1's bit-size
m1 = k*(dimm+dimn)              # set length of vector s1
alpha = BOUND        # set l2-norm bound on vector s1
l = 0                  # set length of vector m


n = [k*(dimm+dimn)]                # set lengths of vectors bounded in l2 norm
B = [BOUND]                  # set l2 norm bounds


nprime = k*dimn             # set length of vector bounded in linf norm
Bprime = (P+dimm*P*dprime*S+E)/p             # set linf norm bound
