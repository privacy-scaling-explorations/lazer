# Create a header file with LNP proof system parameters for
# proving knowledge of a witness w in Rp^n (Rp = Zp[X]/(X^d + 1))
# such that
#
#   1. w satisfies a linear relation over Rp: Aw + t = 0
#   2. each element in a partition of w either ..
#      2.1 has binary coefficients only
#      2.2 satisfies an l2-norm bound
from math import sqrt
vname = "param"           # variable name

deg   = 256               # ring Rp degree d
mod   = 2**214-255        # ring Rp modulus p
m,n   = 32,64
dim   = (m,n)             # dimensions of A in Rp^(m,n)

wpart = [ list(range(n)) ]  # partition of w : [w]
wl2   = [ sqrt(deg*n)    ]  # l2-norm bounds : l2(w)^2 <= 16384
wbin  = [ 0              ]  # binary coeffs  : n/a
wrej  = [ 1              ]  # rej. sampling  : on m w

# Optional: some linf-norm bound on w.
# Tighter bounds result in smaller proofs.
# If not specified, the default is the naive bound max(1,floor(max(wl2))).
wlinf = 1 # optional linf: some linf-norm bound on w.
