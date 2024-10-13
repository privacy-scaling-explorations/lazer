# Create a header file with proof system parameters for
# proving knowledge of a witness w in Rp^n (Rp = Zp[X]/(X^d + 1))
# such that
#
#   1. w satisfies a linear relation over Rp: Aw + u = 0
#   2. each element in a partition of w either ..
#      2.1 has binary coefficients only
#      2.2 satisfies an l2-norm bound

vname = ""  # variable name

deg   =     # Rp degree
mod   =     # Rp modulus
dim   =     # dimensions of A in Rp^(m,n)

wpart = []  # partition of w
wl2   = []  # l2-norm bounds on parts of w
wbin  = []  # binary coeffs condition on parts of w

# wlinf = # optional: linf-norm bound on w
