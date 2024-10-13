# Create a header file with LNP proof system parameters for
# proving knowledge of a witness w in Rp^n (Rp = Zp[X]/(X^d + 1))
# such that
#
#   1. w satisfies a linear relation over Rp: Aw + t = 0
#   2. each element in a partition of w either ..
#      2.1 has binary coefficients only
#      2.2 satisfies an l2-norm bound

vname = "p2_params"                     # variable name

deg   = 512                             # ring Rp degree d
mod   = 12289                           # ring Rp modulus p
dim   = (1,5)                           # dimensions of A in Rp^(m,n)

wpart = [ [0,1], [2],          [3,4] ]  # partition of w : [r1,r2], [tau], [s1,s2]
wl2   = [   109,   0, sqrt(34034726) ]  # l2-norm bounds    : l2(r1,r2) <= 109, l2(s1,s2) <= sqrt(34034726)
wbin  = [     0,   1,              0 ]  # binary coeffs     : tau is binary
wrej  = [     0,   1,              0 ]  # rejection sampling: on tau only

# Optional: some linf-norm bound on w.
# Tighter bounds result in smaller proofs.
# If not specified, the default is the naive bound max(1,floor(max(wl2))).
wlinf = 5833  # optional linf: some linf-norm bound on w.
