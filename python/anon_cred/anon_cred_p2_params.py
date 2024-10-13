from math import sqrt

vname = "p2_param"                     # variable name

deg   = 64                           # ring Rp degree d
mod   = 12289                           # ring Rp modulus p
dim   = (8,48)                           # dimensions of A in Rp^(m,n)

wpart = [ list(range(0,16)), list(range(16,24)), list(range(24,40)), list(range(40,48)) ]  # partition of w : [r1,r2], [tau], [s1,s2], [cred]
wl2   = [   109,   0, sqrt(34034726), 0 ]  # l2-norm bounds    : l2(r1,r2) <= 109, l2(s1,s2) <= sqrt(34034726)
wbin  = [     0,   1,              0, 1 ]  # binary coeffs     : tau and msg is binary
#wrej  = [     0,   1,              1, 1 ]  # rejection sampling: on tau, cred, and s_i

# Optional: some linf-norm bound on w.
# Tighter bounds result in smaller proofs.
# If not specified, the default is the naive bound max(1,floor(max(wl2))).
wlinf = 5833  # optional linf: some linf-norm bound on w.
