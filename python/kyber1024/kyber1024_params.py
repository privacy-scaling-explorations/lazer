from math import sqrt
vname = "param"               # variable name

deg   = 256                   # ring Rp degree d
mod   = 3329                  # ring Rp modulus p
m,n   = 4,8
dim   = (m,n)                 # dimensions of A in Rp^(m,n)

wpart = [ list(range(n))   ]  # partition of w : [w]
wl2   = [ 1.2*sqrt(deg*n)  ]  # l2-norm bounds : l2(w) <= 1.2 * sqrt(deg*n) 
wbin  = [ 0                ]  # binary coeffs  : n/a
#rej  = [ 1                ]  # rej. sampling  : on w

# Optional: some linf-norm bound on w.
# Tighter bounds result in smaller proofs.
# If not specified, the default is the naive bound max(1,floor(max(wl2))).
wlinf = 2 # optional linf: some linf-norm bound on w.