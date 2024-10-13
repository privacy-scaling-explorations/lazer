import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
import hashlib          # for SHAKE128
import secrets          # for RNG

# public randomness
shake128 = hashlib.shake_128(bytes.fromhex("00"))
KYBERPP = shake128.digest(32)   # kyber public randomness
shake128 = hashlib.shake_128(bytes.fromhex("01"))
P1PP = shake128.digest(32)      # proof system public randomness

from kyber1024_params import mod, deg, m, n     # import kyber parameters

from _kyber1024_params_cffi import lib          # import proof system parameters
prover = lin_prover_state_t(P1PP, lib.get_params("param"))
verifier = lin_verifier_state_t(P1PP, lib.get_params("param"))

# define kyber1024 ring (R), expand public parameters (A) and generate kyber key pair (sk,pk)
R = polyring_t(deg, mod)
A1 = polymat_t.urandom_static(R, m, m, mod, KYBERPP, 0)
A2 = polymat_t.identity(R, m)
A = polymat_t(R, m, n, [A1, A2])
sk = polyvec_t.brandom_static(R, n, 2, secrets.token_bytes(32), 0)
pk = A*sk

# prover

prover.set_statement(A, -pk)
prover.set_witness(sk)

print("generate proof ...")
proof = prover.prove()
print_stopwatch_lnp_prover_prove(0)

# verifier

verifier.set_statement(A, -pk)

print(f"len {len(proof)} ")
print("verify proof ... ")
try:
    verifier.verify(proof)
except VerificationError:
    print("reject")
else:
    print("accept")
print_stopwatch_lnp_verifier_verify(0)
