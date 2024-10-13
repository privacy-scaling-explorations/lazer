import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module

def main():
    # setup

    from demo_params import mod, deg, dim
    d, p, m, n = deg, mod, dim[0], dim[1]

    seed = b'\0' * 32
    from _demo_params_cffi import lib
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))

    Rp = polyring_t(d, p)

    A = polymat_t(Rp, m, n)
    A.urandom(p, seed, 0)


    # prover
    #s=polyvec_t.brandom_static(Rp,n,1,seed,0)
    s = polyvec_t(Rp, n)

    s.brandom(1, seed, 0)
    #s.print()

    t = -A*s

    prover.set_statement(A, t)
    prover.set_witness(s)

    print("generate proof ...")
    proof = prover.prove()
    print_stopwatch_lnp_prover_prove(0)



    # verifier

    verifier.set_statement(A, t)

    print("verify proof ... ")
    try:
        verifier.verify(proof)
    except VerificationError:
        print("reject")
    else:
        print("accept")
    print_stopwatch_lnp_verifier_verify(0)

if __name__ == "__main__":
    main()

