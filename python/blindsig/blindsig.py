import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
import hashlib          # for SHAKE128
import secrets          # for internal coins

from _blindsig_params_cffi import lib
from blindsig_p1_params import mod, deg, wl2


# public randomness
shake128 = hashlib.shake_128(bytes.fromhex("00"))
BLINDSIGPP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("01"))
P1PP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("02"))
P2PP = shake128.digest(32)

# expand blindsig public parameters from seeds
RING = polyring_t(deg, mod) # falcon ring
BND = int((mod-1)/2)
AR1, AR2, AM, ATAU = poly_t(RING), poly_t(RING), poly_t(RING), poly_t(RING)
AR1.urandom_bnd(-BND, BND, BLINDSIGPP, 0)
AR2.urandom_bnd(-BND, BND, BLINDSIGPP, 1)
AM.urandom_bnd(-BND, BND, BLINDSIGPP, 2)
ATAU.urandom_bnd(-BND, BND, BLINDSIGPP, 3)
B1 = poly_t(RING, {0: 1})


class InvalidMaskedMsg(Exception):
    pass


class InvalidBlindSig(Exception):
    pass


class InvalidSignature(Exception):
    pass


class user_t:
    def __init__(self, pk: falcon_pkenc):
        #self.B2 = poly_t(RING, falcon_decode_pk(pk))
        self.B2 = falcon_decode_pk(pk)
        self.p1_prover = lin_prover_state_t(P1PP, lib.get_params("p1_param"))
        self.p2_prover = lin_prover_state_t(P2PP, lib.get_params("p2_param"))

    def maskmsg(self, msg: bytes):
        if len(msg) != 64:
            raise ValueError("msg must be 32 bytes.")

        # sample (r1,r1)
        r1, r2 = poly_t(RING), poly_t(RING)
        m = poly_t(RING, msg)
        seed = secrets.token_bytes(32)  # internal coins
        logsigma = 1                    # sigma = 1.55*2^logsigma
        l2sqr_bnd = wl2[0] * wl2[0]     # l2(r1,r2)^2
        ctr = 0
        while (True):
            r1.grandom(logsigma, seed, ctr)
            ctr += 1
            r2.grandom(logsigma, seed, ctr)
            ctr += 1

            l2sqr = r1.l2sq() + r2.l2sq()
            if (l2sqr <= l2sqr_bnd):
                break

        t = AR1*r1 + AR2*r2 + AM*m

        # prove [Ar1,Ar2,Am]*(r1,r2,m) + (-t) = 0
        # as PoK(w): Aw + u = 0
        A = polymat_t(RING, 1, 3, [AR1, AR2, AM])
        u = polyvec_t(RING, 1, [-t])
        w = polyvec_t(RING, 3, [r1, r2, m])

        self.p1_prover.set_statement(A, u)
        self.p1_prover.set_witness(w)
        proof = self.p1_prover.prove()

        # encode t
        coder = coder_t()
        coder.enc_begin(22000)
        coder.enc_urandom(mod, t)
        tenc = coder.enc_end()

        # save randomness and message
        self.r1 = r1
        self.r2 = r2
        self.m = m

        # return masked message t,p1enc
        return tenc + proof

    def sign(self, blindsig: bytes):
        # decode blindsig tau,s1,s2

        tau, s1, s2 = bytes(64), poly_t(RING), poly_t(RING)
    
        try:
            coder = coder_t()
            coder.dec_begin(blindsig)
            coder.dec_bytes(tau)
            coder.dec_grandom(165, s1)
            coder.dec_grandom(165, s2)
            coder.dec_end()
        except DecodingError:
            raise InvalidMaskedMsg

        tau_ = poly_t(RING, tau)

        # prove [Ar1,Ar2,Atau,-B1,-B2]*(r1,r2,tau,s1,s2) + Am*m = 0
        # as PoK(w): Aw + u = 0
        A = polymat_t(RING, 1, 5, [AR1, AR2, ATAU, -B1, -self.B2])
        u = polyvec_t(RING, 1, [AM * self.m])
        w = polyvec_t(RING, 5, [self.r1, self.r2, tau_, s1, s2])

        self.p2_prover.set_statement(A, u)
        self.p2_prover.set_witness(w)
        proof = self.p2_prover.prove()

        return proof


class signer_t:
    def __init__(self, sk: falcon_skenc):
        self.sk = sk
        self.p1_verifier = lin_verifier_state_t(P1PP, lib.get_params("p1_param"))

    def sign(self, masked_msg: bytes):
        # decode masked message
        t = poly_t(RING)
        try:
            coder = coder_t()
            coder.dec_begin(masked_msg)
            coder.dec_urandom(mod, t)
            tlen = coder.dec_end()
        except DecodingError:
            raise InvalidMaskedMsg
        
        # verify proof for PoK(w): Aw + u = 0
        A = polymat_t(RING, 1, 3, [AR1, AR2, AM])
        u = polyvec_t(RING, 1, [-t])
        proof = masked_msg[tlen:]

        self.p1_verifier.set_statement(A, u)
        try:
            self.p1_verifier.verify(proof)
        except VerificationError:
            raise InvalidMaskedMsg("Masked message invalid.")

        # internal coins
        tau_ = secrets.token_bytes(64)
        tau = poly_t(RING, tau_)

        # preimage sampling
        s1, s2 = falcon_preimage_sample(self.sk, ATAU * tau + t)

        # encode blindsig tau,s1,s2
        coder = coder_t()
        coder.enc_begin(2000)
        coder.enc_bytes(tau_)
        coder.enc_grandom(165, s1)
        coder.enc_grandom(165, s2)
        blindsig = coder.enc_end()

        return blindsig



class verifier_t:
    def __init__(self, pk: falcon_pkenc):
        #self.B2 = poly_t(RING, falcon_decode_pk(pk), )
        self.B2 = falcon_decode_pk(pk)        
        self.p2_verifier = lin_verifier_state_t(P2PP, lib.get_params("p2_param"))

    def verify(self, msg: bytes, sig: bytes):
        m = poly_t(RING, msg)

        # verify proof for PoK(w): Aw + u = 0
        A = polymat_t(RING, 1, 5, [AR1, AR2, ATAU, -B1, -self.B2])
        u = polyvec_t(RING, 1, [AM * m])

        self.p2_verifier.set_statement(A, u)
        try:
            self.p2_verifier.verify(sig)
        except VerificationError:
            raise InvalidSignature("Signature invalid.")


# demo
def main():
    print("lazer blind-signature demo")
    print("--------------------------\n")

    msg = bytes.fromhex(
        "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef")
    sk, pk, _ = falcon_keygen()

    print("Initialize user with public key ... ", end='')
    user = user_t(pk)
    print("[OK]\n")

    print("Initialize signer with public and private key ... ", end='')
    signer = signer_t(sk)
    print("[OK]\n")

    print("Initialize verifier with public key ... ", end='')
    verifier = verifier_t(pk)
    print("[OK]\n")

    print("User outputs masked message (including a proof of its well-formedness) ... ", end='')
    masked_msg = user.maskmsg(msg)
    print("[OK]\n")

    print(f"masked message (t,P1): {len(masked_msg)} bytes\n")

    print("Signer checks the proof and if it verifies outputs a blind signature ... ", end='')
    try:
        blindsig = signer.sign(masked_msg)
    except InvalidMaskedMsg:
        print("masked message is invalid.")
        sys.exit(1)
    print("[OK]\n")

    print(f"blind signature (tau,s1,s2): {len(blindsig)} bytes\n")

    print("User outputs a signature on the message ... ", end='')
    try:
        sig = user.sign(blindsig)
    except InvalidBlindSig:
        print("decoding blindsig failed.")
        sys.exit(1)
    print("[OK]\n")

    print(f"signature (P2): {len(sig)} bytes\n")

    print("Verfifier verifies the signature on the message ... ", end='')
    try:
        verifier.verify(msg, sig)
    except InvalidSignature:
        print("signature invalid.")
        sys.exit(1)
    print("[OK]\n")

    print_stopwatch_lnp_prover_prove(0)
    print_stopwatch_lnp_verifier_verify(0)

if __name__ == "__main__":
    main()
