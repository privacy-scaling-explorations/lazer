import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer modules
import hashlib          # for SHAKE128
import secrets          # for internal coins
import sys
import time

from _anon_cred_params_cffi import lib
from anon_cred_p1_params import mod, deg, wl2

# 64-deg ring
RING = polyring_t(deg, mod)

# public randomness
shake128 = hashlib.shake_128(bytes.fromhex("00"))
BLINDSIGPP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("01"))
P1PP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("02"))
P2PP = shake128.digest(32)

# expand blindsig public parameters from seeds
#BND = int((mod-1)/2)
AR=polymat_t.urandom_static(RING,8,16,mod,BLINDSIGPP,1)
AM=polymat_t.urandom_static(RING,8,8,mod,BLINDSIGPP,2)
ATAU=polymat_t.urandom_static(RING,8,8,mod,BLINDSIGPP,3)
B1 = polymat_t.identity(RING,8)


class InvalidMaskedMsg(Exception):
    pass


class InvalidBlindSig(Exception):
    pass


class InvalidSignature(Exception):
    pass


class user_t:
    def __init__(self, pk: falcon_pkenc):
        self.B2 = falcon_decode_pk(pk,RING) #
        self.p1_prover = lin_prover_state_t(P1PP, lib.get_params("p1_param"))
        self.p2_prover = lin_prover_state_t(P2PP, lib.get_params("p2_param"))

    def maskmsg(self, msg: bytes):
        m = polyvec_t(RING,8,msg)
        seed = secrets.token_bytes(32)  # internal coins
        logsigma = 1                    # sigma = 1.55*2^logsigma
        l2sqr_bnd = wl2[0] * wl2[0]     # l2(r1,r2)^2
        r=polyvec_t.grandom_static(RING,16,logsigma,seed,0,0,l2sqr_bnd)
        t = AR*r + AM*m

        # prove [Ar,Am]*(r,m) + (-t) = 0
        # as PoK(w): Aw + u = 0
        A = polymat_t(RING, 8, 24, [AR, AM])
        u = polyvec_t(RING, 8, [-t])
        w = polyvec_t(RING, 24, [r, m])
   
        self.p1_prover.set_statement(A, u)
        self.p1_prover.set_witness(w)
        
        proof = self.p1_prover.prove()

        # encode t
        coder = coder_t()
        coder.enc_begin(22000)
        coder.enc_urandom(mod, t)
        tenc = coder.enc_end()

        # save randomness and message
        self.r = r
        self.m = m
        # return masked message t,p1enc
        return tenc + proof

    def sign(self, blindsig: bytes, pub_mvec: list):
        # decode blindsig tau,s1,s2
        
        tau, s1, s2 = bytes(64), polyvec_t(RING,8), polyvec_t(RING,8)
        try:
            coder = coder_t()
            coder.dec_begin(blindsig)
            coder.dec_bytes(tau)
            coder.dec_grandom(165, s1)
            coder.dec_grandom(165, s2)
            coder.dec_end()
        except DecodingError:
            raise InvalidMaskedMsg

        tau_ = polyvec_t(RING, 8, tau)

        # prove [Ar,Atau,-B1,-B2]*(r,tau,s1,s2) + Am*m = 0
        # as PoK(w): Aw + u = 0
        
        m_priv=self.m.zero_out_pols(pub_mvec)
        #m_priv=self.m.copy()
        AM_priv=AM.zero_out_cols(pub_mvec)
        if len(pub_mvec)==0:
            u=polyvec_t(RING,8)
        else:
            m_pub=self.m.get_pol_list(pub_mvec)
            AM_pub=AM.get_col_list(pub_mvec)
            u=polyvec_t(RING, 8, [AM_pub * m_pub])
      

        A = polymat_t(RING, 8, 48, [AR, ATAU, -B1, -self.B2, AM_priv])        
        w = polyvec_t(RING, 48, [self.r, tau_, s1, s2, m_priv])

        self.p2_prover.set_statement(A, u)
        self.p2_prover.set_witness(w)
        proof = self.p2_prover.prove()

        return proof


class signer_t:
    def __init__(self, pk: falcon_pkenc, sk: falcon_skenc):
        self.sk = sk
        self.p1_verifier = lin_verifier_state_t(P1PP, lib.get_params("p1_param"))

    def sign(self, masked_msg: bytes):
        # decode masked message
        t = polyvec_t(RING,8)
        try:
            coder = coder_t()
            coder.dec_begin(masked_msg)
            coder.dec_urandom(mod, t)
            tlen = coder.dec_end()
        except DecodingError:
            raise InvalidMaskedMsg
        
        # verify proof for PoK(w): Aw + u = 0
        A = polymat_t(RING, 8, 24, [AR, AM])
        u = polyvec_t(RING, 8, [-t])
        proof = masked_msg[tlen:]

        self.p1_verifier.set_statement(A, u)
        try:
            self.p1_verifier.verify(proof)
        except VerificationError:
            raise InvalidMaskedMsg("Masked message invalid.")

        # internal coins
        tau_ = secrets.token_bytes(64)
        tau = polyvec_t(RING,8,tau_)

        # preimage sampling
        s1, s2 = falcon_preimage_sample(self.sk, ATAU * tau + t,RING)
        #s1, s2 = falcon_preimage_sample(self.sk, ATAU * tau + t)

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
        self.B2 = falcon_decode_pk(pk,RING)
        self.p2_verifier = lin_verifier_state_t(P2PP, lib.get_params("p2_param"))

    #the private part should already be zeroed out in pub_msg 
    def verify(self, pub_msg: bytes, pub_mvec:list, sig: bytes):
        
        m = polyvec_t(RING, 8, pub_msg)
        AM_priv=AM.zero_out_cols(pub_mvec)
        if len(pub_mvec)==0:
            u=polyvec_t(RING,8)
        else:
            m_pub=m.get_pol_list(pub_mvec)
            AM_pub=AM.get_col_list(pub_mvec)
            u=polyvec_t(RING, 8, [AM_pub * m_pub])

        # verify proof for PoK(w): Aw + u = 0
        A = polymat_t(RING, 8, 48, [AR, ATAU, -B1, -self.B2, AM_priv])
        #u = polyvec_t(RING, 8, [AM_pub * m_pub])
        proof = sig
        
        self.p2_verifier.set_statement(A, u)
        try:
            self.p2_verifier.verify(proof)
        except VerificationError:
            raise InvalidSignature("Signature invalid.")


# demo
def main():
    print("lazer anonymous credentials demo")
    print("--------------------------\n")

    msg = bytes.fromhex(
        "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef")
    sk, pk, _ = falcon_keygen()

    pub_mvec=[0,4,5]
    priv_mvec=list(set(range(8))-set(pub_mvec))

    print("Initialize user with public key ... ", end='')
    user = user_t(pk)
    print("[OK]\n")

    print("Initialize signer with public and private key ... ", end='')
    signer = signer_t(pk, sk)
    print("[OK]\n")

    print("Initialize verifier with public key ... ", end='')
    verifier = verifier_t(pk)
    print("[OK]\n")

    issue_start=time.perf_counter()
    print("User outputs masked credentials (including a proof of well-formedness) ... ", end='')
    masked_msg = user.maskmsg(msg)
    print("[OK]\n")


    print_stopwatch_lnp_prover_prove(0)
    print(f"masked credentials (t,P1): {len(masked_msg)} bytes\n")

    print("Signer checks the proof and if it verifies outputs blinded credentials ... ", end='')
    try:
        blindsig = signer.sign(masked_msg)
    except InvalidMaskedMsg:
        print("masked credentials are invalid.")
        sys.exit(1)
    print("[OK]\n")
    issue_end=time.perf_counter()
    print(f"blind credentials (tau,s1,s2): {len(blindsig)} bytes\n")

    show_start=time.perf_counter()
    print("User outputs a signature on the hidden credentials ... ", end='')
    try:
        sig = user.sign(blindsig,pub_mvec)
    except InvalidBlindSig:
        print("decoding failed.")
        sys.exit(1)
    print("[OK]\n")

    print(f"signature (P2): {len(sig)} bytes\n")

    print("Verfifier verifies the signature of the blinded credentials ... ", end='')
    try:
        msg_pub=zero_out_bytes(msg,priv_mvec,RING.deg//8)
        verifier.verify(msg_pub,pub_mvec,sig)
    except InvalidSignature:
        print("signature invalid.")
        sys.exit(1)
    print("[OK]\n")
    show_end=time.perf_counter()
    print_stopwatch_lnp_prover_prove(0)
    print("Issue time: ",-issue_start+issue_end)
    print("Show time: ",-show_start+show_end)

if __name__ == "__main__":
    main()
