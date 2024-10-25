import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer modules

import time
import secrets          # for internal coins
import hashlib          # for SHAKE128

from cbdc_p1_params import mod, deg, wl2
from _cbdc_params_cffi import lib

ATTRLEN = 32    # byte-length of an attribute
NATTRS = 5      # number of attributes

# 64-deg ring
RING = polyring_t(deg, mod)

# public randomness
shake128 = hashlib.shake_128(bytes.fromhex("00"))
BLINDSIGPP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("01"))
P1PP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("02"))
P2PP = shake128.digest(32)
shake128 = hashlib.shake_128(bytes.fromhex("03"))
POPENPP = shake128.digest(32)

# expand blindsig public parameters from seeds
# BND = int((mod-1)/2)
AR = polymat_t.urandom_static(RING, 8, 16, mod, BLINDSIGPP, 1)
AM = polymat_t.urandom_static(RING, 8, NATTRS, mod, BLINDSIGPP, 2)
ATAU = polymat_t.urandom_static(RING, 8, 8, mod, BLINDSIGPP, 3)
B1 = polymat_t.identity(RING, 8)
A = polymat_t.urandom_static(RING, 2, 11, mod, POPENPP, 1) # commitment key


class InvalidMaskedMsg(Exception):
    pass


class InvalidBlindSig(Exception):
    pass


class InvalidSignature(Exception):
    pass


def redc16(x):
    assert -7 - 16 <= x and x <= 8 + 16
    if x > 8:
        return x - 16
    elif x < -7:
        return x + 16
    else:
        return x

def encattr(attr):
    # Encode 256 attribute bits as dim 64 coefficient vector:
    # Encode each attribute in base 2^4=16 i.e, each attribute is represented
    # by 256/4 = 64 coefficients in [-7,8].
    coeffvec = [0 for _ in range(64)]
    for j in range(ATTRLEN):
        hi = redc16(attr[j] // 16)  # 4 high bits
        lo = redc16(attr[j] % 16)   # 4 low bits
        coeffvec[j * 2] = hi
        coeffvec[j * 2 + 1] = lo
    return coeffvec

def encattrs(attrs):
    # Encode 5*256 attribute bits as 5 dim 64 coefficient vectors:
    # Encode each attribute in base 2^4=16 i.e, each attribute is represented
    # by 256/4 = 64 coefficients in [-7,8].
    # So the resulting vector has 5*256/4 = 5*64 = 320 coefficients.
    coeffvec = [0 for _ in range(NATTRS * 64)]
    for i in range(NATTRS):
        attr = attrs[i]
        for j in range(ATTRLEN):
            hi = redc16(attr[j] // 16)  # 4 high bits
            lo = redc16(attr[j] % 16)   # 4 low bits
            coeffvec[i * ATTRLEN * 2 + j * 2] = hi
            coeffvec[i * ATTRLEN * 2 + j * 2 + 1] = lo
    return coeffvec


class User:
    def __init__(self, pk: falcon_pkenc):
        self.B2 = falcon_decode_pk(pk, RING)
        self.p1_prover = lin_prover_state_t(P1PP, lib.get_params("p1_param"))
        self.p2_prover = lin_prover_state_t(P2PP, lib.get_params("p2_param"))

    # instead of taking the randomness as an input, this tosses coins internally and
    # returns the result.
    def commit(self, a):
        coeffvec = encattr(a)
        attr = polyvec_t(RING, 1, coeffvec)

        seed = secrets.token_bytes(32)  # internal coins
        r = polyvec_t.urandom_bnd_static(RING, 10, -2, 2, seed, 0)

        s = polyvec_t(RING, 1 + 10, [attr, r])
        nym = A*s

        return nym, r
    
    # assume m is actually the 256 bit hash of a message
    def sokprove(self, nym, m, a, r):
        coeffvec = encattr(a)
        attr = polyvec_t(RING, 1, coeffvec)

        s = polyvec_t(RING, 1 + 10, [attr, r])

        popen_prover = lin_prover_state_t(m, lib.get_params("popen_param"))
        popen_prover.set_statement(A, -nym)
        popen_prover.set_witness(s)
        proof = popen_prover.prove()
        return proof

    def maskattrs(self, attrs: list, idx_attrs_rev: list):
        coeffvec = encattrs(attrs)

        m = polyvec_t(RING, NATTRS, coeffvec)
        seed = secrets.token_bytes(32)  # internal coins
        logsigma = 1                    # sigma = 1.55*2^logsigma
        l2sqr_bnd = wl2[0] * wl2[0]     # l2(r1,r2)^2
        r = polyvec_t.grandom_static(RING, 16, logsigma, seed, 0, 0, l2sqr_bnd)

        t = AR*r + AM*m # XXX move that line down and it fails due to AM_priv XXX

        m_priv = m.zero_out_pols(idx_attrs_rev)
        AM_priv = AM.zero_out_cols(idx_attrs_rev)
        m_pub = m.get_pol_list(idx_attrs_rev)
        AM_pub = AM.get_col_list(idx_attrs_rev)    

        A = polymat_t(RING, 8, 16 + NATTRS, [AR, AM_priv])
        u = -(t - AM_pub * m_pub)
        w = polyvec_t(RING, 16 + NATTRS, [r, m_priv])

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

    def presgen(self, cred: bytes, idx_attrs_rev: list, nyms, nyms_r, idx_attrs_com):
        # decode blindsig tau,s1,s2
        tau, s1, s2 = bytes(64), polyvec_t(RING, 8), polyvec_t(RING, 8)
        try:
            coder = coder_t()
            coder.dec_begin(cred)
            coder.dec_bytes(tau)
            coder.dec_grandom(165, s1)
            coder.dec_grandom(165, s2)
            coder.dec_end()
        except DecodingError:
            raise InvalidMaskedMsg

        tau_ = polyvec_t(RING, 8, tau)

        # prove [Ar,Atau,-B1,-B2]*(r,tau,s1,s2) + Am*m = 0
        # as PoK(w): Aw + u = 0

        m_priv = self.m.zero_out_pols(idx_attrs_rev)
        # m_priv=self.m.copy()
        AM_priv = AM.zero_out_cols(idx_attrs_rev)
        if len(idx_attrs_rev) == 0:
            u = polyvec_t(RING, 8 + 3*2)
        else:
            m_pub = self.m.get_pol_list(idx_attrs_rev)
            AM_pub = AM.get_col_list(idx_attrs_rev)
            pad = polyvec_t(RING, 3*2)
            u = polyvec_t(RING, 8 + 3*2, [AM_pub * m_pub, pad])

        #### create A ####
        A = polymat_t(RING, 8+3*2, 40+5*1+3*10)
        for i in range(0,8):
            for j in range(0,16):
                elem = AR.get_elem(i, j)
                A.set_elem(elem, i, j)
        for i in range(0,8):
            for j in range(8):
                elem = ATAU.get_elem(i, j)
                A.set_elem(elem, i, 16 + j)
        for i in range(0,8):
            for j in range(8):
                elem = B1.get_elem(i, j)
                A.set_elem(-elem, i, 24 + j)
        for i in range(0,8):
            for j in range(8):
                elem = self.B2.get_elem(i, j)
                A.set_elem(-elem, i, 32 + j)
        for i in range(0,8):
            for j in range(5):
                elem = AM_priv.get_elem(i, j)
                A.set_elem(elem, i, 40 + j)
        #### create r ####
        assert len(idx_attrs_com) <= 3
        r = polyvec_t(RING, 3*10)
        for i in range(len(idx_attrs_com)):
            for j in range(10):
                elem = nyms_r[idx_attrs_com[i]].get_elem(j)
                r.set_elem(elem, i * 10 + j)
        ####

        #A = polymat_t(RING, 8, 45, [AR, ATAU, -B1, -self.B2, AM_priv])
        w = polyvec_t(RING, 40 + 5*1 + 3*10, [self.r, tau_, s1, s2, m_priv, r])

        #XXXres=A*w+u
        #res.print()

        self.p2_prover.set_statement(A, u)
        self.p2_prover.set_witness(w)
        proof = self.p2_prover.prove()

        return proof


class Issuer:
    def __init__(self, pk: falcon_pkenc, sk: falcon_skenc):
        self.sk = sk
        self.p1_verifier = lin_verifier_state_t(
            P1PP, lib.get_params("p1_param"))

    def issue(self, attrs_rev: bytes, idx_attrs_rev: list, masked_attrs: bytes):
        coeffvec = encattrs(attrs_rev)

        m = polyvec_t(RING, NATTRS, coeffvec)
        AM_priv = AM.zero_out_cols(idx_attrs_rev)
        m_pub = m.get_pol_list(idx_attrs_rev)
        AM_pub = AM.get_col_list(idx_attrs_rev)
        
        # decode masked attributes
        t = polyvec_t(RING, 8)
        try:
            coder = coder_t()
            coder.dec_begin(masked_attrs)
            coder.dec_urandom(mod, t)
            tlen = coder.dec_end()
        except DecodingError:
            raise InvalidMaskedMsg

        # verify proof for PoK(w): Aw + u = 0
        A = polymat_t(RING, 8, 16 + NATTRS, [AR, AM_priv])
        u = polyvec_t(RING, 8, [-(t - AM_pub * m_pub)])
        proof = masked_attrs[tlen:]

        self.p1_verifier.set_statement(A, u)
        try:
            self.p1_verifier.verify(proof)
        except VerificationError:
            raise InvalidMaskedMsg("Attributes invalid.")

        # internal coins
        tau_ = secrets.token_bytes(64)
        tau = polyvec_t(RING, 8, tau_)

        # preimage sampling
        s1, s2 = falcon_preimage_sample(self.sk, ATAU * tau + t, RING)

        # encode blindsig tau,s1,s2
        coder = coder_t()
        coder.enc_begin(2000)
        coder.enc_bytes(tau_)
        coder.enc_grandom(165, s1)
        coder.enc_grandom(165, s2)
        cred = coder.enc_end()

        return cred


class Verifier:
    def __init__(self, pk: falcon_pkenc):
        self.B2 = falcon_decode_pk(pk, RING)
        self.p2_verifier = lin_verifier_state_t(
            P2PP, lib.get_params("p2_param"))

    # the private part should already be zeroed out in attrs_rev
    def presver(self, attrs_rev: bytes, idx_attrs_rev: list, pres: bytes, nyms, idx_attrs_com):
        coeffvec = encattrs(attrs_rev)

        m = polyvec_t(RING, NATTRS, coeffvec)
        AM_priv = AM.zero_out_cols(idx_attrs_rev)
        if len(idx_attrs_rev) == 0:
            u = polyvec_t(RING, 8 + 3*2)
        else:
            m_pub = m.get_pol_list(idx_attrs_rev)
            AM_pub = AM.get_col_list(idx_attrs_rev)
            pad = polyvec_t(RING, 3*2)
            u = polyvec_t(RING, 8 + 3*2, [AM_pub * m_pub, pad])

        # verify proof for PoK(w): Aw + u = 0
                #### create A ####
        A = polymat_t(RING, 8+3*2, 40+5*1+3*10)
        for i in range(0,8):
            for j in range(0,16):
                elem = AR.get_elem(i, j)
                A.set_elem(elem, i, j)
        for i in range(0,8):
            for j in range(8):
                elem = ATAU.get_elem(i, j)
                A.set_elem(elem, i, 16 + j)
        for i in range(0,8):
            for j in range(8):
                elem = B1.get_elem(i, j)
                A.set_elem(-elem, i, 24 + j)
        for i in range(0,8):
            for j in range(8):
                elem = self.B2.get_elem(i, j)
                A.set_elem(-elem, i, 32 + j)
        for i in range(0,8):
            for j in range(5):
                elem = AM_priv.get_elem(i, j)
                A.set_elem(elem, i, 40 + j)
        #A = polymat_t(RING, 8, 40 + NATTRS, [AR, ATAU, -B1, -self.B2, AM_priv])
        # u = polyvec_t(RING, 8, [AM_pub * m_pub])
        proof = pres

        self.p2_verifier.set_statement(A, u)
        try:
            self.p2_verifier.verify(proof)
        except VerificationError:
            raise InvalidSignature("Presentation invalid.")
        
    # assume m is actually the 256 bit hash of a message
    def sokverify(self, nym, m, proof):
        popen_verifier = lin_verifier_state_t(m, lib.get_params("popen_param"))
        popen_verifier.set_statement(A, -nym)
        try:
            popen_verifier.verify(proof)
        except VerificationError:
            raise InvalidSignature("Proof invalid.")
        return proof


def main():
    print("lazer CBDC demo")
    print("---------------\n")

    # create a keypair for testing
    isk, ipk, _ = falcon_keygen()

    # create a list of 5 random 32-byte attributes for testing
    attrs = [secrets.token_bytes(ATTRLEN) for _ in range(NATTRS)]
    # create a random (hash of a) 256 bit message for testing
    m = secrets.token_bytes(32)

    # attributes revealed to the issuer
    idx_attrs_rev = [1]             # set of indices of attributes to reveal
    assert set(idx_attrs_rev).issubset(set(range(NATTRS))) == True
    attrs_rev = [b"\0" * 32 for i in range(NATTRS)]
    for  i in idx_attrs_rev:
        attrs_rev[i] = attrs[i]

    print("Initialize user with issuer public key ... ", end='')
    user = User(ipk)
    print("[OK]\n")

    print("Initialize issuer with issuer public and private key ... ", end='')
    issuer = Issuer(ipk, isk)
    print("[OK]\n")

    print("Initialize verifier with issuer public key ... ", end='')
    verifier = Verifier(ipk)
    print("[OK]\n")

    print("User outputs masked attributes ... ", end='')
    masked_attrs = user.maskattrs(attrs, idx_attrs_rev)
    print("[OK]\n")

    print(f"masked attributes (t,P1): {len(masked_attrs)} bytes\n")

    print("Issuer checks the masked attributes and outputs credentials ... ", end='')
    try:
        cred = issuer.issue(attrs_rev, idx_attrs_rev, masked_attrs)
    except InvalidMaskedMsg:
        print("masked attributes are invalid.")
        sys.exit(1)
    try:
        cred = issuer.issue(attrs_rev, list(range(5)), masked_attrs)
    except InvalidMaskedMsg:
        pass
    else:
        print("issuing for invalid indices should fail, but succeeded.")
        sys.exit(1)
    try:
        attrs_rev_inv = attrs_rev.copy()
        attrs_rev_inv[idx_attrs_rev[0]] = b"\0" * 32
        cred = issuer.issue(attrs_rev_inv, idx_attrs_rev, masked_attrs)
    except InvalidMaskedMsg:
        pass
    else:
        print("issuing for invalid attributes should fail, but succeeded.")
        sys.exit(1)
    print("[OK]\n")
    print(f"credentials (tau,s1,s2): {len(cred)} bytes\n")

    # attributes revealed to the verifier
    idx_attrs_rev = [0]             # set of indices of attributes to reveal
    idx_attrs_com = [1,2]           # set of indices of attributes to commit to
    assert set(idx_attrs_rev).issubset(set(range(NATTRS))) == True
    assert set(idx_attrs_com).issubset(set(range(NATTRS))) == True
    assert set(idx_attrs_rev).isdisjoint(set(idx_attrs_com)) == True
    attrs_rev = [b"\0" * 32 for i in range(NATTRS)]
    for  i in idx_attrs_rev:
        attrs_rev[i] = attrs[i]

    print("User creates nyms ... ", end='')
    nyms = [0 for _ in range(NATTRS)]
    nyms_r = [0 for _ in range(NATTRS)]
    for i in idx_attrs_com:
        nyms[i], nyms_r[i] = user.commit(attrs[i])
        # test standalone verification
        proof = user.sokprove(nyms[i], m, attrs[i], nyms_r[i])
        verifier.sokverify(nyms[i], m, proof)
    print("[OK]\n")

    print("User generates a presentation ... ", end='')
    try:
        #pres = user.presgen(cred, idx_attrs_rev, [], [], [])
        pres = user.presgen(cred, idx_attrs_rev, nyms, nyms_r, idx_attrs_com)
    except InvalidBlindSig:
        print("decoding failed.")
        sys.exit(1)
    print("[OK]\n")

    print(f"presentation (P2): {len(pres)} bytes\n")

    print("Verfifier verifies the presentation ... ", end='')
    try:
        #verifier.presver(attrs_rev, idx_attrs_rev, pres, [], [])
        verifier.presver(attrs_rev, idx_attrs_rev, pres, nyms, idx_attrs_com)
    except InvalidSignature:
        print("signature invalid.")
        sys.exit(1)
    try:
        verifier.presver(attrs_rev, list(range(5)), pres, [], []) #XXX
    except InvalidSignature:
        pass
    else:
        print("iverifying for invalid indices should fail, but succeeded.")
        sys.exit(1)
    try:
        attrs_rev_inv = attrs_rev.copy()
        attrs_rev_inv[idx_attrs_rev[0]] = b"\0" * 32
        verifier.presver(attrs_rev_inv, idx_attrs_rev, pres, [], []) #XXX
    except InvalidSignature:
        pass
    else:
        print("verifying for invalid attributes should fail, but succeeded.")
        sys.exit(1)

    print("[OK]\n")

if __name__ == "__main__":
    main()
