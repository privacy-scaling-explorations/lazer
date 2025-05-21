#!/usr/bin/python

from lazer import *
from lazer import _invmod, _center_list, _l2sq_list
from labrador import *
import hashlib      # for SHAKE128
import secrets      # for internal coins
import sys
import time


VBIG=0 # counts how many signatures are too big

#Total Signatures
sig_num=int(sys.argv[1])
norms=[17017363,17017363,round(1248245003*.75)]

#Falcon parameters
mod=12289
deg=512

# falcon ring
FALCON_RING=polyring_t(deg,mod)
BIGMOD_RING=polyring_t(deg,LAB_RING_48.mod)
PRIMESIZE=str(math.ceil(math.log2(BIGMOD_RING.mod)))

#use the same sk/pk falcon key to save time on key generation
# makes no difference for the ZK benchmark
SAME_KEY=1

ID=int_to_poly(1,BIGMOD_RING)


# public randomness
shake128 = hashlib.shake_128(bytes.fromhex("00"))
TARGPP = shake128.digest(32)

inv_fal_mod=_invmod(mod,BIGMOD_RING.mod)

deg_list=[deg]*(3*sig_num)
num_pols_list=[1]*(3*sig_num)
norm_list=norms*sig_num
num_constraints=sig_num
PS=proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)

keytime_start=time.perf_counter()
if SAME_KEY:
    skenc,pkenc,pkpol=falcon_keygen()
    l_pk=pkpol.lift(BIGMOD_RING)
else:
    sk_list=[]
    pk_list=[]
    for i in range(sig_num):
        skenc,pkenc,pkpol=falcon_keygen()
        l_pk=pkpol.lift(BIGMOD_RING)
        sk_list+=[skenc]
        pk_list+=[l_pk]
keytime_end=time.perf_counter()

j=0
sig_start=time.perf_counter()
while j<sig_num:

    f_t=poly_t.urandom_static(FALCON_RING,FALCON_RING.mod,TARGPP,0)
    l_t=f_t.lift(BIGMOD_RING)

    if not SAME_KEY:
        skenc=sk_list[j]
        l_pk=pk_list[j]

    l_s1, l_s2 = falcon_preimage_sample(skenc, f_t) # s_1+s_2*pkpol=t, return poly_t

    l_s1=l_s1.lift(BIGMOD_RING)
    l_s2=l_s2.lift(BIGMOD_RING)

    v=poly_t(BIGMOD_RING)
    v=(l_t-l_s1-l_pk*l_s2)*inv_fal_mod

    v.redc()

    if l_s1.l2sq()<norms[0] and l_s2.l2sq()<norms[1] and v.l2sq() < norms[2]:
        assert l_s1.linf()<2**14 and l_s2.linf()<2**14 and v.linf()<2**14
        # the below makes sure that the equation we are trying to prove does not wrap around the labrador modulus - i.e. we're proving it over the integers
        assert math.sqrt(l_pk.l2sq())*math.sqrt(norms[1])+math.sqrt(norms[0])+mod*math.sqrt(norms[2]) < BIGMOD_RING.mod //2
        stat_left=[l_pk,ID,ID*mod]
        wit=[l_s2,l_s1,v]
        PS.fresh_statement(stat_left,wit,l_t)
        j+=1
    else:
        print("Too BIG in ",j)
        print(v.l2sq() < norms[2])
        VBIG+=1

sig_end=time.perf_counter()

stmnt=PS.output_statement()

PS.smpl_verify()
prove_start=time.perf_counter()
proof = PS.pack_prove()
prove_end=time.perf_counter()
ver_start=time.perf_counter()
if proof[0] == 0:
    pack_verify(proof[1:3],stmnt,PRIMESIZE)
ver_end=time.perf_counter()
print("Key creation: ",keytime_end-keytime_start)
print("Signature creation: ",sig_end-sig_start)
print("BENCHMARK: Proof Time: ",prove_end-prove_start)
print("BENCHMARK: Verification Time: ",ver_end-ver_start)
#print(VBIG)
