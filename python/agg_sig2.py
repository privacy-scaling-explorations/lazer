from lazer import *
from lazer import _invmod, _center_list, _l2sq_list
from labrador import *
import hashlib      # for SHAKE128
import secrets      # for internal coins
import sys
import time

sig_num=1024
norms=[17017363,17017363,round(1248245003*.75)]

VBIG=0 # counts how many signatures are too big

#Total Signatures
sig_num=250
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

ZERO=int_to_poly(0,BIGMOD_RING)
ID=int_to_poly(1,BIGMOD_RING)

# public randomness
shake128 = hashlib.shake_128(bytes.fromhex("00"))
TARGPP = shake128.digest(32)

inv_fal_mod=_invmod(mod,BIGMOD_RING.mod)

deg_list=[deg]*(3*sig_num)
num_pols_list=[sig_num]*(3*sig_num)
norm_list=norms*sig_num
num_constraints=sig_num
PS=proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)

# generating falcon pk/sk list
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

# prepare reshaped witnesses
# reference: https://hackmd.io/5R5f2o30SCeEz8U52wiUTQ
rho = round(math.sqrt(sig_num))
vs=[]
# stores s_1 of each signature
y_1=[[ZERO for _ in range(sig_num)] for _ in range(math.ceil(sig_num / rho))]
# stores s_2 of each signature
y_2=[[ZERO for _ in range(sig_num)] for _ in range(math.ceil(sig_num / rho))]
# stores s_1'
y_1_prime = [[ZERO for _ in range(sig_num)] for _ in range(rho)]
# stores s_2'
y_2_prime = [[ZERO for _ in range(sig_num)] for _ in range(rho)]
# delta_i's
delta = [[ID if i == j else ZERO for j in range(sig_num)] for i in range(sig_num)]
# h_i * delta_i's
if SAME_KEY:
    h_times_delta = [[l_pk if i == j else ZERO for j in range(sig_num)] for i in range(sig_num)]
else:
    h_times_delta = [[pk_list[i] if i == j else ZERO for j in range(sig_num)] for i in range(sig_num)]
# q * delta_i's
q_times_delta = [[ID*mod if i == j else ZERO for j in range(sig_num)] for i in range(sig_num)]
# t_i's
ts=[]

def index(i: int) -> int:
    return math.ceil(i / rho)

def index_prime(i: int) -> int:
    return ((i - 1) % rho) + 1

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
        # stat_left=[l_pk,ID,ID*mod]
        # wit=[l_s2,l_s1,v]
        # PS.fresh_statement(stat_left,wit,l_t)
        vs.append(v)
        y_1[index(j)-1][j] = l_s1
        y_2[index(j)-1][j] = l_s2
        ts.append(l_t)
        j+=1
    else:
        print("Too BIG in ",j)
        print(v.l2sq() < norms[2])
        VBIG+=1

for i in range(sig_num):
    delta_i = polyvec_t(BIGMOD_RING, sig_num, delta[i])
    h_times_delta_i = polyvec_t(BIGMOD_RING, sig_num, h_times_delta[i])
    q_times_delta_i = polyvec_t(BIGMOD_RING, sig_num, q_times_delta[i])
    y_1_i = polyvec_t(BIGMOD_RING, sig_num, y_1[index(i) - 1])
    y_2_i = polyvec_t(BIGMOD_RING, sig_num, y_2[index(i) - 1])
    stat_left=[delta_i,h_times_delta_i,q_times_delta_i]
    vs_polyvec = polyvec_t(BIGMOD_RING, sig_num, vs)
    wit=[y_1_i,y_2_i,vs_polyvec]
    PS.fresh_statement(stat_left,wit,ts[i])

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
print("Proof Time: ",prove_end-prove_start)
print("Verification Time: ",ver_end-ver_start)
#print(VBIG)

