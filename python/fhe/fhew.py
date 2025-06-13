import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
from lazer import _invmod
import hashlib      # for SHAKE128
import time
from labrador import *

def bootstrap_init():
    DEG=128
    RING=polyring_t(DEG,LAB_RING_32.mod)
    PRIMESIZE=str(math.ceil(math.log2(RING.mod)))

    # DEG is N in circom-fhe
    q = 32
    # Q = 134215681
    Q = RING.mod
    Qks = 1 << 14  
    Bks = 128
    Bg = 128
    Br = 32

    # f maps [3q/8, 7q/8) -> -Q/8; [-q/8, 3q/8) -> Q/8
    f = []
    Q8 = Q // 8 + 1
    Q8Neg = Q - Q8

    for i in range(q):
        if 3*q <= 8*i < 7*q:
            f.append(Q8Neg)
        else:
            f.append(Q8)


    plist=[1,-1,10,11]
    plist.extend([0]*(DEG-5))
    plist.extend([7])
    qlist=[-6,2,-5,1]
    qlist.extend([0]*(DEG-5))
    qlist.extend([-1])
    p=poly_t(RING,plist)
    q=poly_t(RING,qlist)
    r=polyvec_t(RING,2,[p,q])
    ans=poly_t(RING)
    ans=q*p+r*r+p*q
    r.print()

    #
    zero_poly = poly_t(RING, [0]*DEG)
    one_poly = poly_t(RING, [1] + [0]*(DEG-1))

    acc_init_witness_list = [1, 2, 3, 4]
    acc_init_witness_list.extend([0]*(DEG-4))
    acc_init_witness = poly_t(RING, acc_init_witness_list)

    acc_init_nand_list = [1, 2, 3, 4]
    acc_init_nand_list.extend([0]*(DEG-4))
    acc_init_nand = poly_t(RING, acc_init_nand_list)

    witness_num = 4

    deg_list=[DEG] * witness_num
    num_pols_list = [1,2,1,1]
    # num_pols_list = [1,2,1]
    # num_pols_list = [1]
    norm_list = [100000] * witness_num
    num_constraints = 2
    ps = proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)

    # LaBRADOR (linear) constraints are written as \sum a_i w_i = t.
    # The first fresh_statement argument is a list of polynomials a_i.
    # The second fresh_statement argument is a list of witness vectors.
    # Note: in the first two arguments, you can use polynomial vectors as well (not
    # only polynomials), but you need to specify the length of the vector in num_pols_list.

    ps.fresh_statement([q,r,p,zero_poly], [p,r,q,acc_init_witness], ans)
    ps.fresh_statement([one_poly], [3], acc_init_nand)
    # ps.fresh_statement([q,r,p], [0,1,2], ans)

    print("")
    ans.print()
    print("")

    stmnt=ps.output_statement()
    proof = ps.pack_prove()
    if proof[0]==0:
        pack_verify(proof[1:3],stmnt,PRIMESIZE)
    else:
        print("proof failed with error ",proof[0])
    ps.smpl_verify()

def bootstrap():
    bootstrap_init()

if __name__ == "__main__":
    bootstrap()