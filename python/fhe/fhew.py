import sys
sys.path.append('..')   # path to lazer module
from lazer import *     # import lazer python module
from lazer import _invmod
import hashlib      # for SHAKE128
import time
from labrador import *

BIGMOD_RING=polyring_t(256,LAB_RING_40.mod)
PRIMESIZE=str(math.ceil(math.log2(BIGMOD_RING.mod)))

def bootstrap_init():
    DEG=512
    reps=10
    # RING=polyring_t(DEG,12289)

    RING=polyring_t(DEG,LAB_RING_40.mod)
    PRIMESIZE=str(math.ceil(math.log2(BIGMOD_RING.mod)))


    class ptar:
        def __init__(self,s):
            self.apt=ffi.new("int64_t[]",s)
            self.s=10

        def double_size(self):
            apt2=ffi.new("int64_t[]",2*self.s)
            for i in range(self.s):
                apt2[i]=self.apt[i]
                apt2[10+i]=-1
            self.s=2*self.s
            self.apt=apt2


    #def make_array(vvec:polyvec_t):
    sz=2
    pvec=ffi.new("int64_t * []",sz)
    #    for i in range(vvec.dim):
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

    print("")
    print("")
    print("plist")
    print(plist)
    print("")
    print("")
    print("qlist")
    print(qlist)
    print("")
    print("")
    print(ans)
    print("")

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
    blah = proof_statement(deg_list,num_pols_list,norm_list,num_constraints,PRIMESIZE)
    # start=time.time()
    # for i in range(1):
    #     blah.append_witness(r)
    #     blah.append_witness(q)
    # print(time.time()-start)
    # for i in range(5):
    #     abc=polyvec_to_ar(r)
    #     printi64ar(abc,r.ring.deg*r.dim)
    #for i in range(3):
    #    blah.append_witness(r)
    #    blah.append_statement([q,r,q,p],[0,1,2,3],q)
    #    print("-----------------")

    blah.fresh_statement([q,r,p,zero_poly], [p,r,q,acc_init_witness], ans)
    blah.fresh_statement([one_poly], [3], acc_init_nand)
    # blah.fresh_statement([q,r,p], [0,1,2], ans)

    # blah.fresh_statement([q,r,p,acc_init_witness], [p,r,q,0], ans)
    # blah.fresh_statement([q,r,p,acc_init_witness], [0,1,2,0], ans)


    ans.print()
    for i in range(10000000):
        a=1
    #blah.append_witness(p)

    #ppvec=ffi.new("int64_t ** []",sz)
    #ppvec[0]=pvec
    #print(ppvec[0][1][0])


    a=ptar(10)

    for i in range(100000):
        a=1

    # for i in range(10):
    #     a.apt[i]=i
    # a.double_size()

    # for i in range(20):
    #     print(a.apt[i])
    # b=ptar(20)
    # for i in range(10):
    #     b.apt[i]=10*i
    # b.double_size()

    # for i in range(20):
    #     print(a.apt[i])

    stmnt=blah.output_statement()
    proof = blah.pack_prove()
    if proof[0]==0:
        pack_verify(proof[1:3],stmnt,PRIMESIZE)
    else:
        print("proof failed with error ",proof[0])
    blah.smpl_verify()




    # blah.simple_verify()
    # print("SUCCESS")


if __name__ == "__main__":
    bootstrap_init()