from lazer import *

LAB_RING_24=polyring_t(64,2**24-3)
LAB_RING_32=polyring_t(64,2**32-99)
LAB_RING_40=polyring_t(64,2**40-195)
LAB_RING_48=polyring_t(64,2**48-59)
MAX_C=2**63-1

def printi64ar(ar,ar_size):
    for i in range(ar_size):
        print(ar[i]," ",end="")
    print()

def l2sq_ar(ar,ar_size):
    l2=0
    for i in range(ar_size):
        l2+=ar[i]*ar[i]
    return l2

def int64_to_type(v,v_size,outtype):
    if outtype=="int64":
        return v
    elif outtype=="size_t":
        s_ar=ffi.new("size_t []",v_size)
    elif outtype=="uint64":
        s_ar=ffi.new("uint64_t []",v_size)
    elif outtype=="int16":
        s_ar=ffi.new("int16_t []",v_size)
    else:
        print("output type unknown")
    for i in range(v_size):
        s_ar[i]=v[i]
    return s_ar

def poly_to_ar(v: poly_t,outtype="int64"):
    pv=ffi.new("int64_t []",v.ring.deg)
    lib.poly_get_coeffvec_i64(pv, v.ptr)
    return int64_to_type(pv,v.ring.deg,outtype)


def polyvec_to_ar(v: polyvec_t,outtype="int64"):
    pvec=ffi.new("int64_t []",v.ring.deg*v.dim)
    for i in range(v.dim):
        temp=poly_to_ar(v.get_elem(i))
        for j in range(v.ring.deg):
            pvec[v.ring.deg*i+j]=temp[j]
    return int64_to_type(pvec,v.ring.deg*v.dim,outtype)

class proof_statement:
    """This is the prover class for the succinct proof system. It collects the statement / witness to be proved and creates a proof. 
    
    Attributes:
        cur_witness_num (int): the current witness number, when creating the witness
        cur_constraint_num (int): the current constraint number, when creating the statements
        witness_ptr: the pointer to a C structure holding the witness
        smplstmnt_ptr: the pointer to a C structure holding the statement
        commitment_ptr: a pointer to the C structure to the commitment to the witness
        composite_ptr: a pointer to an internal C structure holding the statement

    """

    def __init__(self, deg_list:list, num_pols_list:list ,norm_list:list ,num_constraints:int, primesize: str):
        """Initializer function. One should think of the witness vector as a list of either polynomials (poly_t),
        or polynomial vectors (polyvec_t). 
        
        Args:
            deg_list ([int]): a list of the degree of witness i
            num_pols_list ([int]): a list of the number polynomials in witness i
            norm_list: ([int]): a list of l2-squared norm bounds for witness i
            num_constraints (int): the number of equations (over a polynomial ring) in the statement
            primesize (str): either "24","32","40", or "48".  The prime used in the proof system is ~ 2**primesize
        """
        assert len(deg_list)==len(num_pols_list)==len(norm_list)
        print("initializing")
        self.cur_witness_num=0 #next witness number, should be len(num_witness_vectors)
        self.cur_constraint_num=0
        self.func_choose_define(primesize)

        self.witness_polys=0 #number of 64-dimensional ring elements in the witness
        # create witness_class
        self.witness_ptr=ffi.new("labrador"+primesize+"_"+"witness *")
        self.smplstmnt_ptr=ffi.new("labrador"+primesize+"_"+"smplstmnt *")

        self.commitment_ptr=ffi.new("labrador"+primesize+"_"+"commitment *")
        self.composite_ptr=ffi.new("labrador"+primesize+"_"+"composite *")


        #print(norm_list)
        #print("number of polys ",end="")
        dim_ar=ffi.new("size_t []",len(num_pols_list))
        for i in range(len(num_pols_list)):
            dim_ar[i]=num_pols_list[i]*deg_list[i]//64
        
        # lib.init_witness_raw(self.witness_ptr,size_of_array,array_dimensions -- how many 64 dim vectors)
        self.init_witness_raw(self.witness_ptr,len(num_pols_list),dim_ar)

        norms_ar=ffi.new("uint64_t []",norm_list)
        #create statement
        self.init_smplstmnt_raw(self.smplstmnt_ptr,len(num_pols_list),dim_ar,norms_ar,num_constraints)    
        # 0 is l_inf norm <= 1

    def __del__(self):
        self.free_commitment(self.commitment_ptr)
        self.free_witness(self.witness_ptr)
        self.free_composite(self.composite_ptr)
        self.free_smplstmnt(self.smplstmnt_ptr)
        self.free_comkey()

    def smpl_verify(self):
        """
        A sanity check to make sure that the input statement and witness actually satisfy the linear 
        statement and the norm bounds. May be useful in debugging. Outputs 0 if everything is correct.
        """
        print("Trying to Simple Verify")
        out=self.simple_verify(self.smplstmnt_ptr,self.witness_ptr)
        print("Simple Verify=",out)

    def pack_prove(self):
        """
        Creates the succinct proof.

        Returns:
            self.composite_ptr (C type): the proof
            self.commitment_ptr (C type): the commitment to the witness
        """
        print("Trying to Pack Prove")
        self.free_composite(self.composite_ptr) # if there was another proof created before, deallocate the memory
        error=self.composite_prove_simple(self.composite_ptr,self.commitment_ptr,self.smplstmnt_ptr,self.witness_ptr)
        # error = 0 means everything is good
        return error,self.composite_ptr,self.commitment_ptr

    # def pack_verify(self):
    #     print("Trying to Pack Verify")
    #     out=self.composite_verify_simple(self.composite_ptr,self.commitment_ptr,self.smplstmnt_ptr)
    #     print("Pack Verify=",out)

    def append_witness(self,v):
        """ Append a new witness
        
        Args:
            v (poly_t/polyvec_t): a witness to be added to the witness set. it can then be accessed
                later using its witness number, which gets consecutively increased every time a witness
                is added
        
        Returns:
            

        """
        
        if type(v) is poly_t:
            pvec=poly_to_ar(v)
            pols=1
            self.witness_polys+=v.ring.deg//64
        elif type(v) is polyvec_t:
            pvec = polyvec_to_ar(v)
            pols=v.dim
            self.witness_polys+=pols*v.ring.deg//64
        output=self.set_witness_vector_raw(self.witness_ptr,self.cur_witness_num,pols,v.ring.deg//64,pvec)    
        assert output==0
        # witness coefficients should all be < 16 bits, make sure they're centered
        self.cur_witness_num+=1
        return self.cur_witness_num-1
        
    def append_statement(self,stat_list,witnum_list,right_pol:poly_t,sec_witnum_list=None):
        """ Adds a linear statement to the proof system.

        Args:
            stat_list (list): a list of poly_t or polyvec_t (or a mix) elements
            witnum_list(list): a list of integers corresponding to the witnesses that were added
                using append_witness() or fresh_statement().  witness i in the witnum_list gets
                multiplied by element i in the stat_list
            right_pol (poly_t): the polynomial t in <ste_list,witnum_list> = t
        
        """


        assert len(stat_list)==len(witnum_list)
        stat_size=0
        ring=right_pol.ring
        len_ar=ffi.new("size_t []",len(stat_list))
      
        witnum_ar=ffi.new("size_t []",witnum_list)
        if type(sec_witnum_list) is list:
            assert len(witnum_list)==len(sec_witnum_list)
            sec_wit_ar=ffi.new("size_t []",witnum_list)
            quadratic=True
        else:
            quadratic=False

        for i in range(len(stat_list)):
            if type(stat_list[i]) is poly_t:
                stat_size+=1 #number of polynomials in stat_list increases by 1
                len_ar[i]=1 #number of polynomials 
            elif type(stat_list[i]) is polyvec_t:
                stat_size+=stat_list[i].dim
                len_ar[i]=stat_list[i].dim
            else:
                print("Error")
        
        stat_vec=polyvec_t(ring,stat_size,stat_list) # concatenation of all poly/polyvec in stat_list into one polyvec
        stat_ar=polyvec_to_ar(stat_vec) # convert polyvec to array
        right_ar=poly_to_ar(right_pol)

        output=self.set_smplstmnt_lincnst_raw(self.smplstmnt_ptr,self.cur_constraint_num,len(stat_list),witnum_ar,len_ar,ring.deg//64,stat_ar,right_ar)
        #print("Set Statement Output",output)
        assert output==0
        # should all be less than q, should centralize everything
        self.cur_constraint_num+=1

    def fresh_statement(self,stat_list,wit_list,right_pol:poly_t):
        """ Input a new constraint of the form <a,w>=t over a polynomial ring to the statement 

        Args:
            stat_list ([poly_t/polyvec_t]): a list of poly_t or polyvec_t comprising a
            wit_list ([poly_t/polyvec_t/int]): the list consists of i elements such that the ith element of
                wit_list is multiplied by the ith element of stat_list (make sure they're both either poly_t
                or polyvec_t!). If the ith element of stat_list is to be multiplied by a witness added in a
                prior fresh_statement, then one should enter the witness number of that witness (it will then
                be necessary to keep track of the order in which the witnesses are being entered)
            right_pol (poly_t): the t part of the constraint
        """
        witnum_list=[0]*len(wit_list)
        for i in range(len(wit_list)):
            if type(wit_list[i]) is poly_t or type(wit_list[i]) is polyvec_t:
                witnum_list[i]=self.cur_witness_num
                self.append_witness(wit_list[i])
            else: # type is integer
                assert wit_list[i]<self.cur_witness_num
                witnum_list[i]=wit_list[i]
        self.append_statement(stat_list,witnum_list,right_pol)

    def output_statement(self):
        """Outputs the statement that will be proved
        """
        return self.smplstmnt_ptr

    # def fresh_mat_vec(self,stat_mat_list,wit_vec_list,right_vec:polyvec_t):
    #     wit_list_num=[0]*len(wit_vec_list)
    #     for i in range(len(wit_vec_list)):
    #         wit_list_num[i]=self.cur_witness_num
    #         self.append_witness(wit_vec_list[i])
    #     tot_stat=stat_mat_list[0].rows
    #     for i in range(tot_stat):
    #         stat_list=[]
    #         right_pol=right_vec.get_elem(i)
    #         for j in range(len(stat_mat_list)):
    #             if type(stat_mat_list[j]) is polymat_t:
    #                 temp_vec=stat_mat_list[j].get_row(i)
    #             elif type(stat_mat_list[j]) is polyvec_t:
    #                 temp_vec=stat_mat_list[j].get_elem(i)
    #             stat_list.append(temp_vec)
    #         self.fresh_statement(stat_list,wit_list_num,right_pol)
    
    #based on the size of the prime, different C functions are used
    def func_choose_define(self,primesize: str):
        if primesize == "24":
                self.init_witness_raw=lib.labrador24_init_witness_raw
                self.init_smplstmnt_raw=lib.labrador24_init_smplstmnt_raw
                self.free_commitment=lib.labrador24_free_commitment
                self.free_witness=lib.labrador24_free_witness
                self.free_composite=lib.labrador24_free_composite
                self.free_smplstmnt=lib.labrador24_free_smplstmnt
                self.simple_verify=lib.labrador24_simple_verify
                self.composite_prove_simple=lib.labrador24_composite_prove_simple
                self.composite_verify_simple=lib.labrador24_composite_verify_simple
                self.set_witness_vector_raw=lib.labrador24_set_witness_vector_raw
                self.set_smplstmnt_lincnst_raw=lib.labrador24_set_smplstmnt_lincnst_raw
                self.init_comkey=lib.labrador24_init_comkey
                self.free_comkey=lib.labrador24_free_comkey
        elif primesize == "32":
                self.init_witness_raw=lib.labrador32_init_witness_raw
                self.init_smplstmnt_raw=lib.labrador32_init_smplstmnt_raw
                self.free_commitment=lib.labrador32_free_commitment
                self.free_witness=lib.labrador32_free_witness
                self.free_composite=lib.labrador32_free_composite
                self.free_smplstmnt=lib.labrador32_free_smplstmnt
                self.simple_verify=lib.labrador32_simple_verify
                self.composite_prove_simple=lib.labrador32_composite_prove_simple
                self.composite_verify_simple=lib.labrador32_composite_verify_simple
                self.set_witness_vector_raw=lib.labrador32_set_witness_vector_raw
                self.set_smplstmnt_lincnst_raw=lib.labrador32_set_smplstmnt_lincnst_raw
                self.init_comkey=lib.labrador32_init_comkey
                self.free_comkey=lib.labrador32_free_comkey
        elif primesize == "40":
                self.init_witness_raw=lib.labrador40_init_witness_raw
                self.init_smplstmnt_raw=lib.labrador40_init_smplstmnt_raw
                self.free_commitment=lib.labrador40_free_commitment
                self.free_witness=lib.labrador40_free_witness
                self.free_composite=lib.labrador40_free_composite
                self.free_smplstmnt=lib.labrador40_free_smplstmnt
                self.simple_verify=lib.labrador40_simple_verify
                self.composite_prove_simple=lib.labrador40_composite_prove_simple
                self.composite_verify_simple=lib.labrador40_composite_verify_simple
                self.set_witness_vector_raw=lib.labrador40_set_witness_vector_raw
                self.set_smplstmnt_lincnst_raw=lib.labrador40_set_smplstmnt_lincnst_raw
                self.init_comkey=lib.labrador40_init_comkey
                self.free_comkey=lib.labrador40_free_comkey
        elif primesize == "48":
                self.init_witness_raw=lib.labrador48_init_witness_raw
                self.init_smplstmnt_raw=lib.labrador48_init_smplstmnt_raw
                self.free_commitment=lib.labrador48_free_commitment
                self.free_witness=lib.labrador48_free_witness
                self.free_composite=lib.labrador48_free_composite
                self.free_smplstmnt=lib.labrador48_free_smplstmnt
                self.simple_verify=lib.labrador48_simple_verify
                self.composite_prove_simple=lib.labrador48_composite_prove_simple
                self.composite_verify_simple=lib.labrador48_composite_verify_simple
                self.set_witness_vector_raw=lib.labrador48_set_witness_vector_raw
                self.set_smplstmnt_lincnst_raw=lib.labrador48_set_smplstmnt_lincnst_raw
                self.init_comkey=lib.labrador48_init_comkey
                self.free_comkey=lib.labrador48_free_comkey
    

def pack_verify(proof,stmnt_ptr,primesize: str):
    """The verification procedure for the succinct proof, which takes in the output of pack_prove and the statement
    
    """
    comp_ptr,comm_ptr = proof[0],proof[1]
    print("Trying to Pack Verify")
    if primesize == "24":
            out=lib.labrador24_composite_verify_simple(comp_ptr,comm_ptr,stmnt_ptr)
    elif primesize == "32":
            out=lib.labrador32_composite_verify_simple(comp_ptr,comm_ptr,stmnt_ptr)
    elif primesize == "40":
            out=lib.labrador40_composite_verify_simple(comp_ptr,comm_ptr,stmnt_ptr)
    elif primesize == "48":
            out=lib.labrador48_composite_verify_simple(comp_ptr,comm_ptr,stmnt_ptr)
    print("Pack Verify=",out)
