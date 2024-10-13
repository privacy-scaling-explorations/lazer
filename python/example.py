from lazer import *     # import everything from the lazer python module
SEED=[0]


print("")
print("working with integers ...")

a = int_t(1)

print("a =",a)

b = int_t(2)
c = int_t(4)
d = int_t(8)
y = int_t((2**64)-1)
z = int_t((2**64))

print("a =", a)
print("b =", b)
print("c =", c)
print("d =", d)
print("y =", y)
print("z =", z)

print("(a + c) * (-b) =", (a + c) * (-b))
print("a =", a)
print("2*a=", 2*a)

d += a
print("d += a, d =", d)

d -= b
print("d -= b, d =", d)

assert a < b
assert -a < a
assert -b < -a
assert a <= b
assert a <= a
assert -a <= b
assert -b <= -a
assert -a <= -a
assert b > a
assert a > -a
assert -a > -b
assert b >= a
assert a >= a
assert b >= -a
assert -a >= -b
assert -a >= -a
assert a == a
assert -a == -a 

print("")
print("creating a ring ...")

Rq = polyring_t(d=64, q=12289)

print("ring degree  =", Rq.deg)
print("ring modulus =", Rq.mod)


print("")
print("working with polynomials ...")

f = poly_t(Rq, 61 * [0] + [1, 0, -2])
g = poly_t(Rq, {1: 3})

print("f =", f)
print("g = ", g)
print("f * 2 =", f*2)
assert f == f

h = f * g
print("f * g =", h)


h = f + g
print("f + g =", h)

h = f - g
print("f - g =", h)

f += g
print("f += g, f =", f)

f -= g
print("f -= g, f =", f)

print("-f =", -f)

f.set_coeffs({1: 8, 31: 3, 62: -9})
print("f =", f)

f.urandom(int_t(7),bytes([1]*32),2)
print("f =",f)

f.urandom(7,bytes([1]*32),2)
print("f =",f)

#v=intvec_t(5,1,[0,1,2,3,4])
#v.print()
#v=v*3
#v.print()

v=polyvec_t(Rq,2)
v.print()
v.set_elem(5,1,3)
v.print()
v1=v.get_elem(1,3)
print(v1)


v=intvec_t(101)
v.urandom_bnd(3,4,SEED,0,0)
v.print()
v.urandom_bnd(3,4,SEED,1,1)
v.print()
v.urandom_bnd(3,4,SEED,0,1)
v.print()
w=v*v
print(w)
v.print()
m=intmat_t(100,101)
m.urandom_bnd(0,5,SEED,0,1)
m.print()
(2*m*v).print()
print(v.l2sqr())
m2=m*(m.transpose())
m2.print()
m3=m2.copy()
m3.print()

p=poly_t(Rq)
p.urandom_bnd(-1,1,SEED,0,1)
print(p)
p.urandom_bnd(-1,1,SEED,0,1)
print("p=")
print(p)
v1=polyvec_t(Rq,3)
v2=polyvec_t(Rq,3)
v1.urandom_bnd(0,3,SEED,0,1)
v2.urandom_bnd(-3,0,SEED,0,1)
m=polymat_t(Rq,2,3)
m.set_row(0,v1)
m.set_row(1,v2)
m.print()
(m*2).print()
(m*v1).print()
v3=polyvec_t(Rq,2)
v3.urandom_bnd(-3,0,SEED,0,1)
(v3*m).print()
(-v3*m).print()
m.print()
m2=m.copy()
m2.set_elem(p,0,0)
m2.print()
print(m==m2)
m.brandom(3,SEED,0,1)
m.print()
m3=polymat_t(Rq,3,3)
m3.brandom(2,SEED,0,1)
m3.print()
(m3*m3).print()


p=falcon_pol([1]*512)
p.print()
p.set_pos(1,12)
p.print()
print(p.get_pos(1))
sump=(p+p+p)
sump.print()
p.set_pos(1,0)
p.print()
#lib.falcon_decode_pubkey(p, (0).to_bytes(32,'little'))

print("")
print("creating a ring ...")
Rp = polyring_t(512,q=12289)

t=falcon_pol([1]*512)
skenc,pkenc,pk=falcon_keygen()
s1,s2=falcon_preimage_sample(skenc,t)

pkmat=pk.to_polymat()
s2vec=s2.to_polyvec()

#pkq,skq=lin_to_isoring(Rq,Rp,pkmat,s2vec)
pkq,skq=falcon_poly_mul_toisoring(Rq,pk,s2)
firmul=pkq*skq
(firmul).print()
fp=from_isoring_tofalconpol(pkq*skq)
(fp+s1).print()

#s1.print()
#s2.print()
s1_t=s1.to_isoring(Rq)
(firmul+s1_t).print()


#keygen
Rq=polyring_t(64,q=12289)
Rp = polyring_t(512,q=12289)
skenc,pkenc,pk=falcon_keygen()
Fseed=SEED[0]
Fmat=polymat_t.urandom_static(Rq,8,4,12289,SEED,0,1)
B1seed=SEED[0]
B1mat=polymat_t.urandom_static(Rq,8,4,12289,SEED,0,1)
B2seed=SEED[0]
B2mat=polymat_t.urandom_static(Rq,8,16,12289,SEED,0,1)

#user move
rvec=polyvec_t.grandom_static(Rq,16,6,SEED,0,1)
mvec=polyvec_t.urandom_bnd_static(Rq,4,0,1,SEED,0,1)
msg1=B1mat*mvec+B2mat*rvec

#user sends msg1
#signer move
xvec=polyvec_t.urandom_bnd_static(Rq,4,0,1,SEED,0,1)
preq=msg1+Fmat*xvec
s1,s2=falcon_preimage_sample(skenc,preq)

#signer sends s1,s2, and xvec
Fmatv=polymat_t.urandom_static(Rq,8,4,12289,Fseed,0)
B1v=polymat_t.urandom_static(Rq,8,4,12289,B1seed,0)
B2v=polymat_t.urandom_static(Rq,8,16,12289,B2seed,0)
(msg1-pk*s2-s1+Fmatv*xvec).print()

#for anonymous credentials, some part of mvec can be revealed
pub_mvec_col=[0,1,2]
priv_mvec_col=list(set(range(4))-set(pub_mvec_col))
mvec_pub=mvec.get_pol_list(pub_mvec_col)
mvec_priv=mvec.get_pol_list(priv_mvec_col)
B1v_pub=B1v.get_col_list(pub_mvec_col)
B1v_priv=B1v.get_col_list(priv_mvec_col)

(pk*s2+s1-B1v*mvec-B2v*rvec-Fmatv*xvec).print()

print("full")
(B1v*mvec).print()
print("sum")
(B1v_pub*mvec_pub+B1v_priv*mvec_priv).print()
print("pub")
(B1v_pub*mvec_pub).print()
print("priv")
(B1v_priv*mvec_priv).print()
# print(s1.linf())
# print(s2.l2sqr())
# print(rvec.linf())

##############
#different rings add and subtract
polp=polyvec_t.urandom_bnd_static(Rp,2,0,3,SEED,0,1)
polq=polyvec_t.urandom_bnd_static(Rq,16,-3,3,SEED,0,1)
polp2=polp+polp
polpq=polp+polq
(-polpq+polp2-polp+polq).print()
#rvec.print()
Rr=polyring_t(256,12289)
Br=polymat_t.urandom_static(Rr,2,2,12,SEED,1)
#Br.print()
vr=polyvec_t.brandom_static(Rr,2,1,SEED,0,1)
(Br*vr-B1v*mvec+s1-Br*vr+B1v*mvec).print()
res=vr+s1-vr
res=res.from_isoring(Rp)
res.print()
print(res.ring.deg)
s1.print()

e0=polyvec_t.urandom_bnd_static(Rr,2,-1,1,SEED,0,1)
e1=poly_t.urandom_bnd_static(Rr,-1,1,SEED,0,1)
e2=polyvec_t.urandom_bnd_static(Rr,2,-1,1,SEED,0,1)

(Br*vr).print()
t=Br*vr+e0
vr2=polyvec_t.brandom_static(Rr,2,1,SEED,0,1)
c1=vr2*Br+e2
c2=vr2*t+e1
#e1.print()
(c1*vr-c2).print()
Br.print()
print("vr=")
print(vr.ring.deg)
vr.print()
print("br=")
(vr*vr).print()
# res=Br*vr
# test_print("result dimension ")
# test_print(res.dim)
# (res).print()
# res=(vr*Br)
# test_print("result2 dimension ")
# test_print(res.dim)
# (res).print()

"""
fvec=polyvec_t(Rq,8)
fvec.urandom_bnd(0,5,SEED,0,1)
(fmat*fvec+firmul).print()





p1=poly_t(Rp,[1]*512)
p1.urandom_bnd(-2,2,SEED,0,1)
p1.print()
tp=p1.to_isoring(Rq)
tp.print()

print("ring degree  =", Rp.deg)
print("ring modulus =", Rp.mod)
"""






"""
v=intvec_t(3,[1,2,3])
v.print()
m.row[1]=intvec_t(m.ncols,v)
m.print()
"""