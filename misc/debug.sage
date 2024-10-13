# Some stuff that may help debugging

Rq = QuotientRing(Integers(q)['x'], Integers(q)['x'].ideal(x ^ d+1), 'x')


def redc(z, p):
    z = int(int(z) % int(p))
    if z > (p-1)/2:
        z = z - p
    if z < -(p-1)/2:
        z = z + p
    return z


def l2sqr(vec, q):
    acc = 0
    for i in range(len(vec)):
        elem = vec[i]
        for j in range(elem.parent().degree()):
            acc += redc(elem[j], q) ^ 2
    return acc


def auto(p):
    global Rq
    global x
    d = Rq.degree()
    r = Rq.zero()
    r += Rq(p[0])
    for i in range(1, d):
        r += Rq(-ZZ(p[d-i])*x ^ i)
    return r


def auto_(v):
    global Rq
    global x
    d = Rq.degree()
    r = vector(Rq, [Rq.zero() for i in range(len(v))])
    for i in range(len(v)):
        r[i] = auto(v[i])
    return r


def auto__(m):
    global Rq
    global x
    d = Rq.degree()
    r = matrix(Rq, m.nrows(), m.ncols())
    for i in range(m.nrows()):
        for j in range(m.ncols()):
            r[i, j] = auto(m[i, j])
    return r


def Ux(v):
    assert len(v) % 2 == 0
    for i in range(len(v) / 2):
        tmp = v[2*i]
        v[2*i] = v[2*i+1]
        v[2*i+1] = tmp
    return v


def UxU(m):
    assert m.nrows() % 2 == 0
    assert m.ncols() % 2 == 0
    for i in range(m.nrows() / 2):
        for j in range(m.ncols() / 2):
            tmp = m[2*i, 2*j]
            m[2*i, 2*j] = m[2*i+1, 2*j+1]
            m[2*i+1, 2*j+1] = tmp

            tmp = m[2*i, 2*j+1]
            m[2*i, 2*j+1] = m[2*i+1, 2*j]
            m[2*i+1, 2*j] = tmp
    return m


def autov(s1, m):
    global Rq
    r = vector(Rq, [Rq.zero() for i in range(2*(len(s1)+len(m)))])
    for i in range(len(s1)):
        r[2*i] = s1[i]
        r[2*i+1] = auto(s1[i])
    for i in range(len(m)):
        j = 2*len(s1)
        r[j+2*i] = m[i]
        r[j+2*i+1] = auto(m[i])
    return r


def tr(poly):
    return (poly+auto(poly))/2


def mkpoly(c0, cdby2):
    global Rq
    global x
    d = Rq.degree()
    return Rq(ZZ(c0) + ZZ(cdby2) * x ^ (ZZ(d/2)))
