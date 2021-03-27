// uses NTL
//   http://www.shoup.net/ntl

#include "EC_p.h"

ZZ EC_p::d;
ZZ_pX EC_p::f;

void EC_p::init(const ZZ_p& c)
// y^2 = x^3 + cx^2 + x
// assume ZZ_p has been initialized
// kill f to avoid "internal error: can't grow this _ntl_gbigint"
{
    f.kill();
    SetCoeff(f,3);
    SetCoeff(f,1);
    SetCoeff(f,2,c);
}

long EC_p::IsSingular()
// return (discriminant 4a^3+27b^2 == 0)
{
    ZZ_p u;
    sqr(u, f[2]);
    u-=4;
    return IsZero(u);
}

long compare_x(const EC_p& a, const EC_p& b)
// return (a.x/a.z == b.x/b.z)
{
    ZZ_p s,t;
    mul(s, a.x, b.z);
    mul(t, a.z, b.x);
    return s==t;
}

std::ostream& operator<<(std::ostream& s, const EC_p& a) {
    s << '[' << a.x << ' ' << a.z << ']';
    return s;
}

void addh(EC_p& p, const EC_p& a, const EC_p& b, const EC_p& q)
// set p=a+b; require q=|a-b|; assume &p!=&q
// reference: Crandall & Pomerance, eq(7.6) with (A,B)=(1,0)
{
    if(IsZero(a)) p=b;
    else if(IsZero(b)) p=a;
    else if(IsZero(q)) doubleh(p,a);
    else if(IsZero(q.x)) {
        ZZ_p s,t,u,v,w;
        mul(s, a.x, b.z);
        mul(t, a.z, b.x);
        add(u,s,t);
        sub(v,s,t);
        mul(w,s,t);
        w *= EC_p::c();
        w *= 2;
        add(s, a.x, a.z);
        add(t, b.x, b.z);
        s*=t;
        s-=u;
        sqr(p.z, v);
        mul(p.x, u, s);
        p.x += w;
        p.x *= 2;
    }
    else {
        ZZ_p s,t,u;
        mul(s, a.x, b.z);
        mul(t, a.z, b.x);
        s-=t;
        sub(t, a.x, a.z);
        add(u, b.x, b.z);
        t*=u;
        t-=s;
        sqr(p.z, s);
        sqr(p.x, t);
        p.z *= q.x;
        p.x *= q.z;
    }
}

void doubleh(EC_p& p, const EC_p& a)
// set p=a+a
// reference: Crandall & Pomerance, eq(7.7) with (A,B)=(1,0)
{
    if(IsZero(a)) { clear(p); }
    else {
        ZZ_p x2,z2,xz;
        sqr(x2, a.x);
        sqr(z2, a.z);
        mul(xz, a.x, a.z);
        sub(p.x, x2, z2);
        sqr(p.x, p.x);
        mul(p.z, xz, EC_p::c());
        p.z += x2;
        p.z += z2;
        p.z *= xz;
        p.z *= 4;
    }
}

void mul(EC_p& p, const EC_p& a, const ZZ& n)
// set p = n*a; assume n>0
// reference: Crandall & Pomerance, Algorithm 7.2.7
{
    if(IsZero(a) || IsZero(n)) clear(p);
    else if(IsOne(n)) p=a;
    else if(n==2) doubleh(p,a);
    else if(&p==&a) {
        EC_p b(a);
        mul(p,b,n);
    }
    else {
        EC_p q;
        p = a;
        doubleh(q,a);
        for(int i=NumBits(n)-2; i>=1; i--) {
            if(bit(n,i)) {
                addh(p,q,p,a);
                doubleh(q,q);
            }
            else {
                addh(q,q,p,a);
                doubleh(p,p);
            }
        }
        if(bit(n,0)) addh(p,q,p,a);
        else doubleh(p,p);
    }
}

void mul(EC_p& p, EC_p& q, const EC_p& a, const ZZ& n)
// set p = n*a, set q = (n+1)*a; assume n>0
{
    if(IsZero(a)) { clear(p); clear(q); }
    else if(IsZero(n)) { clear(p); q=a; } 
    else if(IsOne(n)) { p=a; doubleh(q,a); }
    else if(n==2) { doubleh(p,a); addh(q,p,a,a); }
    else if(&p==&a || &q==&a) {
        EC_p b(a);
        mul(p,q,b,n);
    }
    else {
        p = a;
        doubleh(q,a);
        for(int i=NumBits(n)-2; i>=0; i--) {
            if(bit(n,i)) {
                addh(p,q,p,a);
                doubleh(q,q);
            }
            else {
                addh(q,q,p,a);
                doubleh(p,p);
            }
        }
    }
}

void mul(EC_p& p, const EC_p& a, long n)
// set p = n*a; assume n>0
{
    if(IsZero(a) || n==0) clear(p);
    else if(n==1) p=a;
    else if(n==2) doubleh(p,a);
    else if(&p==&a) {
        EC_p b(a);
        mul(p,b,n);
    }
    else {
        unsigned long k(1L<<(NumBits(n)-1));
        EC_p q;
        p = a;
        doubleh(q,a);
        for(k>>=1; k>1; k>>=1) {
            if(n&k) {
                addh(p,q,p,a);
                doubleh(q,q);
            }
            else {
                addh(q,q,p,a);
                doubleh(p,p);
            }
        }
        if(n&1) addh(p,q,p,a);
        else doubleh(p,p);
    }
}

long x_coord(ZZ_p& x, const EC_p& a)
// if a.z is invertible mod p, set x=a.x/a.z and return 0
// else set d=(divisor of p) and return 1
{
    ZZ s,t;
    XGCD(EC_p::d, s, t, rep(a.z), ZZ_p::modulus());
    if(!IsOne(EC_p::d)) return 1;
    if(sign(s) < 0) s += ZZ_p::modulus();
    conv(x,s);
    x *= a.x;
    return 0;
}

long z_coprime(const EC_p& a)
// if a.z is invertible mod p, return 1
// else set d=(divisor of p) and return 0
{
    GCD(EC_p::d, rep(a.z), ZZ_p::modulus());
    return IsOne(EC_p::d);
}

void random(EC_p& p)
// set p = random point on the curve (except infty)
{
    do {
        random(p.x);
        eval(p.z, EC_p::f, p.x);
    } while(Jacobi(rep(p.z), ZZ_p::modulus()) < 0);
    set(p.z);
}