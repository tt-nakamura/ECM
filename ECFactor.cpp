// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   R. Crandall and C. Pomerance
//    "Prime Numbers: A Computational Perspective"
//     2nd edition, Algorithm 7.4.4

#include "EC_p.h"

#define EC_BOUND  (1<<13)// stage-1 limit for 1st try
#define EC_BOUND2 (1<<7) // (stage-2 limit)/(stage-1 limit)
#define EC_NSTEP  (1<<8) // memory size in stage-2
#define EC_BMAX   (1<<30)// maximum limit

long ECFactor(ZZ& d, const ZZ& n, double timeout)
// factorization by elliptic curve method
// d = proper divisor of n
// n = composite integer, n>=33
// quit factoring after timeout seconds
// return 0 if successful, -1 if failure
{
    long p,i,j,k,B(EC_BOUND), B2, D(EC_NSTEP<<1);
    ZZ_p::init(n);
    ZZ_p s,u,v,c,xz[EC_NSTEP];
    EC_p P,Q,R;
    EC_p S[EC_NSTEP], &T(S[EC_NSTEP-1]);
    double lnB;
    timeout += GetTime();
a:  ;
    if(GetTime() > timeout) return -1;
    // Crandall & Pomerance, Theorem 7.4.3
    do random(s); while(rep(s)<=5);
    sqr(u,s); u-=5;
    mul(v,s,4);
    power(P.x, u, 3);
    power(P.z, v, 3);
    sub(c,v,u);
    power(c,c,3);
    mul(s,u,3);
    s+=v;
    c*=s;
    mul(s, P.x, v);
    s*=4;
    GCD(d,n,rep(s));
    if(d==n) goto a;
    if(!IsOne(d)) return 0;
    c/=s;
    c-=2;
    EC_p::init(c);
    // stage-one
    lnB = log(B);
    PrimeSeq ps;
    while((p = ps.next()) && p<=B) {
        k = long(lnB/log(p));
        for(i=0; i<k; i++) mul(P,P,p);
    }
    GCD(d, n, rep(P.z));
    if(d==n) { B>>=1; goto a; }
    if(!IsOne(d)) return 0;
    // stage-two
    B2 = B*EC_BOUND2;// stage-two limit
    if(B2>EC_BMAX || B2<0) return -1;
    doubleh(S[0], P);
    doubleh(S[1], S[0]);
    for(i=2; i<EC_NSTEP; i++) addh(S[i], S[i-1], S[0], S[i-2]);
    for(i=0; i<EC_NSTEP; i++) mul(xz[i], S[i].x, S[i].z);
    set(c);
    j=B-1;
    mul(Q,P,j-D);
    mul(R,P,j);
    for(; j<=B2; j+=D) {
        mul(s, R.x, R.z);
        k = min(j+D, B2);
        while(p && p<=k) {
            i = ((p-j)>>1)-1;
            sub(u, R.x, S[i].x);
            add(v, R.z, S[i].z);
            u*=v;
            u-=s;
            u+=xz[i];
            c*=u;
            p = ps.next();
        }
        GCD(d, n, rep(c));
        if(!IsOne(d)) break;
        if(GetTime() > timeout) return -1;
        addh(P,R,T,Q);
        Q = R;
        R = P;
    }
//  GCD(d, n, rep(c));
    if(d==n) goto a;
    if(!IsOne(d)) return 0;
    if((B<<=1)>EC_BMAX || B<0) return -1;
    goto a;// retry with doubled limit
}