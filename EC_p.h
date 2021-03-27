// uses NTL
//   http://www.shoup.net/ntl

#ifndef __EC_p_h__
#define __EC_p_h__

#include<NTL/ZZ_pX.h>
using namespace NTL;

struct EC_p {// Elliptic Curve over finite field F_p
    ZZ_p x,z;// montgomery coordinates
    static ZZ d;// d=gcd(z,p) to check d==1
    static ZZ_pX f;// cubic polynomial x^3+cx^2+x
    inline static ZZ_p& c() { return f[2]; }
    static long IsSingular();
    static void init(const ZZ_p&);
};

inline long IsZero(const EC_p& a) { return IsZero(a.z); }
inline void clear(EC_p& a) { clear(a.z); }
inline void set(EC_p& a, const ZZ_p& x1) { a.x=x1; set(a.z); }
long compare_x(const EC_p&, const EC_p&);
void addh(EC_p&, const EC_p&, const EC_p&, const EC_p&);
void doubleh(EC_p&, const EC_p&);
void mul(EC_p&, const EC_p&, const ZZ&);
void mul(EC_p&, const EC_p&, long);
void mul(EC_p&, EC_p&, const EC_p&, const ZZ&);
long x_coord(ZZ&, const EC_p&);
long z_coprime(const EC_p&);
void random(EC_p& a);
std::ostream& operator<<(std::ostream&, const EC_p&);

#endif // __EC_p_h__