#include<NTL/ZZ.h>
using namespace NTL;

long ECFactor(ZZ&, const ZZ&, double);

main() {
    double T(10);
    ZZ p,q,n,d;
    GenPrime(p,50);
    GenPrime(q,50);
    mul(n,p,q);
    std::cout << "factoring " << n << std::endl;
    if(ECFactor(d,n,T)) std::cout << "time out";
    else if(d!=p && d!=q) std::cout << "wrong";
    else std::cout << '=' << d << '*' << (n/d);
    std::cout << std::endl;
}