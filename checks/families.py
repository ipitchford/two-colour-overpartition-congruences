"""Exact corroboration; the universal result is the written theta proof."""
from math import isqrt
import json
from pathlib import Path
from flint import nmod_poly
def need(b,msg):
    if not b: raise ValueError(msg)
def odd(n):
    while n%2==0:n//=2
    return n
def coefficients(ell,L,mod=16):
    def theta(d):
        c=[0]*L;c[0]=1
        for k in range(1,isqrt((L-1)//d)+1):c[d*k*k]=2*(-1)**k%mod
        return nmod_poly(c,mod)
    a=theta(1).inverse_series_trunc(L)
    b=theta(ell).inverse_series_trunc(L)
    return [int(x) for x in (a*b).coeffs()[:L]]
def direct(L):
    p=[1]+[0]*(L-1)
    for k in range(1,L):
        old=p.copy()
        for n in range(k,L):p[n]+=p[n-k]+old[n-k]
    return p
def main():
    p=direct(1001); count8=0; count16=0; support=0
    for ell in range(1,121):
        h=coefficients(ell,4000)
        if ell%4 in (0,1):
            need(all(h[n]%8==0 for n in range(3,4000,4)), 'mod8 sufficient family'); count8+=1000
        else:
            n=ell if ell%4==3 else ell+1
            need(h[n]%8==(2 if ell%4==3 else 4), 'mod8 exclusion witness')
        need(all(h[n]==sum(p[n-ell*j]*p[j] for j in range(n//ell+1))%16 for n in range(181)), 'direct product cross-check')
    targets=[n for n in range(1,20001) if odd(n)%8==7]
    params=sorted({d*d for d in range(1,21)}|{2*d*d for d in range(1,21)}|{ell for ell in range(1,121) if odd(ell)%8==1})
    for ell in params:
        h=coefficients(ell,20001)
        need(all(h[n]==0 for n in targets), 'dyadic family'); count16+=len(targets)
    represented={x*x+c*y*y for c in (1,2) for x in range(32) for y in range(32) if x*x+c*y*y<=1000}
    for ell in sorted({d*d for d in range(1,13)}|{2*d*d for d in range(1,13)}):
        h=coefficients(ell,1001)
        for n in range(1,1001):
            if n not in represented:need(h[n]==0,'nonrepresentation theorem');support+=1
    value=lambda ell,n:sum(p[n-ell*j]*p[j] for j in range(n//ell+1))
    expected={(8,28):86368,(6,17):3648,(6,23):25024,(9,69):594986112,(4,3):8,(9,15):1552,(17,21):9932}
    for pair,v in expected.items():need(value(*pair)==v,'sharpness/boundary witness')
    need(21 not in represented and value(17,21)%16==12,'invalid enlargement must fail')
    print(json.dumps({'mod8_instances':count8,'dyadic_parameters':params,'dyadic_targets_each':len(targets),'dyadic_instances':count16,'support_instances':support,'status':'passed','scope':'finite corroboration, not universal proof'}))
if __name__=='__main__':main()
