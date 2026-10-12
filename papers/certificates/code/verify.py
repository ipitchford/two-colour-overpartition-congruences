"""Rigorous-certificate verifier.
Input: l, M, r, j (claim pbar_l(Mn+r) ≡ 0 mod 2^j), N, eta vector of F (json {delta: exp}).
Checks (exact integer arithmetic): Gordon-Hughes-Newman/Ligozat hypotheses, weight, character data,
Sturm bound; then computes G = U_M F and checks all coefficients up to the Sturm bound are ≡ 0 mod 2^j.
Also cross-checks G = q^e * prod f_{d/M}^{b_d} * sum pbar_l(Mn+r) q^n numerically."""
import sys, json
from math import gcd
from fractions import Fraction as Fr
from flint import fmpz_poly

def require(condition, message):
    if not condition:
        raise ValueError(message)

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def primes(n):
    return [p for p in range(2, n + 1) if n % p == 0 and all((p % q for q in range(2, p)))]

def euler(L):
    c = [0] * L
    k = 0
    while True:
        done = True
        for kk in [k] if k == 0 else [k, -k]:
            e = kk * (3 * kk - 1) // 2
            if e < L:
                c[e] += -1 if kk % 2 else 1
                done = False
        if done and k > 0:
            break
        k += 1
    return c

def inv_series(c, L):
    nz = [(i, v) for i, v in enumerate(c[:L]) if v and i > 0]
    out = [0] * L
    out[0] = 1
    for n in range(1, L):
        s = 0
        for i, v in nz:
            if i > n:
                break
            s -= v * out[n - i]
        out[n] = s
    return out

def stretch(c, d, L):
    o = [0] * L
    for i, v in enumerate(c):
        if i * d < L:
            o[i * d] = v
        else:
            break
    return o

def etaprod(vec, L):
    f = euler(L)
    fi = inv_series(f, L)
    P = fmpz_poly([1])
    for d, e in vec.items():
        if e == 0:
            continue
        base = fmpz_poly(stretch(f if e > 0 else fi, d, L))
        P = P.mul_low(base.pow_trunc(abs(e), L), L)
    return P

def main(l, M, r, j, N, F):
    F = {int(k): v for k, v in F.items()}
    require(l > 0 and M > 0 and (N > 0) and (j > 0) and (0 <= r < M), 'invalid parameter domain')
    require(N % M == 0, 'U_M requires M dividing N')
    require(all((d > 0 and isinstance(e, int) for d, e in F.items())), 'invalid eta vector')
    A = {}
    for d, e in [(1, -2), (2, 1), (l, -2), (2 * l, 1)]:
        A[d] = A.get(d, 0) + e
    B = {d: F.get(d, 0) - A.get(d, 0) for d in F.keys() | A.keys()}
    B = {d: e for d, e in B.items() if e}
    require(all((d % M == 0 for d in B)), 'multiplier must be a function of q^M')
    require(all((N % d == 0 for d in F)), 'eta indices must divide level')
    tot = sum(F.values())
    require(tot % 2 == 0 and tot > 0, 'invalid weight')
    k = tot // 2
    s1 = sum((d * e for d, e in F.items()))
    s2 = sum((N // d * e for d, e in F.items()))
    require(s1 % 24 == 0 and s2 % 24 == 0, 'eta congruence conditions failed')
    c = s1 // 24
    orders = {cc: Fr(N, 24) * sum((Fr(gcd(cc, d) ** 2 * e, d) for d, e in F.items())) / (gcd(cc, N // cc) * cc) for cc in divisors(N)}
    require(min(orders.values()) >= 0, 'negative cusp order')
    require(c >= 0 and (c + r) % M == 0, 'invalid extraction shift')
    e0 = (c + r) // M
    sq = 1
    num = {}
    for d, e in F.items():
        for p in primes(d):
            v = 0
            t = d
            while t % p == 0:
                t //= p
                v += 1
            num[p] = num.get(p, 0) + v * e
    idx = N
    for p in primes(N):
        idx = idx * (p + 1) // p
    sturm = Fr(k * idx, 12)
    Bd = int(sturm)
    print(f'l={l} M={M} r={r} claim mod 2^{j}={2 ** j}')
    print(f'F eta vector {dict(sorted(F.items()))}; multiplier B={dict(sorted(B.items()))}')
    print(f'weight k={k}, level N={N}, sum d r_d={s1} (q-order c={c}), sum (N/d) r_d={s2}')
    print(f'character chi(d)=((-1)^k s / d), s = prod d^r_d = ' + ' * '.join((f'{p}^{v}' for p, v in sorted(num.items()))))
    print('orders at cusps c|N (Ligozat):', {cc: str(v) for cc, v in orders.items()})
    print(f'[SL2:Gamma0(N)]={idx}, Sturm bound k*idx/12={sturm}; G=U_M F = q^{e0}(...); check m=0..{Bd}')
    L = M * Bd - c + 1
    P = etaprod(F, L)
    G = [int(P[M * m - c]) if 0 <= M * m - c < L else 0 for m in range(Bd + 1)]
    bad = [m for m, g in enumerate(G) if g % 2 ** j]
    print('first G coeffs:', G[:8])
    print('all G coeffs m<=Sturm divisible by 2^j:', not bad, '| failures:', bad[:5])
    La = M * Bd + r + 1
    Aser = etaprod(A, La)
    pb = [int(Aser[M * n + r]) for n in range(Bd + 1) if M * n + r < La]
    print('direct pbar_l(Mn+r) n<=Sturm all ≡0 mod 2^j:', all((x % 2 ** j == 0 for x in pb)), 'count', len(pb), '| max 2-adic sharpness: min v2 =', min(((x & -x).bit_length() - 1 for x in pb if x)))
    Lg = Bd + 1 - e0
    Q = etaprod({d // M: e for d, e in B.items()}, Lg).mul_low(fmpz_poly(pb[:Lg]), Lg)
    print('cross-check G = q^e0 prod f_{d/M}^{b_d} * sum pbar(Mn+r)q^n:', [int(Q[i]) for i in range(Lg)] == G[e0:e0 + Lg] and all((g == 0 for g in G[:e0])))
    direct_ok = all((x % 2 ** j == 0 for x in pb))
    identity_ok = [int(Q[i]) for i in range(Lg)] == G[e0:e0 + Lg] and all((g == 0 for g in G[:e0]))
    require(not bad, 'Sturm coefficients failed')
    require(direct_ok, 'direct progression failed')
    require(identity_ok, 'extraction identity failed')
    return True
if __name__ == '__main__':
    a = sys.argv
    ok = main(int(a[1]), int(a[2]), int(a[3]), int(a[4]), int(a[5]), json.loads(a[6]))
    print('RESULT', 'CERTIFICATE PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)
