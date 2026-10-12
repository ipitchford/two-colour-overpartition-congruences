"""Independent exact-integer audit of four BSC certificates.

Only the Python standard library is used. No source code is imported from
the submitted bundle. Coefficients of A_l use the inverse-theta recurrence;
coefficients of the small multiplier g use its logarithmic derivative.
The finite calculation verifies the hypotheses needed by the separately
inspected modular-form proof; it is not by itself an infinite proof.
"""
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json
import platform
import sys
import time
ROOT = Path(__file__).resolve().parent
CASES = [dict(name='A', ell=8, M=32, r=28, j=5, N=128, B={32: 69, 64: -29, 128: -2}, k=18, c=4, S=288), dict(name='B17', ell=6, M=24, r=17, j=6, N=288, B={24: 61, 48: -31, 72: -16, 96: 2, 144: 8}, k=11, c=7, S=528), dict(name='B23', ell=6, M=24, r=23, j=6, N=288, B={24: 59, 48: -25, 72: -16, 96: -2, 144: 8}, k=11, c=1, S=528), dict(name='C', ell=9, M=36, r=33, j=7, N=216, B={36: 88, 72: -43, 108: -26, 216: 13}, k=15, c=3, S=540)]

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def factors(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def overpartition_series(nmax):
    """1/theta(-q): P[n] = -2 sum (-1)^k P[n-k^2]."""
    p = [1] + [0] * nmax
    squares = [(k * k, 2 if k % 2 else -2) for k in range(1, isqrt(nmax) + 1)]
    for n in range(1, nmax + 1):
        value = 0
        for sq, coefficient in squares:
            if sq > n:
                break
            value += coefficient * p[n - sq]
        p[n] = value
    return p

def pair_series(p, ell, nmax):
    """A_l=P(q)/theta(-q^ell), by sparse series inversion."""
    h = [0] * (nmax + 1)
    squares = [(ell * k * k, 2 if k % 2 else -2) for k in range(1, isqrt(nmax // ell) + 1)]
    for n in range(nmax + 1):
        value = p[n]
        for sq, coefficient in squares:
            if sq > n:
                break
            value += coefficient * h[n - sq]
        h[n] = value
    return h

def elementary_overpartition_series(nmax):
    """Small sanity check directly from prod (1+q^m)/(1-q^m)."""
    a = [1] + [0] * nmax
    for m in range(1, nmax + 1):
        for n in range(nmax, m - 1, -1):
            a[n] += a[n - m]
        for n in range(m, nmax + 1):
            a[n] += a[n - m]
    return a

def multiplier_series(b, nmax):
    """g=prod f_d^b_d via n*g[n]=sum L[m]*g[n-m].

    L[m] is the q^m coefficient of q*g'/g, namely
    -sum_{d|m} b_d*d*sigma_1(m/d). Division is exact in Z.
    """
    sigma = [0] * (nmax + 1)
    for divisor in range(1, nmax + 1):
        for multiple in range(divisor, nmax + 1, divisor):
            sigma[multiple] += divisor
    log_derivative = [0] * (nmax + 1)
    for d, exponent in b.items():
        for m in range(d, nmax + 1, d):
            log_derivative[m] -= exponent * d * sigma[m // d]
    g = [1] + [0] * nmax
    for n in range(1, nmax + 1):
        numerator = sum((log_derivative[m] * g[n - m] for m in range(1, n + 1)))
        quotient, remainder = divmod(numerator, n)
        if not remainder == 0:
            raise ValueError((n, remainder))
        g[n] = quotient
    return g

def audit_case(case, p):
    start = time.perf_counter()
    ell, M, r, j, N = [case[key] for key in ('ell', 'M', 'r', 'j', 'N')]
    b = case['B']
    f = defaultdict(int)
    for d, e in [(1, -2), (2, 1), (ell, -2), (2 * ell, 1), *b.items()]:
        f[d] += e
    f = dict(sorted(((d, e) for d, e in f.items() if e)))
    if not all((N % d == 0 for d in f)):
        raise ValueError('verification condition failed')
    if not all((d % M == 0 and N % d == 0 for d in b)):
        raise ValueError('verification condition failed')
    if not N % M == 0:
        raise ValueError('verification condition failed')
    k = Fraction(sum(f.values()), 2)
    s1 = sum((d * e for d, e in f.items()))
    s2 = sum((N // d * e for d, e in f.items()))
    if not (k.denominator == 1 and k > 0):
        raise ValueError('verification condition failed')
    if not (s1 % 24 == 0 and s2 % 24 == 0):
        raise ValueError('verification condition failed')
    c = s1 // 24
    if not (c >= 0 and (c + r) % M == 0):
        raise ValueError('verification condition failed')
    shift = (c + r) // M
    orders = {den: Fraction(N, 24 * gcd(den, N // den) * den) * sum((Fraction(gcd(den, d) ** 2 * e, d) for d, e in f.items())) for den in divisors(N)}
    if not min(orders.values()) >= 0:
        raise ValueError('verification condition failed')
    s_factors = defaultdict(int)
    for d, e in f.items():
        for prime, power in factors(d).items():
            s_factors[prime] += power * e
    squarefree_part = 1
    for prime, power in s_factors.items():
        if power % 2:
            squarefree_part *= prime
    discriminant = (-1 if int(k) % 2 else 1) * squarefree_part
    if not discriminant in (1, -3):
        raise ValueError('verification condition failed')
    conductor = abs(discriminant)
    if not N % conductor == 0:
        raise ValueError('verification condition failed')
    index = N
    for prime in factors(N):
        index = index // prime * (prime + 1)
    bound = k * index / 12
    S = bound.numerator // bound.denominator
    if not (int(k), c, S) == (case['k'], case['c'], case['S']):
        raise ValueError('verification condition failed')
    max_required = M * S - c
    max_direct = M * S + r
    h = pair_series(p, ell, max_direct)
    progression = [h[M * n + r] for n in range(S + 1)]
    modulus = 2 ** j
    if not all((value % modulus == 0 for value in progression)):
        raise ValueError('verification condition failed')
    g = multiplier_series({d // M: e for d, e in b.items()}, S - shift)
    G = [0] * shift + [sum((g[i] * progression[n - i] for i in range(n + 1))) for n in range(S - shift + 1)]
    if not len(G) == S + 1:
        raise ValueError('verification condition failed')
    if not all((value % modulus == 0 for value in G)):
        raise ValueError('verification condition failed')
    valuation = min(((v & -v).bit_length() - 1 for v in progression if v))
    witnesses = [{'progression_index': n, 'argument': M * n + r, 'value': value, 'quotient_by_modulus': value // modulus} for n, value in enumerate(progression) if value % (2 * modulus) != 0]
    result = dict(name=case['name'], ell=ell, M=M, r=r, modulus=modulus, F=f, B=b, k=int(k), sum_dr=s1, sum_N_over_d_r=s2, c=c, shift=shift, cusp_orders={d: str(v) for d, v in orders.items()}, character_s_prime_powers=dict(s_factors), character_discriminant=discriminant, character_conductor=conductor, index=index, sturm_bound=str(bound), S=S, G_count=len(G), progression_required_count=S - shift + 1, progression_actual_checked_count=len(progression), max_required_A_argument=max_required, max_direct_argument=max_direct, all_G_divisible=True, all_progression_divisible=True, minimum_progression_v2=valuation, first_sharpness_witness=witnesses[0], first_eight_G=G[:8], elapsed_seconds=time.perf_counter() - start)
    payload = {'summary': result, 'G': G, 'progression': progression, 'g': g}
    filename = ROOT / (case['name'] + '_exact_coefficients.json')
    filename.write_text(json.dumps(payload, indent=2) + '\n')
    result['coefficient_artifact'] = filename.name
    result['coefficient_artifact_sha256'] = sha256(filename.read_bytes()).hexdigest()
    print(json.dumps(result), flush=True)
    return result

def main():
    start = time.perf_counter()
    started = datetime.now(timezone.utc).isoformat()
    nmax = max((c['M'] * c['S'] + c['r'] for c in CASES))
    p = overpartition_series(nmax)
    if not p[:181] == elementary_overpartition_series(180):
        raise ValueError('verification condition failed')
    results = [audit_case(case, p) for case in CASES]
    receipt = dict(status='passed', started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(), wall_seconds=time.perf_counter() - start, python=sys.version, executable=sys.executable, platform=platform.platform(), dependencies='Python standard library only', script=Path(__file__).name, script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), exact_overpartition_max_argument=nmax, direct_product_sanity_check_through=180, cases=results)
    (ROOT / 'independent_bsc_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('ALL FOUR CERTIFICATES PASSED; total wall seconds', receipt['wall_seconds'], flush=True)
if __name__ == '__main__':
    main()
