import math
from math import log2

def fast_pow(base, exp, mod):
    result = 1
    base %= mod
    while exp:
        if exp & 1: result = (result * base) % mod
        base = (base * base) % mod
        exp >>= 1
    return result

def miller_rabin_fast(n, witnesses):
    if n < 2: return False
    if n in {2, 3}: return True
    s, d = 0, n - 1
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in witnesses:
        if a >= n: continue
        x = fast_pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for r in range(1, s):
            x = (x * x) % mod
            if x == n - 1: break
        else:
            return False
    return True

def lattice_partition(N, depth):
    if N <= 2: return []
    intervals = []
    step = N // (1 << depth)
    if step == 0: step = 1
    for i in range(1 << depth):
        start = i * step
        end = (i + 1) * step if i < (1 << depth) - 1 else N
        intervals.append((start, end))
    return intervals

def cri_pre_filter(n, depth):
    if n < 3: return False
    b = bin(n)[2:].zfill(depth)
    ones = b.count('1')
    if ones in {0, depth}:
        H = 0.0
    else:
        p1 = ones / depth
        p0 = 1 - p1
        H = -(p1 * log2(p1) + p0 * log2(p0)) if 0 < p1 < 1 else 0.0
    pairs = sum(1 for i in range(0, depth - 1, 2) if b[i] == b[i + 1])
    align = pairs / (depth // 2)
    cri = 0.4 * align + 0.3 / (1 + abs(H - 1.0)) + 0.3
    return cri

def build_cri_cache(start, end, depth, threshold=0.63):
    cache = {}
    step = max(1, (end - start) // 256)
    for s in range(start, min(end, start + 256 * step), step):
        n = s | 1  # force odd
        if n < 3 or n >= end: continue
        cri_val = cri_pre_filter(n, depth)
        cache[n:n + step] = (cri_val > threshold, cri_val)
    return cache

def sum_primes_below(N, depth=16, threshold=0.63):
    if N < 2: return 0
    if N == 2: return 2
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    intervals = lattice_partition(N, depth)
    total = 0
    is_first_interval = True

    for start, end in intervals:
        local = 0
        cache = build_cri_cache(start, end, depth, threshold)

        # Start from first odd >= max(3, start)
        n = max(3, start + (start % 2 == 0))  # align to odd
        if is_first_interval and start == 0:
            local += 2
            n = 3
        is_first_interval = False

        while n < end:
            # Boundary safety: skip if n in next interval
            if n >= end: break
            key = n - (n % 256)
            if key in cache and n < end and cache[key][0]:
                if miller_rabin_fast(n, witnesses):
                    local += n
            n += 2
        total += local
    return total
