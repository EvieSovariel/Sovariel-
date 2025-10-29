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
    if n in {2,3}: return True
    s, d = 0, n-1
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in witnesses:
        if a >= n: continue
        x = fast_pow(a, d, n)
        if x == 1 or x == n-1: continue
        for r in range(1, s):
            x = (x * x) % n
            if x == n-1: break
        else:
            return False
    return True

def lattice_partition(N, depth):
    intervals = []
    step = N // (1 << depth)
    for i in range(1 << depth):
        start = i * step
        end = (i + 1) * step if i < (1 << depth) - 1 else N
        intervals.append((start, end))
    return intervals

def cri_pre_filter_v3(n, depth):
    if n < 3: return False
    b = bin(n)[2:].zfill(depth)
    ones = b.count('1')
    if ones in {0, depth}: H = 0
    else:
        p1 = ones / depth
        p0 = 1 - p1
        H = -(p1 * log2(p1) + p0 * log2(p0)) if p1 not in {0,1} else 0
    pairs = sum(b[i]==b[i+1] for i in range(0, depth-1, 2))
    align = pairs / (depth // 2)
    cri = 0.4 * align + 0.3 / (1 + abs(H - 1.0)) + 0.3
    return cri

def cri_v5(n, depth, local_gap_avg, global_avg_log):
    base = cri_pre_filter_v3(n, depth)
    boost = min(0.08, 0.05 * (local_gap_avg / global_avg_log))
    return base + boost

def build_cri_cache(start, end, depth, threshold=0.63):
    cache = {}
    step = max(1, (end - start) // 256)
    for s in range(start, end, step):
        n = s | 1
        if n < 3: continue
        cri = cri_pre_filter_v3(n, depth)
        cache[n:n+step] = (cri > threshold, cri)
    return cache

def sum_primes_below(N, depth=2048, base_threshold=0.63):
    if N < 2: return 0
    witnesses = [2,3,5,7,11,13,17,19,23]
    intervals = lattice_partition(N, depth)
    total = 2
    global_avg_log = math.log(N) if N > 1 else 1
    local_gaps = []
    last_prime = 2

    for start, end in intervals:
        local = 0
        cache = build_cri_cache(start, end, depth, base_threshold)
        n = max(3, (start + 1)//2*2 + 1)
        local_gap_sum = 0
        gap_count = 0
        while n < end:
            key = n - (n % 256)
            if key in cache and cache[key][0]:
                if miller_rabin_fast(n, witnesses):
                    local += n
                    if last_prime:
                        gap = n - last_prime
                        local_gap_sum += gap
                        gap_count += 1
                    last_prime = n
            n += 2
        if gap_count:
            local_gap_avg = local_gap_sum / gap_count
            local_gaps.append(local_gap_avg)
        total += local
    return total
