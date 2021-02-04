#!/usr/bin/env python
# coding: utf-8

import random
from time import process_time

def decompose(x, sboxes):
    t = len(sboxes)
    dec = [0]*t
    for k in range(t-1, -1, -1):
        dec[k] = x % sboxes[k]
        x = (x - dec[k]) // sboxes[k]
    return dec

def small_s_box(x, v):
    if x < v:
        return (x**(v-2)) % v
    return x

def compose(dec, sboxes):
    t = len(sboxes)
    x = 0
    for i in range(t):
        x *= sboxes[i]
        x += dec[i]
    return x

def bar_function(x, order, sboxes, v, s_box=small_s_box):
    x = decompose(x, sboxes)
    x = [s_box(y, v) for y in x]
    x = compose(x, sboxes)
    return x % order

def bar(state, order, sboxes, v, s_box=small_s_box):
    new_state = [bar_function(x, order, sboxes, v, s_box=s_box) for x in state]
    return new_state

def test_compose_decompose():
    for sboxes in [[10,10,10], [4,5,6], [11,3,17,53]]:
        for i in range(reduce(operator.mul, sboxes, 1)):
            assert i == compose(decompose(i, sboxes), sboxes)
test_compose_decompose()

def test_bar():
    assert bar([12, 123, 456, 789], 1009, [10, 10, 10], 7) == [14, 145, 236, 789]
    assert bar([36, 37, 38, 41, 64001, 18, 1618], 70657, [40, 40, 40], 37) == [36, 37, 38, 41, 1, 35, 1635]
# test_bar()

def miller_rabin(n, k):
    # If number is even, it's a composite number
    if n == 2: return True
    if n % 2 == 0: return False
    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def mult_minus_print_list(pair_list):
    min_v, max_box = float('inf'), -float('inf')
    out = None
    acc = 1
    length = len(pair_list)
    for j in range(0, length//2):
        mul = pair_list[length-1-2*j]
        dif = pair_list[length-1-2*j-1]
        out = f"{out}*{mul}" if out else f"{mul}"
        acc *= mul
        if dif: out = f"({out}-{dif})"
        acc -= dif
        max_box = max(max_box, mul)
        min_v = min(min_v, mul - dif)
    print(f"{out} = {acc}")
    print(f"max box size = {max_box} min v value = {min_v}")
    return acc

def find_decomposition_opt(number, v, max_d, pre_list):
    # print(f"n = {number} v = {v} max_d = {max_d}")
    if v <= number and number <= v + max_d:
        pre_list = pre_list + [0, number]
        # mult_minus_print_list(pre_list)
        return pre_list;
    if number < v:
        return None
    for k in range(max_d+1):
        if miller_rabin(number+k, 40): # primality test
            continue
        for q in range(v+max_d, v+k-1, -1):
            if q and q != 1 and (number+k) % q == 0: # q == 0: invalid modulus. q == 1: infinite recursion.
                new_n = (number+k)//q
                maybe_res = find_decomposition_opt(new_n, v, max_d, pre_list+[k,q])
                if maybe_res and len(maybe_res):
                    return maybe_res
    return None

def findDecBest(number, num_sboxes, max_d, desired_v):
    avg_sbox_size = math.floor(math.pow(number, 1/num_sboxes))
    min_v = max(avg_sbox_size - 33, 4) # magic numbers from original reference code
    max_v = max(avg_sbox_size - 25, 6)
    for v in range(min_v, max_v):
        # print(f"number: {number} num sboxes: {num_sboxes} approx size: {avg_sbox_size} v: {v} max_d: {max_d}")
        res = find_decomposition_opt(number, v, max_d, []);
        if res and len(res):
            res = res[1::2]
            min_dec = min(decompose(number, res))
            if min_dec > desired_v:
                v = previous_prime(min_dec + 1)
                return v, res
    return None

def decomposeField(fieldSize, maxd, num_sboxes=27, min_v=2):
    return findDecBest(fieldSize, num_sboxes, maxd, min_v)

def test_decomposeField():
    BLSr = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    BN254 = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    assert decomposeField(BLSr,42) == [32, 693, 35, 696, 20, 694, 4, 668, 4, 679, 12, 695, 3, 691, 9, 693, 30, 700, 3, 688, 27, 700, 27, 694, 20, 701, 31, 694, 12, 699, 32, 701, 39, 701, 22, 701, 2, 695, 11, 698, 7, 697, 42, 703, 11, 702, 3, 691, 4, 688, 28, 703, 0, 679]
    assert decomposeField(BN254,41) == [8, 651, 7, 658, 7, 656, 14, 666, 10, 663, 8, 654, 11, 668, 17, 677, 13, 681, 26, 683, 21, 669, 21, 681, 36, 680, 34, 677, 31, 675, 0, 668, 32, 675, 33, 683, 29, 681, 32, 683, 20, 683, 4, 655, 19, 680, 27, 683, 13, 667, 4, 678, 0, 673]
# test_decomposeField()

def generate_prime_list(base=10, lower_limit=5, upper_limit=16, min_v=20):
    prime_list = []
    for prime in [random_prime(base**(i+1), lbound=base**i) for i in range(lower_limit, upper_limit)]:
        max_d = 42
        num_sboxes = math.floor(math.log(prime, 10)) // 2
        num_sboxes = max(num_sboxes, 2)
        res = decomposeField(prime, max_d, num_sboxes, min_v)
        v, dec, symbol = "– ", [], "-"
        if res:
            v, dec, symbol = *res, "+"
            prime_list += [(prime, v, dec)]
        if get_verbose() >= 1:
            print(f"[{symbol}] p = {prime:>15} num_boxes = {len(dec)}({num_sboxes})", end=" ")
            print(f"v = {v:>3} s_boxes = {dec if dec else '– '}")
    return prime_list

def modulo_polynomial(ring, modulus):
    field_size = len(ring.base_ring())
    points = [(i, i % modulus) for i in range(field_size)]
    poly = ring.lagrange_polynomial(points)
    return poly

# returned polynomial vanishes for lower < var ≤ upper
def interval_polynomial(var, lower, upper):
    return reduce(operator.mul, [var - j for j in range(lower+1, upper + 1)], 1)

def test_interval_polynomial(p=17):
    field.<x> = GF(p)[]
    for u in range(p-1):
        for l in range(-1, u-1):
            poly = interval_polynomial(x, l, u)
            for j in range(l+1):
                assert poly(j)
            for j in range(l+1, u+1):
                assert not poly(j)
            for j in range(u+1, p):
                assert poly(j)
# test_interval_polynomial()

def invert_by_v_poly(ring, v):
    assert v < len(ring.base_ring())
    points = [(i, (i^(v-2)) % v) for i in range(v + 1)]
    poly = ring.lagrange_polynomial(points)
    return poly

def test_invert_by_v_poly(p=17):
    ring.<x> = GF(p)[]
    for v in primes(2,p):
        poly = invert_by_v_poly(ring, v)
        for j in [0, v]:
            assert int(poly(j))*int(j) % v == 0
        for j in range(1, v):
            assert int(poly(j))*int(j) % v == 1
# test_invert_by_v_poly()

def maybe_invert_by_v_poly(ring, v, sbox):
    points = [(i, (i^(v-2)) % v) for i in range(v)]
    points += [(i, i) for i in range(v, sbox)]
    poly = ring.lagrange_polynomial(points)
    return poly

def test_maybe_invert_by_v_poly(p=17):
    ring.<x> = GF(p)[]
    for v in primes(2,p):
        poly = maybe_invert_by_v_poly(ring, v, p)
        for j in [0, v]:
            assert int(poly(j))*int(j) % v == 0
        for j in range(1, v):
            assert int(poly(j))*int(j) % v == 1
        for j in range(v, p):
            assert int(poly(j)) == j
# test_maybe_invert_by_v_poly()

def decomposition_poly(variables, decomposition):
    assert len(variables) == len(decomposition) + 1
    poly = -variables[0]
    for i in range(1, len(decomposition)+1):
        poly += reduce(operator.mul, decomposition[i:], variables[i])
    return poly

def test_decomposition_poly():
    for sboxes in [[8,9,10], [11,5,6], [3,3,3], [693, 696]]:
        size = reduce(operator.mul, sboxes, 1)
        ring = PolynomialRing(GF(next_prime(size)), len(sboxes)+1, 'x', order='degrevlex')
        variables = ring.gens()
        poly = decomposition_poly(variables, sboxes)
        for i in range(size):
            dec = decompose(i, sboxes)
            assert not poly(i, *dec)
# test_decomposition_poly()

def bar_poly_system(order, decomposition, v, s_box=small_s_box):
    num_s_boxes = len(decomposition)
    num_vars = num_s_boxes*2 + 2
    ring = PolynomialRing(GF(order), 'x', num_vars)
    variables = ring.gens() # x, x_0, …, x_{n-1}, y, y_0, …, y_{n-1}
    x = variables[0]
    system =  [decomposition_poly(variables[:num_vars//2], decomposition)] # x and x_i correspond
    system += [decomposition_poly(variables[num_vars//2:], decomposition)] # y and y_i correspond
    s_min, s_max = min(sboxes), max(sboxes)
    q_i_min = interval_polynomial(x, -1, s_min - 1)
    q_i_max = interval_polynomial(x, -1, s_max - 1)
    uni_ring = GF(order)[x] # lagrange interpolation only works in univariate rings in sage
    s_box_points = [(i, s_box(i, v)) for i in range(s_max)]
    for i in range(num_s_boxes):
        x_i = variables[i + 1]
        y_i = variables[i + 1 + num_s_boxes + 1]
        s_i = decomposition[i]
        uni_ring = GF(order)[x_i]
        system += [q_i_min.subs({x:x_i}) * interval_polynomial(x_i, s_min-1, s_i-1)] # ensure x_i < s_i
        system += [ring(uni_ring.lagrange_polynomial(s_box_points[:s_i])) - y_i] # interpolate the S-Box
    # Additional constraints – you might need to set variable 'testing' to False down below!
    # system += [variables[0] - variables[num_vars//2]] # Bar(x) == x
    return system

def bar_pow_bar_poly_system(order, decomposition, v, exponent=2, s_box=small_s_box):
    num_s_boxes = len(decomposition)
    num_vars = num_s_boxes*4 + 4 # twice of bar_poly_system: left half for first bar, right half for second bar
    ring = PolynomialRing(GF(order), 'x', num_vars)
    var = ring.gens()
    shift_dict = {var[i] : var[i+num_vars//2] for i in range(num_vars//2)} # substitution shifts to right half of variables
    bar_sys = bar_poly_system(order, decomposition, v, s_box=s_box)
    system = [ring(poly) for poly in bar_sys]
    system += [var[num_vars//4]^exponent - var[num_vars//2]] # Output of first Bar to the exp is input to second Bar
    system += [ring(poly).subs(shift_dict) for poly in bar_sys]
    return system


from sage.rings.polynomial.toy_buchberger import spol
def is_groebner_basis(gb):
    for f in gb:
        for g in gb:
            try:
                redu = spol(f, g).reduce(gb)
            except NotImplementedError as nie:
                print(f"[!] No Gröbner Basis consistency check performed")
                print(f"[!] {nie}")
                return True
            if redu:
                print(f"f:              {f}")
                print(f"g:              {g}")
                print(f"s_poly:         {spol(f,g)}")
                print(f"reduced s_poly: {spol(f,g).reduce(gb)}")
                assert False
    return True

def random_s_box(field_size, degree = 5, terms = 15):
    ring.<x> = GF(field_size)[]
    f = ring.random_element(degree, terms)
    s_box = lambda x, v: int(f(x)) if x < v else x
    return s_box, f

def test_sboxes_too_tight_collision():
    for prime, v, sboxes in prime_dec_list:
        assert reduce(operator.mul, sboxes, 1) > prime
        assert compose([v]*len(sboxes), sboxes) < prime
        lower = compose([s - 1 for s in sboxes[:-1]] + [0], sboxes)
        #col_candidates = compose([s - 1 for s in sboxes[:-1]] + [], sboxes)
        print(f" ————————————")
        print(f"prime = {prime} v = {v} sboxes = {sboxes}")


        my_list = map(operator.sub, decompose(prime, sboxes), [s - 1 for s in sboxes])
        i = next((i for i, x in enumerate(my_list) if x), None)
        if i and decompose(prime, sboxes)[i] < v:
            print(f"[-] predicting collisions")
        if lower < prime:
            print(f"[!] danger zone: {lower}")
        for i in range(prime - floor(log(prime,10)), prime):
            permute = bar([i], prime, sboxes, v)[0]
            double_permute = bar([permute], prime, sboxes, v )[0]
            if i != double_permute:
                print(f"[!!] collision: bar({i:>9}) == bar({double_permute:>4}) == {permute:>4}", end=" | ")
                x = decompose(i, sboxes)
                print(f"{i:>5} = {x}", end=" ")
                x = [small_s_box(y, v) for y in x]
                print(f"→ {x}", end=" ")
                x = compose(x, sboxes)
                print(f"→ {x}")
        print(f" ————————————")
        print(f"")
# test_sboxes_too_tight_collision()


if __name__ == "__main__":
    set_verbose(3)
    testing = False
    time_it = True
    box_type = ['default', 'random', 'iden'][1]

    # Proposed by Dmitry 2021-02-04
    prime_dec_list = [
        # 2 small S-Boxes
        (1030297, 101, [101, 101, 101]), # has 3 boxes. only 2 desired?
        (12541, 107, [112, 112]),
        # 3 small S-Boxes
        (1295027, 109, [109, 109, 109]),
        (1191013, 101, [106, 106, 106]),
        # 4 small S-Boxes
        (2248087, 131, [131, 131, 131]), # has 3 boxes. 4 desired?
        (1003875853, 173, [178, 178, 178, 178]),
        # 5 small S-Boxes
        (51888844697, 139, [139, 139, 139, 139, 139]),
        (210906087421, 179, [184, 184, 184, 184, 184]),
    ]
    # Random generation of primes and their decompositions
    prime_dec_list = generate_prime_list(base=10, lower_limit=3, upper_limit=5, min_v=0)


    for prime, v, sboxes in prime_dec_list:
        print(f"————————————————————————————")
        print(f"p = {prime}, v = {v}, sboxes = {sboxes}")
        s_box, f = small_s_box, f"x^(v-2) % v"
        if box_type == 'random': s_box, f = random_s_box(v, degree=v, terms=2*v)
        elif box_type == 'iden': s_box, f = (lambda x, v: x, 1)
        if get_verbose() >= 2:
            print(f"f in sbox = {f}")
        time_sys_start = process_time()
        system = bar_pow_bar_poly_system(prime, sboxes, v, s_box=s_box)
        time_sys_stop = process_time()
        if get_verbose() >= 3:
            print(f"——————————————")
            [print(f"{poly}") for poly in system]
            print(f"——————————————")
        if testing:
            tmp = reduce(operator.mul, sboxes, 1)
            assert tmp >= prime, f"[!] S-Boxes too restrictive: {tmp} < {prime}"
            tmp = compose([v]*len(sboxes), sboxes)
            assert tmp < prime, f"[!] [v,…,v] is no field element (potential collisions): {tmp} >= {prime}"
            assert all([x >= v for x in decompose(prime, sboxes)])
            for inpu in [0, 1, prime-1] + [randint(1,prime-2) for _ in range(1000)]:
                outp = bar([inpu], prime, sboxes, v, s_box=s_box)[0]
                inpu_dec = decompose(inpu, sboxes)
                outp_dec = decompose(outp, sboxes)
                check = [poly(inpu, *inpu_dec, outp, *outp_dec) for poly in system]
                if any(check):
                    print(f"prime: {prime} v: {v} sboxes: {sboxes}")
                    print(f"input:  {inpu}")
                    x = decompose(inpu, sboxes)
                    print(f"  decomp: {x}")
                    x = [small_s_box(y, v) for y in x]
                    print(f"     inv: {x}")
                    x = compose(x, sboxes)
                    print(f"    comp: {x}")
                    print(f"output: {outp}")
                    print(f"  decomp:  {decompose(outp, sboxes)}")
                    print(f"polys:  {check}")
                    assert False
        if time_it:
            print(f"time system:    {float(n(time_sys_stop - time_sys_start, digits=5)):>8.5} sec")
        print(f"polys deg's:    {[poly.degree() for poly in system]}")
        print(f"macaulay bound: {1 + sum([poly.degree() - 1 for poly in system])}")
        print(f"estimated dreg: {floor(2*len(sboxes)*prime**(1/len(sboxes)))}")

        time_gb_start = process_time()
        gb = Ideal(system).groebner_basis()
        time_gb_stop = process_time()
        if testing:
            assert is_groebner_basis(gb)
        print(f"time gb:        {float(n(time_gb_stop - time_gb_start, digits=5)):>8.5} sec")
        gb = sorted(gb, key=lambda poly : poly.degree())
        print(f"actual dreg:    {gb[-1].degree()}")
        print(f"#elts in gb:    {len(gb)}")
        if get_verbose() >= 1:
            print(f"gb's deg's:     {[poly.degree() for poly in gb]}")
        if get_verbose() >= 3:
            print(f"——————————————")
            [print(f"{poly}") for poly in gb]
            print(f"——————————————")
        print(f"————————————————————————————")
        print(f"")
