#!/usr/bin/env python
# coding: utf-8

import random
from time import process_time

class PlookupHash():
    def __init__(self, order, constants, mult_matrix, sboxes, v, s_box_f=None):
        self.order = order
        self.constants = matrix(constants)
        self.matrix = matrix(mult_matrix)
        self.sboxes = sboxes
        self.v = v
        self.s_box_f = s_box_f
        if not s_box_f:
            self.s_box_f = lambda x : (x**(v-2)) % v

    def __call__(self, state):
        brick, concrete, bar = self.brick, self.concrete, self.bar
        state = concrete(state, 0)
        state = brick(state)
        state = concrete(state, 1)
        state = brick(state)
        state = concrete(state, 2)
        state = bar(state)
        state = concrete(state, 3)
        state = brick(state)
        state = concrete(state, 4)
        state = brick(state)
        state = concrete(state, 5)
        return state

    def _compose(self, dec):
        sboxes = self.sboxes
        t = len(sboxes)
        x = 0
        for i in range(t):
            x *= sboxes[i]
            x += dec[i]
        return x

    def _decompose(self, x):
        sboxes = self.sboxes
        t = len(sboxes)
        dec = [0]*t
        for k in range(t-1, -1, -1):
            dec[k] = x % sboxes[k]
            x = (x - dec[k]) // sboxes[k]
        return dec

    def _small_s_box(self, x):
        if x < self.v:
            return self.s_box_f(x)
        return x

    def _bar_function(self, x):
        order = self.order
        sboxes = self.sboxes
        compose = self._compose
        decompose = self._decompose
        s_box = self._small_s_box
        v = self.v
        x = decompose(x)
        x = [s_box(y) for y in x]
        x = compose(x)
        return x % order

    def bar(self, state):
        order = self.order
        sboxes = self.sboxes
        v = self.v
        new_state = [self._bar_function(x) for x in state]
        return new_state

    def brick(self, state):
        order = self.order
        x, y, z = state
        a = z**5 % order
        b = x*(z**2 + 1*z + 2) % order
        c = y*(x**2 + 3*x + 4) % order
        return [a, b, c]

    def concrete(self, state, cnst_idx):
        order = self.order
        matrix = self.matrix
        new_state = vector(self.constants[cnst_idx])
        new_state += matrix * vector(state)
        new_state = [s % order for s in new_state]
        return new_state

class TestPlookupHash():
    def __init__(self):
        self._test_compose_decompose()
        self._test_bar()
        if get_verbose() >= 2: print("Testing of PlookupHash completed.")

    def _test_compose_decompose(self):
        _ = None
        for sboxes in [[10,10,10], [4,5,6], [11,3,17,53]]:
            ph = PlookupHash(_, _, _, sboxes, _)
            for i in range(reduce(operator.mul, sboxes, 1)):
                assert i == ph._compose(ph._decompose(i))

    def _test_bar(self):
        _ = None
        ph = PlookupHash(1009, _, _, [10, 10, 10], 7)
        assert ph.bar([12, 123, 456, 789]) == [14, 145, 236, 789]
        ph = PlookupHash(70657, _, _, [40, 40, 40], 37)
        assert ph.bar([36, 37, 38, 41, 64001, 18, 1618]) == [36, 37, 38, 41, 1, 35, 1635]

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

def bar_poly_system(order, decomposition, s_box):
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
    s_box_points = [(i, s_box(i)) for i in range(s_max)]
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

def bar_pow_bar_poly_system(order, decomposition, s_box, exponent=5):
    num_s_boxes = len(decomposition)
    num_vars = num_s_boxes*4 + 4 # twice of bar_poly_system: left half for first bar, right half for second bar
    ring = PolynomialRing(GF(order), 'x', num_vars)
    var = ring.gens()
    shift_dict = {var[i] : var[i+num_vars//2] for i in range(num_vars//2)} # substitution shifts to right half of variables
    bar_sys = bar_poly_system(order, decomposition, s_box)
    system = [ring(poly) for poly in bar_sys]
    system += [var[num_vars//4]^exponent - var[num_vars//2]] # Output of first Bar to the exp is input to second Bar
    system += [ring(poly).subs(shift_dict) for poly in bar_sys]
    return system

def lbar_pow_rbar_poly_system(order, decomposition, s_box, exponent=5):
    num_s_boxes = len(decomposition)
    num_l_sboxes = num_s_boxes//2
    num_r_sboxes = num_s_boxes - num_l_sboxes
    num_l_vars = (num_l_sboxes + 2) * 2
    num_r_vars = (num_r_sboxes + 2) * 2
    num_vars = num_l_vars + num_r_vars
    ring = PolynomialRing(GF(order), 'x', num_vars)
    var = ring.gens()
    shift_dict = {var[i] : var[i+num_r_vars] for i in range(num_r_vars)} # substitution shifts to right “half” of variables
    print(shift_dict)
    r_big_box = reduce(operator.mul, decomposition[num_l_sboxes:], 1)
    l_big_box = reduce(operator.mul, decomposition[:num_l_sboxes], 1)
    up_sys = bar_poly_system(order, decomposition[:num_l_sboxes] + [r_big_box], s_box)
    lo_sys = bar_poly_system(order, [l_big_box] + decomposition[num_l_sboxes:], s_box)
    system = [ring(poly) for poly in up_sys]
    system += [var[num_l_vars//2]^exponent - var[num_l_vars]] # Output of first Bar to the exp is input to second Bar
    system += [ring(poly).subs(shift_dict) for poly in lo_sys]
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

def random_s_box_f(field_size, degree=5, terms=15):
    ring.<x> = GF(field_size)[]
    f = ring.random_element(degree, terms)
    s_box_f = lambda x: int(f(x))
    return s_box_f, f

if __name__ == "__main__":
    set_verbose(2)
    testing = False
    time_it = True
    box_type = ['default', 'random', 'iden'][0]

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
    # prime_dec_list = generate_prime_list(base=10, lower_limit=2, upper_limit=5, min_v=0)

    # parameters that kill Ferdinand's machine when using bar_pow_bar
    prime_dec_list = [ (5701, 53, [84, 68]) ]

    constants = [
        [3**100, 2**100, 5**50],
        [3**110, 2**110, 5**60],
        [3**120, 2**120, 5**70],
        [3**130, 2**130, 5**80],
        [3**140, 2**140, 5**90],
        [3**150, 2**150, 5**100],
    ]
    mult_matrix = [
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ]

    for prime, v, sboxes in prime_dec_list:
        print(f"————————————————————————————")
        print(f"p = {prime}, v = {v}, sboxes = {sboxes}")
        s_box_f, f = None, f"x^(v-2) % v"
        if box_type == 'random': s_box_f, f = random_s_box_f(v, degree=v, terms=2*v)
        elif box_type == 'iden': s_box_f, f = (lambda x: x, 'ID')
        if get_verbose() >= 2:
            print(f"f in sbox = {f}")
        ph = PlookupHash(prime, constants, mult_matrix, sboxes, v, s_box_f=s_box_f)
        time_sys_start = process_time()
        system = bar_poly_system(prime, sboxes, ph._small_s_box)
        time_sys_stop = process_time()
        if get_verbose() >= 3:
            print(f"——————————————")
            [print(f"{poly}") for poly in system]
            print(f"——————————————")
        if testing:
            TestPlookupHash()
            tmp = reduce(operator.mul, sboxes, 1)
            assert tmp >= prime, f"[!] S-Boxes too restrictive: {tmp} < {prime}"
            tmp = ph._compose([v]*len(sboxes))
            assert tmp < prime, f"[!] [v,…,v] is no field element (potential collisions): {tmp} >= {prime}"
            assert all([x >= v for x in ph._decompose(prime)])
            for inpu in [0, 1, prime-1] + [randint(1,prime-2) for _ in range(1000)]:
                outp = ph.bar([inpu])[0]
                inpu_dec = ph._decompose(inpu)
                outp_dec = ph._decompose(outp)
                check = [poly(inpu, *inpu_dec, outp, *outp_dec) for poly in system]
                if any(check):
                    print(f"prime: {prime} v: {v} sboxes: {sboxes}")
                    print(f"input:  {inpu}")
                    x = ph._decompose(inpu)
                    print(f"  decomp: {x}")
                    x = [ph._small_s_box(y) for y in x]
                    print(f"     inv: {x}")
                    x = ph._compose(x)
                    print(f"    comp: {x}")
                    print(f"output: {outp}")
                    print(f"  decomp:  {ph._decompose(outp)}")
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
