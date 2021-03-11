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
        if get_verbose() >= 2: print("Testing of PlookupHash complete.")

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
    _ = None
    avg_sbox_size = math.floor(math.pow(number, 1/num_sboxes))
    min_v = max(avg_sbox_size - 33, 4) # magic numbers from original reference code
    max_v = max(avg_sbox_size - 25, 6)
    for v in range(min_v, max_v):
        # print(f"number: {number} num sboxes: {num_sboxes} approx size: {avg_sbox_size} v: {v} max_d: {max_d}")
        res = find_decomposition_opt(number, v, max_d, []);
        if res and len(res):
            res = res[1::2]
            ph = PlookupHash(_, _, _, res, _)
            min_dec = min(ph._decompose(number))
            if min_dec > desired_v:
                v = previous_prime(min_dec + 1)
                return v, res
    return None

def decomposeField(fieldSize, maxd, num_sboxes=27, min_v=2):
    return findDecBest(fieldSize, num_sboxes, maxd, min_v)

def test_decomposeField():
    BLSr = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    BN254 = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    assert decomposeField(BLSr,42)[1] == [693, 696, 694, 668, 679, 695, 691, 693, 700, 688, 700, 694, 701, 694, 699, 701, 701, 701, 695, 698, 697, 703, 702, 691, 688, 703, 679]
    assert decomposeField(BN254,41)[1] == [651, 658, 656, 666, 663, 654, 668, 677, 681, 683, 669, 681, 680, 677, 675, 668, 675, 683, 681, 683, 683, 655, 680, 683, 667, 678, 673]
    if get_verbose() >= 2: print(f"Testing of field decomposition complete.")

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
    if get_verbose() >= 2: print(f"Testing of interval polynomial complete.")

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
    if get_verbose() >= 2: print(f"Testing of invert-by-v polynomial complete.")

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
    if get_verbose() >= 2: print(f"Testing of maybe-invert-by-v polynomial complete.")

def decomposition_poly(variables, decomposition):
    assert len(variables) == len(decomposition) + 1
    poly = -variables[0]
    for i in range(1, len(decomposition)+1):
        poly += reduce(operator.mul, decomposition[i:], variables[i])
    return poly

def test_decomposition_poly():
    _ = None
    for sboxes in [[8,9,10], [11,5,6], [3,3,3], [693, 696]]:
        size = reduce(operator.mul, sboxes, 1)
        ring = PolynomialRing(GF(next_prime(size)), len(sboxes)+1, 'x', order='degrevlex')
        variables = ring.gens()
        poly = decomposition_poly(variables, sboxes)
        ph = PlookupHash(_, _, _, sboxes, _)
        for i in range(size):
            dec = ph._decompose(i)
            assert not poly(i, *dec), f"Buggy decomposition polynomial."
    if get_verbose() >= 2: print(f"Testing of decomposition polynomial complete.")

def bar_poly_system(order, decomposition, s_box):
    num_s_boxes = len(decomposition)
    num_vars = num_s_boxes*2 + 2
    ring = PolynomialRing(GF(order), 'x', num_vars)
    variables = ring.gens() # x, x_0, …, x_{n-1}, y, y_0, …, y_{n-1}
    x = variables[0]
    system =  [decomposition_poly(variables[:num_vars//2], decomposition)] # x and x_i correspond
    system += [decomposition_poly(variables[num_vars//2:], decomposition)] # y and y_i correspond
    s_min, s_max = min(decomposition), max(decomposition)
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

def test_bar_poly_system(prime=5701, sboxes=[84, 68], v=53):
    _ = None
    ph = PlookupHash(prime, _, _, sboxes, v)
    system = bar_poly_system(prime, sboxes, ph._small_s_box)
    for inpu in [0, 1, prime-1] + [randint(1,prime-2) for _ in range(100)]:
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
            assert False, f"Polynomials of Bar's polynomial system do not correspond to PlookupHash."
    if get_verbose() >= 2: print(f"Testing of Bar's poly system complete.")

def conc_poly_system(order, constants, mult_matrix):
    assert len(constants) == mult_matrix.nrows(), f"Dimensions of constants and matrix mismatch: {len(constants)} vs {mult_matrix.nrows()}"
    state_size = len(constants)
    ring = PolynomialRing(GF(order), 'x', 2*state_size)
    var = ring.gens() # in_0, …, in_{s-1}, out_0, …, out_{s-1}
    invars = vector(var[:state_size])
    polys = [ constants[i] + mult_matrix[i]*invars - var[i+state_size] for i in range(state_size)]
    return polys

def test_conc_poly_system(prime=5701,
                          constants=[[3**100, 2**100, 5**50], [3**110, 2**110, 5**60]],
                          mult_matrix=matrix([[2, 1, 1], [1, 2, 1], [1, 1, 2]])):
    assert all([len(constants[0]) == len(constants[i]) for i in range(len(constants))]), f"All constants' vectors need to be of the same length."
    assert len(constants[0]) == mult_matrix.nrows(), f"Dimensions of constants and matrix mismatch: {len(constants[0])} vs {mult_matrix.nrows()}"
    _ = None
    constants = [[c % prime for c in cnts] for cnts in constants]
    state_size = len(constants[0])
    ph = PlookupHash(prime, constants, mult_matrix, _, _)
    for i in range(len(constants)):
        polys = conc_poly_system(prime, constants[i], mult_matrix)
        for _ in range(100):
            vals = [randint(0, prime) for _ in range(state_size)]
            vals += ph.concrete(vals, i)
            assert not any([p(vals) for p in polys])
    if get_verbose() >= 2: print(f"Testing of Concrete's poly system complete.")

def brick_poly_system(order, state_size=3):
    assert state_size == 3, f"It's not obvious how to generalize Brick from the specification of PlookupHash."
    ring = PolynomialRing(GF(order), 'x', 2*state_size)
    x = ring.gens() # in_0, in_1, in_2, out_0, out_1, out_2
    polys = [x[2]**5 - x[3],
             x[0]*x[2]**2 + x[0]*x[2] + 2*x[0] - x[4],
             x[0]**2*x[1] + 3*x[0]*x[1] + 4*x[1] - x[5]]
    return polys

def test_brick_poly_system(prime=5701):
    _ = None
    ph = PlookupHash(prime, _, _, _, _)
    polys = brick_poly_system(prime)
    for _ in range(100):
        vals = [randint(0, prime) for _ in range(3)]
        vals += ph.brick(vals)
        assert not any([p(vals) for p in polys])
    if get_verbose() >= 2: print(f"Testing of Brick's poly system complete.")

def conc_bar_conc_poly_system(order, constants, mult_matrix, decomposition, s_box):
    assert len(constants) >= 2, f"Multiple 'concrete' require multiple lists of constants"
    assert len(constants[0]) == len(constants[1]), f"The lists of constants have to have the same length"
    assert len(constants[0]) == mult_matrix.nrows(), f"Dimensions of constants and matrix mismatch: {len(constants[0])} vs {mult_matrix.nrows()}"
    num_s_boxes = len(decomposition)
    state_size = len(constants[0])
    num_vars = 2*state_size + state_size*(num_s_boxes*2 + 1) + state_size
    ring = PolynomialRing(GF(order), 'x', num_vars)
    var = ring.gens()
    all_shift_dict_bar = []
    for s in range(state_size):
        shift_dict_bar = {var[0] : var[state_size + s]} # output of concrete is input of bar
        shift_dict_bar.update( {var[i] : var[i + 2*state_size + s*(num_s_boxes*2 + 1) - 1] for i in range(1, num_s_boxes*2 + 2)} ) # take next free variables
        all_shift_dict_bar += [shift_dict_bar]
    shift_dict_conc = {var[i] : var[2*state_size + i*(num_s_boxes*2 + 1) + num_s_boxes] for i in range(state_size)} # collect the y's from the bars
    shift_dict_conc.update( {var[i] : var[i + 2*state_size + state_size*(num_s_boxes*2 + 1) - state_size] for i in range(state_size, 2*state_size)} )
    conc_sys_0 = conc_poly_system(order, constants[0], mult_matrix)
    bar_sys = bar_poly_system(order, decomposition, s_box)
    conc_sys_1 = conc_poly_system(order, constants[1], mult_matrix)
    system = [ring(p) for p in conc_sys_0]
    for shift_dict_bar in all_shift_dict_bar:
        system += [ring(p).subs(shift_dict_bar) for p in bar_sys]
    system += [ring(p).subs(shift_dict_conc) for p in conc_sys_1]
    return system

def test_conc_bar_conc_poly_system(prime=5701,
                                   constants=[[3**100, 2**100, 5**50], [3**110, 2**110, 5**60]],
                                   mult_matrix=matrix([[2, 1, 1], [1, 2, 1], [1, 1, 2]]),
                                   s_boxes=[84, 68],
                                   v=53):
    constants = [[c % prime for c in cnts] for cnts in constants]
    state_size = len(constants[0])
    ph = PlookupHash(prime, constants, mult_matrix, s_boxes, v)
    polys = conc_bar_conc_poly_system(prime, constants, mult_matrix, s_boxes, ph._small_s_box)
    for _ in range(100):
        state_0 = [randint(0, prime) for _ in range(state_size)]
        state_1 = ph.concrete(state_0, 0)
        state_2 = ph.bar(state_1)
        state_3 = ph.concrete(state_2, 1)
        vals = state_0 + state_1
        for i in range(state_size): # all the bars
            vals += ph._decompose(state_1[i])
            vals += [state_2[i]]
            vals += ph._decompose(state_2[i])
        vals += state_3
        assert not any([p(vals) for p in polys])
    if get_verbose() >= 2: print(f"Testing of Conc-Bar-Conc's poly system complete.")


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
                assert False, f"Some S-Polynomial did not reduce to 0: see output above."
    return True

def random_s_box_f(field_size, degree=5, terms=15):
    ring.<x> = GF(field_size)[]
    f = ring.random_element(degree, terms)
    s_box_f = lambda x: int(f(x))
    return s_box_f, f

if __name__ == "__main__":
    set_verbose(2)
    testing = 2
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

    if testing >= 2:
        TestPlookupHash()
        test_decomposeField()
        test_interval_polynomial()
        test_maybe_invert_by_v_poly()
        test_decomposition_poly()
        test_bar_poly_system()
        test_conc_poly_system()
        test_brick_poly_system()
        test_conc_bar_conc_poly_system()

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
            # Do sboxes and prime correspond?
            tmp = reduce(operator.mul, sboxes, 1)
            assert tmp >= prime, f"[!] S-Boxes too restrictive: {tmp} < {prime}"
            tmp = ph._compose([v]*len(sboxes))
            assert tmp < prime, f"[!] [v,…,v] is no field element (potential collisions): {tmp} >= {prime}"
            assert all([x >= v for x in ph._decompose(prime)])
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