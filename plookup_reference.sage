from random import randrange

primeBN = 21888242871839275222246405745257275088548364400416034343698204186575808495617
decompositionBN = [693, 696, 694, 668, 679, 695, 691, 693, 700, 688, 700, 694, 701, 694, 699, 701, 701, 701, 695, 698, 697, 703, 702, 691, 688, 703, 679]
vBN = 661 # 661
matrixBN = [
    [2, 1, 1],
    [1, 2, 1],
    [1, 1, 2]
]
constantsBN = [
    [3**100, 2**100, 5**50],
    [3**110, 2**110, 5**60],
    [3**120, 2**120, 5**70],
    [3**130, 2**130, 5**80],
    [3**140, 2**140, 5**90],
    [3**150, 2**150, 5**100]
]

# inverse if smaller than v
# def smallSbox(x, v):
#     if x < v:
#         return x**(v-2) % v
#     return x

# Lookup-table-based Sbox
# def Bar(state, order, sbox_sizes, v):
#     t = len(sbox_sizes)
#     new_state=[0]*len(state)
#     nibbles = [0]*t #temporary nibble array
#     for i in range(len(state)):
#         tmp_val = state[i]
#         #1. Decomposition
#         for k in range(t): #nibble loop
#             nibbles[k] = tmp_val % sbox_sizes[k]
#             tmp_val = (tmp_val - nibbles[k]) // sbox_sizes[k] #reduce state value
#         #2. Sbox
#             nibbles[k] = smallSbox(nibbles[k], v)

#         #3. Composition
#         for k in range(t-1, -1, -1):
#             new_state[i] *= sbox_sizes[k]
#             new_state[i] += nibbles[k]
#         new_state[i] = new_state[i] % order
#     return new_state

def decompose(x, sboxes):
    t = len(sboxes)
    dec = [0]*t
    for k in range(t-1, -1, -1):
        dec[k] = x % sboxes[k]
        x = (x - dec[k]) // sboxes[k]
    return dec

def small_s_box(x, ν):
    if x < ν:
        return (x**(ν-2)) % ν
    return x

def compose(dec, sboxes):
    t = len(sboxes)
    x = 0
    for i in range(t):
        x *= sboxes[i]
        x += dec[i]
    return x

def bar_function(x, order, sboxes, ν):
    x = decompose(x, sboxes)
    x = [small_s_box(y, ν) for y in x]
    x = compose(x, sboxes)
    return x % order

def bar(state, order, sboxes, ν):
    new_state = [bar_function(x, order, sboxes, ν) for x in state]
    return new_state

def Brick(state, order):
    x, y, z = state
    c = z**5 % order
    b = x*(z**2 + 1*z + 2) % order
    a = y*(x**2 + 3*x + 4) % order
    return [a, b, c]

# Apply affine transformation to state
def Concrete(state, order, matrix, constants):
    new_state = constants   # initialize with constant vector
    for i in range(len(state)):  # matrix multiplication
        for j in range(len(state)):
            new_state[i] += matrix[i][j]*state[j]
    return [x % order for x in new_state]

def fullHash(state, decomposition, v, matrix, constants, order):
    new_state = state
    new_state = Concrete(new_state, order, matrix, constants[0])
    new_state = Brick(new_state, order)
    new_state = Concrete(new_state, order, matrix, constants[1])
    new_state = Brick(new_state, order)
    new_state = Concrete(new_state, order, matrix, constants[2])
    new_state = Bar(new_state, order, decomposition, v)
    new_state = Concrete(new_state, order, matrix, constants[3])
    new_state = Brick(new_state, order)
    new_state = Concrete(new_state, order, matrix, constants[4])
    new_state = Brick(new_state, order)
    new_state = Concrete(new_state, order, matrix, constants[5])
    return new_state

def test(decomposition, v, matrix, constants, order):
    #1. -1 -> -1
    input_1 = [order - 1, order - 2, order - 3]
    output_1 = Bar(input_1, order, decomposition, v)
    if input_1 != output_1:
        print("Test 1 not passed")
        print("input:")
        print(f"\t{input_1[0]}\n\t{input_1[1]}\n\t{input_1[2]}")
        print("output:")
        print(f"\t{output_1[0]}\n\t{output_1[1]}\n\t{output_1[2]}")

    #3. 0 -> 0
    input_2 = [0, 0, 0]
    output_2 = Bar(input_2, order, decomposition, v)
    if input_2 != output_2:
        print("Test 2 not passed")

    #4 1 ->1
    input_3 = [1, 1, 1]
    output_3 = Bar(input_3, order, decomposition, v)
    if input_3 != output_3:
        print("Test 3 not passed")

    #5. 5->167
    input_4 = [v//2]*3
    output_4 = Bar(input_4, order, decomposition, v)
    if output_4 != [(v//2)**(v-2) % v]*3:
        print("Test 4 not passed")

def approximate_identity_fraction_of_bar(decomposition, v, order):
    # print([s_i - v for s_i in decomposition]) # if these values are all small, things should be good.
    expected_ratio = [1 - ((v-2)/s_i) for s_i in decomposition] # S-Box is identity for inputs 0, 1, and anything greater than s_i
    expected_ratio = reduce(operator.mul, expected_ratio, 1)
    expected_ratio = expected_ratio.n(50)
    batch_size = 1000
    cnt_ttl = 0
    cnt_id = 0
    while True:
        input_value = [randrange(0, order) for _ in range(batch_size)]
        output_value = Bar(input_value, order, decomposition, v)
        cnt_ttl += batch_size
        cnt_id += sum([input_value[i] == output_value[i] for i in range(batch_size)])
        ratio = (cnt_id/cnt_ttl).n(50)
        print(f"expected: {expected_ratio} tests: {cnt_ttl:>10} ids: {cnt_id:>5} ratio: {ratio}", end="\r")

#MAIN
if __name__ == "__main__":
    # test(decompositionBN, vBN, matrixBN, constantsBN, primeBN)
    # approximate_identity_fraction_of_bar(decompositionBN, vBN, primeBN)

    prime = 263
    sboxes = [7, 7, 7]
    ν = 5

    new_state = [65, 104, 0]
    new_state = Concrete(new_state, prime, matrixBN, constantsBN[0])
    new_state = Brick(new_state, prime)
    new_state = Concrete(new_state, prime, matrixBN, constantsBN[1])
    new_state = Brick(new_state, prime)
    new_state = Concrete(new_state, prime, matrixBN, constantsBN[2])
    y_0 = bar(new_state, prime, sboxes, ν)

    new_state = [153, 95, 0]
    new_state = Concrete(new_state, prime, matrixBN, constantsBN[0])
    new_state = Brick(new_state, prime)
    new_state = Concrete(new_state, prime, matrixBN, constantsBN[1])
    new_state = Brick(new_state, prime)
    new_state = Concrete(new_state, prime, matrixBN, constantsBN[2])
    y_1 = bar(new_state, prime, sboxes, ν)

    print(y_0, y_1)
    exit()

    for prime in [263]:
        print(prime)
        sboxes = [7, 7, 7]
        ν = 5
        list_of_coll_inputs = []
        list_of_coll_candidates = []
        for x in range(prime):
            for y in range(prime):
                for z in [0]:
                    new_state = [x, y, z]
                    new_state = Concrete(new_state, prime, matrixBN, constantsBN[0])
                    new_state = Brick(new_state, prime)
                    new_state = Concrete(new_state, prime, matrixBN, constantsBN[1])
                    new_state = Brick(new_state, prime)
                    new_state = Concrete(new_state, prime, matrixBN, constantsBN[2])
                    if new_state != bar(bar(new_state, prime, sboxes, ν), prime, sboxes, ν):
                        list_of_coll_inputs += [[x, y, z]]
                        list_of_coll_candidates += [bar(new_state, prime, sboxes, ν)]
        for x in range(prime):
            for y in range(prime):
                for z in [0]:
                    new_state = [x, y, z]
                    new_state = Concrete(new_state, prime, matrixBN, constantsBN[0])
                    new_state = Brick(new_state, prime)
                    new_state = Concrete(new_state, prime, matrixBN, constantsBN[1])
                    new_state = Brick(new_state, prime)
                    new_state = Concrete(new_state, prime, matrixBN, constantsBN[2])
                    for i in range(len(list_of_coll_inputs)):
                        if new_state == list_of_coll_candidates[i]:
                            print([x, y, z])
                            print(list_of_coll_inputs[i])
                            print()
