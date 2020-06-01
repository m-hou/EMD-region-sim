from scipy.optimize import linprog
import numpy as np
import random
from fractions import Fraction

def load_V_from_cache(H, eta, G):
    # Todo
    return [H]

def construct_V(H, eta, G, num_histogram_to_sample, log=False, use_cache=True):
    V = load_V_from_cache(H, eta, G) if use_cache else [H]
    for i in range(num_histogram_to_sample):
        if log:
            print(i)
        sample_point = sample_histogram_within_EMD(H, eta, G)
        if not point_within_affine_V(V, sample_point):
            V.append(sample_point)
            V = remove_redundent_points(V)
    return V

def construct_V_until_stable(H, eta, G, iters_to_stable, log=False, use_cache=True):
    V = load_V_from_cache(H, eta, G) if use_cache else [H]
    stable_iters = 0
    while stable_iters < iters_to_stable:
        if log:
            print(stable_iters)
        sample_point = sample_histogram_within_EMD(H, eta, G)
        if not point_within_affine_V(V, sample_point):
            V.append(sample_point)
            V = remove_redundent_points(V)
            stable_iters = 0
        stable_iters += 1
    return V

def sample_histogram_within_EMD(H, eta, G):

    LEFT, RIGHT = 1, -1

    def reverse_move(move_index, move_direction):
        return (move_index + move_direction , -move_direction)

    def generate_possible_moves(len_H):
        left_moves = set((i, LEFT) for i in range(len_H - 1))
        right_moves = set(reverse_move(*move) for move in left_moves)
        return left_moves | right_moves

    available_move_index_direction = generate_possible_moves(len(H))

    new_H = list(H)
    budget = round(eta / G)

    keep_going = True
    while keep_going:
        move_index, move_direction = random.choice(list(available_move_index_direction))
        granules_to_move = min(random.randint(1, budget), round(new_H[move_index] / G))
        apply_move_to_histogram(new_H, move_index, move_direction, granules_to_move * G)

        available_move_index_direction.remove((move_index, move_direction))
        available_move_index_direction.remove(reverse_move(move_index, move_direction))
        budget -= granules_to_move
        keep_going = available_move_index_direction and budget > 0 and random.choice([True, False])
    
    return new_H

def apply_move_to_histogram(histogram, move_index, move_direction, mass_to_move):
    histogram[move_index] -= mass_to_move
    histogram[move_index + move_direction] += mass_to_move

def point_within_affine_V(V, P):
    c = np.repeat(1, len(V))
    A_eq = np.transpose(np.array([np.array(v) for v in V]))
    b_eq = np.array(P)
    A_ub = np.transpose(np.array([np.array([1]) for _ in V]))
    b_ub = np.array([1])

    solution = linprog(c, A_ub, b_ub, A_eq, b_eq)
    return solution.success

def remove_redundent_points(V):
    new_V = []
    for i in range(len(V)):
        test_point = V[i]
        V_prime = V[:i] + V[i+1:]
        if not point_within_affine_V(V_prime, test_point):
            new_V.append(test_point)
    return new_V

def test_V(V, points):
    return all([point_within_affine_V(V, point) for point in points])
