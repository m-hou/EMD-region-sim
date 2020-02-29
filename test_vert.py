from scipy.optimize import linprog
import numpy as np
import random

H = [0.1, 0, 0, 0.9]
eta = 0.2
G = 0.05
num_histogram_to_sample = 1000

def construct_V(H, eta, G):
    V = [H]
    for _ in range(num_histogram_to_sample):
        sample_point = sample_histogram_within_EMD(H, eta, G)
        if not point_within_affine_V(V, sample_point):
            V.append(sample_point)
            V = remove_redundent_points(V)
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
    budget = eta // G

    keep_going = True
    while keep_going:
        move_index, move_direction = random.choice(list(available_move_index_direction))
        granules_to_move = min(random.randint(1, budget), new_H[move_index] // G)
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
    A_eq = np.transpose(np.array([np.array(v + [1]) for v in V]))

    b_eq = np.array(P + [1])

    solution = linprog(c, A_eq=A_eq, b_eq=b_eq)
    return solution.success

def remove_redundent_points(V):
    new_V = []
    for i in range(len(V)):
        test_point = V[i]
        V_prime = V[:i] + V[i+1:]
        if not point_within_affine_V(V_prime, test_point):
            new_V.append(test_point)
    return new_V

print(construct_V(H, eta, G))