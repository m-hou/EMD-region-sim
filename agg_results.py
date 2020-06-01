from collections import defaultdict
from test_vert import construct_V_until_stable
from fractions import Fraction
import matplotlib.pyplot as plt

H = [Fraction(1, 4)] * 4
eta = Fraction(1, 2)
Gs = [Fraction(1, 4), Fraction(1, 8), Fraction(1, 12), Fraction(1, 16), Fraction(1, 20), Fraction(1, 24)]
samples = [1000, 2000, 4000, 8000, 16000, 32000]

def count_zeros(histogram):
    count = 0
    for i in histogram:
        if i == 0:
            count += 1
    return count

def plot_results():
    fig, ax = plt.subplots(
        nrows=1,
        ncols=len(Gs),
        sharex='row',
        sharey='row')

    H_str = ", ".join([f"{h}" for h in H])
    fig.suptitle(f"Number of 0s in vertices (H = [{H_str}], n = {len(H)}, η = {eta})")

    counts = {}
    verts = {}
    for i, G in enumerate(Gs):
        print(f"Plot {i+1}: G = {G}")
        sim_V = construct_V_until_stable(H, eta, G, samples[i], True)
        zero_counts = defaultdict(int)
        for v in sim_V:
            zero_counts[count_zeros(v)] += 1
        counts[G] = zero_counts
        verts[G] = sim_V

    print(counts)
    print(verts)

    for i, key in enumerate(counts):
        col = ax[i]
        xs = counts[key].keys()
        ys = counts[key].values()
        col.bar(xs, ys, width=1)
        col.set_title(f"G = {key}")
        col.set_xticks(range(0, len(H) // 2 + 1))
    plt.show()

plot_results()