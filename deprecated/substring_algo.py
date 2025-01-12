import random
import numpy as np


def generate_random_nucleotide_sequence(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))

def mutate_sequence(sequence, m):
    nucleotides = ['A', 'T', 'C', 'G']
    sequence = list(sequence)  # Convert to a list for easy mutation
    for _ in range(m):
        idx = random.randint(0, len(sequence) - 1)  # Random index to mutate
        current_nucleotide = sequence[idx]
        new_nucleotide = random.choice([n for n in nucleotides if n != current_nucleotide])  # Random nucleotide but not the current one
        sequence[idx] = new_nucleotide
    return ''.join(sequence)

s1 = generate_random_nucleotide_sequence(10)
s2 = mutate_sequence(s1, 3)

print(s1)
print(s2)


def prop_identical_blk(s1, s2, m):
    matches = 0
    n = len(s1)
    num_blocks = n-m+1

    i=0
    prev_match = False
    while i <= n-m:
        if (s1[i : i+m] == s2[i : i+m]).all():
            matches += 1
            i += 1
            prev_match = True
        elif prev_match:
            i += m  
        else:
            i += 1

    return matches/num_blocks

m = 2
result = prop_identical_blk(s1, s2, m)
print(f"Proportion of identical overlapping blocks: {result}")
