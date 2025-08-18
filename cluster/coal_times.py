import pickle
import argparse
import tskit

parser = argparse.ArgumentParser(
                    prog='coal_times')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

coal_times = []

i = 0
for tree in mts.trees():

    if i == 100:
        break
    else:
        i+=1

    for u in tree.nodes():
        if tree.is_leaf(u):
            continue
        
        k = len(tree.children(u))
        coal_times.extend([tree.time(u)] * k)


# dump to file
with open(args.input + "_coal_times", 'wb') as file:
    pickle.dump(coal_times, file)

print("---coal_times computed---")
