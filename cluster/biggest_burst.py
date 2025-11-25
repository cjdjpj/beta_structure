import tskit
import argparse

parser = argparse.ArgumentParser(
                    prog='biggestburst')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

### T OF LARGEST BURST
# For populations which have a coalescent event with arity > MIN_ARITY
# When is the event with largest arity?
trees_to_check = 100
min_arity = 10
highest_arity = 0
highest_arity_T = None
c=0
for tree in mts.trees():
    if c > trees_to_check:
        break
    for node in tree.nodes():
        if tree.num_children(node) >= min_arity and tree.num_children(node) > highest_arity:
            highest_arity = tree.num_children(node)
            highest_arity_T = tree.time(node)
    c+=1

with open(args.input + "_biggestburst", "w") as file:
    file.write(str(highest_arity) + "," + str(highest_arity_T))

print("---biggest burst computed---")
