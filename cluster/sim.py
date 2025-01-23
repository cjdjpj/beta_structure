import re
import json, argparse
import msprime

parser = argparse.ArgumentParser(
                    prog='sim')
parser.add_argument('--output', type=str, default="output")
parser.add_argument('--Ne', type=int, default=25000000)
parser.add_argument('--length', type=int, default=5000000)
parser.add_argument('--track_length', type=int, default=5000)
parser.add_argument('--nsample', type=int, default=50)
parser.add_argument('--mu', type=float, default=0.0000000006)
parser.add_argument('--r_m', type=float, default=0.0)
parser.add_argument('--model', type=str, default="kingman")

args = parser.parse_args()

l = args.length  # number of genes
t = args.track_length  # tract length
nsample = args.nsample  # the number of genomes sampled
mu = args.mu  # mutation rate
r_m = args.r_m # r/m

r = r_m * mu # per base recomb rate
print("r="+str(r))

match = re.match(r"^beta([0-9]*\.[0-9]+)$", args.model)
if match:
    a = float(match.group(1))
    model = msprime.BetaCoalescent(alpha=a)
elif  args.model == "kingman":
    model = None
else:
    raise ValueError(f"Invalid model argument: {args.model}")

ts = msprime.sim_ancestry(nsample,
                          model=model,
                          population_size=args.Ne,
                          ploidy=1,
                          sequence_length=l,
                          gene_conversion_rate=r,
                          gene_conversion_tract_length=t,
                          record_provenance = True,
                          )

print("---ancestry simulation done---")

mts = msprime.sim_mutations(ts, rate=mu)

print("---mutation simulation done---")

# save tree_sequence
mts.dump(args.output)

# save params to json
params_dict = vars(args)
with open(args.output + ".json", 'w') as metadata_file:
    json.dump(params_dict, metadata_file, indent=4)
