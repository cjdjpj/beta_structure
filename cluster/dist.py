import tskit
import pickle
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

pairs = list(combinations(range(mts.num_samples), 2))

# compute
dist = mts.diversity(pairs, mode='site')

# dump to file
with open(args.input + "_dist", 'wb') as file:
    pickle.dump(dist, file)

print("---distance matrix computed---")

### Paralellized ###
# import time
# import tskit
# import pickle
# import argparse
# import math
# from itertools import combinations
# from multiprocessing import Pool
#
# _global_ts = None
#
# def _init_worker(ts_path):
#     global _global_ts
#     _global_ts = tskit.load(ts_path)
#
# def compute_diversity_chunk(args):
#     start, end, pair_list = args
#     return start, _global_ts.diversity(pair_list, mode="site")
#
# def main():
#     p = argparse.ArgumentParser()
#     p.add_argument("--input", default="output")
#     p.add_argument("--threads", type=int, default=1)
#     args = p.parse_args()
#
#     ts = tskit.load(args.input)
#     pairs = list(combinations(range(ts.num_samples), 2))
#     total = len(pairs)
#
#     start_time = time.time()
#     nproc = min(args.threads, total)
#     chunk = math.ceil(total / nproc)
#     tasks = []
#     for i in range(nproc):
#         s = i * chunk
#         e = min(s + chunk, total)
#         tasks.append((s, e, pairs[s:e]))
#
#     with Pool(processes=nproc, initializer=_init_worker, initargs=(args.input,)) as pool:
#         distances = [None] * total
#         for start, dist_chunk in pool.imap_unordered(compute_diversity_chunk, tasks):
#             distances[start : start + len(dist_chunk)] = dist_chunk.tolist()
#
#     end_time = time.time()
#     print(end_time-start_time)
#     with open(args.input + "_dist", "wb") as f:
#         pickle.dump(distances, f)
#
#
#     print(f"{total} diversities across {nproc} processes")
#
#
# if __name__ == "__main__":
#     main()
#
