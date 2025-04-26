import os
import json
import csv

os.chdir(os.path.dirname(os.path.abspath(__file__)))
runs_json = [f for f in os.listdir() if f.endswith(".json")]
runs_json.sort()

json_data = []

for filename in runs_json:
    with open(filename, "r") as f:
        data = json.load(f)

        # filename
        if "output" in data:
            data["output"] = data["output"].removeprefix("runs/")

        run_index = filename.removesuffix(".json")

        # check for stats
        data["_dist"] = os.path.exists(run_index + "_dist")
        data["_iden_blk"] = os.path.exists(run_index + "_frac_iden_blk")
        data["_clonal"] = os.path.exists(run_index + "_frac_clonal")
        data["_trueclonal"] = os.path.exists(run_index + "_frac_trueclonal")
        
        # add dict
        json_data.append(data)

keys = ['output', 'model', 'alpha', 'r_m', "_dist", "_iden_blk", "_trueclonal", "_clonal", 'pi', 'mu', 'length', 'track_length', 'nsample']

with open("catalogue.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=keys)
    writer.writeheader()
    for row in json_data:
        writer.writerow({key: row.get(key, "") for key in keys})
