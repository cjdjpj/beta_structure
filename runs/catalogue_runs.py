import os
import json
import csv

os.chdir(os.path.dirname(os.path.abspath(__file__)))
runs_json = [f for f in os.listdir() if f.endswith(".json")]
runs_json.sort()

json_data = []

for file in runs_json:
    with open(file, "r") as f:
        data = json.load(f)
        if "output" in data:
            data["output"] = data["output"].removeprefix("runs/")
        json_data.append(data)

keys = ['output', 'model', 'alpha', 'r_m', 'Ne', 'pi', 'mu', 'length', 'track_length', 'nsample']

with open("catalogue.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=keys)
    writer.writeheader()
    for row in json_data:
        writer.writerow({key: row.get(key, "") for key in keys})

