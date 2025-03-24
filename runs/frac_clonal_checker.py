import pickle
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
runs_json = [f for f in os.listdir() if f.endswith(".json")]
runs_json.sort()

json_data = []

for filename in runs_json:
    with open(filename, "r") as f:
        run_index = filename.removesuffix(".json")

        if os.path.exists(run_index + "_frac_clonal"):
            with open(run_index + "_frac_clonal", "rb") as file:
                clonal_tmrca = pickle.load(file)
                print(run_index + ": " +  str(len(clonal_tmrca)))
        else:
            print(run_index + ": " + "lacking")

