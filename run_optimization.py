from config import *
from optimize_seq import run_optimization_seq
import os


path = os.path.dirname(__file__)
os.chdir(path)

# id can be useful to test different sets of selected farms
id = 1
temp_path = "temp_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3] + '_' + str(id)
out_path = "output_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3] + '_seq' + '_' + str(id)

print('preparing data to run optimization....')
print('selected farms: ', selected_farm_ids)
run_optimization_seq(selected_farm_ids, temp_path, out_path)