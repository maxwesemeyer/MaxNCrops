########################################################################################################################
# by Maximilian Wesemeyer
# for questions contact wesemema@hu-berlin.de
########################################################################################################################


# the agg_length controls the size of the landscape pixels; agg_length of 100 and 10m pixel equals a km²
# Zero is the no data value for farms; farm_id == 0 will be ignored.
agg_length = 100
count_pixel_per_block = agg_length ** 2

# Tolerance value in percent; Controls the allowed deviation per crop and farm
tolerance = 5

# rasterization necessary? can be set to False to speed up the process if run a second time
rasterize = True

# state here the column names in the Shapefile
crop_type_column = 'ID_KTYP_2'
farm_id_column = 'new_farm_i' #farm_id

# What diversity to calculate? chose either 'attainable' or 'potential'
# If 'potential' is selected, there are no farm acreage constraints
diversity_type = 'attainable'

# set this to True to print information about the data preprocessing;
verbatim = True

# No data value
nd_value = -999

# Should the crop sequences be considered or not?

seq = True
temp_path = "temp_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3]

if seq:
    # creates the paths according to the parameters stated above, so we get a new folder when we change the parameters
    out_path = "output_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3] + '_seq'
else:
    # oy = one year
    out_path = "output_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3] + '_oy'


# crop dictionary, croptype name as key, croptype integer id as value
crop_names_dict = {'no_data': 0, 'maize': 1, 'winter_cereals': 2, 'beets': 3, 'rapeseed': 4, 'potato': 5,
                   'spring_cereals': 6, 'legumes': 12, 'arable_grass': 13, 'sunflowers': 60, 'unknown': 80}
crop_names_dict_reversed = {value: key for key, value in crop_names_dict.items()}
# here you can select all farms of which the sequences will be optimized;
import numpy as np
import pandas as pd
np.random.seed(13) # 42
farm_ids_uckermark = np.unique(pd.read_csv('./delete/uckermark_farmids.csv')[farm_id_column])
selected_farm_ids = np.random.choice(farm_ids_uckermark, size=10, replace=False)
selected_farm_ids = [5392, 5517, 5462, 5322, 5461] # uckermark
selected_farm_ids = [1983, 2203, 1836, 1980, 1977] # oder
selected_farm_ids = [5140, 5085, 5165, 5176, 5161] # teltow
print(selected_farm_ids)