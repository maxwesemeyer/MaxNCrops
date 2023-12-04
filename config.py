########################################################################################################################
# by Maximilian Wesemeyer
# for questions contact wesemema@hu-berlin.de
########################################################################################################################


# the agg_length controls the size of the landscape pixels; agg_length of 100 and 10m pixel equals a km²
# Zero is the no data value for farms; farm_id == 0 will be ignored.
agg_length = 100
count_pixel_per_block = agg_length ** 2

# Tolerance value in percent; Controls the allowed deviation per crop and farm
tolerance = 10

# rasterization necessary? can be set to False to speed up the process if run a second time
rasterize = False

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

# creates the paths according to the parameters stated above so we get a new folder when we change the parameters
temp_path = "temp_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3]
out_path = "output_" + str(agg_length) + '_' + str(tolerance) + '_' + diversity_type[0:3]

# crop dictionary, croptype name as key, croptype integer id as value
crop_names_dict = {'no_data': 0, 'maize': 1, 'winter_cereals': 2, 'beets': 3, 'rapeseed': 4, 'potato': 5,
                   'spring_cereals': 6, 'legumes': 12, 'arable_grass': 13, 'sunflowers': 60, 'unknown': 80}
crop_names_dict_reversed = {value: key for key, value in crop_names_dict.items()}

# Crop rotation rules
# not the same crop as last year, except when the farmer didn't care
# if legumes are assigned we don't change that
# cultivation breaks for potatos, beets and rapeseed (rape =4; potato = 5, beets = 3)

selected_farm_ids = [5008]
