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

# Model infeasible? If your model is infeasible set this to True and it will generate a "my_iis.ilp" file
# The Irreducible Inconsistent Subsystem (iis) helps to understand why a model is infeasible
m_infeas = True

# set this to True to print information about the data preprocessing;
verbatim = True

# No data value
nd_value = -999

# creates the paths according to the parameters stated above so we get a new folder when we change the parameters
temp_path = "temp_" + str(agg_length) + '_' + str(tolerance)
out_path = "output_" + str(agg_length) + '_' + str(tolerance)

# Crop rotation rules
# not the same crop as last year, except when the farmer didn't care
# if legumes are assigned we don't change that
# cultivation breaks for potatos, beets and rapeseed (rape =4; potato = 5, beets = 3)
# 1365, 1056, 1062, 1263, 1065, 1052, 1011, 1010, 1115