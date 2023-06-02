# the agg_length controls the size of the landscape pixels; agg_length of 100 and 10m pixel equals a km²
# Zero is the no data value for farms; farm_id == 0 will be ignored.
agg_length = 100
count_pixel_per_block = agg_length ** 2

# Tolerance value in percent; Controls the allowed deviation per crop and farm
tolerance = 10

# rasterization necessary? can be set to False to speed up the process if run a second time
rasterize = True

# state here the column names in the Shapefile
crop_type_column = 'ID_KTYP'
farm_id_column = 'farm_id'

# What diversity to calculate? chose either 'attainable' or 'potential'
# If 'potential' is selected, there are no farm acreage constraints
diversity_type = 'attainable'

# Model infeasible? If your model is infeasible set this to True and it will generate a "my_iis.ilp" file
# The Irreducible Inconsistent Subsystem (iis) helps to understand why a model is infeasible
m_infeas = False

# set this to True to print information about the data preprocessing;
verbatim = False

# No data value
nd_value = -999