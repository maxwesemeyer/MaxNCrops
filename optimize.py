from create_dicts_script import *
from rasterize_script import *
from analyse_solution import *

"""
by Maximilian Wesemeyer

"""
############################################################################
# the agg_length controls the size of the landscape pixels; agg_length of 100 and 10m pixel equals a km²
# Zero is the no data value for farms; farm_id == 0 will be ignored.
agg_length = 100
count_pixel_per_block = agg_length ** 2

# Tolerance value in percent; Controls the allowed deviation per crop and farm
tolerance = 1

# rasterization necessary? can be set to False to speed up the process if run a second time
rasterize = True

# state here the column names in the Shapefile
crop_type_column = 'ID_KTYP'
farm_id_column = 'farm_id'

# What diversity to calculate? chose either 'attainable' or 'potential'
# If 'potential' is selected, there are no farm acreage constraints
diversity_type = 'attainable'

# Model infeasible? If your model is infeasible set this to True and it will generate a "my
m_infeas = False
############################################################################


def run_optimization():
    if not os.path.exists("temp"):
        # Create the temp directory if it does not exist
        os.makedirs("temp")
    if not os.path.exists("output"):
        # Create the output directory if it does not exist
        os.makedirs("output")
    if rasterize:
        print('rasterizing...')
        rasterize_input_shp(crop_type_column=crop_type_column, farm_id_column=farm_id_column)
    crop_arr, field_id_arr, farmid_arr, sparse_idx, \
    unique_crops, unique_field_ids, iacs_gp, unique_farms, field_id_arr_1d, \
    size_of_test_data, len_raster, num_blocks, num_block_y, start_vals = prepare_data(agg_length=agg_length, crop_type_column=crop_type_column, farm_id_column=farm_id_column)

    print(len(unique_field_ids), 'number of decision units')

    shares_croptypes = pd.read_csv('./temp/shares_iacs.csv')

    with open('./temp/historic_croptypes_dict.pkl', 'rb') as f:
        historic_croptypes_dict = pickle.load(f)

    ############################################################################
    with open('./temp/farm_field_dict.pkl', 'rb') as f:
        farm_field_dict = pickle.load(f)
    # the shape is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
    ############################################################################
    # block id: [share of fieldid1 in this block, share of fieldid2 in this block,
    # share of fieldid3 in this block, (=number of decision units)...]

    with open('./temp/block_dict.pkl', 'rb') as f:
        block_dict = pickle.load(f)
    ####################################################################################################################
    # Create the Gurobi model
    m = gp.Model("Maximum Landscape Entropy")

    vars = {}
    print(len(unique_crops), 'unique crops check all crop types', unique_crops)
    for i, crop in enumerate(unique_crops):
        crop = str(crop)
        vars["{0}".format(crop)] = m.addVars(len(unique_field_ids), vtype=GRB.BINARY, name=crop)
        if crop != str(0):
            print('setting start', crop)
            vars["{0}".format(crop)].Start = start_vals[i-1]

    # The helper will be used as variable that is 1 if a crop type exists in a landscape and else is 0
    for crop in unique_crops:
        crop = str(crop)
        vars["{0}".format('helper_' + crop)] = m.addVars(num_blocks, vtype=GRB.BINARY, name='helper_' + crop)

    ############################################################################
    # https://math.stackexchange.com/questions/2500415/how-to-write-if-else-statement-in-linear-programming
    # big M constraints
    b = 1
    small_number = 0.0001
    M = 1e5

    for i in range(num_blocks):
        print(i, '/', num_blocks, 'adding if else constraints')
        tempvars = {}

        indices_block_i = block_dict[i].indices
        for i_crop, crop in enumerate(unique_crops):
            crop = str(crop)

            # Tempvar = Sum per block, which is the area of a crop per block; This could be used as input for an
            # Information Entropy function
            tempvars["{0}".format('tempvar_' + crop)] = m.addVar(vtype=GRB.INTEGER, lb=0, ub=count_pixel_per_block)

            m.addConstr(tempvars['tempvar_' + crop] == gp.quicksum(
                [vars[crop][id_] * block_dict[i][0, id_] for id_ in indices_block_i]))
            # if var sum > 1 helper should be one else 0
            m.addConstr(tempvars['tempvar_' + crop] >= b + small_number - M * (1 - vars['helper_' + crop][i]))
            m.addConstr(tempvars['tempvar_' + crop] <= b + M * vars['helper_' + crop][i])

    ############################################################################
    # this is the Objective function; We maximize the helper variable, which is 1 if a crop exists in a landscape
    obj = gp.quicksum(gp.quicksum([vars['helper_' + str(crop)][i] for crop in unique_crops]) for i in range(num_blocks))

    m.addConstrs((gp.quicksum([vars[str(str(crop))][i] for crop in unique_crops]) == 1 for i in range(len(unique_field_ids))),
                 name='Only_one_LU')
    del block_dict

    # this constraint ensures that the no data variable is 1 where it was no data before
    m.addConstr(vars[str(0)][0] == 1, 'nodata_fixed')

    for i, id in enumerate(unique_field_ids):
        print(i, ' of ', len(unique_field_ids), 'multiple constraints')
        if i == 0:
            print('skipping')
        else:
            # this constraint ensures that the nodata class cannot be applied to agricultural parcels
            m.addConstr(vars[str(0)][i] == 0, 'nodata_fixed_2' + str(i))
            try:
                taboo_crops = historic_croptypes_dict[id]
                for crop in unique_crops:
                    if crop in taboo_crops and crop < 99:
                        m.addConstr(vars[str(crop)][i] == 0, 'taboo_crop_' + str(crop) + '_' + str(id))
            except:
                None

    del historic_croptypes_dict

    ############################################################################
    # A constraint for each crop type and each farm using the block_dict_farms
    # This constraint ensures that the acreage per crop and farm remains similar with a tolerance, which is stated above
    # This constraint is only enforced if diversity_type = 'attainable'
    if diversity_type == 'attainable':
        for ct, farm in enumerate(unique_farms):
            print(ct, 'of', len(unique_farms), 'crop proportion per farm constraints')
            if farm == 0 or farm == -999.0:
                print('skipping farm', farm)
                continue
            indices_farm_i = farm_field_dict[farm].indices
            for i_crop, crop in enumerate(unique_crops):
                if crop == 0:
                    continue
                # try to find the farm in the shares_croptypes dataframe if not possible set thrs to 0
                try:
                    thrs = (shares_croptypes.loc[(shares_croptypes[crop_type_column] == crop) & (
                            shares_croptypes[farm_id_column] == farm), 'area_m2'].values[0])

                except:
                    thrs = 0

                m.addConstr(gp.quicksum(
                    [(farm_field_dict[farm][0, id_]) * vars[str(crop)][id_] for id_ in indices_farm_i]) >= thrs - thrs * (
                                    tolerance / 100), '{0}'.format(crop) + '_' + str(farm) + '_1')
                m.addConstr(gp.quicksum(
                    [(farm_field_dict[farm][0, id_]) * vars[str(crop)][id_] for id_ in indices_farm_i]) <= thrs + thrs * (
                                    tolerance / 100), '{0}'.format(crop) + '_' + str(farm) + '_2' )

    # in case the model is not feasible try this:
    if m_infeas:
        iis = m.computeIIS()
        m.write('my_iis.ilp')

    # default is minimize
    m.setObjective(obj, GRB.MAXIMIZE)
    #m.write('maxent_lp.lp')
    m.optimize()


    ############################################################################
    # extract the binary solution for each crop from the model m
    out_imgs = []
    for crop in unique_crops:
        sol = m.getAttr("X", vars["{0}".format(str(crop))]).values()
        out_imgs.append(sol)
    all_vars = m.getVars()
    #values = m.getAttr("X", all_vars)
    #names = m.getAttr("VarName", all_vars)
    #print(names)

    fids_list = []
    crop_type_list = []
    for id, indices_field in enumerate(sparse_idx):
        # id = field id
        # indices_field the indices of field field_id
        for original_crop_class, croptype in zip(unique_crops, out_imgs):
            # only larger than 0.1 will be interpreted as 1;
            if croptype[id] > 0.1:
                field_id_arr[indices_field] = original_crop_class
                fids_list.append(id)
                crop_type_list.append(original_crop_class)
    field_id_arr = field_id_arr.astype(int)
    opt_frame = pd.DataFrame({'field_id': fids_list, 'OPT_KTYP': crop_type_list})
    iacs_gp = iacs_gp.merge(opt_frame, on='field_id')
    iacs_gp.to_file('./output/iacs_opt.shp')
    write_array_disk_universal(np.expand_dims(field_id_arr, axis=0), './temp/reference_raster.tif', outPath='./output/maxent_croptypes_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0)

    analyse_solution(tolerance=tolerance)

    diss_init = iacs_gp.dissolve(by=[crop_type_column], as_index=False)
    diss_init['area_init'] = diss_init.area * 0.0001

    diss_opt = iacs_gp.dissolve(by=['OPT_KTYP'], as_index=False)
    diss_init['area_opt'] = diss_opt.area * 0.0001
    # diss_init = diss_init.drop('geometry')
    diss_init = diss_init[[crop_type_column, 'area_init', 'area_opt']]
    diss_init.to_csv('output/shares_bb_iacs.csv')

    diss_init = iacs_gp.dissolve(by=[farm_id_column, crop_type_column], as_index=False).copy()
    diss_init['area_init'] = diss_init.area * 0.0001
    print(diss_init)
    diss_opt = iacs_gp.dissolve(by=[farm_id_column, 'OPT_KTYP'], as_index=False).copy()
    diss_opt['area_opt'] = diss_opt.area * 0.0001
    merged = pd.merge(diss_init, diss_opt, left_on=['farm_id', 'ID_KTYP'], right_on=['farm_id', 'OPT_KTYP'])

    if diversity_type == 'attainable':
        # Check if farm crop acreage constraint was violated
        error = 0
        for i, row in merged.iterrows():
            # print(row['area_opt'], row['area_init'])
            if (row['area_opt'] > row['area_init'] + row['area_init'] * (tolerance/100)) or (
                    row['area_opt'] < row['area_init'] - row['area_init'] * (tolerance/100)):
                #print(error)
                error += 1
        print('errors: ', error)


if __name__ == '__main__':
    run_optimization()


