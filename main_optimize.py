from create_dicts import *
from rasterize import *
from analyse_solution import *
from config import *

########################################################################################################################
# by Maximilian Wesemeyer
# for questions contact wesemema@hu-berlin.de
########################################################################################################################


def run_optimization():
    if not os.path.exists(temp_path):
        # Create the temp directory if it does not exist
        os.makedirs(temp_path)
    if not os.path.exists(out_path):
        # Create the output directory if it does not exist
        os.makedirs(out_path)
    if rasterize:
        print('rasterizing...')
        rasterize_input_shp(crop_type_column=crop_type_column, farm_id_column=farm_id_column)

    field_id_arr, farmid_arr, sparse_idx, unique_crops, unique_field_ids, iacs_gp, unique_farms, num_blocks, \
    start_vals = prepare_data(agg_length=agg_length, crop_type_column=crop_type_column, farm_id_column=farm_id_column)

    if verbatim:
        print(len(unique_field_ids), 'number of decision units')
    ####################################################################################################################
    shares_croptypes = pd.read_csv('./' + temp_path + '/' + 'shares_iacs.csv')

    with open('./' + temp_path + '/' + 'taboo_croptypes_dict.pkl', 'rb') as f:
        taboo_croptypes_dict = pickle.load(f)
    # the structure of farm_field_dict is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
    with open('./' + temp_path + '/' + 'farm_field_dict.pkl', 'rb') as f:
        farm_field_dict = pickle.load(f)

    ####################################################################################################################
    # structure of block_dict
    # block id: [share of fieldid1 in this block, share of fieldid2 in this block,
    # share of fieldid3 in this block, (=number of decision units)...]

    with open('./' + temp_path + '/' + 'block_dict.pkl', 'rb') as f:
        block_dict = pickle.load(f)
    if verbatim:
        print(len(unique_crops), 'unique crops check all crop types', unique_crops)
    ####################################################################################################################
    # Create the Gurobi model
    m = gp.Model("Maximum Landscape Entropy")

    vars = {}
    for i, crop in enumerate(unique_crops):
        crop = str(crop)
        vars["{0}".format(crop)] = m.addVars(len(unique_field_ids), vtype=GRB.BINARY, name=crop)
        if crop != str(0):
            if verbatim:
                print('setting start', crop)
            vars["{0}".format(crop)].Start = start_vals[i-1]

    #crop_translator = dict(zip(crop_names_dict_reversed.keys(), range(len(crop_names_dict_reversed.keys()))))
    #print('translator', crop_translator)
    for i, row in iacs_gp.iterrows():
        print('iterating i: ', i, row)
        if row[farm_id_column] not in selected_farm_ids:
            print('adding forced constraint', i)
            m.addConstr(vars["{0}".format(row[crop_type_column])][i+1] == 1)

    # The helper will be used as variable that is 1 if a crop type exists in a landscape and else is 0
    for crop in unique_crops:
        crop = str(crop)
        vars["{0}".format('helper_' + crop)] = m.addVars(num_blocks, vtype=GRB.BINARY, name='helper_' + crop)

    ####################################################################################################################
    # https://math.stackexchange.com/questions/2500415/how-to-write-if-else-statement-in-linear-programming
    # big M constraints to model if/else
    # tempvars['tempvar_' + crop] is 1 if a crop type exists in a landscape and else is 0
    b = 1
    small_number = 0.0001
    M = 1e5

    for i in range(num_blocks):
        if verbatim:
            print(i+1, '/', num_blocks, 'adding if else constraints')
        tempvars = {}

        indices_block_i = block_dict[i].indices
        block_empty = False
        if len(indices_block_i) <= 1 and indices_block_i[0] == 0:
            print('continuing')
            block_empty = True

        for i_crop, crop in enumerate(unique_crops):
            crop = str(crop)
            if block_empty:
                m.addConstr(vars['helper_' + crop][i] == 0)
                continue
            # Tempvar = Sum per block, which is the area of a crop per block; This could be used as input for an
            # Information Entropy function
            tempvars["{0}".format('tempvar_' + crop + str(i))] = m.addVar(vtype=GRB.INTEGER, lb=0, ub=count_pixel_per_block)

            m.addConstr(tempvars['tempvar_' + crop + str(i)] == gp.quicksum(
                [vars[crop][id_] * block_dict[i][0, id_] for id_ in indices_block_i]))
            # EXAMPLE:
            # tempvar = 5 ; b = 10
            # 5 >= 10 + 0.1 - 1e5 -> helper must be 0
            # 5 <= 10 + 0.1 + 1e5 / 0 -> helper can be 1 or 0
            # tempvar = 20 ; b = 10
            # 20 >= 10 + 0.1 - 0 -> helper can be 0 or 1
            # 20 <= 10 + 0.1 + 100000 -> helper must be 1

            m.addConstr(tempvars['tempvar_' + crop + str(i)] >= b + small_number - M * (1 - vars['helper_' + crop][i]))
            m.addConstr(tempvars['tempvar_' + crop + str(i)] <= b + M * vars['helper_' + crop][i])

    ####################################################################################################################
    # Objective function; We maximize the sum of helper variables for each block (landscape),
    # where each is 1 if a crop exists in a landscape
    obj = gp.quicksum(gp.quicksum([vars['helper_' + str(crop)][i] for crop in unique_crops]) for i in range(num_blocks))

    # This constraint ensures that each field has only one crop assigned
    m.addConstrs((gp.quicksum([vars[str(str(crop))][i] for crop in unique_crops]) == 1 for i in range(len(unique_field_ids))),
                 name='Only_one_LU')
    del block_dict

    # this constraint ensures that the no data variable is 1 where it was no data before
    m.addConstr(vars[str(0)][0] == 1, 'nodata_fixed')

    for i, id in enumerate(unique_field_ids):
        if verbatim:
            print(i+1, ' of ', len(unique_field_ids), 'taboo constraints')
        if i == 0:
            if verbatim:
                print('skipping')
        else:
            # this constraint ensures that the nodata class cannot be applied to agricultural parcels
            m.addConstr(vars[str(0)][i] == 0, 'nodata_fixed_2_' + str(i))
            try:
                taboo_crops = taboo_croptypes_dict[id]
                for crop in unique_crops:
                    if crop in taboo_crops:
                        # this constraint ensures that the field specific taboo crop cannot be applied to the field
                        m.addConstr(vars[str(crop)][i] == 0, 'taboo_crop_' + str(crop) + '_' + str(id))
            except:
                None

    del taboo_croptypes_dict

    ####################################################################################################################
    # A constraint for each crop type and each farm using the block_dict_farms
    # This constraint ensures that the acreage per crop and farm remains similar with a tolerance, which is stated above
    # This constraint is only enforced if diversity_type = 'attainable'
    if diversity_type == 'attainable':
        for ct, farm in enumerate(unique_farms):
            if verbatim:
                print(ct+1, 'of', len(unique_farms), 'crop composition per farm constraints')
            if farm == 0 or farm == nd_value:
                if verbatim:
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
                    # This happens when a farm does not cultivate a crop
                    thrs = 0

                m.addConstr(gp.quicksum(
                    [(farm_field_dict[farm][0, id_]) * vars[str(crop)][id_] for id_ in indices_farm_i]) >= thrs - thrs * (
                                    tolerance / 100), '{0}'.format(crop) + '_' + str(farm) + '_1')
                m.addConstr(gp.quicksum(
                    [(farm_field_dict[farm][0, id_]) * vars[str(crop)][id_] for id_ in indices_farm_i]) <= thrs + thrs * (
                                    tolerance / 100), '{0}'.format(crop) + '_' + str(farm) + '_2' )

    ####################################################################################################################
    # default is minimize
    m.setObjective(obj, GRB.MAXIMIZE)

    m._vars = vars
    #m.write('maxent_lp.lp')
    m.optimize()

    if m.status == gp.GRB.Status.INFEASIBLE:
        # in case the model is not feasible
        # it will generate a "my_iis.ilp" file
        # The Irreducible Inconsistent Subsystem (iis) helps to understand why a model is infeasible
        # this can easily happen when farmers the crop rotation rules are too strict (because farmers violated them)
        iis = m.computeIIS()
        m.write('my_iis.ilp')
    ####################################################################################################################
    # extract the binary solution for each crop from the model m
    out_imgs = []
    for crop in unique_crops:
        #sol = np.array(m.getAttr("X", vars["{0}".format(str(crop))]).values()[0])
        sol = np.array(list(m.getAttr("X", vars["{0}".format(str(crop))]).values()))
        out_imgs.append(sol)

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
    iacs_gp.to_file('./' + out_path + '/' + 'iacs_opt.shp')
    write_array_disk_universal(np.expand_dims(field_id_arr, axis=0), './' + temp_path + '/' + 'reference_raster.tif', outPath='./' + out_path + '/' + 'opt_crop_allocation_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0)
    ####################################################################################################################
    analyse_solution()
    get_change_map()
    diss_init = iacs_gp.dissolve(by=[crop_type_column], as_index=False)
    diss_init['area_init'] = diss_init.area * 0.0001

    diss_opt = iacs_gp.dissolve(by=['OPT_KTYP'], as_index=False)
    diss_init['area_opt'] = diss_opt.area * 0.0001
    # diss_init = diss_init.drop('geometry')
    diss_init = diss_init[[crop_type_column, 'area_init', 'area_opt']]
    diss_init.to_csv('' + out_path + '/' + 'shares_bb_iacs.csv')

    diss_init = iacs_gp.dissolve(by=[farm_id_column, crop_type_column], as_index=False).copy()
    diss_init['area_init'] = diss_init.area * 0.0001

    diss_opt = iacs_gp.dissolve(by=[farm_id_column, 'OPT_KTYP'], as_index=False).copy()
    diss_opt['area_opt'] = diss_opt.area * 0.0001
    merged = pd.merge(diss_init, diss_opt, left_on=[farm_id_column, 'ID_KTYP'], right_on=[farm_id_column, 'OPT_KTYP'])

    if diversity_type == 'attainable':
        # Check if farm crop acreage constraint was violated; There should be zero violations;
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


