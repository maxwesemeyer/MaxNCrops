from create_dicts import *
from rasterize import *
from analyse_solution import *
from config import *


########################################################################################################################
# by Maximilian Wesemeyer
# for questions contact wesemema@hu-berlin.de
########################################################################################################################
# THIS IS NOT WORKING RELIABLY YET (especially for new datasets)
# some crop types are hard coded
########################################################################################################################


def run_optimization_seq():
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
    shares_croptypes = pd.read_csv('./' + temp_path + '/shares_iacs_seq.csv')

    with open('./' + temp_path + '/historic_croptypes_dict.pkl', 'rb') as f:
        historic_croptypes_dict = pickle.load(f)

    ####################################################################################################################
    # A bit of a messy workaround if not all crop types are cultivated in a given year; This is the case
    # especially for very small study areas
    all_crops = []
    for i, id in enumerate(unique_field_ids):
        all_crops.append(list(np.unique(np.array(historic_croptypes_dict[id]))))
    unique_crops = np.unique(list(chain.from_iterable(all_crops)))
    unique_crops = np.setdiff1d(unique_crops, [255, 99])
    unique_crops = np.insert(unique_crops, 0, 0)
    ####################################################################################################################

    n_years = len(historic_croptypes_dict[1])

    # the structure of farm_field_dict is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
    with open('./' + temp_path + '/farm_field_dict.pkl', 'rb') as f:
        farm_field_dict = pickle.load(f)
    ####################################################################################################################
    # structure of block_dict
    # block id: [share of fieldid1 in this block, share of fieldid2 in this block,
    # share of fieldid3 in this block, (=number of decision units)...]

    with open('./' + temp_path + '/block_dict.pkl', 'rb') as f:
        block_dict = pickle.load(f)
    ####################################################################################################################
    # Create the Gurobi model
    m = gp.Model("Maximum Landscape Entropy")

    vars = {}
    if verbatim:
        print(len(unique_crops), 'unique crops check all crop types', unique_crops)

    for i, crop in enumerate(unique_crops):
        crop = str(crop)
        vars["{0}".format(crop)] = m.addMVar((len(unique_field_ids), n_years), vtype=GRB.BINARY, name=crop)

    # The helper will be used as variable that is 1 if a crop type exists in a landscape and else is 0
    for crop in unique_crops:
        crop = str(crop)
        vars["{0}".format('helper_' + crop)] = m.addMVar((num_blocks, n_years), vtype=GRB.BINARY, name='helper_' + crop)

    ####################################################################################################################
    # https://math.stackexchange.com/questions/2500415/how-to-write-if-else-statement-in-linear-programming
    # big M constraints to model if/else
    # tempvars['tempvar_' + crop] is 1 if a crop type exists in a landscape and else is 0
    b = 1
    small_number = 0.0001
    M = 1e5

    for i in range(num_blocks):
        if verbatim:
            print(i, '/', num_blocks, 'adding if else constraints')
        tempvars = {}

        indices_block_i = block_dict[i].indices
        for i_crop, crop in enumerate(unique_crops):
            for year in range(n_years):
                crop = str(crop)

                # Tempvar = Sum per block, which is the area of a crop per block; This could be used as input for an
                # Information Entropy function
                tempvars["{0}".format('tempvar_' + crop + str(year))] = m.addVar(vtype=GRB.INTEGER, lb=0, ub=count_pixel_per_block)

                m.addConstr(tempvars['tempvar_' + crop + str(year)] == gp.quicksum(
                    [vars[crop][id_, year] * block_dict[i][0, id_] for id_ in indices_block_i]))
                # if var sum > 1 helper should be one else 0
                m.addConstr(tempvars['tempvar_' + crop + str(year)] >= b + small_number - M * (1 - vars['helper_' + crop][i, year]))
                m.addConstr(tempvars['tempvar_' + crop + str(year)] <= b + M * vars['helper_' + crop][i, year])

    ####################################################################################################################
    # Objective function; We maximize the sum of helper variables for each block (landscape),
    # where each is 1 if a crop exists in a landscape
    obj = gp.quicksum(gp.quicksum(gp.quicksum([vars['helper_' + str(crop)][i, year] for crop in unique_crops]) for i in range(num_blocks)) for year in range(n_years))

    m.addConstrs((sum(sum([vars[str(crop)][i, year] for crop in unique_crops]) for year in range(n_years)) == n_years for i in range(len(unique_field_ids))),
                 name='Only_one_per_year')
    # This constraint ensures that each field has only one crop per year assigned
    m.addConstrs((gp.quicksum([vars[str(crop)][i, year] for crop in unique_crops]) == 1 for i in range(len(unique_field_ids)) for year in range(n_years)),
                 name='Only_one_LU')
    del block_dict

    # this constraint ensures that the no data variable is 1 where it was no data before
    m.addConstr(vars[str(0)][0, :].sum() == n_years, 'nodata_fixed')
    for i, id in enumerate(unique_field_ids):
        if verbatim:
            print(i, ' of ', len(unique_field_ids), 'multiple constraints')
        if i == 0:
            if verbatim:
                print('skipping')
        else:
            try:
                historic_crops = np.array(historic_croptypes_dict[id])
                for crop in unique_crops:
                    # TODO this assumes that 255 is the no data value in the historic crop types raster
                    if crop == 0:
                        # for no data value we fix the position in the sequence
                        crop_na_position = np.where((historic_crops == 255) | (historic_crops == 99), 1, 0)
                        #print(historic_crops, crop, 'sumCrop:', crop_na_position, 'i:', i)

                        m.addConstr(sum(vars[str(crop)][i, :] * crop_na_position) == crop_na_position.sum(),
                                    name='fix_na_position_1_' + str(crop) + '_' + str(id))
                        m.addConstr(sum(vars[str(crop)][i, :] * crop_na_position) == vars[str(crop)][i, :].sum(),
                                    name='fix_na_position_2_' + str(crop) + '_' + str(id))
                    else:
                        # we set the number of occurences per crop in the sequence of each field
                        sum_crop = np.where(historic_crops == crop, True, False).sum()
                        #print(historic_crops, crop, 'sumCrop:', sum_crop, 'i:', i)
                        m.addConstr(vars[str(crop)][i, :].sum() == sum_crop, name='taboo_crop_' + str(crop) + '_' + str(id))
            except:
                None

    del historic_croptypes_dict

    ####################################################################################################################
    # A constraint for each crop type and each farm using the block_dict_farms
    # This constraint ensures that the acreage per crop and farm and year remains similar with a tolerance, which is stated above
    # This constraint is only enforced if diversity_type = 'attainable'
    if diversity_type == 'attainable':
        for ct, farm in enumerate(unique_farms):
            if verbatim:
                print(ct, 'of', len(unique_farms), 'crop proportion per farm constraints')
            if farm == 0 or farm == nd_value:
                if verbatim:
                    print('skipping farm', farm)
                continue
            indices_farm_i = farm_field_dict[farm].indices
            for year in range(n_years):
                for i_crop, crop in enumerate(unique_crops):
                    if crop == 0:
                        continue
                    # try to find the farm in the shares_croptypes dataframe if not possible set thrs to 0
                    try:
                        thrs = (shares_croptypes.loc[(shares_croptypes[crop_type_column] == crop) & (
                                shares_croptypes[farm_id_column] == farm) & (
                                shares_croptypes['year'] == year), 'area_m2'].values[0])
                        # if thrsh is nan the crop was not cultivated by that farm in that year
                        if np.isnan(thrs):
                            thrs = 0
                    except:
                        # This happens when a farm does not cultivate a crop
                        thrs = 0

                    m.addConstr(gp.quicksum(
                        [(farm_field_dict[farm][0, id_]) * vars[str(crop)][id_, year] for id_ in indices_farm_i]) >= thrs - thrs * (
                                        tolerance / 100),
                                name='acreage_' + str(crop) + '_' + str(farm) + '_' + str(year) + '_1')
                    m.addConstr(gp.quicksum(
                        [(farm_field_dict[farm][0, id_]) * vars[str(crop)][id_, year] for id_ in indices_farm_i]) <= thrs + thrs * (
                                        tolerance / 100),
                                name='acreage_' + str(crop) + '_' + str(farm) + '_' + str(year) + '_2' )

    def beet_pot_RotCnstr(model_, vals_, year_, crop_, i_):
        if year_ >= 3:
            if crop_ == 3 or crop_ == 5:
                # we only enfore this specific constraint for (rape =4; potato = 5, beets = 3)
                # Enforce potato and beets constraint only if potato and beets in year AND in any of
                # the preceding two years
                if vals_[str(crop_)][i_, year_] > 0.5 and any(vals_[str(crop_)][i_, year_ - n] > 0.5 for n in range(1, 4)):
                    model_.cbLazy(model_._vars[str(crop)][i_, year_].item() +
                                 model_._vars[str(crop)][i_, year_ - 1].item() +
                                 model_._vars[str(crop)][i_, year_ - 2].item() +
                                 model_._vars[str(crop)][i_, year_ - 3].item() <= 1)

    def rape_RotCnstr(model_, vals_, year_, crop_, i_):
        if year_ >= 2:
            if crop_ == 4:
                # Enforce rapeseed constraint only if rapeseed in year AND in any of the preceding
                # two years
                if vals_[str(crop_)][i_, year_] > 0.5 and any(vals_[str(crop_)][i_, year_ - n] > 0.5 for n in range(1, 3)):
                    model_.cbLazy(model_._vars[str(crop_)][i_, year_].item() +
                                 model_._vars[str(crop_)][i_, year_ - 1].item() +
                                 model_._vars[str(crop_)][i_, year_ - 2].item() <= 1)

    def legume_RotCnstr(model_, vals_, year_, crop_, i_):
        if year_ >= 0:
            if crop_ == 12:
                # Enforce legume constraint only if legumes in two consecutive years
                if vals_[str(crop_)][i_, year_] > 0.5 and vals_[str(crop_)][i_, year_ - 1] > 0.5:
                    # legume constraints
                    model_.cbLazy(model_._vars[str(crop_)][i_, year_].item() +
                                 model_._vars[str(crop_)][i_, year_ - 1].item() <= 1)

    def CropRotRules(model, where):
        if where == gp.GRB.Callback.MIPSOL:
            print('enforcing crop rot rules')
            vals = model.cbGetSolution(model._vars)
            for i_crop, crop in enumerate(unique_crops):
                # implementing the sugar beet constraint first
                for i, id in enumerate(unique_field_ids):
                    for year in range(n_years):
                        #legume_RotCnstr(model, vals, year, crop, i)
                        beet_pot_RotCnstr(model, vals, year, crop, i)
                        rape_RotCnstr(model, vals, year, crop, i)

    ####################################################################################################################
    # default is minimize
    m.setObjective(obj, GRB.MAXIMIZE)
    # in case the model is not feasible try this:

    if m_infeas:
        iis = m.computeIIS()
        m.write('my_iis.ilp')

    m.write('maxent_lp.lp')
    m.params.LazyConstraints = 1
    m.params.Heuristics = 0.9
    print('starting optimization')
    m._vars = vars  # Store variables for use in the callback
    m.optimize(CropRotRules)
    # CropRotRules
    ####################################################################################################################
    # extract the binary solution for each crop from the model m
    out_imgs_allyears = []
    for year in range(n_years):
        out_imgs = []
        for crop in unique_crops:
            solution_arr = vars["{0}".format(crop)].getAttr("x")
            out_imgs.append(solution_arr[:, year])
        out_imgs_allyears.append(out_imgs)

    field_id_arr_allyears = []
    for year in range(n_years):
        fids_list = []
        crop_type_list = []
        for id, indices_field in enumerate(sparse_idx):
            # id = field id
            # indices_field the indices of field field_id
            for original_crop_class, croptype in zip(unique_crops, out_imgs_allyears[year]):
                # only larger than 0.1 will be interpreted as 1;
                if croptype[id] > 0.1:
                    field_id_arr[indices_field] = original_crop_class
                    fids_list.append(id)
                    crop_type_list.append(original_crop_class)
        field_id_arr = field_id_arr.astype(int)
        field_id_arr_allyears.append(field_id_arr)
        opt_frame = pd.DataFrame({'field_id': fids_list, 'OPT_KTYP_' + str(year): crop_type_list})
        iacs_gp = iacs_gp.merge(opt_frame, on='field_id')
    iacs_gp.to_file('./' + out_path + '/iacs_opt.shp')

    write_array_disk_universal(np.stack(field_id_arr_allyears, axis=0), './' + temp_path + '/reference_raster.tif', outPath='./' + out_path + '/maxent_croptypes_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0)
    ####################################################################################################################
    analyse_solution_seq(tolerance=tolerance)
    get_change_map_seq(n_years)
    for year in range(n_years):
        diss_init = iacs_gp.dissolve(by=['crp_yr_' + str(year)], as_index=False)
        diss_init['area_init'] = diss_init.area * 0.0001

        diss_opt = iacs_gp.dissolve(by=['OPT_KTYP_' + str(year)], as_index=False)
        diss_init['area_opt'] = diss_opt.area * 0.0001
        # diss_init = diss_init.drop('geometry')
        diss_init = diss_init[['crp_yr_' + str(year), 'area_init', 'area_opt']]
        diss_init.to_csv('./' + out_path + '/shares_bb_iacs' + str(year) + '.csv')

        diss_init = iacs_gp.dissolve(by=[farm_id_column, 'crp_yr_' + str(year)], as_index=False).copy()
        diss_init['area_init'] = diss_init.area * 0.0001
        diss_opt = iacs_gp.dissolve(by=[farm_id_column, 'OPT_KTYP_' + str(year)], as_index=False).copy()
        diss_opt['area_opt'] = diss_opt.area * 0.0001
        merged = pd.merge(diss_init, diss_opt, left_on=[farm_id_column, 'crp_yr_' + str(year)], right_on=[farm_id_column, 'OPT_KTYP_' + str(year)])

        if diversity_type == 'attainable':
            # Check if farm crop acreage constraint was violated; There should be zero violations;
            error = 0
            for i, row in merged.iterrows():
                # print(row['area_opt'], row['area_init'])
                if (row['area_opt'] > row['area_init'] + row['area_init'] * (tolerance/100)) or (
                        row['area_opt'] < row['area_init'] - row['area_init'] * (tolerance/100)):
                    error += 1
            print('errors: ', error)


if __name__ == '__main__':
    run_optimization_seq()


