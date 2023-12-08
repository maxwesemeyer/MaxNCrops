from create_dicts import *
from rasterize import *
from analyse_solution import *
from config import *
from CropRotRules import *

########################################################################################################################
# by Maximilian Wesemeyer
# for questions contact wesemema@hu-berlin.de
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
    # A bit of a workaround if not all crop types are cultivated in a given year; This is the case
    # mostly for very small study areas
    all_crops = []
    for i, id in enumerate(unique_field_ids):
        all_crops.append(list(np.unique(np.array(historic_croptypes_dict[id]))))
    unique_crops = np.unique(list(chain.from_iterable(all_crops)))
    unique_crops = np.setdiff1d(unique_crops, [255, 99])
    #unique_crops = np.insert(unique_crops, 0, 0)

    ####################################################################################################################
    n_years = len(historic_croptypes_dict[1])
    # this dictionary contains 1 for each Crop rotation constraint the farmers violated themselves
    CropRotViolation_dict, longest_seq_dict = check_CropRotRules(historic_croptypes_dict)
    crop_rot_freq = get_rotations(historic_croptypes_dict)

    relevant_fields_list = []
    for i, id in enumerate(unique_field_ids):
        if i == 0:
            continue
        # we skip all fields of farms that are not selected because the rotation on these
        # fields of these farms are fixed
        selected_row = iacs_gp[iacs_gp['field_id'] == id]
        if selected_row[farm_id_column].iloc[0] in selected_farm_ids:
            relevant_fields_list.append(id)
    crop_rot_freq = crop_rot_freq[crop_rot_freq['fid'].isin(relevant_fields_list)]

    crop_rot_freq.to_csv('./' + out_path + '/crop_rot_freq_init.csv')
    del crop_rot_freq
    print(longest_seq_dict)
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
        vars["{0}".format(crop)] = m.addMVar((len(unique_field_ids), n_years), vtype=GRB.BINARY, name='crop_bin_' + crop)

    # setting those farms that are not selected as static (same values as in historic crop types)
    selected_farms_field_ids = []
    for i, row in iacs_gp.iterrows():
        if row[farm_id_column] not in selected_farm_ids:
            for year in range(n_years):
                crop_type_column_yr = 'crp_yr_' + str(year)
                # setting those farms that are not selected as static (same values as in historic crop types)
                m.addConstr(vars["{0}".format(row[crop_type_column_yr])][i+1, year] == 1)
        else:
            selected_farms_field_ids.append(row['field_id'])

    # The helper will be used as variable that is 1 if a crop type exists in a landscape and else is 0
    for crop in unique_crops:
        crop = str(crop)
        vars["{0}".format('helper_' + crop)] = m.addMVar((num_blocks, n_years), vtype=GRB.BINARY, name='helper_' + crop)
    ####################################################################################################################
    # select all fields that are in the same blocks as the selected farms
    # the relevant fields is a list of all field ids that are located within a landscape, where one of the selected farms
    # has a field
    relevant_block_list = []
    for i in range(num_blocks):
        indices_block_i = block_dict[i].indices

        if any(value in selected_farms_field_ids for value in indices_block_i):
            relevant_block_list.append(i)

    ####################################################################################################################
    # https://math.stackexchange.com/questions/2500415/how-to-write-if-else-statement-in-linear-programming
    # big M constraints to model if/else
    # tempvars['tempvar_' + crop] is 1 if a crop type exists in a landscape and else is 0
    b = 1
    small_number = 0.0001
    M = 1e5

    for i in range(num_blocks):
        if verbatim and i % 1000 == 0:
            print(i, '/', num_blocks, 'adding landscape-level constraints')
        tempvars = {}

        indices_block_i = block_dict[i].indices
        block_irrelevant = False
        if len(indices_block_i) <= 1 and indices_block_i[0] == 0:
            # if no fields are located in the block
            block_irrelevant = True
        if i not in relevant_block_list:
            # if no relevant fields of the selected farms are in the block
            block_irrelevant = True
        for i_crop, crop in enumerate(unique_crops):
            crop = str(crop)
            if block_irrelevant:
                m.addConstr(gp.quicksum([vars['helper_' + crop][i, year] for year in range(n_years)]) == 0)
                continue
            for year in range(n_years):
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
        if i == 0:
            continue
        selected_row = iacs_gp[iacs_gp['field_id'] == id]
        if selected_row[farm_id_column].iloc[0] in selected_farm_ids:
            if verbatim:
                print(i, ' of ', len(unique_field_ids), 'field-level constraints')
            if i == 0:
                if verbatim:
                    print('skipping')
            else:
                try:
                    historic_crops = np.array(historic_croptypes_dict[id])
                    for crop in unique_crops:
                        if crop == 0:
                            # for no data value we fix the position in the sequence
                            #crop_na_position = np.where((historic_crops == 255) | (historic_crops == 99), 1, 0)
                            crop_na_position = np.where(historic_crops == 0, 1, 0)

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

    #del historic_croptypes_dict

    ####################################################################################################################
    # A constraint for each crop type and each farm using the block_dict_farms
    # This constraint ensures that the acreage per crop and farm and year remains similar with a tolerance, which is stated above
    # This constraint is only enforced if diversity_type = 'attainable'
    if diversity_type == 'attainable':
        for ct, farm in enumerate(unique_farms):
            # we skip all farms that are not selected because the fields of these farms are fixed
            if farm in selected_farm_ids:
                if verbatim:
                    print(ct, 'of', len(unique_farms), 'adding farm-level constraints')
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

    ####################################################################################################################

    def max_seq_x_RotCnstr(model_, vals_, crop_, i_, x):
        # Currently works only for x == maize or x == winter cereals
        if crop_ == x:
            # index 1 == cereals, index 0 == maize
            if x == crop_names_dict['winter_cereals']:
                longest_farmer_seq_i = longest_seq_dict[i_][1]
            elif x == crop_names_dict['maize']:
                longest_farmer_seq_i = longest_seq_dict[i_][0]

            bin_seq = [1 if vals_[str(crop_)][i_, iter] > 0.5 else 0 for iter in range(n_years)]
            max_optimized_seq_length, max_start_index = longest_sequence(bin_seq)
            if max_optimized_seq_length > longest_farmer_seq_i:
                # the sum of the longest sequence needs to be smaller or
                # equal to the max_optimized_seq_length - 1, which means we introduce a 0 in it
                # EXAMPLE: [1, 1, 1, 1] and longest sequence is supposed to be 2 in this case we cannot set the sum
                # of it to 2 because this case would be ok too: [1, 1, 0, 1] (sum = 3); so we set a constraint that
                # introduces a zero by setting the sum to max_optimized_seq_length - 1; we do this until constraint is not violated
                #print('field: ', i_)
                #print('4:', bin_seq, max_optimized_seq_length, max_start_index, longest_farmer_seq_i)
                model_.cbLazy(gp.quicksum([model_._vars[str(crop_)][i_, max_start_index + indexo].item() for indexo
                                                   in range(max_optimized_seq_length)]) <= max_optimized_seq_length - 1)

    def no_x_after_y_RotCnstr(model_, vals_, year_, crop_, i_, x, y):
        # EXAMPLE: we can reformulate any no x before y to no x after y:
        # no sunflowers after legumes [12, 60] x = sunflowers , y = legumes
        # no sunflowers before legumes = [60, 12] == no legumes after sunflowers x = legumes, y = sunflowers
        # can also be used when x == y; (e.g. no sunflowers after sunflowers)
        if year_ <= n_years-2:
            if crop_ == y:
                if x == y:
                    # if farmers cultivated y more than 3 times in 7 years we cannot enfore the constraint
                    if sum(vals_[str(y)][i_, iter] > 0.5 for iter in range(n_years)) > 3:
                        return

                if vals_[str(y)][i_, year_] > 0.5 and vals_[str(x)][i_, year_ + 1] > 0.5:
                    model_.cbLazy(model_._vars[str(y)][i_, year_].item() +
                                  model_._vars[str(x)][i_, year_ + 1].item() <= 1)


    def min_return_t_for_x(model_, vals_, year_, crop_, i_, t, x):
        # we do t-1 because then setting t is then more intuitive. Setting t = 2 means that a crop can be grown every
        # second year [1, 0, 1]
        t = t -1
        # adapt t to n_years, so we don't index more than the list length
        t_adapted_to_n_years = t
        while t_adapted_to_n_years + year_ > n_years-1:
            t_adapted_to_n_years = t_adapted_to_n_years - 1

        if crop_ == x:
            if vals_[str(x)][i_, year_] > 0.5 and any(vals_[str(x)][i_, year_ + add] > 0.5 for add in range(1, t_adapted_to_n_years + 1)):
                model_.cbLazy(gp.quicksum([model_._vars[str(crop_)][i_, year_ + add].item()
                                               for add in range(t_adapted_to_n_years + 1)]) <= 1)


    def CropRotRules_lazy(model, where):
        if where == gp.GRB.Callback.MIPSOL:
            print('enforcing crop rot rules')
            vals = model.cbGetSolution(model._vars)
            for i_crop, crop in enumerate(unique_crops):
                # implementing the sugar beet constraint first
                for i, id in enumerate(unique_field_ids):
                    if i == 0:
                        continue
                    # we skip all fields of farms that are not selected because the rotation on these
                    # fields of these farms are fixed
                    selected_row = iacs_gp[iacs_gp['field_id'] == id]
                    if selected_row[farm_id_column].iloc[0] in selected_farm_ids:

                        # the cereal constraint needs the entire sequence as input
                        max_seq_x_RotCnstr(model, vals, crop, i, x=crop_names_dict['winter_cereals'])
                        for year in range(n_years):
                            # no legumes after legumes
                            if CropRotViolation_dict[i][3] == 0:
                                # we enforce this constraint only it was respected in the initial solution on that field
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                  x=crop_names_dict['legumes'], y=crop_names_dict['legumes'])
                            if CropRotViolation_dict[i][8] == 0:
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                  x=crop_names_dict['sunflowers'], y=crop_names_dict['legumes'])
                            if CropRotViolation_dict[i][4] == 0:
                                # no rapeseed after sunflowers
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                  x=crop_names_dict['rapeseed'], y=crop_names_dict['sunflowers'])
                                # no sunflowers after rapeseed
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                  x=crop_names_dict['sunflowers'], y=crop_names_dict['rapeseed'])

                            if CropRotViolation_dict[i][9] == 0:
                                # no rapeseed after beets
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                      x=crop_names_dict['rapeseed'], y=crop_names_dict['beets'])
                                # no beets after rapeseed
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                      x=crop_names_dict['beets'], y=crop_names_dict['rapeseed'])
                            if CropRotViolation_dict[i][7] == 0:
                                # no legumes after rapeseed
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                      x=crop_names_dict['legumes'], y=crop_names_dict['rapeseed'])
                            if CropRotViolation_dict[i][6] == 0:
                                # no rapeseed after maize
                                no_x_after_y_RotCnstr(model, vals, year, crop, i,
                                                      x=crop_names_dict['rapeseed'], y=crop_names_dict['maize'])


                            # max return time for rapeseed = 3 years; beets and potato = 4
                            # max return time of 3 means  this is ok: [1, 0, 0, 1, 0, 0, 1]
                            if CropRotViolation_dict[i][0] == 0:
                                min_return_t_for_x(model, vals, year, crop, i, t=3, x=crop_names_dict['rapeseed'])
                            else:
                                # set crop sequence to original sequence
                                if historic_croptypes_dict[i][year] == crop:
                                    model.cbLazy(gp.quicksum([model._vars[str(crop)][i, year].item()]) == 1)
                            if CropRotViolation_dict[i][1] == 0:
                                min_return_t_for_x(model, vals, year, crop, i, t=4, x=crop_names_dict['potato'])
                            if CropRotViolation_dict[i][2] == 0:
                                min_return_t_for_x(model, vals, year, crop, i, t=4, x=crop_names_dict['beets'])
                            if CropRotViolation_dict[i][5] == 0:
                                min_return_t_for_x(model, vals, year, crop, i, t=4, x=crop_names_dict['sunflowers'])

    ####################################################################################################################
    # default is minimize
    m.setObjective(obj, GRB.MAXIMIZE)

    #m.write('maxent_lp.lp')
    m.params.LazyConstraints = 1
    m.params.Heuristics = 0.3
    print('starting optimization')
    m._vars = vars  # Store variables for use in the callback
    # CropRotRules
    m.optimize(CropRotRules_lazy)
    m.write('maxent_lp.lp')

    if m.status == gp.GRB.Status.INFEASIBLE:
        # in case the model is not feasible
        # it will generate a "my_iis.ilp" file
        # The Irreducible Inconsistent Subsystem (iis) helps to understand why a model is infeasible
        # this can easily happen when farmers the crop rotation rules are too strict (because farmers violated them)
        iis = m.computeIIS()
        m.write('my_iis.ilp')
        #constraint_to_check_relaxed = m.getConstrByName(constraint_to_check)
    ####################################################################################################################
    # extract the binary solution for each crop from the model m
    field_id_arr_allyears = []
    for year in range(n_years):
        fids_list = []
        crop_type_list = []
        field_id_arr_copy = field_id_arr.copy()
        for id, indices_field in enumerate(sparse_idx):
            for crop_counter, original_crop_class in enumerate(unique_crops):
                    # id = field id
                    # indices_field the indices of field field_id
                    # only larger than 0.1 will be interpreted as 1;
                    if m._vars[str(original_crop_class)].X[id, year].item() > 0.1:
                        field_id_arr_copy[indices_field] = original_crop_class
                        fids_list.append(id)
                        crop_type_list.append(original_crop_class)
        field_id_arr_copy = field_id_arr_copy.astype(int)
        field_id_arr_allyears.append(field_id_arr_copy)
        opt_frame = pd.DataFrame({'field_id': fids_list, 'OPT_KTYP_' + str(year): crop_type_list})
        iacs_gp = iacs_gp.merge(opt_frame, on='field_id')

    iacs_gp['sel_farm'] = iacs_gp[farm_id_column].apply(lambda val: 1 if val in selected_farm_ids else 0)
    iacs_gp.to_file('./' + out_path + '/iacs_opt.shp')

    write_array_disk_universal(np.stack(field_id_arr_allyears, axis=0), './' + temp_path + '/reference_raster.tif', outPath='./' + out_path + '/maxent_croptypes_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0)

    ####################################################################################################################
    analyse_solution_seq()
    get_change_map_seq(n_years)
    get_shares_seq(iacs_gp, n_years)

    ####################################################################################################################
    # calculate historic_croptypes_dict and check for violations of the rules again
    taboo_croptypes_dict, historic_croptypes_dict = get_historic_croptypes(field_id_array=field_id_arr.copy(), historic_croptypes_array=np.stack(field_id_arr_allyears, axis=0), unique_croptypes=unique_crops)
    CropRotViolation_dict, longest_seq_dict_opt = check_CropRotRules(historic_croptypes_dict)
    crop_rot_freq = get_rotations(historic_croptypes_dict)

    crop_rot_freq = crop_rot_freq[crop_rot_freq['fid'].isin(relevant_fields_list)]
    crop_rot_freq.to_csv('./' + out_path + '/crop_rot_freq_' + str(tolerance) + '.csv')
    crop_rot_figures()
    for i in range(len(unique_field_ids)):
        if longest_seq_dict_opt[i][1] > longest_seq_dict[i][1]:
            print('error', longest_seq_dict[i], longest_seq_dict_opt[i], i)
    print(longest_seq_dict)
    print(longest_seq_dict_opt)


if __name__ == '__main__':
    run_optimization_seq()


